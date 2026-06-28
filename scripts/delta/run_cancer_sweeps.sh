#!/usr/bin/env bash
# Cancer parameter sweeps (report §4 item "#2") — drive the cancer pipelines across
# ONE axis at a time and write a per-axis manifest for collect_cancer_sweeps.sh.
#
# Axes (set AXIS=...):
#   purity    Somatic recall vs tumor purity, with the matched normal PINNED at a
#             fixed coverage (NORMAL_FIX) via total = NORMAL_FIX/(1-purity). This is
#             the clean re-run of the §3.9 finding: the earlier collapse at purity
#             0.9 came from splitting a FIXED budget, which starved the normal to
#             ~6×. Pinning the normal removes that confound, so if recall no longer
#             collapses, the collapse was experimental design, not rneat. (Tumor-
#             sample depth then rises with purity, which only reinforces this.)
#   coverage  Somatic recall vs total tumor depth, at fixed purity.
#   tmr       Somatic recall/precision + truth count vs tumor somatic mutation rate.
#   tissue    Per-tissue SV spectra (BRCA/skin/lung) via sv_pipeline.sbatch + the
#             bundled per-tissue SV models — do tissue-specific spectra (BRCA DUP-,
#             skin BND-, lung DEL-leaning) survive downstream to Manta/Delly calls?
#
# cancer_pipeline.sbatch / sv_pipeline.sbatch are idempotent, so re-running an axis
# only does the unfinished points. Each point is its own SLURM job.
#
# Usage (from repo root):
#   AXIS=purity   bash scripts/delta/run_cancer_sweeps.sh        # start here
#   AXIS=coverage bash scripts/delta/run_cancer_sweeps.sh
#   AXIS=tmr      bash scripts/delta/run_cancer_sweeps.sh
#   AXIS=tissue   bash scripts/delta/run_cancer_sweeps.sh
#   # overridable: PURITIES, COVERAGES, TMRS, TISSUES, NORMAL_FIX, REFERENCE, PURITY, TOTAL_COVERAGE
# When an axis's jobs finish:
#   bash scripts/delta/collect_cancer_sweeps.sh <printed manifest path>
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
source "$REPO_ROOT/scripts/delta/lib_report.sh"   # $SCRATCH + RESULTS_DIR resolution

AXIS="${AXIS:?set AXIS=purity|coverage|tmr|tissue}"
REFERENCE="${REFERENCE:-$SCRATCH/neat_data/chr22.fa}"
MANIFEST="${MANIFEST:-${RESULTS_DIR:-$HOME}/cancer_sweep_${AXIS}.manifest}"
CANCER="$REPO_ROOT/scripts/delta/cancer_pipeline.sbatch"
SVPIPE="$REPO_ROOT/scripts/delta/sv_pipeline.sbatch"
mkdir -p "$(dirname "$MANIFEST")"
echo "# axis point jobid outdir total_cov purity tmr model" > "$MANIFEST"

# REPS>1 submits each point REPS times with a DISTINCT per-rep seed into <out>_rep<r>,
# for mean±sd error bars (see docs/replication_audit.md). REPS=1 = original behavior.
REPS="${REPS:-1}"

# submit_point <extra manifest cols> — caller sets PIPELINE, POINT, OUTBASE and the
# pipeline env array PIPE_ENV=(KEY=val ...). Writes one manifest row per rep; uses a
# distinct RNG_SEED_ROOT per rep when REPS>1 (both pipelines honor it).
submit_point() {
    local extra="$1" rr o lbl jid
    for rr in $(seq 1 "$REPS"); do
        if [[ "$REPS" -gt 1 ]]; then
            o="${OUTBASE}_rep${rr}"; lbl="${POINT}_rep${rr}"
            jid=$(env "${PIPE_ENV[@]}" OUTDIR="$o" PRUNE="${PRUNE:-0}" RNG_SEED_ROOT="cancer $AXIS $POINT rep$rr" \
                  sbatch --parsable "$PIPELINE")
        else
            o="$OUTBASE"; lbl="$POINT"
            jid=$(env "${PIPE_ENV[@]}" OUTDIR="$o" PRUNE="${PRUNE:-0}" sbatch --parsable "$PIPELINE")
        fi
        echo "$AXIS $lbl $jid $o $extra" | tee -a "$MANIFEST"
    done
}

case "$AXIS" in
  purity)
    NORMAL_FIX="${NORMAL_FIX:-30}"               # matched-normal coverage, held constant (pinned mode)
    FIXED_TOTAL="${FIXED_TOTAL:-}"               # if set: §3.7-style FIXED budget (normal shrinks)
    PURITIES="${PURITIES:-0.3 0.5 0.7 0.9}"
    if [[ -n "$FIXED_TOTAL" ]]; then             # distinct manifest so it never clobbers the pinned run
        MANIFEST="${MANIFEST%.manifest}_fixed.manifest"
        echo "# axis point jobid outdir total_cov purity tmr model" > "$MANIFEST"
    fi
    PIPELINE="$CANCER"
    for p in $PURITIES; do
        if [[ -n "$FIXED_TOTAL" ]]; then
            total="$FIXED_TOTAL"; OUTBASE="$SCRATCH/sweep_purity_p${p}_fixed${FIXED_TOTAL}"
        else
            total=$(awk -v n="$NORMAL_FIX" -v p="$p" 'BEGIN{printf "%d", (n/(1-p))+0.5}')
            OUTBASE="$SCRATCH/sweep_purity_p${p}_n${NORMAL_FIX}"
        fi
        POINT="$p"; PIPE_ENV=(REFERENCE="$REFERENCE" TOTAL_COVERAGE="$total" PURITY="$p")
        submit_point "$total $p - pancancer"
    done
    if [[ -n "$FIXED_TOTAL" ]]; then
        echo "NOTE: FIXED budget ${FIXED_TOTAL}× (normal=(1-purity)·total, shrinks) — the §3.7-style 'before' arm."
    else
        echo "NOTE: matched normal pinned at ${NORMAL_FIX}× (total=NORMAL_FIX/(1-purity)) — the 'after' arm."
    fi
    ;;
  coverage)
    PURITY="${PURITY:-0.6}"
    COVERAGES="${COVERAGES:-30 60 100 200}"
    PIPELINE="$CANCER"
    for c in $COVERAGES; do
        OUTBASE="$SCRATCH/sweep_cov_c${c}_p${PURITY}"; POINT="$c"
        PIPE_ENV=(REFERENCE="$REFERENCE" TOTAL_COVERAGE="$c" PURITY="$PURITY")
        submit_point "$c $PURITY - pancancer"
    done
    ;;
  tmr)
    PURITY="${PURITY:-0.6}"; TOTAL_COVERAGE="${TOTAL_COVERAGE:-60}"
    TMRS="${TMRS:-1e-6 1e-5 5e-5 1e-4}"
    PIPELINE="$CANCER"
    for rate in $TMRS; do
        OUTBASE="$SCRATCH/sweep_tmr_r${rate}_c${TOTAL_COVERAGE}_p${PURITY}"; POINT="$rate"
        PIPE_ENV=(REFERENCE="$REFERENCE" TOTAL_COVERAGE="$TOTAL_COVERAGE" PURITY="$PURITY" TUMOR_MUTATION_RATE="$rate")
        submit_point "$TOTAL_COVERAGE $PURITY $rate pancancer"
    done
    ;;
  tissue)
    PURITY="${PURITY:-0.8}"; TOTAL_COVERAGE="${TOTAL_COVERAGE:-60}"
    SV_RATE_SCALE="${SV_RATE_SCALE:-1.5}"
    TISSUES="${TISSUES:-BRCA skin lung}"
    PIPELINE="$SVPIPE"
    for t in $TISSUES; do
        model="$REPO_ROOT/tools/cosmic_pancancer_sv_${t}.json.gz"
        [[ -f "$model" ]] || { echo "missing model: $model — skipping $t" >&2; continue; }
        OUTBASE="$SCRATCH/sweep_tissue_${t}"; POINT="$t"
        PIPE_ENV=(REFERENCE="$REFERENCE" TOTAL_COVERAGE="$TOTAL_COVERAGE" PURITY="$PURITY" \
                  SV_RATE_SCALE="$SV_RATE_SCALE" TUMOR_MODEL="$model" RUN_DELLY=1 RUN_CNV=1)
        submit_point "$TOTAL_COVERAGE $PURITY - $model"
    done
    ;;
  *) echo "unknown AXIS=$AXIS (expected purity|coverage|tmr|tissue)" >&2; exit 1 ;;
esac

echo
echo "Manifest: $MANIFEST   (watch: squeue --me)"
echo "When jobs finish: bash scripts/delta/collect_cancer_sweeps.sh $MANIFEST"

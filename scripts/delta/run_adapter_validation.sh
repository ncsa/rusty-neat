#!/usr/bin/env bash
# Adapter-readthrough validation (#125): does enabling 3' sequencing-adapter
# readthrough keep the data realistic AND keep it callable by a standard germline
# pipeline? Submits the germline_e2e variant-caller pipeline (rneat → BWA-MEM2 →
# GATK HaplotypeCaller → hap.py) under three conditions, each replicated REPS
# times for mean±sd, then collect_adapter_validation.sh tabulates and flags any
# statistical difference vs the no-adapter baseline.
#
#   off       adapters disabled (baseline = the path the regression already gates)
#   on_raw    adapters enabled, reads aligned RAW (exposes soft-clipped adapter tails)
#   on_trim   adapters enabled, fastp-trimmed before alignment (realistic Illumina)
#
# Paired by seed: rep r uses the SAME rng_seed across all three conditions, so the
# truth variant set is identical (adapters don't touch variant generation — proven
# in-tree) and differences are attributable to adapter handling, not seed draw.
#
# Usage (from repo root, on Delta after setup.sh):
#   bash scripts/delta/run_adapter_validation.sh
#   REPS=5 COVERAGE=50 ADAPTER_PRESET=nextera \
#     REFERENCE=$SCRATCH/neat_data/chr22.fa \
#     bash scripts/delta/run_adapter_validation.sh
#
# Re-runnable: germline_e2e is idempotent per outdir, so re-submitting only
# recomputes missing stages. When all jobs finish:
#   bash scripts/delta/collect_adapter_validation.sh $MANIFEST
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
source "$REPO_ROOT/scripts/delta/lib_report.sh"

REFERENCE="${REFERENCE:-$SCRATCH/neat_data/chr22.fa}"
COVERAGE="${COVERAGE:-30}"
REPS="${REPS:-3}"
ADAPTER_PRESET="${ADAPTER_PRESET:-truseq}"
# SHORT-INSERT library — REQUIRED for this test to mean anything. 3' adapter
# readthrough only occurs when the insert is shorter than the read length
# (generate_fragments.rs:76 keeps frag < read_len only when adapters=1), so with
# the germline default FRAG_MEAN=350 >> read_len=151 NO readthrough ever happens
# and off/on_raw/on_trim produce identical reads (the 2026-07-01 null run). At
# 180±60 with 151 bp reads ~30% of inserts fall below the read length → real
# readthrough for the `on` arms, while ~63% stay >=161 so the `off` baseline can
# still fill its fragment pool.
FRAG_MEAN="${FRAG_MEAN:-180}"
FRAG_SD="${FRAG_SD:-60}"
MANIFEST="${MANIFEST:-${RESULTS_DIR:-$HOME}/adapter_validation.manifest}"

[[ -f "$REFERENCE" ]] || { echo "REFERENCE not found: $REFERENCE" >&2; exit 1; }
mkdir -p "$(dirname "$MANIFEST")"
echo "# condition rep jobid outdir reference preset" > "$MANIFEST"

# condition -> "ADAPTERS TRIM" env values
submit() {
    local cond="$1" adapters="$2" trim="$3" rep="$4"
    local out="$SCRATCH/adapter_${cond}_rep${rep}"
    local jid
    jid=$(REFERENCE="$REFERENCE" TOOLS=rneat COVERAGE="$COVERAGE" OUTDIR="$out" \
          SEED="adapter-validation rep$rep" \
          ADAPTERS="$adapters" TRIM="$trim" ADAPTER_PRESET="$ADAPTER_PRESET" \
          FRAG_MEAN="$FRAG_MEAN" FRAG_SD="$FRAG_SD" \
          MEASURE_REALISM=1 PRUNE="${PRUNE:-1}" LOCAL_STAGE="${LOCAL_STAGE:-1}" \
          sbatch --parsable "$REPO_ROOT/scripts/delta/germline_e2e.sbatch")
    echo "$cond $rep $jid $out $REFERENCE $ADAPTER_PRESET" | tee -a "$MANIFEST"
}

echo "Adapter validation: ref=$(basename "$REFERENCE") cov=${COVERAGE}x reps=$REPS preset=$ADAPTER_PRESET frag=${FRAG_MEAN}±${FRAG_SD} (short-insert → readthrough active)"
for r in $(seq 1 "$REPS"); do
    submit off     0 0 "$r"
    submit on_raw  1 0 "$r"
    submit on_trim 1 1 "$r"
done

echo
echo "Manifest: $MANIFEST   (watch: squeue --me)"
echo "When all finish: bash scripts/delta/collect_adapter_validation.sh $MANIFEST"

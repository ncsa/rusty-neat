#!/usr/bin/env bash
# Input-variety validation (report §4 item 4): rneat germline fidelity across a
# size + composition range of genomes, to confirm SNP/indel recall and Ts/Tv are
# NOT overtuned to human/chr22. For each FASTA in $GENOMES it submits germline_e2e
# (rneat only — fast; the NEAT arm isn't needed here) into variety_<name>/;
# collect_variety.sh then tabulates the per-genome metrics for a chart.
#
# Stage the genomes in $SCRATCH first (the job builds the BWA index + faidx).
# Usage (from repo root):
#   GENOMES="$SCRATCH/neat_data/ecoli.fa $SCRATCH/neat_data/yeast.fa \
#            $SCRATCH/neat_data/chr22.fa" \
#     bash scripts/delta/run_input_variety.sh
#   # add more once scp'd:  ... c_elegans.fa rice.fa
#
# Re-runnable: a genome whose variety_<name>/rneat_scored.summary.csv already
# exists is skipped by germline_e2e's own idempotency, so expanding the list and
# re-running only does the new genomes.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
source "$REPO_ROOT/scripts/delta/lib_report.sh"

GENOMES="${GENOMES:?set GENOMES to a space-separated list of reference FASTAs}"
COVERAGE="${COVERAGE:-30}"
# REPS>1 submits each genome REPS times with a DISTINCT seed into variety_<g>_rep<r>,
# for mean±sd error bars (see docs/archive/replication_audit.md). REPS=1 = original behavior.
REPS="${REPS:-1}"
MANIFEST="${MANIFEST:-${RESULTS_DIR:-$HOME}/variety.manifest}"
mkdir -p "$(dirname "$MANIFEST")"
echo "# name jobid outdir reference" > "$MANIFEST"

for ref in $GENOMES; do
    if [[ ! -f "$ref" ]]; then echo "MISSING: $ref — skipping" >&2; continue; fi
    base="$(basename "$ref")"; base="${base%.*}"
    for r in $(seq 1 "$REPS"); do
        if [[ "$REPS" -gt 1 ]]; then
            name="${base}_rep${r}"; out="$SCRATCH/variety_${base}_rep${r}"
            jid=$(REFERENCE="$ref" TOOLS=rneat COVERAGE="$COVERAGE" OUTDIR="$out" \
                  SEED="variety $base rep$r" \
                  sbatch --parsable "$REPO_ROOT/scripts/delta/germline_e2e.sbatch")
        else
            name="$base"; out="$SCRATCH/variety_$base"   # unchanged default seed
            jid=$(REFERENCE="$ref" TOOLS=rneat COVERAGE="$COVERAGE" OUTDIR="$out" \
                  sbatch --parsable "$REPO_ROOT/scripts/delta/germline_e2e.sbatch")
        fi
        echo "$name $jid $out $ref" | tee -a "$MANIFEST"
    done
done

echo
echo "Manifest: $MANIFEST   (watch: squeue --me)"
echo "When all finish: bash scripts/delta/collect_variety.sh $MANIFEST"

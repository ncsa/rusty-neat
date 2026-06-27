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
MANIFEST="${MANIFEST:-${RESULTS_DIR:-$HOME}/variety.manifest}"
mkdir -p "$(dirname "$MANIFEST")"
echo "# name jobid outdir reference" > "$MANIFEST"

for ref in $GENOMES; do
    if [[ ! -f "$ref" ]]; then echo "MISSING: $ref — skipping" >&2; continue; fi
    name="$(basename "$ref")"; name="${name%.*}"
    out="$SCRATCH/variety_$name"
    jid=$(REFERENCE="$ref" TOOLS=rneat COVERAGE="$COVERAGE" OUTDIR="$out" \
          sbatch --parsable "$REPO_ROOT/scripts/delta/germline_e2e.sbatch")
    echo "$name $jid $out $ref" | tee -a "$MANIFEST"
done

echo
echo "Manifest: $MANIFEST   (watch: squeue --me)"
echo "When all finish: bash scripts/delta/collect_variety.sh $MANIFEST"

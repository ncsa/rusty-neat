#!/usr/bin/env bash
# Merge sharded whole-genome simulation outputs into one dataset.
#
# After the genome_array.sbatch array finishes, concatenates the per-shard
# outputs under SHARD_OUTROOT into whole-genome FASTQ + truth VCF:
#   - FASTQ: raw concat (gzip members concatenate; read names carry contig+coords
#            and shards are disjoint, so names are globally unique — no rename).
#   - VCF:   bcftools concat (disjoint, ordered shards) then sort, bgzip, index.
#
# Run on a login node or as a dependent job (sbatch --dependency=afterok:<arrayjobid>).
# Needs bcftools (bioinf env). Verifies every expected shard is present first.
#
# Usage:
#   SHARD_OUTROOT=$SCRATCH/wg_<arrayjobid> bash scripts/delta/merge_shards.sh
#   (optional) OUT_PREFIX=$SCRATCH/wg_<id>/merged  EXPECT_SHARDS=<N>

set -euo pipefail

REPO_ROOT="${RNEAT_REPO:-$(cd "$(dirname "$0")/../.." && pwd)}"
source "$REPO_ROOT/scripts/delta/lib_report.sh"   # $SCRATCH + setup_conda/conda_activate

SHARD_OUTROOT="${SHARD_OUTROOT:?set SHARD_OUTROOT to the array OUTROOT (e.g. \$SCRATCH/wg_<jobid>)}"
OUT_PREFIX="${OUT_PREFIX:-$SHARD_OUTROOT/merged}"

setup_conda
conda_activate bioinf   # bcftools

mapfile -t SHARDS < <(find "$SHARD_OUTROOT" -maxdepth 1 -type d -name 'shard_*' | sort)
[[ "${#SHARDS[@]}" -gt 0 ]] || { echo "no shard_* dirs under $SHARD_OUTROOT" >&2; exit 1; }
echo "Found ${#SHARDS[@]} shard dirs under $SHARD_OUTROOT"

# Completeness check: every shard must have both outputs, and the count must
# match EXPECT_SHARDS if given (catches a partially-failed array before merge).
missing=0
r1_list=(); vcf_list=()
for d in "${SHARDS[@]}"; do
    if [[ -f "$d/s_r1.fastq.gz" && -f "$d/s_r2.fastq.gz" && -f "$d/s.vcf.gz" ]]; then
        r1_list+=("$d/s_r1.fastq.gz"); vcf_list+=("$d/s.vcf.gz")
    else
        echo "  INCOMPLETE: $d" >&2; missing=$(( missing + 1 ))
    fi
done
if (( missing > 0 )); then
    echo "ERROR: $missing shard(s) incomplete — rerun those array tasks before merging." >&2
    exit 1
fi
if [[ -n "${EXPECT_SHARDS:-}" && "${#SHARDS[@]}" -ne "$EXPECT_SHARDS" ]]; then
    echo "ERROR: found ${#SHARDS[@]} shards, expected $EXPECT_SHARDS." >&2
    exit 1
fi

echo "[1/2] Concatenating FASTQ (${#r1_list[@]} shards)..."
cat "${r1_list[@]}"                         > "${OUT_PREFIX}_r1.fastq.gz"
cat "${SHARDS[@]/%//s_r2.fastq.gz}"         > "${OUT_PREFIX}_r2.fastq.gz"

echo "[2/2] Concatenating + sorting truth VCF..."
for v in "${vcf_list[@]}"; do [[ -f "$v.tbi" || -f "$v.csi" ]] || bcftools index -f -t "$v"; done
# concat (disjoint regions) → sort → bgzip → index. -a allows any order.
bcftools concat -a -Ou "${vcf_list[@]}" \
    | bcftools sort -Oz -o "${OUT_PREFIX}.vcf.gz"
bcftools index -f -t "${OUT_PREFIX}.vcf.gz"

reads=$(( $(zcat "${OUT_PREFIX}_r1.fastq.gz" | wc -l) / 4 ))
vars=$(bcftools view -H "${OUT_PREFIX}.vcf.gz" | wc -l)
echo
echo "════════════════════════════════════════════════════════════════"
echo "Merged whole-genome dataset:"
echo "  ${OUT_PREFIX}_r1.fastq.gz / _r2.fastq.gz   ($reads read pairs)"
echo "  ${OUT_PREFIX}.vcf.gz                        ($vars truth variants)"
echo "════════════════════════════════════════════════════════════════"

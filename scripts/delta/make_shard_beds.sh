#!/usr/bin/env bash
# Partition a reference genome into anchor-windows for sharded simulation.
#
# Writes one BED per shard (shard_0000.bed, shard_0001.bed, …) into OUTDIR, each
# a contiguous window of the genome (large contigs split at SHARD_BP boundaries).
# Prints the shard COUNT to stdout — feed it to `sbatch --array=0-$((N-1))`.
#
# rneat's target_bed restricts read GENERATION to these regions in ABSOLUTE
# reference coordinates and anchors each fragment in exactly one window (reads
# may extend past the window edge; ownership is by anchor), so the union of all
# shards reconstructs a whole-genome run with no double-counting or gaps — and
# outputs merge by simple concatenation (read names carry contig+coords).
#
# Usage:
#   bash scripts/delta/make_shard_beds.sh REF.fa.fai [SHARD_BP] OUTDIR
#   N=$(bash scripts/delta/make_shard_beds.sh $SCRATCH/neat_data/GRCh38.fa.fai 25000000 $SCRATCH/shards)
#   sbatch --array=0-$((N-1)) ... scripts/delta/genome_array.sbatch
#
# SHARD_BP default 25,000,000 (25 Mb) → ~124 shards for GRCh38 primary assembly.

set -euo pipefail

FAI="${1:?usage: make_shard_beds.sh REF.fa.fai [SHARD_BP] OUTDIR}"
SHARD_BP="${2:-25000000}"
OUTDIR="${3:?need OUTDIR}"
[[ -f "$FAI" ]] || { echo "fai not found: $FAI (run: samtools faidx REF.fa)" >&2; exit 1; }

mkdir -p "$OUTDIR"
rm -f "$OUTDIR"/shard_*.bed

idx=0
# fai columns: name<TAB>length<TAB>offset<TAB>linebases<TAB>linewidth
while IFS=$'\t' read -r name len _rest; do
    [[ -n "$name" && "$len" =~ ^[0-9]+$ ]] || continue
    start=0
    while (( start < len )); do
        end=$(( start + SHARD_BP )); (( end > len )) && end="$len"
        printf '%s\t%d\t%d\n' "$name" "$start" "$end" \
            > "$(printf '%s/shard_%04d.bed' "$OUTDIR" "$idx")"   # BED: 0-based, half-open
        idx=$(( idx + 1 ))
        start="$end"
    done
done < "$FAI"

echo "$idx"   # shard count → use as `--array=0-$((idx-1))`

#!/usr/bin/env bash
# Drive replicated (and K-swept) whole-genome simulation runs on EXCLUSIVE Delta
# nodes to measure the OPTIMAL, reproducible wall-clock — the clean counterpart
# to a naive shared-partition run (which we measured at 3h41m, contention-bound).
#
# For each shards-per-node value K in $KS and each replicate r in 1..$REPS it
# submits one genome_array.sbatch array (exclusive nodes, K windows/node) into
# its own OUTROOT with a distinct seed, and records (K, rep, jobid, nodes,
# outroot) to a manifest. The distinct seeds mean the reps are independent draws
# — so they bound wall-clock variance AND exercise data variety (anti-overtuning)
# at the same time. After the jobs finish, aggregate_reps.sh turns the manifest
# into a wall-clock table (mean +/- sd per K) for report section 3.6.
#
# Usage (from repo root, after setup.sh + make_shard_beds.sh):
#   SHARD_DIR=$SCRATCH/shards_hg38 REFERENCE=$SCRATCH/neat_data/GRCh38.fa \
#     KS="4 8 16" REPS=3 MAXNODES=20 TAG=hg38 \
#     bash scripts/delta/run_genome_reps.sh
#
# Then, once `squeue --me` is clear:
#   bash scripts/delta/aggregate_reps.sh $SCRATCH/wg_sweep_hg38.manifest
#
# To produce the final merged whole-genome dataset, pick the fastest clean rep
# and run merge_shards.sh on its OUTROOT (timing reps don't need merging).
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
source "$REPO_ROOT/scripts/delta/lib_report.sh"

REFERENCE="${REFERENCE:-$SCRATCH/neat_data/GRCh38.fa}"
SHARD_DIR="${SHARD_DIR:-$SCRATCH/shards_hg38}"
KS="${KS:-16 32 64}"         # shards-per-node values to sweep (fixed build packs ~99% to 64)
REPS="${REPS:-3}"            # replicates per K
MAXNODES="${MAXNODES:-20}"   # concurrent exclusive nodes (the %throttle)
TAG="${TAG:-hg38}"
COVERAGE="${COVERAGE:-30}"
SV_RATE_SCALE="${SV_RATE_SCALE:-0}"

[[ -f "$REFERENCE" ]] || { echo "reference not found: $REFERENCE" >&2; exit 1; }
N=$(ls "$SHARD_DIR"/shard_*.bed 2>/dev/null | wc -l)
(( N > 0 )) || { echo "no shard_*.bed under $SHARD_DIR — run make_shard_beds.sh first" >&2; exit 1; }

# Manifest lives on DURABLE storage, never scratch — scratch can throw transient
# EIO / hit quota mid-run, and the driver must not half-submit because a tiny
# bookkeeping write failed. OUTROOTs (the bulky FASTQ) still go to scratch.
MANIFEST="${MANIFEST:-${RESULTS_DIR:-$HOME}/rneat-sweeps/wg_sweep_${TAG}.manifest}"
mkdir -p "$(dirname "$MANIFEST")"
{ echo "# K rep jobid nodes outroot"; } > "$MANIFEST"
echo "shards=$N  KS=[$KS]  reps=$REPS  maxnodes=$MAXNODES  tag=$TAG"
# Rough footprint warning: ~35-40 GB FASTQ per whole-genome 30x run.
echo "NOTE: ~$(( $(echo "$KS" | wc -w) * REPS * 38 )) GB of FASTQ total across all (K,rep) runs — check scratch quota."

for K in $KS; do
    T=$(( (N + K - 1) / K ))               # node-tasks = ceil(N / K)
    # Per-task walltime. On the fixed build (#340) packing is ~99% efficient to
    # 64/node — windows run in ~1-2 min even packed deep, no superlinear slowdown —
    # so a flat 1 h is generous headroom (the old 30+20K ask scaled multi-hour on a
    # false premise and just hurt scheduling). Override with WALL_MINUTES if needed.
    MINUTES="${WALL_MINUTES:-60}"
    WALL=$(printf '%02d:%02d:00' $(( MINUTES / 60 )) $(( MINUTES % 60 )))
    for r in $(seq 1 "$REPS"); do
        OUTROOT="$SCRATCH/wg_${TAG}_k${K}_rep${r}"
        jid=$(SHARD_DIR="$SHARD_DIR" OUTROOT="$OUTROOT" REFERENCE="$REFERENCE" \
              SHARDS_PER_NODE="$K" COVERAGE="$COVERAGE" SV_RATE_SCALE="$SV_RATE_SCALE" \
              SEED_ROOT="rneat wg ${TAG} k${K} rep${r} shard" \
              sbatch --parsable --time="$WALL" --array=0-$((T-1))%${MAXNODES} \
                     "$REPO_ROOT/scripts/delta/genome_array.sbatch")
        echo "$K $r $jid $T $OUTROOT" >> "$MANIFEST"
        echo "  K=$K rep=$r -> job $jid  ($T node-tasks, walltime $WALL)  $OUTROOT"
    done
done

echo
echo "Manifest: $MANIFEST"
echo "Watch:    squeue --me"
echo "Aggregate when clear: bash scripts/delta/aggregate_reps.sh $MANIFEST"

#!/bin/bash
# SLURM array job: sharded whole-genome simulation on Delta.
#
# Each array task simulates ONE genome window (shard_<taskid>.bed from
# make_shard_beds.sh) with a SINGLE-THREADED rneat process. This is rneat's fast
# path on Delta: in-process threading regresses on the dual-socket EPYC (memory-
# bandwidth/allocator bound), but the work is embarrassingly parallel ACROSS
# processes, and SLURM spreads array tasks over many nodes so each gets
# independent memory bandwidth. Reads/variants come out in absolute genome
# coordinates (target_bed is a generation-time filter), so merge_shards.sh
# concatenates the per-shard outputs into a whole-genome dataset.
#
# Workflow:
#   samtools faidx REF.fa
#   N=$(bash scripts/delta/make_shard_beds.sh REF.fa.fai 25000000 $SCRATCH/shards)
#   SHARD_DIR=$SCRATCH/shards OUTROOT=$SCRATCH/wg_$(date +%s) REFERENCE=REF.fa \
#     sbatch --array=0-$((N-1))%64 scripts/delta/genome_array.sbatch
#   # then, after the array finishes:
#   SHARD_OUTROOT=<OUTROOT> bash scripts/delta/merge_shards.sh
#
# The `%64` throttle caps concurrently-running tasks (tune to balance node
# spread vs. per-node memory-bandwidth contention — see the procs-per-node note).
#
# cpus-per-task here is the packing knob, NOT parallelism (rneat runs 1 thread):
# a few cores + the mem request limit how many shards co-reside on a node, which
# governs per-node bandwidth contention. 4 cpus / 16 GB → ~15 shards/node.
# The bandwidth-optimal value is what the procs-per-node sweep determines.

#SBATCH --job-name=rneat-wgshard
#SBATCH --partition=cpu
#SBATCH --account=bhrd-delta-cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

set -euo pipefail

REPO_ROOT="${RNEAT_REPO:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/../.." && pwd)}}"
source "$REPO_ROOT/scripts/delta/lib_report.sh"   # $SCRATCH resolution

REFERENCE="${REFERENCE:-$SCRATCH/neat_data/GRCh38.fa}"
SHARD_DIR="${SHARD_DIR:-$SCRATCH/shards}"
OUTROOT="${OUTROOT:-$SCRATCH/wg_$SLURM_ARRAY_JOB_ID}"
COVERAGE="${COVERAGE:-30}"
READ_LEN="${READ_LEN:-151}"
FRAG_MEAN="${FRAG_MEAN:-350}"
FRAG_SD="${FRAG_SD:-50}"
PLOIDY="${PLOIDY:-2}"
SV_RATE_SCALE="${SV_RATE_SCALE:-0}"        # germline-style by default; >0 for SVs
SEED_ROOT="${SEED_ROOT:-rneat wg shard}"
CARGO_TARGET_DIR="${CARGO_TARGET_DIR:-$SCRATCH/cargo-target/rusty-neat}"
RNEAT_BIN="$CARGO_TARGET_DIR/release/rneat"

source "$HOME/.cargo/env" 2>/dev/null || true

TASK="${SLURM_ARRAY_TASK_ID:?this script must be submitted with --array}"
BED="$(printf '%s/shard_%04d.bed' "$SHARD_DIR" "$TASK")"
[[ -f "$BED" ]] || { echo "shard BED not found: $BED" >&2; exit 1; }
[[ -f "$REFERENCE" ]] || { echo "reference not found: $REFERENCE" >&2; exit 1; }

OUTDIR="$OUTROOT/shard_$(printf '%04d' "$TASK")"
mkdir -p "$OUTDIR"

# Idempotent: skip a shard whose outputs already exist (cheap re-submit of a
# partially-completed array).
if [[ -f "$OUTDIR/s_r1.fastq.gz" && -f "$OUTDIR/s.vcf.gz" ]]; then
    echo "shard $TASK already done — skipping."
    exit 0
fi

echo "=== wg shard $TASK : $(cat "$BED") ==="
echo "ref=$REFERENCE  cov=${COVERAGE}x  sv_rate_scale=$SV_RATE_SCALE  out=$OUTDIR"

# Single-thread (fastest per-process on Delta); per-shard distinct seed so the
# shards are independent draws. target_bed restricts generation to this window.
cat > "$OUTDIR/sim.yml" <<YAML
reference: $REFERENCE
read_len: $READ_LEN
coverage: $COVERAGE
ploidy: $PLOIDY
paired_ended: true
fragment_mean: $FRAG_MEAN
fragment_st_dev: $FRAG_SD
produce_fastq: true
produce_vcf: true
produce_bam: false
overwrite_output: true
output_dir: $OUTDIR
output_filename: s
rng_seed: $SEED_ROOT $TASK
sv_rate_scale: $SV_RATE_SCALE
target_bed: $BED
num_threads: 1
YAML

/usr/bin/time -v -o "$OUTDIR/time.txt" "$RNEAT_BIN" gen-reads -c "$OUTDIR/sim.yml"
echo "shard $TASK complete: $(grep -m1 'Elapsed' "$OUTDIR/time.txt" 2>/dev/null || true)"

#!/bin/bash
# SLURM array job: sharded whole-genome simulation on Delta — EXCLUSIVE-node packing.
#
# Each array task owns ONE whole node (--exclusive) and runs a BATCH of
# SHARDS_PER_NODE windows on it concurrently, one single-threaded rneat process
# per window. This is the key to predictable performance on Delta:
#
#   rneat is memory-bandwidth-bound. On the SHARED `cpu` partition, a per-shard
#   16-cpu slice lands our process on a node alongside other tenants' jobs, and
#   they steal the memory bandwidth rneat needs — we measured 8-min windows
#   blowing out to 100-170 min (10-20x) purely from co-tenant contention, with
#   a bimodal time distribution (median ~2.5 min, p90 ~131 min) that made the
#   whole-genome wall-clock unreproducible (one run: 3h41m, almost all of it
#   stragglers on contended nodes).
#
#   --exclusive removes the strangers: we hold the entire node, so the ONLY
#   contention is our own SHARDS_PER_NODE processes competing for the node's
#   bandwidth — and that is controllable and reproducible. The procs-per-node
#   sweep put the efficient zone at <=8 procs/node (>=29% per-proc efficiency);
#   raise SHARDS_PER_NODE for more per-node throughput (fewer nodes, cheaper),
#   lower it for faster per-shard (more nodes, faster wall-clock).
#
# Reads/variants come out in absolute genome coordinates (target_bed is a
# generation-time filter), so merge_shards.sh concatenates the per-shard outputs
# into a whole-genome dataset.
#
# Workflow:
#   samtools faidx REF.fa
#   N=$(bash scripts/delta/make_shard_beds.sh REF.fa.fai 25000000 $SCRATCH/shards)
#   K=8; T=$(( (N + K - 1) / K ))          # node-tasks = ceil(N / SHARDS_PER_NODE)
#   SHARD_DIR=$SCRATCH/shards OUTROOT=$SCRATCH/wg_$(date +%s) REFERENCE=REF.fa \
#     SHARDS_PER_NODE=$K \
#     sbatch --array=0-$((T-1))%20 scripts/delta/genome_array.sbatch
#   # then, after the array finishes:
#   SHARD_OUTROOT=<OUTROOT> bash scripts/delta/merge_shards.sh
#
# For replicated / swept runs, drive it with run_genome_reps.sh instead.
#
# The `%20` throttle caps concurrently-running NODES (each task = 1 exclusive
# node), tune to your node availability. Straggler-proofing: --time is tight
# (a healthy K<=8 batch finishes well inside it) and --requeue is set, so a task
# stuck on a sick node is killed and rescheduled elsewhere — finished shards in
# the batch are skipped on the retry (idempotent), only the unfinished ones rerun.

#SBATCH --job-name=rneat-wgshard
#SBATCH --partition=cpu
#SBATCH --account=bhrd-delta-cpu
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --requeue
#SBATCH --time=01:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

# NOTE: not `-e` — we manage per-shard failures explicitly so one bad window
# doesn't abort the whole node-batch before the others finish.
set -uo pipefail

REPO_ROOT="${RNEAT_REPO:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/../.." && pwd)}}"
source "$REPO_ROOT/scripts/delta/lib_report.sh"   # $SCRATCH resolution

REFERENCE="${REFERENCE:-$SCRATCH/neat_data/GRCh38.fa}"
SHARD_DIR="${SHARD_DIR:-$SCRATCH/shards}"
OUTROOT="${OUTROOT:-$SCRATCH/wg_$SLURM_ARRAY_JOB_ID}"
SHARDS_PER_NODE="${SHARDS_PER_NODE:-8}"    # windows packed onto each exclusive node
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
[[ -f "$REFERENCE" ]] || { echo "reference not found: $REFERENCE" >&2; exit 1; }

# This node-task owns shard ids [FIRST, LAST].
FIRST=$(( TASK * SHARDS_PER_NODE ))
LAST=$(( FIRST + SHARDS_PER_NODE - 1 ))
echo "=== wg node-task $TASK : shards $FIRST..$LAST on $(hostname) (exclusive, ${SHARDS_PER_NODE}/node) ==="
echo "ref=$REFERENCE  cov=${COVERAGE}x  sv_rate_scale=$SV_RATE_SCALE  outroot=$OUTROOT"

# Simulate one window. Single-threaded (fastest per-process on Delta); per-shard
# distinct seed so the shards are independent draws; target_bed restricts
# generation to this window. Idempotent: skips a window already complete (so a
# requeued task only redoes the unfinished windows in its batch).
run_one_shard() {
    local id="$1"
    local bed; bed="$(printf '%s/shard_%04d.bed' "$SHARD_DIR" "$id")"
    [[ -f "$bed" ]] || return 0                 # past the end of the shard list
    local outdir="$OUTROOT/shard_$(printf '%04d' "$id")"
    mkdir -p "$outdir"
    if [[ -f "$outdir/s_r1.fastq.gz" && -f "$outdir/s.vcf.gz" ]]; then
        echo "  shard $id already done — skipping."; return 0
    fi
    cat > "$outdir/sim.yml" <<YAML
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
output_dir: $outdir
output_filename: s
rng_seed: $SEED_ROOT $id
sv_rate_scale: $SV_RATE_SCALE
target_bed: $bed
num_threads: 1
YAML
    echo "  shard $id start: $(cat "$bed")"
    /usr/bin/time -v -o "$outdir/time.txt" "$RNEAT_BIN" --log-level warn gen-reads -c "$outdir/sim.yml" \
        > "$outdir/run.log" 2>&1
}

# Launch the whole batch concurrently on the exclusive node, then reap.
pids=(); ids=()
for (( id=FIRST; id<=LAST; id++ )); do
    run_one_shard "$id" &
    pids+=("$!"); ids+=("$id")
done

rc=0
for i in "${!pids[@]}"; do
    if wait "${pids[$i]}"; then
        e="$(grep -m1 'Elapsed' "$OUTROOT/shard_$(printf '%04d' "${ids[$i]}")/time.txt" 2>/dev/null || true)"
        echo "  shard ${ids[$i]} complete: $e"
    else
        echo "  shard ${ids[$i]} FAILED (see shard_$(printf '%04d' "${ids[$i]}")/run.log)" >&2
        rc=1
    fi
done

# Non-zero exit lets SLURM --requeue retry this node-batch on another node;
# the idempotent skip above means only the unfinished windows rerun.
exit $rc

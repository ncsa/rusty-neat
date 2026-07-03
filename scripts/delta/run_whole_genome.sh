#!/usr/bin/env bash
# ONE-COMMAND staged whole-genome simulation on Delta — fail-fast, fire-and-forget.
#
# Chains the whole pipeline with SLURM dependencies so a bad build/config/reference
# dies on a tiny smoke job in minutes instead of wasting a big exclusive-node grab:
#
#   [1] smoke   — one small window + FASTQ integrity check (wg_smoke.sh).
#                 A non-exclusive ~few-minute job. If it FAILS, every downstream
#                 job is auto-cancelled (afterok) — you never burn nodes on a
#                 config that produces malformed output (cf. the #125 collapse).
#   [2] array   — the sharded whole-genome run (genome_array.sbatch), exclusive
#                 nodes, K windows/node, depends on afterok:[1].
#   [3] timing  — per-window wall-clock summary (wg_timing_summary.sh), afterany:[2].
#                 This is the "simulation-only timing" deliverable — throughput
#                 numbers users can plan around, with NO downstream align/call.
#   [4] merge   — concatenate shards into one whole-genome FASTQ + truth VCF
#                 (merge_shards.sh), afterok:[2]. Skip with DO_MERGE=0 for a
#                 pure timing run.
#
# Everything is idempotent + --requeue (genome_array), so a task killed on a sick
# node reruns elsewhere and only the unfinished windows redo.
#
# USAGE (from repo root, after setup.sh has built the binary + staged the reference):
#   REFERENCE=$SCRATCH/neat_data/GRCh38.fa TAG=hg38 \
#     bash scripts/delta/run_whole_genome.sh
#
#   # bigger/messier genome, smaller shards, timing-only (no merge):
#   REFERENCE=$SCRATCH/neat_data/soy.fa TAG=soy SHARD_BP=10000000 K=4 DO_MERGE=0 \
#     bash scripts/delta/run_whole_genome.sh
#
# KNOBS (env): REFERENCE TAG SHARD_BP(25M) K(4=shards/node) MAXNODES(20) COVERAGE(30)
#   READ_LEN(151) FRAG_MEAN(350) FRAG_SD(50) PLOIDY(2) SV_RATE_SCALE(0)
#   WALL_MINUTES(auto=30+20K) DO_MERGE(1) SKIP_SMOKE(0) OUTROOT(auto)
#
# PACKING (K) — see docs/hpc_guide.md §5. rneat is memory-bandwidth-bound, so the
#   node saturates by K≈2; low K = fastest per-window but wants many whole nodes
#   (hard to schedule on a busy cluster); high K = slower per-window but schedules
#   sooner. K=4 is a schedulable default with no deadline; use K=1–2 only with a
#   reservation / idle cluster. For big/messy genomes lower SHARD_BP so no single
#   window's wall-clock blows past the per-task walltime.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
source "$REPO_ROOT/scripts/delta/lib_report.sh"   # $SCRATCH + RESULTS_DIR
D="$REPO_ROOT/scripts/delta"

REFERENCE="${REFERENCE:-$SCRATCH/neat_data/GRCh38.fa}"
TAG="${TAG:-wg}"
SHARD_BP="${SHARD_BP:-25000000}"
K="${K:-4}"                       # shards per exclusive node
MAXNODES="${MAXNODES:-20}"        # concurrent exclusive nodes (%throttle)
COVERAGE="${COVERAGE:-30}"
READ_LEN="${READ_LEN:-151}"
FRAG_MEAN="${FRAG_MEAN:-350}"
FRAG_SD="${FRAG_SD:-50}"
PLOIDY="${PLOIDY:-2}"
SV_RATE_SCALE="${SV_RATE_SCALE:-0}"
DO_MERGE="${DO_MERGE:-1}"
SKIP_SMOKE="${SKIP_SMOKE:-0}"
ACCT="${ACCESS_ACCOUNT:-bhrd-delta-cpu}"
PART="${PARTITION:-cpu}"
RD="${RESULTS_DIR:-$HOME}"
OUTROOT="${OUTROOT:-$SCRATCH/wg_${TAG}}"
SHARD_DIR="${SHARD_DIR:-$SCRATCH/shards_${TAG}}"

CARGO_TARGET_DIR="${CARGO_TARGET_DIR:-$SCRATCH/cargo-target/rusty-neat}"
RNEAT_BIN="${RNEAT_BIN:-$CARGO_TARGET_DIR/release/rneat}"
source "$HOME/.cargo/env" 2>/dev/null || true
module load samtools/1.22-cce19.0.0 2>/dev/null || true

[[ -x "$RNEAT_BIN" ]] || { echo "rneat binary missing: $RNEAT_BIN (run setup.sh)" >&2; exit 1; }
[[ -f "$REFERENCE" ]] || { echo "reference not found: $REFERENCE" >&2; exit 1; }

# ── prep (login node — light: index + write small BED files) ────────────────────
[[ -f "$REFERENCE.fai" ]] || { echo "indexing reference..."; samtools faidx "$REFERENCE"; }
N=$(bash "$D/make_shard_beds.sh" "$REFERENCE.fai" "$SHARD_BP" "$SHARD_DIR")
[[ "$N" =~ ^[0-9]+$ && "$N" -gt 0 ]] || { echo "make_shard_beds produced no shards" >&2; exit 1; }
T=$(( (N + K - 1) / K ))                              # node-tasks
MINUTES="${WALL_MINUTES:-$(( 30 + K * 20 ))}"
WALL=$(printf '%02d:%02d:00' $(( MINUTES / 60 )) $(( MINUTES % 60 )))

echo "════════════════════════════════════════════════════════════════"
echo "Whole-genome run '$TAG'"
echo "  reference : $REFERENCE ($("$RNEAT_BIN" --version 2>/dev/null || echo '?'))"
echo "  shards    : $N  (${SHARD_BP} bp windows)   packing K=$K -> $T node-tasks"
echo "  nodes     : up to $MAXNODES concurrent (exclusive)   per-task walltime $WALL"
echo "  coverage  : ${COVERAGE}x  sv_rate_scale=$SV_RATE_SCALE  ploidy=$PLOIDY"
echo "  outroot   : $OUTROOT     merge=$DO_MERGE  smoke=$([ "$SKIP_SMOKE" = 1 ] && echo skip || echo yes)"
echo "════════════════════════════════════════════════════════════════"

# ── [1] smoke — fail-fast gate ──────────────────────────────────────────────────
DEP=""
if [[ "$SKIP_SMOKE" != "1" ]]; then
    smoke_jid=$(sbatch --parsable --job-name="wg-smoke-${TAG}" --account="$ACCT" --partition="$PART" \
        --nodes=1 --ntasks=1 --cpus-per-task=4 --mem=16G --time=00:30:00 \
        --output="$RD/wg-smoke_${TAG}_%j.out" \
        --wrap="cd '$REPO_ROOT' && REFERENCE='$REFERENCE' COVERAGE=$COVERAGE READ_LEN=$READ_LEN \
                FRAG_MEAN=$FRAG_MEAN FRAG_SD=$FRAG_SD PLOIDY=$PLOIDY SV_RATE_SCALE=$SV_RATE_SCALE \
                bash scripts/delta/wg_smoke.sh")
    echo "[1] smoke   -> job $smoke_jid  (gate: downstream runs only afterok)"
    DEP="--dependency=afterok:$smoke_jid"
fi

# ── [2] array — the whole-genome sharded run ────────────────────────────────────
array_jid=$(SHARD_DIR="$SHARD_DIR" OUTROOT="$OUTROOT" REFERENCE="$REFERENCE" \
    SHARDS_PER_NODE="$K" COVERAGE="$COVERAGE" READ_LEN="$READ_LEN" \
    FRAG_MEAN="$FRAG_MEAN" FRAG_SD="$FRAG_SD" PLOIDY="$PLOIDY" \
    SV_RATE_SCALE="$SV_RATE_SCALE" SEED_ROOT="rneat wg ${TAG} shard" \
    sbatch --parsable $DEP --job-name="wg-array-${TAG}" \
           --time="$WALL" --array=0-$((T-1))%${MAXNODES} \
           "$D/genome_array.sbatch")
echo "[2] array   -> job $array_jid  ($T tasks%$MAXNODES)$([ -n "$DEP" ] && echo '  afterok smoke')"

# ── [3] timing summary — simulation-only throughput (deliverable #5) ────────────
timing_jid=$(sbatch --parsable --dependency=afterany:"$array_jid" \
    --job-name="wg-timing-${TAG}" --account="$ACCT" --partition="$PART" \
    --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=4G --time=00:20:00 \
    --output="$RD/wg-timing_${TAG}_%j.out" \
    --wrap="bash '$D/wg_timing_summary.sh' '$OUTROOT' $N")
echo "[3] timing  -> job $timing_jid  (afterany array)   summary: $RD/wg-timing_${TAG}_${timing_jid}.out"

# ── [4] merge — one whole-genome dataset ────────────────────────────────────────
if [[ "$DO_MERGE" == "1" ]]; then
    merge_jid=$(SHARD_OUTROOT="$OUTROOT" EXPECT_SHARDS="$N" \
        sbatch --parsable --dependency=afterok:"$array_jid" \
               --job-name="wg-merge-${TAG}" "$D/merge_shards.sh")
    echo "[4] merge   -> job $merge_jid  (afterok array)   -> $OUTROOT/merged_*"
else
    echo "[4] merge   -> skipped (DO_MERGE=0; timing-only run)"
fi

echo
echo "Fire-and-forget. Watch: squeue --me"
echo "If smoke FAILS, the array/merge auto-cancel (DependencyNeverSatisfied) — nothing wasted."

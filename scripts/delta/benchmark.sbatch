#!/bin/bash
# SLURM job: NEAT 4 (Python) vs rneat (Rust) throughput benchmark on Delta.
#
# Two phases (see config below):
#   A. Head-to-head rneat vs NEAT 4 at 1 thread across genome sizes -> the
#      "rneat is faster + uses less memory" claim (wall-clock + peak RSS).
#   B. rneat-only thread-scaling sweep on a MULTI-CONTIG genome -> rneat's
#      parallel speedup curve (NEAT's 1-thread time from A is the baseline).
#
# This replaces the old full (genomes x 8 threads x reps x 2 tools) cross
# product, which redundantly swept threads on single-contig genomes (where
# rneat can't parallelize) and ran NEAT at every thread count — that overran
# the walltime. Results land in $OUTDIR/median.tsv (one row per genome/threads/
# tool) and per_rep.tsv (raw).
#
# Usage:
#   sbatch scripts/delta/benchmark.sbatch
# Overrides (no edits needed):
#   SCALING_GENOME=chrset:/path/multi.fa SCALING_THREAD_MODES="1 8 64" \
#   REPS=2 sbatch scripts/delta/benchmark.sbatch

#SBATCH --job-name=rneat-benchmark
#SBATCH --partition=cpu
#SBATCH --account=bhrd-delta-cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=240G
# Lean A/B design: ~3 genomes x 2 tools @ 1 thread + ~4 rneat scaling points.
# Far smaller than the old cross-product; 6h is ample (NEAT at 1 thread is the
# slowest single run, ~13 min/rep on chr22 at 10x).
#SBATCH --time=06:00:00
# Exclusive node: this is the headline rneat-vs-NEAT4 performance comparison, so
# we want the whole node to ourselves — co-located jobs add noise to wall-clock
# and memory-bandwidth measurements. (Justified only for the benchmark; the
# other jobs are I/O- or single-contig-bound and don't need it.)
#SBATCH --exclusive
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -euo pipefail

# ── Paths (edit these) ──────────────────────────────────────────────────
# Under sbatch, $0 is a spooled copy — derive the repo from the submit dir
# (run `sbatch` from the repo root) or override with RNEAT_REPO.
REPO_ROOT="${RNEAT_REPO:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/../.." && pwd)}}"
source "$REPO_ROOT/scripts/delta/lib_report.sh"   # archive_run() + Delta path resolution
DATA_DIR="${DATA_DIR:-$SCRATCH/neat_data}"        # reference FASTAs live here
OUTDIR="${OUTDIR:-$SCRATCH/benchmark_$SLURM_JOB_ID}"

CARGO_TARGET_DIR="${CARGO_TARGET_DIR:-$SCRATCH/cargo-target/rusty-neat}"
RNEAT_BIN="$CARGO_TARGET_DIR/release/rneat"
CONDA_ENV_NAME="neat4"

# PHASE A — head-to-head rneat vs NEAT 4 at a single thread count (HEADTOHEAD_THREADS)
# across genome sizes. This is the core "rneat is faster + leaner" claim. We do
# NOT sweep threads with both tools: most of these genomes are few-contig, and
# rneat parallelizes per-contig (sub-contig chunking off by default), so a thread
# sweep on a single-contig genome (e.g. chr22) is pure redundancy — every mode
# returns the same number — and NEAT-at-every-thread-count is what blew the
# 8h walltime. Missing files are skipped.
GENOMES=(
    "ecoli:$DATA_DIR/ecoli.fa"        # single-contig, tiny
    "yeast:$DATA_DIR/yeast.fa"        # MULTI-contig (~16) — also the Phase B genome
    "chr22:$DATA_DIR/chr22.fa"        # single-contig, ~50 Mb
    "human:$DATA_DIR/GRCh38.fa"       # multi-contig, large — staged only for a big run
)
HEADTOHEAD_THREADS="${HEADTOHEAD_THREADS:-1}"

# PHASE B — rneat thread-scaling curve. Only meaningful on a MULTI-CONTIG genome
# (rneat's unit of parallelism is the contig), so chr22 alone can't show scaling.
# rneat only: NEAT 4's single-thread time from Phase A is the reference baseline;
# sweeping NEAT across threads is the expensive, low-value part we're cutting.
SCALING_GENOME="${SCALING_GENOME:-yeast:$DATA_DIR/yeast.fa}"
SCALING_THREAD_MODES=(${SCALING_THREAD_MODES:-1 4 16 64 128})

# Simulation parameters — must be identical between NEAT 4 and rneat for
# a fair comparison.
COVERAGE=10
READ_LEN=151
FRAG_MEAN=350
FRAG_SD=50
REPS=3       # repetitions per config; median is reported

# ── Environment ─────────────────────────────────────────────────────────
source "$HOME/.cargo/env"
setup_conda
conda_activate "$CONDA_ENV_NAME"
NEAT_BIN="$(conda run -n "$CONDA_ENV_NAME" which neat 2>/dev/null || echo '')"

NCORES=$SLURM_CPUS_PER_TASK
mkdir -p "$OUTDIR"

echo "=== rneat benchmark: $SLURM_JOB_ID ==="
echo "Cores available: $NCORES"
echo "rneat:  $RNEAT_BIN"
echo "NEAT 4: ${NEAT_BIN:-NOT FOUND}"
echo "Output: $OUTDIR"
echo "Genomes: ${GENOMES[*]}"

PERREP="$OUTDIR/per_rep.tsv"
MEDIAN="$OUTDIR/median.tsv"
printf "genome\tsize_mb\tthreads\ttool\trep\twall_s\tpeak_rss_mb\texit\n" > "$PERREP"
printf "genome\tsize_mb\tthreads\ttool\twall_s_med\tpeak_rss_mb_med\tn_ok\n" > "$MEDIAN"

# ── Helpers ─────────────────────────────────────────────────────────────
parse_time() {
    awk '
        /Elapsed \(wall clock\)/ {
            n=split($NF,a,":");
            if(n==3) wall=a[1]*3600+a[2]*60+a[3];
            else if(n==2) wall=a[1]*60+a[2];
            else wall=a[1];
        }
        /Maximum resident set size/ { rss=$NF/1024 }
        END { printf "%.2f %.1f", wall, rss }
    ' "$1"
}

median_of() {
    printf "%s\n" "$@" | sort -n | awk '{v[NR]=$1}
        END{ if(NR==0){print "NA"; exit}
             if(NR%2){print v[(NR+1)/2]}
             else{printf "%.2f",(v[NR/2]+v[NR/2+1])/2} }'
}

run_one() {
    local label="$1" fasta="$2" nthreads="$3" tool="$4" cfg="$5" run_dir="$6"
    local size_mb walls=() rsss=() ok=0

    size_mb=$(awk -v b="$(stat -c%s "$fasta")" 'BEGIN{printf "%.1f", b/1048576}')

    for rep in $(seq 1 "$REPS"); do
        local tf="$OUTDIR/${tool}_${label}_t${nthreads}_r${rep}.time"
        local log="$OUTDIR/${tool}_${label}_t${nthreads}_r${rep}.log"
        rm -rf "$run_dir"; mkdir -p "$run_dir"
        echo "  [$tool] $label threads=$nthreads rep=$rep/$REPS"

        # Capture the tool's real exit code without tripping `set -e`: the
        # `|| ec=$?` keeps a failed run from aborting the whole sweep, so one
        # bad config is recorded (exit!=0) rather than killing the job.
        local ec=0
        if [[ "$tool" == "neat" ]]; then
            /usr/bin/time -v -o "$tf" "$NEAT_BIN" read-simulator -c "$cfg" \
                -o "$run_dir" -p bench > "$log" 2>&1 || ec=$?
        else
            /usr/bin/time -v -o "$tf" "$RNEAT_BIN" --log-level warn gen-reads -c "$cfg" > "$log" 2>&1 || ec=$?
        fi

        local parsed wall rss
        parsed=$(parse_time "$tf"); wall=${parsed% *}; rss=${parsed#* }
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
            "$label" "$size_mb" "$nthreads" "$tool" "$rep" "$wall" "$rss" "$ec" >> "$PERREP"
        echo "    wall=${wall}s rss=${rss}MB exit=$ec"
        [[ "$ec" -eq 0 ]] && { walls+=("$wall"); rsss+=("$rss"); ok=$((ok+1)); }
    done

    local wmed rmed
    wmed=$(median_of "${walls[@]:-}")
    rmed=$(median_of "${rsss[@]:-}")
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$label" "$size_mb" "$nthreads" "$tool" "$wmed" "$rmed" "$ok" >> "$MEDIAN"
}

# run_rneat / run_neat <label> <fasta> <nthreads>: write the tool's config and
# time it. Factored out so Phase A and Phase B share one definition.
run_rneat() {
    local label="$1" fasta="$2" nthreads="$3"
    local out="$OUTDIR/rneat_${label}_t${nthreads}" cfg="$OUTDIR/rneat_${label}_t${nthreads}.yml"
    cat > "$cfg" <<YAML
reference: $fasta
read_len: $READ_LEN
coverage: $COVERAGE
paired_ended: true
fragment_mean: $FRAG_MEAN
fragment_st_dev: $FRAG_SD
produce_fastq: true
produce_vcf: false
produce_bam: false
overwrite_output: true
rng_seed: delta benchmark seed
num_threads: $nthreads
output_dir: $out
output_filename: bench
YAML
    run_one "$label" "$fasta" "$nthreads" "rneat" "$cfg" "$out"
}
run_neat() {
    local label="$1" fasta="$2" nthreads="$3"
    [[ -n "${NEAT_BIN:-}" ]] || { echo "  NEAT binary not found — skipping NEAT for $label"; return; }
    local cfg="$OUTDIR/neat_${label}_t${nthreads}.yml"
    cat > "$cfg" <<YAML
reference: $fasta
read_len: $READ_LEN
coverage: $COVERAGE
paired_ended: true
fragment_mean: $FRAG_MEAN
fragment_st_dev: $FRAG_SD
produce_fastq: true
produce_bam: false
produce_vcf: false
rng_seed: 42
threads: $nthreads
overwrite_output: true
YAML
    run_one "$label" "$fasta" "$nthreads" "neat" "$cfg" "$OUTDIR/neat_${label}_t${nthreads}"
}

# ── Phase A: head-to-head (rneat vs NEAT) at HEADTOHEAD_THREADS ───────────
echo
echo "=== Phase A: head-to-head @ ${HEADTOHEAD_THREADS} thread(s) ==="
for entry in "${GENOMES[@]}"; do
    label="${entry%%:*}"; fasta="${entry#*:}"
    [[ -f "$fasta" ]] || { echo "MISSING genome: $fasta — skipping"; continue; }
    run_rneat "$label" "$fasta" "$HEADTOHEAD_THREADS"
    run_neat  "$label" "$fasta" "$HEADTOHEAD_THREADS"
done

# ── Phase B: rneat thread-scaling on a multi-contig genome ───────────────
sc_label="${SCALING_GENOME%%:*}"; sc_fasta="${SCALING_GENOME#*:}"
if [[ -f "$sc_fasta" ]]; then
    echo
    echo "=== Phase B: rneat thread-scaling on $sc_label (multi-contig) ==="
    ncontigs=$(grep -c '^>' "$sc_fasta" 2>/dev/null || echo '?')
    echo "    ($sc_label has $ncontigs contigs; rneat parallelizes per contig)"
    for nthreads in "${SCALING_THREAD_MODES[@]}"; do
        # Skip the threads value already measured in Phase A (no need to repeat).
        [[ "$nthreads" == "$HEADTOHEAD_THREADS" ]] && continue
        run_rneat "$sc_label" "$sc_fasta" "$nthreads"
    done
else
    echo "Phase B skipped: scaling genome not staged ($sc_fasta)."
fi

echo
echo "=== MEDIAN RESULTS (REPS=$REPS) ==="
column -t -s$'\t' "$MEDIAN"
echo
echo "Full per-rep data: $PERREP"

# Persist report artifacts off scratch + record provenance/resources.
archive_run benchmark "$OUTDIR" median.tsv per_rep.tsv
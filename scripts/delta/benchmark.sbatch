#!/bin/bash
# SLURM job: NEAT 4 (Python) vs rneat (Rust) throughput benchmark on Delta.
#
# Tests wall-clock time and peak RSS across thread counts 1→64 on multiple
# genome sizes. Results land in $OUTDIR/median.tsv — one row per
# (genome, threads, tool) combination.
#
# Mirrors tools/neat_vs_rneat_benchmark.sh but runs under SLURM on a full
# 64-core Delta CPU node with proper module loading.
#
# Usage:
#   sbatch scripts/delta/benchmark.sbatch
#
# To override defaults without editing this file:
#   OUTDIR=/path/to/out sbatch scripts/delta/benchmark.sbatch

#SBATCH --job-name=rneat-benchmark
#SBATCH --partition=cpu
#SBATCH --account=<YOUR_ALLOCATION>
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=240G
#SBATCH --time=08:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -euo pipefail

# ── Paths (edit these) ──────────────────────────────────────────────────
REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
DATA_DIR="${DATA_DIR:-$SCRATCH/neat_data}"        # reference FASTAs live here
OUTDIR="${OUTDIR:-$SCRATCH/benchmark_$SLURM_JOB_ID}"

RNEAT_BIN="$REPO_ROOT/target/release/rneat"
CONDA_ENV_NAME="neat4"

# Genomes to benchmark — label:fasta_path pairs. Smaller genomes run fast
# and exercise single-thread performance; larger ones expose scaling behaviour.
# The paths below assume you've staged the files into DATA_DIR.
GENOMES=(
    "ecoli:$DATA_DIR/ecoli.fa"
    "chr22:$DATA_DIR/chr22.fa"
    "human:$DATA_DIR/GRCh38.fa"       # comment out if not staged
)

# Thread counts to test. 64 = full node; include 1 as the serial baseline.
THREAD_MODES=(1 2 4 8 16 32 64)

# Simulation parameters — must be identical between NEAT 4 and rneat for
# a fair comparison.
COVERAGE=10
READ_LEN=151
FRAG_MEAN=350
FRAG_SD=50
REPS=3       # repetitions per config; median is reported

# ── Environment ─────────────────────────────────────────────────────────
source "$HOME/.cargo/env"
module load anaconda3_cpu 2>/dev/null || module load miniforge 2>/dev/null || true
conda activate "$CONDA_ENV_NAME"
NEAT_BIN="$(conda run -n $CONDA_ENV_NAME which neat 2>/dev/null || echo '')"

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

        if [[ "$tool" == "neat" ]]; then
            /usr/bin/time -v -o "$tf" "$NEAT_BIN" read-simulator -c "$cfg" \
                -o "$run_dir" -p bench > "$log" 2>&1
        else
            /usr/bin/time -v -o "$tf" "$RNEAT_BIN" gen-reads -c "$cfg" > "$log" 2>&1
        fi

        local ec=$? parsed wall rss
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

# ── Main loop ────────────────────────────────────────────────────────────
for entry in "${GENOMES[@]}"; do
    label="${entry%%:*}"
    fasta="${entry#*:}"
    [[ -f "$fasta" ]] || { echo "MISSING genome: $fasta — skipping"; continue; }

    for nthreads in "${THREAD_MODES[@]}"; do
        # rneat config
        rneat_out="$OUTDIR/rneat_${label}_t${nthreads}"
        rneat_cfg="$OUTDIR/rneat_${label}_t${nthreads}.yml"
        cat > "$rneat_cfg" <<YAML
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
output_dir: $rneat_out
output_filename: bench
YAML
        run_one "$label" "$fasta" "$nthreads" "rneat" "$rneat_cfg" "$rneat_out"

        # NEAT 4 config (skip if binary not found)
        if [[ -n "${NEAT_BIN:-}" ]]; then
            neat_cfg="$OUTDIR/neat_${label}_t${nthreads}.yml"
            cat > "$neat_cfg" <<YAML
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
            run_one "$label" "$fasta" "$nthreads" "neat" "$neat_cfg" \
                "$OUTDIR/neat_${label}_t${nthreads}"
        fi
    done
done

echo
echo "=== MEDIAN RESULTS (REPS=$REPS) ==="
column -t -s$'\t' "$MEDIAN"
echo
echo "Full per-rep data: $PERREP"
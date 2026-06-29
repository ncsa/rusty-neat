#!/bin/bash
# SLURM job: germline recall comparison — rneat vs NEAT 4 — on Delta.
#
# For EACH simulator (rneat and NEAT 4), simulates SNP/indel-only germline reads
# and runs the SAME downstream pipeline (BWA-MEM2 → GATK HaplotypeCaller → hap.py),
# scoring the caller's recall/precision/F1 against THAT tool's own truth VCF. The
# two summaries are reported side-by-side: "how well does a standard germline
# caller recover each simulator's variants" — the fidelity half of goal 2
# (performance is benchmark.sbatch).
#
# Both tools use their DEFAULT mutation model (realistic Ts/Tv ~2.3 for each) and
# identical coverage/read-length/fragment params, so the comparison is fair.
#
# Stages per tool (all idempotent — existing outputs skipped on re-run):
#   1. simulate (SNP/indel only) -> reads + truth VCF
#   2. BWA-MEM2 align -> SAMtools sort/index
#   3. GATK HaplotypeCaller (gVCF -> genotyped VCF)
#   4. hap.py: called vs that tool's truth VCF
#
# Prerequisites: run scripts/delta/setup.sh once (builds rneat, neat4 + bioinf
# conda envs, hap.py SIF). Submit from the repo root.
#
# Usage:
#   sbatch scripts/delta/germline_e2e.sbatch
#   REFERENCE=/path/ref.fa COVERAGE=30 TOOLS="rneat neat" sbatch scripts/delta/germline_e2e.sbatch
#   TOOLS=rneat sbatch ...        # rneat only (skip the NEAT 4 arm)
#   SHUFFLE=1 TOOLS=rneat sbatch ...   # randomize read order before alignment
#       (off by default; run once shuffled vs unshuffled to confirm the caller
#        is read-order-independent — same recall expected)

#SBATCH --job-name=rneat-germline
#SBATCH --partition=cpu
#SBATCH --account=bhrd-delta-cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=240G
#SBATCH --time=24:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -euo pipefail

# ── Parameters ───────────────────────────────────────────────────────────
# Under sbatch, $0 is a spooled copy — derive the repo from the submit dir
# (run `sbatch` from the repo root) or override with RNEAT_REPO.
REPO_ROOT="${RNEAT_REPO:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/../.." && pwd)}}"
source "$REPO_ROOT/scripts/delta/lib_report.sh"   # archive_run() + Delta path resolution

REFERENCE="${REFERENCE:-$SCRATCH/neat_data/chr22.fa}"
OUTDIR="${OUTDIR:-$SCRATCH/germline_$SLURM_JOB_ID}"
TOOLS="${TOOLS:-rneat neat}"          # which simulators to run + score
# Seeds. Under `sbatch --array=1-3` (triplicate), each replicate gets a distinct
# seed (independent draw) and SLURM spreads tasks across nodes; plain submit uses
# the base seeds.
SEED="${SEED:-germline fidelity delta}${SLURM_ARRAY_TASK_ID:+ rep$SLURM_ARRAY_TASK_ID}"
NEAT_SEED="${NEAT_SEED:-42}"
[[ -n "${SLURM_ARRAY_TASK_ID:-}" ]] && NEAT_SEED=$((NEAT_SEED + SLURM_ARRAY_TASK_ID))

COVERAGE="${COVERAGE:-30}"
READ_LEN="${READ_LEN:-151}"
FRAG_MEAN="${FRAG_MEAN:-350}"
FRAG_SD="${FRAG_SD:-50}"
# Read-order randomization before alignment (off by default). rneat/NEAT emit
# reads in genome order; SHUFFLE=1 randomizes (seqkit, paired) to mimic real
# sequencer output and to demonstrate the caller is order-independent. Needs
# seqkit in the bioinf env (setup.sh installs it).
SHUFFLE="${SHUFFLE:-0}"
SHUFFLE_SEED="${SHUFFLE_SEED:-11}"

HAPPY_SIF="${HAPPY_SIF:-$SCRATCH/sif/happy_v0.3.12.sif}"
NEAT_ENV="${NEAT_ENV:-neat4}"         # conda env with NEAT 4
THREADS=$SLURM_CPUS_PER_TASK
CARGO_TARGET_DIR="${CARGO_TARGET_DIR:-$SCRATCH/cargo-target/rusty-neat}"
RNEAT_BIN="$CARGO_TARGET_DIR/release/rneat"

# ── Environment ─────────────────────────────────────────────────────────
source "$HOME/.cargo/env"
module load samtools/1.22-cce19.0.0
module load htslib/1.22-gcc13.3.1
setup_conda
conda_activate bioinf   # bwa-mem2, gatk4, bcftools (NEAT runs via `conda run -n $NEAT_ENV`)

mkdir -p "$OUTDIR"
REF_DIR="$(dirname "$REFERENCE")"
REF_NAME="$(basename "$REFERENCE")"

echo "=== germline fidelity comparison: $SLURM_JOB_ID ==="
echo "Reference: $REFERENCE   Coverage: ${COVERAGE}x   Read: ${READ_LEN}bp   Threads: $THREADS"
echo "Tools:     $TOOLS"
echo "Output:    $OUTDIR"

# ── Shared one-time reference indexing ───────────────────────────────────
# Guard on -s (exists AND non-empty), not -f: an index build killed mid-write
# leaves a 0-byte .bwt.2bit.64 that -f treats as "present", so every later run
# skips indexing and bwa-mem2 dies with "Unexpected end of file". -s rebuilds it.
[[ -f "${REFERENCE}.fai" ]] || samtools faidx "$REFERENCE"
if [[ ! -s "${REFERENCE}.bwt.2bit.64" && ! -s "${REFERENCE}.bwt" ]]; then
    echo "Indexing reference for BWA-MEM2 (one-time)..."
    bwa-mem2 index "$REFERENCE"
fi
REF_DICT="${REFERENCE%.fa}.dict"
[[ -f "$REF_DICT" ]] || gatk CreateSequenceDictionary -R "$REFERENCE" 2>&1 | tail -2

# ── Simulators (each writes <wd>/sample_r1/r2.fastq.gz + a truth VCF) ─────
# rneat: SNP/indel only via sv_rate_scale=0; truth = sample.vcf.gz
simulate_rneat() {
    local wd="$1"
    cat > "$wd/sim.yml" <<YAML
reference: $REFERENCE
read_len: $READ_LEN
coverage: $COVERAGE
ploidy: 2
paired_ended: true
fragment_mean: $FRAG_MEAN
fragment_st_dev: $FRAG_SD
produce_fastq: true
produce_vcf: true
produce_bam: false
overwrite_output: true
output_dir: $wd
output_filename: sample
rng_seed: $SEED
sv_rate_scale: 0
num_threads: $THREADS
YAML
    "$RNEAT_BIN" gen-reads -c "$wd/sim.yml"
}

# NEAT 4: default variant model is already SNP/indel only (no symbolic SVs);
# truth = sample_golden.vcf.gz. Runs in its own conda env.
simulate_neat() {
    local wd="$1"
    cat > "$wd/sim.yml" <<YAML
reference: $REFERENCE
read_len: $READ_LEN
coverage: $COVERAGE
ploidy: 2
paired_ended: true
fragment_mean: $FRAG_MEAN
fragment_st_dev: $FRAG_SD
produce_fastq: true
produce_vcf: true
produce_bam: false
rng_seed: $NEAT_SEED
overwrite_output: true
YAML
    conda run -n "$NEAT_ENV" neat read-simulator -c "$wd/sim.yml" -o "$wd" -p sample
}

# ── Shared align → call → score (scores against that tool's own truth) ────
# align_call_score <tool> <wd> <r1> <r2> <truth_vcf>
align_call_score() {
    local tool="$1" wd="$2" r1="$3" r2="$4" truth="$5"
    local bam="$wd/aligned.bam" gvcf="$wd/raw.g.vcf.gz" called="$wd/called.vcf.gz"

    if [[ -f "$bam" && -f "$bam.bai" ]]; then
        echo "  [$tool 2/3] BAM present — skipping alignment."
    else
        echo "  [$tool 2/3] BWA-MEM2 align..."
        local rg="@RG\tID:$tool\tSM:$tool\tLB:$tool\tPL:ILLUMINA"
        bwa-mem2 mem -t "$THREADS" -R "$rg" "$REFERENCE" "$r1" "$r2" 2>"$wd/bwa.log" \
            | samtools sort -@ "$THREADS" -o "$bam"
        samtools index "$bam"
    fi

    if [[ -f "$called" ]]; then
        echo "  [$tool 3/3] Called VCF present — skipping HaplotypeCaller."
    else
        echo "  [$tool 3/3] GATK HaplotypeCaller (gVCF -> genotyped)..."
        [[ -f "$gvcf" ]] || gatk HaplotypeCaller -R "$REFERENCE" -I "$bam" -O "$gvcf" \
            -ERC GVCF --native-pair-hmm-threads "$THREADS" 2>&1 | tail -5
        gatk GenotypeGVCFs -R "$REFERENCE" -V "$gvcf" -O "$called" 2>&1 | tail -3
    fi

    if [[ -f "$wd/scored.summary.csv" ]]; then
        echo "  [$tool score] hap.py output present — skipping."
    else
        echo "  [$tool score] hap.py (truth=$(basename "$truth"))..."
        [[ -f "$HAPPY_SIF" ]] || { echo "hap.py SIF not found: $HAPPY_SIF — run setup.sh" >&2; exit 1; }
        [[ -f "$truth.tbi" ]] || bcftools index -f -t "$truth"
        # Truth + called both live in $wd, bind it to /work.
        apptainer exec --bind "$REF_DIR:/refs:ro" --bind "$wd:/work" "$HAPPY_SIF" \
            /opt/hap.py/bin/hap.py \
                "/work/$(basename "$truth")" \
                "/work/$(basename "$called")" \
                -r "/refs/$REF_NAME" \
                -o "/work/scored" \
                --engine xcmp --no-roc
    fi
}

# ── Per-tool pipeline ─────────────────────────────────────────────────────
for tool in $TOOLS; do
    wd="$OUTDIR/$tool"; mkdir -p "$wd"
    r1="$wd/sample_r1.fastq.gz"; r2="$wd/sample_r2.fastq.gz"
    echo
    echo "──────── $tool ────────"

    case "$tool" in
        rneat)
            truth="$wd/sample.vcf.gz"
            if [[ -f "$r1" && -f "$truth" ]]; then echo "  [rneat 1/3] sim present — skipping."
            else echo "  [rneat 1/3] simulating (SNP/indel only)..."; simulate_rneat "$wd"; fi
            ;;
        neat)
            truth="$wd/sample_golden.vcf.gz"
            if ! conda run -n "$NEAT_ENV" which neat >/dev/null 2>&1; then
                echo "  [neat] NEAT 4 not found in env '$NEAT_ENV' — skipping NEAT arm." >&2
                continue
            fi
            if [[ -f "$r1" && -f "$truth" ]]; then echo "  [neat 1/3] sim present — skipping."
            else echo "  [neat 1/3] simulating (NEAT 4 default = SNP/indel)..."; simulate_neat "$wd"; fi
            ;;
        *) echo "  unknown tool '$tool' — skipping." >&2; continue ;;
    esac

    # Optional read-order randomization (SHUFFLE=1; off by default). rneat (and
    # NEAT) write reads in genome-coordinate order, whereas a real sequencer's
    # FASTQ is in arbitrary order. The caller is order-independent (the BAM is
    # coordinate-sorted before calling, no order-based dup-marking), so this does
    # not change recall in principle — the flag exists to *demonstrate* that
    # (shuffled vs unshuffled diff) and to better mimic real input. seqkit's
    # per-index permutation with a fixed seed shuffles R1/R2 identically, keeping
    # mates in register.
    if [[ "$SHUFFLE" == "1" ]]; then
        sh1="$wd/shuf_r1.fastq.gz"; sh2="$wd/shuf_r2.fastq.gz"
        if [[ -f "$sh1" && -f "$sh2" ]]; then
            echo "  [$tool shuffle] shuffled FASTQ present — skipping."
        else
            echo "  [$tool shuffle] randomizing read order (seqkit, paired)..."
            seqkit shuffle -s "$SHUFFLE_SEED" -o "$sh1" "$r1"
            seqkit shuffle -s "$SHUFFLE_SEED" -o "$sh2" "$r2"
        fi
        r1="$sh1"; r2="$sh2"
    fi

    align_call_score "$tool" "$wd" "$r1" "$r2" "$truth"
    # Stage the summary at the top level with a tool tag for the report.
    [[ -f "$wd/scored.summary.csv" ]] && cp -f "$wd/scored.summary.csv" "$OUTDIR/${tool}_scored.summary.csv"
done

# ── Side-by-side summary ──────────────────────────────────────────────────
echo
echo "════════════════════════════════════════════════════════════════"
echo "Germline fidelity comparison — Job $SLURM_JOB_ID  (${COVERAGE}x, $REF_NAME)"
echo "Recall/precision/F1 of GATK recovering each simulator's truth set."
for tool in $TOOLS; do
    f="$OUTDIR/${tool}_scored.summary.csv"
    echo
    echo "### $tool"
    [[ -f "$f" ]] && column -t -s',' "$f" || echo "  (no summary — arm skipped or failed)"
done

# Persist report artifacts off scratch + record provenance/resources.
archive_run germline "$OUTDIR" rneat_scored.summary.csv neat_scored.summary.csv

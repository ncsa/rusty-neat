#!/bin/bash
# SLURM job: rneat germline end-to-end recall benchmark on Delta.
#
# Simulates SNPs and indels only (no SVs) then runs the full downstream
# germline variant-calling pipeline and scores precision/recall/F1 against
# the rneat truth VCF. This establishes the baseline recall floor before
# layering in SVs or cancer modeling.
#
# Stages:
#   1. rneat gen-reads — SNP/indel only (sv_rate_scale: 0), paired-end, 30x
#   2. BWA-MEM2 align → SAMtools sort/index
#   3. GATK HaplotypeCaller (diploid germline, gVCF → VCF)
#   4. Score with hap.py against rneat truth VCF (germline mode, not som.py)
#
# Each stage is idempotent — existing outputs are skipped on re-run.
#
# Prerequisites: run scripts/delta/setup.sh once before submitting.
#
# Usage:
#   sbatch scripts/delta/germline_e2e.sbatch
#
# Overrides (no file edits needed):
#   REFERENCE=/path/to/chr22.fa COVERAGE=30 sbatch scripts/delta/germline_e2e.sbatch

#SBATCH --job-name=rneat-germline
#SBATCH --partition=cpu
#SBATCH --account=bhrd-delta-cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=240G
#SBATCH --time=12:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -euo pipefail

# ── Parameters (edit these) ─────────────────────────────────────────────
# Under sbatch, $0 is a spooled copy — derive the repo from the submit dir
# (run `sbatch` from the repo root) or override with RNEAT_REPO.
REPO_ROOT="${RNEAT_REPO:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/../.." && pwd)}}"
source "$REPO_ROOT/scripts/delta/lib_report.sh"   # archive_run() + Delta path resolution

# Start with chr22 to validate the pipeline cheaply, then swap to GRCh38.
REFERENCE="${REFERENCE:-$SCRATCH/neat_data/chr22.fa}"
OUTDIR="${OUTDIR:-$SCRATCH/germline_$SLURM_JOB_ID}"

COVERAGE="${COVERAGE:-30}"
READ_LEN="${READ_LEN:-151}"
FRAG_MEAN="${FRAG_MEAN:-350}"
FRAG_SD="${FRAG_SD:-50}"
# sv_rate_scale: 0 disables all symbolic SVs — SNP/indel only.
# This is the key parameter that separates this job from cancer_pipeline.sbatch.
SV_RATE_SCALE=0

HAPPY_SIF="${HAPPY_SIF:-$SCRATCH/sif/happy_v0.3.12.sif}"
THREADS=$SLURM_CPUS_PER_TASK
CARGO_TARGET_DIR="${CARGO_TARGET_DIR:-$SCRATCH/cargo-target/rusty-neat}"
RNEAT_BIN="$CARGO_TARGET_DIR/release/rneat"
PREFIX="sample"

# ── Environment ─────────────────────────────────────────────────────────
source "$HOME/.cargo/env"
module load samtools/1.22-cce19.0.0
module load htslib/1.22-gcc13.3.1
module load anaconda3_cpu 2>/dev/null || module load miniforge 2>/dev/null || true
conda activate bioinf   # provides bwa-mem2, gatk4, bcftools (not available as modules)

mkdir -p "$OUTDIR"
REF_DIR="$(dirname "$REFERENCE")"
REF_NAME="$(basename "$REFERENCE")"

echo "=== rneat germline e2e: $SLURM_JOB_ID ==="
echo "Reference:  $REFERENCE"
echo "Output:     $OUTDIR"
echo "Coverage:   ${COVERAGE}x  Read length: ${READ_LEN}bp  SV scale: $SV_RATE_SCALE"
echo "Threads:    $THREADS"

SIM_R1="$OUTDIR/${PREFIX}_r1.fastq.gz"
SIM_R2="$OUTDIR/${PREFIX}_r2.fastq.gz"
TRUTH_VCF="$OUTDIR/${PREFIX}.vcf.gz"
ALIGNED_BAM="$OUTDIR/${PREFIX}.bam"
RAW_GVCF="$OUTDIR/${PREFIX}.g.vcf.gz"
CALLED_VCF="$OUTDIR/${PREFIX}.vcf.genotyped.gz"
SCORE_PREFIX="$OUTDIR/scored"

# ── Stage 1: Simulate (SNP/indel only) ──────────────────────────────────
if [[ -f "$SIM_R1" && -f "$TRUTH_VCF" ]]; then
    echo "[1/4] Simulation outputs already present — skipping."
else
    echo "[1/4] Simulating reads (SNP/indel only, sv_rate_scale=0)..."
    SIM_CFG="$OUTDIR/sim.yml"
    cat > "$SIM_CFG" <<YAML
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
output_dir: $OUTDIR
output_filename: $PREFIX
rng_seed: germline e2e delta
sv_rate_scale: $SV_RATE_SCALE
num_threads: $THREADS
YAML
    "$RNEAT_BIN" gen-reads -c "$SIM_CFG"
    echo "[1/4] Simulation complete."
fi

# ── Stage 2: Align ───────────────────────────────────────────────────────
if [[ -f "$ALIGNED_BAM" && -f "${ALIGNED_BAM}.bai" ]]; then
    echo "[2/4] BAM already present — skipping alignment."
else
    echo "[2/4] Aligning reads..."
    if [[ ! -f "${REFERENCE}.bwt.2bit.64" ]] && [[ ! -f "${REFERENCE}.bwt" ]]; then
        echo "  Indexing reference (one-time, ~2 min)..."
        bwa-mem2 index "$REFERENCE"
    fi
    [[ -f "${REFERENCE}.fai" ]] || samtools faidx "$REFERENCE"

    RG="@RG\tID:SAMPLE\tSM:SAMPLE\tLB:SAMPLE\tPL:ILLUMINA"
    bwa-mem2 mem -t "$THREADS" -R "$RG" "$REFERENCE" "$SIM_R1" "$SIM_R2" \
        2>"$OUTDIR/bwa.log" \
        | samtools sort -@ "$THREADS" -o "$ALIGNED_BAM"
    samtools index "$ALIGNED_BAM"
    echo "[2/4] Alignment complete."
fi

# ── Stage 3: GATK HaplotypeCaller ────────────────────────────────────────
# gVCF → genotyped VCF: two-step is more expensive but gives clean per-site
# depth annotation used by hap.py for ROC-style analysis.
if [[ -f "$CALLED_VCF" ]]; then
    echo "[3/4] Genotyped VCF already present — skipping HaplotypeCaller."
else
    REF_DICT="${REFERENCE%.fa}.dict"
    [[ -f "$REF_DICT" ]] || gatk CreateSequenceDictionary -R "$REFERENCE"

    if [[ ! -f "$RAW_GVCF" ]]; then
        echo "[3/4] Running HaplotypeCaller (gVCF mode)..."
        gatk HaplotypeCaller \
            -R "$REFERENCE" \
            -I "$ALIGNED_BAM" \
            -O "$RAW_GVCF" \
            -ERC GVCF \
            --native-pair-hmm-threads "$THREADS" \
            2>&1 | tail -10
    else
        echo "[3/4] gVCF already present — skipping to genotyping."
    fi

    echo "[3/4] Genotyping gVCF..."
    gatk GenotypeGVCFs \
        -R "$REFERENCE" \
        -V "$RAW_GVCF" \
        -O "$CALLED_VCF" \
        2>&1 | tail -5
    echo "[3/4] HaplotypeCaller complete."
fi

# ── Stage 4: Score with hap.py ───────────────────────────────────────────
# hap.py in germline mode (not som.py) — evaluates diploid genotype calls.
# No confident-regions BED needed for synthetic truth; hap.py evaluates the
# entire reference. SNP and INDEL rows in the output are the numbers to watch.
if [[ -f "${SCORE_PREFIX}.summary.csv" ]]; then
    echo "[4/4] Scoring output already present — skipping."
else
    echo "[4/4] Scoring with hap.py (germline mode)..."
    [[ -f "$HAPPY_SIF" ]] || { echo "hap.py SIF not found: $HAPPY_SIF — run setup.sh first" >&2; exit 1; }

    # Index the truth VCF if needed (hap.py requires a tabix index).
    [[ -f "${TRUTH_VCF}.tbi" ]] || bcftools index -f -t "$TRUTH_VCF"

    apptainer exec \
        --bind "$REF_DIR:/refs:ro" \
        --bind "$OUTDIR:/work" \
        "$HAPPY_SIF" \
        /opt/hap.py/bin/hap.py \
            "/work/$(basename "$TRUTH_VCF")" \
            "/work/$(basename "$CALLED_VCF")" \
            -r "/refs/$REF_NAME" \
            -o "/work/scored" \
            --engine xcmp \
            --no-roc
    echo "[4/4] Scoring complete."
fi

# ── Summary ──────────────────────────────────────────────────────────────
cat <<EOF

════════════════════════════════════════════════════════════════
Germline e2e complete — Job $SLURM_JOB_ID

Settings:
  Reference:    $REFERENCE
  Coverage:     ${COVERAGE}x  Paired-end: ${FRAG_MEAN}±${FRAG_SD}bp fragments
  SV rate:      $SV_RATE_SCALE (0 = SNP/indel only)

Key outputs in $OUTDIR:
  Simulated reads:  ${PREFIX}_r1/r2.fastq.gz
  Truth VCF:        ${PREFIX}.vcf.gz
  Called VCF:       ${PREFIX}.vcf.genotyped.gz
  hap.py summary:   scored.summary.csv   ← SNP / INDEL recall & precision

Top-line results (SNP and INDEL rows are the signal):
EOF
[[ -f "${SCORE_PREFIX}.summary.csv" ]] \
    && column -t -s',' "${SCORE_PREFIX}.summary.csv" \
    || echo "  (summary file not found)"

# Persist report artifacts off scratch + record provenance/resources.
archive_run germline "$OUTDIR" scored.summary.csv scored.extended.csv sim.yml

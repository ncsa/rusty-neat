#!/bin/bash
# SLURM job: end-to-end cancer pipeline on Delta.
#
# Stages:
#   1. Simulate tumor/normal FASTQ + origin-tagged truth VCF via
#      tools/cancer_simulate.sh (wraps rneat gen-reads twice).
#   2. Align normal + merged-tumor FASTQs to reference with BWA-MEM2.
#   3. Sort, index BAMs with SAMtools.
#   4. Call somatic variants with GATK Mutect2 (tumor/normal mode).
#   5. Filter calls with GATK FilterMutectCalls.
#   6. Score filtered calls against the rneat truth VCF using hap.py
#      som.py inside an Apptainer container (Docker not available on HPC).
#
# Each stage is idempotent — outputs that already exist are skipped,
# so re-running after a failure resumes from the last successful stage.
#
# Prerequisites:
#   - Run scripts/delta/setup.sh once before submitting this job.
#   - Stage your reference FASTA into $SCRATCH (see REFERENCE below).
#   - Modules: bwa-mem2, samtools, gatk, bcftools must be available.
#     Adjust `module load` lines to match Delta's current module names.
#
# Usage:
#   sbatch scripts/delta/cancer_pipeline.sbatch
#
# To override key parameters without editing:
#   REFERENCE=/path/to/ref.fa PURITY=0.4 sbatch scripts/delta/cancer_pipeline.sbatch

#SBATCH --job-name=rneat-cancer
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

# ── Parameters (edit these) ─────────────────────────────────────────────
# Under sbatch, $0 is a spooled copy — derive the repo from the submit dir
# (run `sbatch` from the repo root) or override with RNEAT_REPO.
REPO_ROOT="${RNEAT_REPO:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/../.." && pwd)}}"
source "$REPO_ROOT/scripts/delta/lib_report.sh"   # archive_run() + Delta path resolution

# Reference genome. For a smoke test use chr22; for a full-scale run swap
# in GRCh38. Must be uncompressed (.fa, not .fa.gz) for BWA-MEM2 indexing.
REFERENCE="${REFERENCE:-$SCRATCH/neat_data/chr22.fa}"

# Output root — all intermediate files and final scores go here.
OUTDIR="${OUTDIR:-$SCRATCH/cancer_$SLURM_JOB_ID}"

# Simulation parameters.
TOTAL_COVERAGE="${TOTAL_COVERAGE:-30}"
PURITY="${PURITY:-0.6}"          # tumor cell fraction
READ_LEN="${READ_LEN:-151}"
PAIRED="${PAIRED:-true}"
FRAG_MEAN="${FRAG_MEAN:-350}"
FRAG_SD="${FRAG_SD:-50}"
SV_RATE_SCALE="${SV_RATE_SCALE:-1.0}"   # >0 enables SV simulation in tumor pass

# Tumor mutation model (COSMIC pan-cancer). Bundled in the repo.
TUMOR_MODEL="$REPO_ROOT/tools/cosmic_v104_pancancer_model.json.gz"

# Path to the hap.py Apptainer SIF (built by setup.sh).
HAPPY_SIF="${HAPPY_SIF:-$SCRATCH/sif/happy_v0.3.12.sif}"

# Run Strelka2 as a 2nd somatic SNV/indel caller alongside Mutect2 (cross-validation).
RUN_STRELKA="${RUN_STRELKA:-1}"
STRELKA_ENV="${STRELKA_ENV:-strelka}"

THREADS=$SLURM_CPUS_PER_TASK
CARGO_TARGET_DIR="${CARGO_TARGET_DIR:-$SCRATCH/cargo-target/rusty-neat}"
RNEAT_BIN="$CARGO_TARGET_DIR/release/rneat"
PREFIX="$OUTDIR/sample"

# ── Environment ─────────────────────────────────────────────────────────
source "$HOME/.cargo/env"
module load samtools/1.22-cce19.0.0
module load htslib/1.22-gcc13.3.1
setup_conda
conda_activate bioinf   # provides bwa-mem2, gatk4, bcftools (not available as modules)

mkdir -p "$OUTDIR"
REF_DIR="$(dirname "$REFERENCE")"
REF_NAME="$(basename "$REFERENCE")"

echo "=== rneat cancer pipeline: $SLURM_JOB_ID ==="
echo "Reference:  $REFERENCE"
echo "Output:     $OUTDIR"
echo "Coverage:   ${TOTAL_COVERAGE}x  Purity: $PURITY"
echo "Threads:    $THREADS"
echo "SV scale:   $SV_RATE_SCALE"

# ── Stage 1: Simulate ────────────────────────────────────────────────────
MERGED_R1="${PREFIX}_merged_r1.fastq.gz"
NORMAL_R1="${PREFIX}_normal_r1.fastq.gz"
TRUTH_VCF="${PREFIX}_merged_truth.vcf.gz"

if [[ -f "$MERGED_R1" && -f "$NORMAL_R1" && -f "$TRUTH_VCF" ]]; then
    echo "[1/6] Simulation outputs already present — skipping."
else
    echo "[1/6] Simulating tumor/normal reads..."
    PAIRED_FLAG=""
    [[ "$PAIRED" == "true" ]] && PAIRED_FLAG="--paired-ended --fragment-mean $FRAG_MEAN --fragment-st-dev $FRAG_SD"
    # Under `sbatch --array=1-3` (triplicate), each replicate uses a distinct
    # seed root (independent draw); SLURM spreads tasks across nodes. Plain
    # submit uses the default root.
    RNG_SEED_ROOT="${RNG_SEED_ROOT:-cancer-simulate}${SLURM_ARRAY_TASK_ID:+-rep$SLURM_ARRAY_TASK_ID}"
    bash "$REPO_ROOT/tools/cancer_simulate.sh" \
        --reference "$REFERENCE" \
        --output-prefix "$PREFIX" \
        --total-coverage "$TOTAL_COVERAGE" \
        --purity "$PURITY" \
        --read-len "$READ_LEN" \
        --tumor-model "$TUMOR_MODEL" \
        --sv-rate-scale "$SV_RATE_SCALE" \
        --rng-seed "$RNG_SEED_ROOT" \
        --rneat-bin "$RNEAT_BIN" \
        $PAIRED_FLAG
    echo "[1/6] Simulation complete."
fi

# ── Stage 2 & 3: Index reference, align FASTQs ──────────────────────────
NORMAL_BAM="$OUTDIR/normal.bam"
TUMOR_BAM="$OUTDIR/tumor.bam"

index_and_align() {
    local in_fq="$1" out_bam="$2" sample="$3"

    if [[ -f "$out_bam" && -f "${out_bam}.bai" ]]; then
        echo "  $sample BAM already present — skipping."
        return 0
    fi

    # Index reference once (idempotent — bwa-mem2 skips if already present).
    if [[ ! -f "${REFERENCE}.bwt.2bit.64" ]] && [[ ! -f "${REFERENCE}.bwt" ]]; then
        echo "  Indexing reference for BWA-MEM2 (one-time, ~2 min)..."
        bwa-mem2 index "$REFERENCE" 2>&1 | tail -3
    fi
    if [[ ! -f "${REFERENCE}.fai" ]]; then
        samtools faidx "$REFERENCE"
    fi

    echo "  Aligning $sample..."
    local rg="@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA"
    if [[ "$PAIRED" == "true" ]]; then
        local in_r2="${in_fq/_r1./_r2.}"
        bwa-mem2 mem -t "$THREADS" -R "$rg" "$REFERENCE" "$in_fq" "$in_r2" 2>"$OUTDIR/${sample}_bwa.log" \
            | samtools sort -@ "$THREADS" -o "$out_bam"
    else
        bwa-mem2 mem -t "$THREADS" -R "$rg" "$REFERENCE" "$in_fq" 2>"$OUTDIR/${sample}_bwa.log" \
            | samtools sort -@ "$THREADS" -o "$out_bam"
    fi
    samtools index "$out_bam"
    echo "  $sample BAM ready: $out_bam"
}

echo "[2/6] Aligning normal reads..."
index_and_align "$NORMAL_R1" "$NORMAL_BAM" "NORMAL"

echo "[3/6] Aligning tumor (merged) reads..."
# Merged FASTQ is the purity-mixed tumor biopsy.
TUMOR_MERGED_R1="${PREFIX}_merged_r1.fastq.gz"
index_and_align "$TUMOR_MERGED_R1" "$TUMOR_BAM" "TUMOR"

# ── Stage 4: Mutect2 ─────────────────────────────────────────────────────
RAW_VCF="$OUTDIR/mutect2.unfiltered.vcf.gz"
if [[ -f "$RAW_VCF" ]]; then
    echo "[4/6] Mutect2 raw VCF already present — skipping."
else
    echo "[4/6] Running Mutect2 (tumor/normal mode)..."
    # Sequence dictionary required by GATK.
    REF_DICT="${REFERENCE%.fa}.dict"
    [[ -f "$REF_DICT" ]] || gatk CreateSequenceDictionary -R "$REFERENCE"

    gatk Mutect2 \
        -R "$REFERENCE" \
        -I "$TUMOR_BAM"  --tumor-sample  TUMOR \
        -I "$NORMAL_BAM" --normal-sample NORMAL \
        --native-pair-hmm-threads "$THREADS" \
        -O "$RAW_VCF" 2>&1 | tail -20
    echo "[4/6] Mutect2 complete."
fi

# ── Stage 5: FilterMutectCalls ───────────────────────────────────────────
FILTERED_VCF="$OUTDIR/mutect2.filtered.vcf.gz"
if [[ -f "$FILTERED_VCF" ]]; then
    echo "[5/6] FilterMutectCalls output already present — skipping."
else
    echo "[5/6] Filtering Mutect2 calls..."
    gatk FilterMutectCalls \
        -R "$REFERENCE" \
        -V "$RAW_VCF" \
        -O "$FILTERED_VCF"
    echo "[5/6] Filtering complete."
fi

# ── Stage 5b: Strelka2 somatic SNV/indel (2nd caller vs Mutect2) ─────────
# Cross-validation: a second, independently-modeled somatic caller confirms
# rneat's somatic SNV/indel data isn't recovered only by Mutect2. Reuses the same
# tumor/normal BAMs. RUN_STRELKA=0 to skip.
STRELKA_VCF=""
if [[ "$RUN_STRELKA" == "1" ]]; then
    STRELKA_DIR="$OUTDIR/strelka"
    STRELKA_VCF="$OUTDIR/strelka.somatic.vcf.gz"
    if [[ -f "$STRELKA_VCF" ]]; then
        echo "[5b] Strelka2 somatic VCF already present — skipping."
    else
        echo "[5b] Running Strelka2 (tumor/normal somatic SNV/indel)..."
        # Strelka2 is a SECONDARY cross-check and must NEVER abort the run: the
        # bioconda strelka 2.9.10 binary SIGSEGVs on some HPC glibc/CPU combos.
        # Run it guarded — on failure, warn and continue so the primary Mutect2
        # result is still scored (Stage 6).
        run_strelka() {
            conda_activate "$STRELKA_ENV"
            rm -rf "$STRELKA_DIR"
            configureStrelkaSomaticWorkflow.py \
                --normalBam "$NORMAL_BAM" --tumorBam "$TUMOR_BAM" \
                --referenceFasta "$REFERENCE" --runDir "$STRELKA_DIR" \
                > "$OUTDIR/strelka_config.log" 2>&1 || return 1
            "$STRELKA_DIR/runWorkflow.py" -m local -j "$THREADS" \
                > "$OUTDIR/strelka_run.log" 2>&1 || return 1
        }
        if run_strelka; then
            conda_activate bioinf   # Strelka emits SNVs + indels separately; merge.
            bcftools concat -a \
                "$STRELKA_DIR/results/variants/somatic.snvs.vcf.gz" \
                "$STRELKA_DIR/results/variants/somatic.indels.vcf.gz" \
                | bcftools sort -Oz -o "$STRELKA_VCF"
            bcftools index -f -t "$STRELKA_VCF"
            echo "[5b] Strelka2 complete: $(zcat "$STRELKA_VCF" | grep -vc '^#' || true) somatic calls."
        else
            echo "[5b] WARNING: Strelka2 failed (see $OUTDIR/strelka_run.log) — continuing without it; Mutect2 unaffected." >&2
            STRELKA_VCF=""
        fi
    fi
fi

# ── Stage 6: Score each caller with hap.py som.py ───────────────────────
# Somatic-only truth (NEAT_ORIGIN=somatic); germline carry-through (shared) is
# correctly ignored by T/N callers, so scoring against it would be meaningless.
SCORE_TRUTH="$OUTDIR/scoring_truth.vcf.gz"
if [[ ! -f "$SCORE_TRUTH" ]]; then
    echo "[6/6] Filtering truth to somatic records..."
    bcftools view -i 'INFO/NEAT_ORIGIN="somatic"' -O z -o "$SCORE_TRUTH" "$TRUTH_VCF"
    bcftools index -f -t "$SCORE_TRUTH"
fi
somatic_count=$(bcftools view -H "$SCORE_TRUTH" | wc -l)
echo "      Somatic truth records: $somatic_count"
[[ "$somatic_count" -eq 0 ]] && echo "WARNING: no somatic truth records — check NEAT_ORIGIN." >&2
[[ -f "$HAPPY_SIF" ]] || { echo "hap.py SIF not found: $HAPPY_SIF — run setup.sh first" >&2; exit 1; }

# hap.py/som.py hardcodes UCSC chr-prefixing of input VCFs (perl s/^([0-9XYM])/chr/),
# then bcftools-norm's against the reference. With an Ensembl-named reference (1/2/3)
# the renamed VCFs (chr1) no longer match the reference (1) and som.py dies
# ("sequence chr1 not found"). If the reference isn't chr-prefixed, build a
# chr-prefixed copy for scoring so the renamed VCFs match. (chr22.fa etc. pass
# through unchanged.) Alignment/Mutect2 already used $REFERENCE as-is — unaffected.
SCORE_REF_DIR="$REF_DIR"; SCORE_REF_NAME="$REF_NAME"
if ! head -1 "${REFERENCE}.fai" 2>/dev/null | grep -q '^chr'; then
    SCORE_REF="${REFERENCE%.fa}.chr.fa"
    if [[ ! -f "$SCORE_REF" ]]; then
        echo "[6/6] Reference is Ensembl-named; building chr-prefixed copy for som.py..."
        sed '/^>/ s/^>/>chr/' "$REFERENCE" > "$SCORE_REF"
        samtools faidx "$SCORE_REF"
    fi
    SCORE_REF_DIR="$(dirname "$SCORE_REF")"; SCORE_REF_NAME="$(basename "$SCORE_REF")"
fi

# som_score <query.vcf.gz> <out-label>. -N (--normalize-all) runs bcftools norm on
# BOTH truth and query — essential for indels (representation/left-align mismatch
# otherwise double-counts FN+FP). HGREF resolves som.py's reference lookup.
som_score() {
    local query="$1" label="$2"
    [[ -f "$query" ]] || { echo "  ($label: no query VCF to score)"; return 0; }
    if [[ -f "$OUTDIR/${label}.stats.csv" ]]; then
        echo "[6/6] $label score already present — skipping."; return 0
    fi
    echo "[6/6] Scoring $label with som.py (Apptainer)..."
    apptainer exec \
        --env "HGREF=/refs/$SCORE_REF_NAME" \
        --bind "$SCORE_REF_DIR:/refs:ro" \
        --bind "$OUTDIR:/work" \
        "$HAPPY_SIF" \
        /opt/hap.py/bin/som.py \
            "/work/$(basename "$SCORE_TRUTH")" \
            "/work/$(basename "$query")" \
            -r "/refs/$SCORE_REF_NAME" \
            -o "/work/$label" \
            -N \
            --feature-table generic
}

som_score "$FILTERED_VCF" scored                                    # Mutect2
[[ -n "$STRELKA_VCF" ]] && som_score "$STRELKA_VCF" scored_strelka   # Strelka2
echo "[6/6] Scoring complete."

# ── Summary ──────────────────────────────────────────────────────────────
cat <<EOF

════════════════════════════════════════════════════════════════
Cancer pipeline complete — Job $SLURM_JOB_ID

Inputs:
  Reference:      $REFERENCE
  Coverage:       ${TOTAL_COVERAGE}x  Purity: $PURITY
  SV rate scale:  $SV_RATE_SCALE

Key outputs in $OUTDIR:
  Simulation:     sample_merged_{r1,r2}.fastq.gz
  Truth VCF:      sample_merged_truth.vcf.gz  (NEAT_ORIGIN-tagged)
  Somatic truth:  scoring_truth.vcf.gz         (somatic only)
  Mutect2 calls:  mutect2.filtered.vcf.gz   -> scored.stats.csv
  Strelka2 calls: strelka.somatic.vcf.gz    -> scored_strelka.stats.csv  (RUN_STRELKA=1)

Top-line results (Mutect2):
EOF
[[ -f "$OUTDIR/scored.stats.csv" ]] && column -t -s',' "$OUTDIR/scored.stats.csv" || echo "  (Mutect2 score not found)"
if [[ -n "$STRELKA_VCF" && -f "$OUTDIR/scored_strelka.stats.csv" ]]; then
    echo; echo "Top-line results (Strelka2):"
    column -t -s',' "$OUTDIR/scored_strelka.stats.csv"
fi

# Persist report artifacts off scratch + record provenance/resources.
archive_run cancer "$OUTDIR" scored.stats.csv

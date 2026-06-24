#!/bin/bash
# SLURM job: structural-variant (SV) validation pipeline on Delta.
#
# WHY THIS EXISTS: gen-cancer-reads emits structural variants (DEL/DUP/INV/INS/BND
# and CNVs) into the tumor truth VCF, but every validation so far ran through
# Mutect2 + som.py, which score SNVs/indels ONLY. The SV machinery has therefore
# never been checked against a real SV caller. This pipeline closes that gap: it
# simulates an SV-rich tumor/normal pair, calls somatic SVs with Manta, and scores
# recall/precision against rneat's own SV truth with truvari — the same
# "align + call + compare" loop that surfaced the strand bugs (#285/#287) for SNVs.
#
# Stages:
#   1. Simulate tumor/normal FASTQ + SV-tagged truth VCF (cancer_simulate.sh,
#      SV_RATE_SCALE > 0).
#   2-3. Align normal + merged-tumor to reference with BWA-MEM2; sort + index.
#   3b. PRE-FLIGHT: inspect rneat's SV truth records (SVTYPE/SVLEN/END/symbolic
#      ALT) and extract a somatic-SV-only truth set; abort early if it's empty or
#      malformed so we don't waste a Manta run. (First-run sanity for VCF
#      compatibility — adjust the SV filter below to match the schema it prints.)
#   4. Call somatic SVs with Manta (tumor/normal mode).
#   5. Score Manta calls vs the SV truth with truvari (overall + per SV type).
#
# Each stage is idempotent — existing outputs are skipped, so a re-run resumes.
#
# PREREQUISITES (one-time; Manta needs python2, truvari python3 — separate envs):
#   bash scripts/delta/setup.sh              # reference build + bioinf env
#   conda create -y -n manta   -c bioconda -c conda-forge manta
#   conda create -y -n truvari -c bioconda -c conda-forge truvari
#   (override env names with MANTA_ENV / TRUVARI_ENV if you named them differently)
#
# Usage:
#   sbatch scripts/delta/sv_pipeline.sbatch
#   REFERENCE=/path/ref.fa SV_RATE_SCALE=2.0 TUMOR_MODEL=/path/model.json.gz \
#     sbatch scripts/delta/sv_pipeline.sbatch
#   # per-tissue SV spectrum check (BRCA DUP-dominant / skin BND / lung DEL):
#   TUMOR_MODEL=$REPO_ROOT/tools/cosmic_pancancer_sv_BRCA.json.gz sbatch ...

#SBATCH --job-name=rneat-sv
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

# ── Parameters ───────────────────────────────────────────────────────────
REPO_ROOT="${RNEAT_REPO:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/../.." && pwd)}}"
source "$REPO_ROOT/scripts/delta/lib_report.sh"   # archive_run() + Delta path resolution

REFERENCE="${REFERENCE:-$SCRATCH/neat_data/chr22.fa}"
OUTDIR="${OUTDIR:-$SCRATCH/sv_$SLURM_JOB_ID}"

TOTAL_COVERAGE="${TOTAL_COVERAGE:-30}"
PURITY="${PURITY:-0.6}"
READ_LEN="${READ_LEN:-151}"
PAIRED="${PAIRED:-true}"
FRAG_MEAN="${FRAG_MEAN:-350}"
FRAG_SD="${FRAG_SD:-50}"
# SVs are the whole point — must be > 0. Scale up to enrich the truth set.
SV_RATE_SCALE="${SV_RATE_SCALE:-1.0}"

# Tumor model: pan-cancer SV model by default; swap a per-tissue SV model
# (cosmic_pancancer_sv_{BRCA,skin,lung}.json.gz) to check tissue-specific spectra.
TUMOR_MODEL="${TUMOR_MODEL:-$REPO_ROOT/tools/cosmic_v104_pancancer_model.json.gz}"

MANTA_ENV="${MANTA_ENV:-manta}"
TRUVARI_ENV="${TRUVARI_ENV:-truvari}"

THREADS=$SLURM_CPUS_PER_TASK
CARGO_TARGET_DIR="${CARGO_TARGET_DIR:-$SCRATCH/cargo-target/rusty-neat}"
RNEAT_BIN="$CARGO_TARGET_DIR/release/rneat"
PREFIX="$OUTDIR/sample"

if awk "BEGIN{exit !($SV_RATE_SCALE <= 0)}"; then
    echo "ERROR: SV_RATE_SCALE must be > 0 for the SV pipeline (got $SV_RATE_SCALE)." >&2
    exit 1
fi

# ── Environment ─────────────────────────────────────────────────────────
source "$HOME/.cargo/env"
module load samtools/1.22-cce19.0.0
module load htslib/1.22-gcc13.3.1
setup_conda

mkdir -p "$OUTDIR"
REF_NAME="$(basename "$REFERENCE")"

echo "=== rneat SV validation pipeline: $SLURM_JOB_ID ==="
echo "Reference:  $REFERENCE"
echo "Output:     $OUTDIR"
echo "Coverage:   ${TOTAL_COVERAGE}x  Purity: $PURITY  SV scale: $SV_RATE_SCALE"
echo "Tumor model:$TUMOR_MODEL"
echo "Threads:    $THREADS"

# ── Stage 1: Simulate (SV-rich tumor/normal) ─────────────────────────────
conda_activate bioinf   # bwa-mem2, bcftools (samtools via module)
MERGED_R1="${PREFIX}_merged_r1.fastq.gz"
NORMAL_R1="${PREFIX}_normal_r1.fastq.gz"
TRUTH_VCF="${PREFIX}_merged_truth.vcf.gz"

if [[ -f "$MERGED_R1" && -f "$NORMAL_R1" && -f "$TRUTH_VCF" ]]; then
    echo "[1/5] Simulation outputs already present — skipping."
else
    echo "[1/5] Simulating SV-rich tumor/normal reads..."
    PAIRED_FLAG=""
    [[ "$PAIRED" == "true" ]] && PAIRED_FLAG="--paired-ended --fragment-mean $FRAG_MEAN --fragment-st-dev $FRAG_SD"
    RNG_SEED_ROOT="${RNG_SEED_ROOT:-sv-validate}${SLURM_ARRAY_TASK_ID:+-rep$SLURM_ARRAY_TASK_ID}"
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
    echo "[1/5] Simulation complete."
fi

# ── Stage 2 & 3: Index reference, align ──────────────────────────────────
NORMAL_BAM="$OUTDIR/normal.bam"
TUMOR_BAM="$OUTDIR/tumor.bam"

index_and_align() {
    local in_fq="$1" out_bam="$2" sample="$3"
    if [[ -f "$out_bam" && -f "${out_bam}.bai" ]]; then
        echo "  $sample BAM already present — skipping."; return 0
    fi
    if [[ ! -f "${REFERENCE}.bwt.2bit.64" ]] && [[ ! -f "${REFERENCE}.bwt" ]]; then
        echo "  Indexing reference for BWA-MEM2 (one-time)..."
        bwa-mem2 index "$REFERENCE" 2>&1 | tail -3
    fi
    [[ -f "${REFERENCE}.fai" ]] || samtools faidx "$REFERENCE"
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

echo "[2/5] Aligning normal reads..."
index_and_align "$NORMAL_R1" "$NORMAL_BAM" "NORMAL"
echo "[3/5] Aligning tumor (merged) reads..."
index_and_align "${PREFIX}_merged_r1.fastq.gz" "$TUMOR_BAM" "TUMOR"

# ── Stage 3b: PRE-FLIGHT — inspect + extract the somatic SV truth ────────
# rneat's SVs are symbolic records (ALT like <DEL>/<DUP>/<INV>, INFO SVTYPE/
# SVLEN/END). Print a sample so we can confirm the schema, then extract the
# somatic SV subset truvari will score. NEAT_ORIGIN=somatic drops germline
# carry-through; the symbolic-ALT test keeps only SVs (not SNVs/indels).
TRUTH_SV="$OUTDIR/truth_sv.vcf.gz"
echo "[3b/5] Inspecting rneat SV truth..."
echo "  --- SV INFO header lines ---"
bcftools view -h "$TRUTH_VCF" | grep -iE 'INFO=<ID=(SVTYPE|SVLEN|END|NEAT_ORIGIN)' || \
    echo "  (no SVTYPE/SVLEN/END/NEAT_ORIGIN INFO headers found — schema may differ!)"
echo "  --- first symbolic-ALT (SV) records ---"
bcftools view -H "$TRUTH_VCF" | awk '$5 ~ /^</' | head -5 || true
echo "  --- SVTYPE counts in truth (all origins) ---"
bcftools view -H "$TRUTH_VCF" | grep -oE 'SVTYPE=[A-Z]+' | sort | uniq -c || true

if [[ -f "$TRUTH_SV" ]]; then
    echo "  somatic SV truth already extracted — skipping."
else
    # Permissive SV selector: symbolic ALT OR an SVTYPE INFO tag, AND somatic.
    # If the pre-flight above shows a different schema, adjust this -i expression.
    bcftools view -i '(ALT[0]~"<" || INFO/SVTYPE!=".") && INFO/NEAT_ORIGIN="somatic"' \
        -O z -o "$TRUTH_SV" "$TRUTH_VCF"
    bcftools sort -O z -o "$TRUTH_SV" "$TRUTH_SV" 2>/dev/null || true
    bcftools index -f -t "$TRUTH_SV"
fi
SV_TRUTH_N=$(bcftools view -H "$TRUTH_SV" 2>/dev/null | wc -l)
echo "  somatic SV truth records: $SV_TRUTH_N"
if [[ "$SV_TRUTH_N" -eq 0 ]]; then
    echo "ERROR: no somatic SV records extracted. Check the pre-flight schema dump" >&2
    echo "       above and adjust the SV selector. Raising SV_RATE_SCALE also helps." >&2
    exit 1
fi

# Pre-generate per-type truth subsets NOW, while bcftools (bioinf) is active, so
# the truvari scoring stage needs only truvari on its PATH (no bcftools in that
# env). Empty types are removed so the scoring loop skips them.
for svt in DEL DUP INV INS BND; do
    typed="$OUTDIR/truth_sv_${svt}.vcf.gz"
    bcftools view -i "INFO/SVTYPE=\"$svt\"" -O z -o "$typed" "$TRUTH_SV" 2>/dev/null || continue
    bcftools index -f -t "$typed" 2>/dev/null || true
    n=$(bcftools view -H "$typed" 2>/dev/null | wc -l)
    echo "  truth $svt: $n"
    [[ "$n" -eq 0 ]] && rm -f "$typed" "${typed}.tbi"
done

# ── Stage 4: Manta somatic SV calling ────────────────────────────────────
MANTA_DIR="$OUTDIR/manta"
MANTA_SV="$MANTA_DIR/results/variants/somaticSV.vcf.gz"
if [[ -f "$MANTA_SV" ]]; then
    echo "[4/5] Manta somaticSV already present — skipping."
else
    echo "[4/5] Running Manta (tumor/normal somatic SV)..."
    conda_activate "$MANTA_ENV"
    rm -rf "$MANTA_DIR"
    # Manta wants a bgzipped+faidx'd reference and indexed BAMs (we have both).
    configManta.py \
        --normalBam "$NORMAL_BAM" \
        --tumorBam  "$TUMOR_BAM" \
        --referenceFasta "$REFERENCE" \
        --runDir "$MANTA_DIR" 2>&1 | tail -5
    "$MANTA_DIR/runWorkflow.py" -j "$THREADS" 2>&1 | tail -10
    [[ -f "$MANTA_SV" ]] || { echo "ERROR: Manta produced no somaticSV.vcf.gz" >&2; exit 1; }
    echo "[4/5] Manta complete: $(bcftools view -H "$MANTA_SV" | wc -l) somatic SV calls."
fi

# ── Stage 5: Score with truvari (overall + per SV type) ──────────────────
echo "[5/5] Scoring SVs with truvari..."
conda_activate "$TRUVARI_ENV"

run_truvari() {  # <truth.vcf.gz> <out-subdir-name>
    local truth="$1" name="$2" out="$OUTDIR/truvari_$name"
    rm -rf "$out"   # truvari requires the output dir not pre-exist
    truvari bench -b "$truth" -c "$MANTA_SV" -f "$REFERENCE" -o "$out" --passonly 2>&1 | tail -3 || true
    if [[ -f "$out/summary.json" ]]; then
        python - "$out/summary.json" "$name" <<'PY'
import json, sys
s = json.load(open(sys.argv[1])); name = sys.argv[2]
print(f"  {name:<8} recall={s.get('recall',0):.3f}  precision={s.get('precision',0):.3f}  "
      f"f1={s.get('f1',0):.3f}  TP={s.get('TP-base',0)}  FN={s.get('FN',0)}  FP={s.get('FP',0)}")
PY
    else
        echo "  $name: truvari produced no summary.json (see $out)"
    fi
}

echo "  --- overall ---"
run_truvari "$TRUTH_SV" overall

echo "  --- per SV type (recall: did Manta recover each type rneat emitted?) ---"
# Score the per-type truth subsets pre-generated in stage 3b (bcftools not needed
# here — the truvari env may not have it).
for typed in "$OUTDIR"/truth_sv_*.vcf.gz; do
    [[ -f "$typed" ]] || continue            # no typed files matched
    svt=$(basename "$typed" .vcf.gz); svt=${svt#truth_sv_}
    run_truvari "$typed" "$svt"
done

# ── Summary ──────────────────────────────────────────────────────────────
cat <<EOF

════════════════════════════════════════════════════════════════
SV validation pipeline complete — Job $SLURM_JOB_ID

Inputs:
  Reference:    $REFERENCE
  Coverage:     ${TOTAL_COVERAGE}x  Purity: $PURITY  SV scale: $SV_RATE_SCALE
  Tumor model:  $TUMOR_MODEL

Key outputs in $OUTDIR:
  SV truth:     truth_sv.vcf.gz            ($SV_TRUTH_N somatic SVs)
  Manta calls:  manta/results/variants/somaticSV.vcf.gz
  Scores:       truvari_overall/summary.json + truvari_<TYPE>/

Interpretation: low overall recall with high per-type recall for some types =
a representation/format mismatch (truvari SVTYPE/SVLEN matching), not necessarily
a simulator bug — inspect the pre-flight schema dump. Uniformly low recall for a
type rneat emitted heavily = candidate SV-generation bug (the Phase-2 target).
EOF

# Persist report artifacts off scratch + record provenance/resources.
archive_run sv "$OUTDIR" truvari_overall/summary.json

#!/bin/bash
# SLURM job: stage REAL soybean data for the model-builder harness, self-consistently.
#
# Soy has no public pre-aligned BAM (unlike GIAB HG002), so we align a resequencing
# run ourselves — and CALL the VCF from that same BAM. That makes reference + BAM +
# VCF agree by construction, sidestepping the SoyBase-`Gm` vs NCBI-GenBank contig
# naming trap that would otherwise silently empty gen-mut-model. Everything is scoped
# to the reference's largest contig (chromosome 1) to keep it light, mirroring the
# HG002 chr20 approach.
#
# Pipeline: download reads (ENA, no sra-tools) -> bwa-mem2 align to soy.fa -> subset
# BAM + reference to the largest contig -> bcftools-call a VCF on it -> derive a FASTQ.
# Outputs (in $DATA_DIR): soy.chr.fa, soy.chr.bam, soy.chr.fastq.gz, soy.chr.vcf.gz.
# It prints the ready model_builders.sbatch command at the end.
#
# Prereqs: soy.fa staged (fetch_validation_data.sh DATA=soy); bwa-mem2 + bcftools
# (bioinf conda env), samtools/htslib modules. rneat NOT needed here.
#
# Usage:
#   sbatch scripts/delta/stage_soy.sh
#   SRR=SRR2163317 MAX_PAIRS=40000000 sbatch scripts/delta/stage_soy.sh   # MAX_PAIRS=0 -> all
#
# HEAVY: ~10 GB FASTQ download, ~28 GB bwa-mem2 index build, multi-hour alignment.

#SBATCH --job-name=rneat-stagesoy
#SBATCH --partition=cpu
#SBATCH --account=bhrd-delta-cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=80G
#SBATCH --time=12:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -euo pipefail

REPO_ROOT="${RNEAT_REPO:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/../.." && pwd)}}"
source "$REPO_ROOT/scripts/delta/lib_report.sh"   # $SCRATCH + setup_conda/conda_activate

D="${DATA_DIR:-$SCRATCH/neat_data/soy}"
REF="${REFERENCE:-$D/soy.fa}"
SRR="${SRR:-SRR2163317}"                 # SRP062245, paired, ~15x of the 1 Gb genome
MAX_PAIRS="${MAX_PAIRS:-40000000}"       # subsample for speed; 0 = align all 81M pairs
THREADS="${SLURM_CPUS_PER_TASK:-16}"
# ENA direct FASTQ (no sra-tools needed). Path layout: SRR216/007/SRR2163317/...
R1_URL="https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR216/007/${SRR}/${SRR}_1.fastq.gz"
R2_URL="https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR216/007/${SRR}/${SRR}_2.fastq.gz"

mkdir -p "$D"
module load samtools/1.22-cce19.0.0
module load htslib/1.22-gcc13.3.1
setup_conda
conda_activate bioinf                    # bwa-mem2 + bcftools

[[ -f "$REF" ]] || { echo "soy reference not staged: $REF (run fetch_validation_data.sh DATA=soy)" >&2; exit 1; }
echo "=== stage soy: SRR=$SRR  ref=$REF  max_pairs=$MAX_PAIRS  threads=$THREADS ==="

# ── 1. reads (ENA, resumable) ───────────────────────────────────────────────────
[[ -s "$D/${SRR}_1.fastq.gz" ]] || wget -c -O "$D/${SRR}_1.fastq.gz" "$R1_URL"
[[ -s "$D/${SRR}_2.fastq.gz" ]] || wget -c -O "$D/${SRR}_2.fastq.gz" "$R2_URL"

# Optional subsample (head keeps R1/R2 in the same order → pairs stay matched).
R1="$D/${SRR}_1.fastq.gz"; R2="$D/${SRR}_2.fastq.gz"
if [[ "$MAX_PAIRS" -gt 0 ]]; then
    if [[ ! -s "$D/sub_1.fastq.gz" ]]; then
        echo "subsampling first $MAX_PAIRS pairs..."
        zcat "$R1" | head -n $(( MAX_PAIRS * 4 )) | gzip > "$D/sub_1.fastq.gz"
        zcat "$R2" | head -n $(( MAX_PAIRS * 4 )) | gzip > "$D/sub_2.fastq.gz"
    fi
    R1="$D/sub_1.fastq.gz"; R2="$D/sub_2.fastq.gz"
fi

# ── 2. index reference (bwa-mem2 + faidx); -s guard against a 0-byte stub ────────
[[ -f "$REF.fai" ]] || samtools faidx "$REF"
if [[ ! -s "$REF.bwt.2bit.64" ]]; then
    echo "building bwa-mem2 index (one-time, ~28 GB RAM)..."
    bwa-mem2 index "$REF"
fi

# ── 3. align → coordinate-sorted BAM ─────────────────────────────────────────────
if [[ ! -s "$D/soy.full.bam" ]]; then
    echo "aligning..."
    bwa-mem2 mem -t "$THREADS" "$REF" "$R1" "$R2" 2> "$D/bwa.log" \
        | samtools sort -@ "$THREADS" -o "$D/soy.full.bam" -
    samtools index "$D/soy.full.bam"
fi

# ── 4. subset BAM + reference to the largest contig (a real chromosome) ─────────
CHR=$(sort -t$'\t' -k2,2 -nr "$REF.fai" | head -1 | cut -f1)
echo "largest contig: $CHR"
samtools faidx "$REF" "$CHR" > "$D/soy.chr.fa"; samtools faidx "$D/soy.chr.fa"
samtools view -b "$D/soy.full.bam" "$CHR" -o "$D/soy.chr.bam"; samtools index "$D/soy.chr.bam"

# ── 5. call a VCF FROM that BAM (matches soy.chr.fa by construction) ─────────────
bcftools mpileup -f "$D/soy.chr.fa" "$D/soy.chr.bam" 2>/dev/null \
    | bcftools call -mv -Oz -o "$D/soy.chr.vcf.gz"
bcftools index -t "$D/soy.chr.vcf.gz"

# ── 6. FASTQ for the seq-error model (reads on this contig) ─────────────────────
samtools fastq -n "$D/soy.chr.bam" 2>/dev/null | gzip > "$D/soy.chr.fastq.gz"

# ── reclaim the bulky whole-genome intermediates (chr subset is all we need) ────
rm -f "$D/soy.full.bam" "$D/soy.full.bam.bai" "$D/sub_1.fastq.gz" "$D/sub_2.fastq.gz"

nvar=$(bcftools view -H "$D/soy.chr.vcf.gz" | wc -l)
nreads=$(( $(zcat "$D/soy.chr.fastq.gz" | wc -l) / 4 ))
echo
echo "════════════════════════════════════════════════════════════════"
echo "Soy staged (contig $CHR): $nreads reads, $nvar called variants"
echo "Run the model-builders:"
echo "  REFERENCE=$D/soy.chr.fa INPUT_BAM=$D/soy.chr.bam \\"
echo "  INPUT_FASTQ=$D/soy.chr.fastq.gz INPUT_VCF=$D/soy.chr.vcf.gz \\"
echo "    sbatch scripts/delta/model_builders.sbatch"
echo "════════════════════════════════════════════════════════════════"

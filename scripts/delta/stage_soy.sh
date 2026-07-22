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
# FULL_GENOME=1 instead stages the WHOLE messy assembly: keeps the genome-wide
# soy.full.bam, calls a genome-wide VCF (soy.full.vcf.gz, fanned out per contig),
# and points the builders at soy.fa + soy.full.bam + the raw reads. That exposes the
# multi-contig / big-BAM / peak-RSS / OOM behavior the single-chromosome slice hides
# — the actual stress test. NOTE: the chr run reclaims soy.full.bam, so a full run
# re-aligns from scratch.
#
# Prereqs: soy.fa staged (fetch_validation_data.sh DATA=soy); bwa-mem2 + bcftools
# (bioinf conda env), samtools/htslib modules. eidolon NOT needed here.
#
# Usage:
#   sbatch scripts/delta/stage_soy.sh
#   SRR=SRR2163317 MAX_PAIRS=40000000 sbatch scripts/delta/stage_soy.sh   # MAX_PAIRS=0 -> all
#   FULL_GENOME=1 sbatch scripts/delta/stage_soy.sh                        # whole-genome stress
#
# HEAVY: ~10 GB FASTQ download, ~28 GB bwa-mem2 index build, multi-hour alignment.

#SBATCH --job-name=eidolon-stagesoy
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

REPO_ROOT="${EIDOLON_REPO:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/../.." && pwd)}}"
source "$REPO_ROOT/scripts/delta/lib_report.sh"   # $SCRATCH + setup_conda/conda_activate

D="${DATA_DIR:-$SCRATCH/neat_data/soy}"
REF="${REFERENCE:-$D/soy.fa}"
SRR="${SRR:-SRR2163317}"                 # SRP062245, paired, ~15x of the 1 Gb genome
MAX_PAIRS="${MAX_PAIRS:-40000000}"       # subsample for speed; 0 = align all 81M pairs
FULL_GENOME="${FULL_GENOME:-}"           # non-empty -> whole-genome stress inputs (no chr subset)
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
    # Require BOTH mates: a prior run that died mid-subsample (e.g. the pipefail/
    # SIGPIPE abort) can leave sub_1 without sub_2 — redo it rather than march on
    # to bwa-mem2 with a missing R2.
    if [[ ! -s "$D/sub_1.fastq.gz" || ! -s "$D/sub_2.fastq.gz" ]]; then
        echo "subsampling first $MAX_PAIRS pairs..."
        # `head` closes the pipe after N lines → `zcat` gets SIGPIPE (exit 141).
        # That is expected here, not an error, but `set -o pipefail` would promote
        # it to a pipeline failure and `set -e` would abort mid-subsample (leaving
        # sub_1 written but sub_2 never created). Disable pipefail for just these.
        set +o pipefail
        zcat "$R1" | head -n $(( MAX_PAIRS * 4 )) | gzip > "$D/sub_1.fastq.gz"
        zcat "$R2" | head -n $(( MAX_PAIRS * 4 )) | gzip > "$D/sub_2.fastq.gz"
        set -o pipefail
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

# ── FULL-GENOME mode: stage the whole messy assembly, skip the chr subset ───────
# The default (below) scopes to one clean chromosome for a fast correctness check.
# FULL_GENOME keeps the genome-wide BAM and calls a genome-wide VCF so the builders
# see every contig/scaffold, the full read volume, and real peak RSS.
if [[ -n "$FULL_GENOME" ]]; then
    ncontig=$(wc -l < "$REF.fai")
    echo "FULL_GENOME: whole-genome inputs across $ncontig contigs (no chr subset)"

    # Genome-wide VCF. bcftools mpileup is single-threaded per region, so fan out
    # one call per contig (in $REF.fai / genome order) and concat — this keeps a
    # ~1 Gb genome inside the wall clock. Each part is resumable ([-s] guard).
    if [[ ! -s "$D/soy.full.vcf.gz" ]]; then
        parts="$D/vcf_parts"; mkdir -p "$parts"
        echo "calling genome-wide VCF ($ncontig contigs, ${THREADS}-way)..."
        cut -f1 "$REF.fai" | xargs -P "$THREADS" -I CTG bash -c '
            set -o pipefail
            ctg="$1"; out="'"$parts"'/$ctg.vcf.gz"
            [[ -s "$out" ]] && exit 0
            if ! bcftools mpileup -r "$ctg" -f "'"$REF"'" "'"$D"'/soy.full.bam" 2>"$out.log" \
                 | bcftools call -mv -Oz -o "$out.tmp" 2>>"$out.log"; then
                echo "ERROR: per-contig VCF call failed for $ctg (see $out.log)" >&2
                rm -f "$out.tmp"; exit 1
            fi
            mv -f "$out.tmp" "$out"
        ' _ CTG || { echo "ERROR: a per-contig VCF call failed; aborting before concat (logs in $parts/*.log)" >&2; exit 1; }
        # A silently-dropped contig would vanish from the genome-wide VCF; require every
        # fai contig to have produced a part (variant-free → valid header-only file; only
        # a MISSING part is a failure).
        for c in $(cut -f1 "$REF.fai"); do
            [[ -s "$parts/$c.vcf.gz" ]] || { echo "ERROR: missing per-contig VCF part for '$c'" >&2; exit 1; }
        done
        # concat parts in genome (fai) order; same reference header → no -a needed
        cut -f1 "$REF.fai" | sed "s#^#$parts/#; s#\$#.vcf.gz#" > "$D/vcf_parts.list"
        bcftools concat -Oz -f "$D/vcf_parts.list" -o "$D/soy.full.vcf.gz"
        bcftools index -t "$D/soy.full.vcf.gz"
        rm -rf "$parts" "$D/vcf_parts.list"
    fi

    nvar=$(bcftools view -H "$D/soy.full.vcf.gz" | wc -l)
    nreads=$(( $(zcat "$R1" | wc -l) / 4 ))
    echo
    echo "════════════════════════════════════════════════════════════════"
    echo "Soy staged FULL GENOME: $ncontig contigs, ~$nreads read pairs, $nvar variants"
    echo "Run the model-builders (whole-genome stress):"
    echo "  REFERENCE=$REF INPUT_BAM=$D/soy.full.bam \\"
    echo "  INPUT_FASTQ=$R1 INPUT_VCF=$D/soy.full.vcf.gz \\"
    echo "    sbatch scripts/delta/model_builders.sbatch"
    echo "════════════════════════════════════════════════════════════════"
    exit 0
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

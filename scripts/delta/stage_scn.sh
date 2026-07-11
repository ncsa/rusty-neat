#!/bin/bash
# SLURM job: stage REAL soybean-cyst-nematode (SCN, Heterodera glycines) data for the
# model-builder harness, self-consistently — the SCN analogue of stage_soy.sh.
#
# Data is from João Gomes Viana (Matt Hudson lab). SCN strains differ in virulence to
# soybean, so — to mitigate reference bias — each population is aligned to its BEST-FIT
# strain reference (Walden et al. assemblies): MM26 (virulent) and PA3 (avirulent).
#   Reference (Genome) assemblies (NCBI datasets):
#     MM26 = GCA_040805935.1 (virulent)     PA3 = GCA_040805705.1 (avirulent)
#   Population Pool-seq reads (Illumina NextSeq 2000, paired; PRJNA1055977):
#     MM26A  = SRR27329602 (pool of 150 MM26 females)   ← default, matches MM26 ref
#     MM-BD3A= SRR27329600 (pool of 150 MM-BD3 females)  ← contrasting line
#
# NOTE — Pool-seq: the population reads are POOLED (150 individuals → allele frequencies,
# not diploid genotypes). This staging supports **Phase 1** only: build rneat models
# (mut-rate / seq-error / fragment / GC / bam) from real SCN Illumina data mapped to the
# fitting reference, and check standard fidelity. Reproducing the population ALLELE
# FREQUENCIES (João's Phase 2 ask) does NOT map onto rneat's single-sample model yet and
# is a separate design task — do not expect it from this script.
#
# Pipeline (mirrors stage_soy.sh): fetch assembly (NCBI datasets) -> fetch reads (ENA, no
# sra-tools) -> bwa-mem2 align -> genome-wide VCF called FROM that BAM (per-contig fan-out,
# so reference+BAM+VCF agree by construction) -> derive FASTQ. Prints the ready
# model_builders.sbatch command. SCN genome is small (~150 Mb) so no chromosome subset.
#
# Prereqs: bwa-mem2 + bcftools (bioinf conda env), samtools/htslib modules, and NCBI
# `datasets` on PATH (or pre-stage the reference FASTA yourself and pass REFERENCE=).
#
# Usage:
#   sbatch scripts/delta/stage_scn.sh                                   # MM26 ref + MM26A reads
#   ACC=GCA_040805705.1 SRR=SRR27329600 sbatch scripts/delta/stage_scn.sh   # PA3 ref + MM-BD3A reads
#   MAX_PAIRS=0 sbatch scripts/delta/stage_scn.sh                       # align ALL pairs (no subsample)
#
# HEAVY: ~45 GB FASTQ per run (subsampled by default); bwa-mem2 index build; multi-hour.

#SBATCH --job-name=rneat-stagescn
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

D="${DATA_DIR:-$SCRATCH/neat_data/scn}"
ACC="${ACC:-GCA_040805935.1}"            # MM26 (virulent) assembly; PA3 = GCA_040805705.1
SRR="${SRR:-SRR27329602}"                # MM26A pool reads; MM-BD3A = SRR27329600
REF="${REFERENCE:-$D/${ACC}.fa}"
MAX_PAIRS="${MAX_PAIRS:-15000000}"       # ~30x of the ~150 Mb genome; 0 = align all pairs
THREADS="${SLURM_CPUS_PER_TASK:-16}"

mkdir -p "$D"
module load samtools/1.22-cce19.0.0
module load htslib/1.22-gcc13.3.1
setup_conda
conda_activate bioinf                    # bwa-mem2 + bcftools

# ENA direct-FASTQ URL from an SRR accession (no sra-tools). ENA nests runs under a
# subdir derived from the accession length: 7 digits -> 00<last1>, 8 -> 0<last2>,
# 9 -> <last3>; a 6-digit (9-char) accession has no subdir.
ena_url() {
    local acc="$1" mate="$2"
    # Split from the line above: in a single `local`, all RHS are expanded before any
    # assignment, so ${acc:0:6}/${#acc} would see the OLD (empty) acc.
    local dir6="${acc:0:6}" n=${#acc} sub=""
    if   (( n == 10 )); then sub="00${acc: -1}"
    elif (( n == 11 )); then sub="0${acc: -2}"
    elif (( n == 12 )); then sub="${acc: -3}"
    fi
    echo "https://ftp.sra.ebi.ac.uk/vol1/fastq/${dir6}/${sub:+$sub/}${acc}/${acc}_${mate}.fastq.gz"
}

# Robustly fetch one mate ($1 = 1|2) of $SRR into $D. ENA truncates large HTTPS
# transfers mid-stream (TLS drop / transient 403 on resume), so: resume with retry
# backoff, then verify the result against the ENA API's size + md5 and re-resume /
# re-download until it matches. $ENA_META (queried once below) is a single TSV line
# run_accession<TAB>fastq_ftp<TAB>fastq_md5<TAB>fastq_bytes — note the API ALWAYS prepends
# run_accession as col 1 even though we don't request it, so the fields we want are $2/$3/$4;
# each of those is ';'-joined per mate. If the API is unreachable ENA_META is empty and we
# fall back to the computed URL + wget's exit code (unverified). Guards against a truncated
# transfer AND against silently accepting one.
fetch_read() {
    local mate="$1" dst="$D/${SRR}_${mate}.fastq.gz" url md5 bytes rc tries=0 got
    url="$(  awk -F'\t' 'NR==1{print $2}' <<<"$ENA_META" | cut -d';' -f"$mate")"
    md5="$(  awk -F'\t' 'NR==1{print $3}' <<<"$ENA_META" | cut -d';' -f"$mate")"
    bytes="$(awk -F'\t' 'NR==1{print $4}' <<<"$ENA_META" | cut -d';' -f"$mate")"
    if [[ -n "$url" ]]; then url="https://$url"; else url="$(ena_url "$SRR" "$mate")"; fi
    while :; do
        tries=$((tries + 1))
        if wget -c --tries=20 --waitretry=30 --retry-connrefused --read-timeout=60 --timeout=60 \
                -O "$dst" "$url"; then rc=0; else rc=$?; fi
        if [[ -n "$bytes" ]]; then
            got="$(stat -c%s "$dst" 2>/dev/null || echo 0)"
            if [[ "$got" != "$bytes" ]]; then
                echo "  $dst size $got != $bytes (attempt $tries) — resuming" >&2
            elif [[ -n "$md5" ]] && [[ "$(md5sum "$dst" | cut -d' ' -f1)" != "$md5" ]]; then
                echo "  $dst md5 mismatch (attempt $tries) — re-downloading from scratch" >&2
                rm -f "$dst"
            else
                echo "  $dst verified (size${md5:+ + md5})"
                return 0
            fi
        else
            if [[ $rc -eq 0 && -s "$dst" ]]; then
                echo "  $dst downloaded (ENA API gave no size/md5 → UNVERIFIED)"
                return 0
            fi
            echo "  $dst wget rc=$rc (attempt $tries)" >&2
        fi
        if (( tries >= 6 )); then
            echo "ERROR: $dst failed to fetch/verify after $tries attempts" >&2
            return 1
        fi
        sleep 30
    done
}

echo "=== stage SCN: ref=$ACC  reads=$SRR  max_pairs=$MAX_PAIRS  threads=$THREADS ==="

# ── 1. reference assembly (NCBI datasets → single FASTA) ─────────────────────────
if [[ ! -s "$REF" ]]; then
    command -v datasets >/dev/null 2>&1 || {
        echo "NCBI 'datasets' not on PATH — install it, or pre-stage the FASTA and pass REFERENCE=$REF" >&2
        exit 1; }
    echo "downloading assembly $ACC via NCBI datasets..."
    datasets download genome accession "$ACC" --include genome --filename "$D/${ACC}.zip"
    rm -rf "$D/${ACC}_dl"; unzip -q -o -d "$D/${ACC}_dl" "$D/${ACC}.zip"
    # concat the assembly's .fna(s) into one reference FASTA
    find "$D/${ACC}_dl" -name '*.fna' -print0 | xargs -0 cat > "$REF"
    rm -rf "$D/${ACC}_dl" "$D/${ACC}.zip"
fi
[[ -f "$REF.fai" ]] || samtools faidx "$REF"

# ── 2. reads (ENA, resumable + verified against the ENA API) ─────────────────────
# One API call → canonical URLs + md5 + byte sizes for both mates; fetch_read resumes
# on drops and re-checks size/md5 so a broken transfer self-heals instead of aborting
# the job or being silently accepted as complete (the old bare `-s` guard did the latter).
ENA_META="$(curl -s --max-time 60 \
    "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${SRR}&result=read_run&fields=fastq_ftp,fastq_md5,fastq_bytes" \
    2>/dev/null | tail -n +2 | head -1 || true)"
[[ -n "$ENA_META" ]] || echo "WARN: ENA API gave nothing for $SRR — using computed URL, no size/md5 verification" >&2
fetch_read 1
fetch_read 2

R1="$D/${SRR}_1.fastq.gz"; R2="$D/${SRR}_2.fastq.gz"
if [[ "$MAX_PAIRS" -gt 0 ]]; then
    if [[ ! -s "$D/sub_1.fastq.gz" || ! -s "$D/sub_2.fastq.gz" ]]; then
        echo "subsampling first $MAX_PAIRS pairs..."
        # head closes the pipe → zcat gets SIGPIPE (141); expected, so drop pipefail for
        # just these two so `set -e` doesn't abort a step that actually succeeded.
        set +o pipefail
        zcat "$R1" | head -n $(( MAX_PAIRS * 4 )) | gzip > "$D/sub_1.fastq.gz"
        zcat "$R2" | head -n $(( MAX_PAIRS * 4 )) | gzip > "$D/sub_2.fastq.gz"
        set -o pipefail
    fi
    R1="$D/sub_1.fastq.gz"; R2="$D/sub_2.fastq.gz"
fi

# ── 3. index reference (bwa-mem2); -s guard against a 0-byte stub ─────────────────
if [[ ! -s "$REF.bwt.2bit.64" ]]; then
    echo "building bwa-mem2 index (one-time)..."
    bwa-mem2 index "$REF"
fi

# ── 4. align → coordinate-sorted BAM ─────────────────────────────────────────────
BAM="$D/${SRR}.bam"
if [[ ! -s "$BAM" ]]; then
    echo "aligning..."
    bwa-mem2 mem -t "$THREADS" "$REF" "$R1" "$R2" 2> "$D/${SRR}_bwa.log" \
        | samtools sort -@ "$THREADS" -o "$BAM" -
    samtools index "$BAM"
fi

# ── 5. call a genome-wide VCF FROM that BAM (agrees with $REF by construction) ────
# bcftools mpileup is single-threaded per region → fan out one call per contig (fai
# order) and concat. Each part is resumable ([-s] guard).
VCF="$D/${SRR}.vcf.gz"
if [[ ! -s "$VCF" ]]; then
    parts="$D/${SRR}_vcf_parts"; mkdir -p "$parts"
    ncontig=$(wc -l < "$REF.fai")
    echo "calling genome-wide VCF ($ncontig contigs, ${THREADS}-way)..."
    cut -f1 "$REF.fai" | xargs -P "$THREADS" -I CTG bash -c '
        ctg="$1"; out="'"$parts"'/$ctg.vcf.gz"
        [[ -s "$out" ]] && exit 0
        bcftools mpileup -r "$ctg" -f "'"$REF"'" "'"$BAM"'" 2>/dev/null \
            | bcftools call -mv -Oz -o "$out" 2>/dev/null
    ' _ CTG
    cut -f1 "$REF.fai" | sed "s#^#$parts/#; s#\$#.vcf.gz#" > "$D/${SRR}_parts.list"
    bcftools concat -Oz -f "$D/${SRR}_parts.list" -o "$VCF"
    bcftools index -t "$VCF"
    rm -rf "$parts" "$D/${SRR}_parts.list"
fi

# ── 6. FASTQ for the seq-error model (aligned reads) ─────────────────────────────
FQ="$D/${SRR}.fastq.gz"
[[ -s "$FQ" ]] || samtools fastq -n "$BAM" 2>/dev/null | gzip > "$FQ"

# ── reclaim the bulky subsample intermediates ────────────────────────────────────
rm -f "$D/sub_1.fastq.gz" "$D/sub_2.fastq.gz"

nvar=$(bcftools view -H "$VCF" | wc -l)
nreads=$(( $(zcat "$FQ" | wc -l) / 4 ))
echo
echo "════════════════════════════════════════════════════════════════"
echo "SCN staged: ref $ACC, reads $SRR — $nreads reads, $nvar called variants"
echo "Run the model-builders:"
echo "  REFERENCE=$REF INPUT_BAM=$BAM \\"
echo "  INPUT_FASTQ=$FQ INPUT_VCF=$VCF \\"
echo "    sbatch scripts/delta/model_builders.sbatch"
echo "════════════════════════════════════════════════════════════════"

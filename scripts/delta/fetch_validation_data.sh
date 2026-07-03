#!/usr/bin/env bash
# Stage validation datasets for the model-builder + gen-reads harnesses into
# $NEAT_DATA. Three tracks, run whichever you need:
#
#   DATA=hg002  bash scripts/delta/fetch_validation_data.sh   # real human (open)
#   DATA=soy    bash scripts/delta/fetch_validation_data.sh   # messy polyploid (open)
#   DATA=scn    bash scripts/delta/fetch_validation_data.sh   # weird-genetics (colleague)
#
# Design: auto-download only what has a STABLE, verified URL (truth VCFs, reference
# assemblies). For raw reads (BAM/FASTQ) it points at the canonical index instead of
# hardcoding a filename that could 404 — you paste the exact URL (or use the
# align-your-own recipe), keeping everything reference-consistent. Nothing here is
# committed to a big node request, so a wrong URL just fails a wget, not a run.
set -uo pipefail

REPO_ROOT="${RNEAT_REPO:-$(cd "$(dirname "$0")/../.." && pwd)}"
source "$REPO_ROOT/scripts/delta/lib_report.sh" 2>/dev/null || true
NEAT_DATA="${NEAT_DATA:-${SCRATCH:-$HOME}/neat_data}"
DATA="${DATA:-}"
module load samtools/1.22-cce19.0.0 2>/dev/null || true

note() { printf '\n\033[1m%s\033[0m\n' "$*"; }
have() { command -v "$1" >/dev/null 2>&1; }
get()  { # get URL DEST  — resumable, skips if present and non-empty
    local url="$1" dst="$2"
    if [[ -s "$dst" ]]; then echo "  have $(basename "$dst")"; return 0; fi
    mkdir -p "$(dirname "$dst")"
    echo "  fetching $(basename "$dst") ..."
    wget -q --show-progress -c -O "$dst" "$url" || { echo "  FAILED: $url" >&2; return 1; }
}

case "$DATA" in

hg002) # ── Track A: GIAB HG002, GRCh38 (all 5 builders, no gating) ──────────────
    D="$NEAT_DATA/hg002"; base="https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38"
    note "[hg002] curated truth VCF (v4.2.1, GRCh38) — feeds gen-mut-model"
    get "$base/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"            "$D/HG002_GRCh38_v4.2.1_benchmark.vcf.gz"
    get "$base/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi"        "$D/HG002_GRCh38_v4.2.1_benchmark.vcf.gz.tbi"
    get "$base/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed" "$D/HG002_GRCh38_v4.2.1_highconf.bed"
    note "[hg002] reference: must be GRCh38 to match the VCF above"
    echo "  point REFERENCE=\$NEAT_DATA/GRCh38.fa (the same GRCh38 you already sim against)."
    note "[hg002] reads (BAM for frag/gc/bam-models; FASTQ for seq-error) — pick from the GIAB index:"
    echo "  Canonical index (exact, GRCh38-aligned URLs): https://github.com/genome-in-a-bottle/giab_data_indexes"
    echo "    -> AshkenazimTrio/alignment.index.*  (BAM/CRAM)   +  sequence.index.*  (FASTQ)"
    echo "  Then either set HG002_BAM=<url> HG002_FASTQ_R1=<url> HG002_FASTQ_R2=<url> and re-run,"
    echo "  or align a FASTQ subset yourself for guaranteed GRCh38 consistency:"
    echo "     bwa-mem2 mem -t \$THREADS \$NEAT_DATA/GRCh38.fa R1.fq.gz R2.fq.gz | samtools sort -o \$D/hg002.grch38.bam"
    echo "     samtools index \$D/hg002.grch38.bam"
    [[ -n "${HG002_BAM:-}" ]]      && get "$HG002_BAM"      "$D/hg002.bam"
    [[ -n "${HG002_FASTQ_R1:-}" ]] && get "$HG002_FASTQ_R1" "$D/hg002_R1.fastq.gz"
    [[ -n "${HG002_FASTQ_R2:-}" ]] && get "$HG002_FASTQ_R2" "$D/hg002_R2.fastq.gz"
    ;;

soy) # ── Track B: soybean Wm82 (scale/sharding + messy polyploid ref) ───────────
    D="$NEAT_DATA/soy"; SOY_ACC="${SOY_ACC:-GCA_000004515.5}"   # Wm82.a4; a6 T2T = ask Phytozome
    note "[soy] reference assembly $SOY_ACC via NCBI datasets"
    if have datasets; then
        ( cd "$D" 2>/dev/null || { mkdir -p "$D"; cd "$D"; }
          datasets download genome accession "$SOY_ACC" --include genome --filename soy.zip \
            && unzip -o soy.zip -d soy_unzip >/dev/null \
            && find soy_unzip -name '*.fna' -exec mv {} "$D/soy.fa" \; \
            && echo "  -> $D/soy.fa" )
    else
        echo "  NCBI 'datasets' CLI not found. Options:"
        echo "    conda install -c conda-forge ncbi-datasets-cli   # then re-run"
        echo "    or Phytozome v13 (Wm82.a4/a6, needs free JGI login): https://phytozome-next.jgi.doe.gov"
        echo "    or SoyBase: https://www.soybase.org/collections/"
    fi
    note "[soy] resequencing reads (align -> BAM for frag/gc/bam-models; FASTQ for seq-error)"
    echo "  SRA study SRP062245 (106 lines ~17x). Pick one run, e.g. SOY_SRR=SRRxxxxxxx, then:"
    echo "    prefetch \$SOY_SRR && fasterq-dump --split-files -O \$D \$SOY_SRR   # sra-tools"
    echo "    bwa-mem2 mem -t \$THREADS \$D/soy.fa \$D/\${SOY_SRR}_1.fastq \$D/\${SOY_SRR}_2.fastq | samtools sort -o \$D/soy.bam"
    note "[soy] variants (gen-mut-model): SoyBase resequencing/pangenome VCFs — https://www.soybase.org/collections/"
    echo "  IMPORTANT: match the VCF's assembly version (a1/a2/a4) to soy.fa, or contigs won't line up."
    ;;

scn) # ── Track C: soybean cyst nematode (weird genetics; mostly via colleague) ──
    D="$NEAT_DATA/scn"; mkdir -p "$D"
    note "[scn] Heterodera glycines — ~158 Mb, 9 chr, HIGH heterozygosity, 34% repeats"
    echo "  Public start (reference FASTA): SCNBase TN10 chromosome-level assembly (157,982,452 bp)"
    echo "    portal: https://scnbase.org/   (also search NCBI Assembly 'Heterodera glycines')"
    echo "    once you have a GCA accession:  datasets download genome accession GCA_XXXXXXXXX.X --include genome"
    echo
    echo "  Best data = from João P. Gomes Viana (NCSA/UIUC) — pangenome paper, BMC Genomics 2026"
    echo "    doi 10.1186/s12864-025-12493-x ; PMC12892471 ; preprint rs-7716171"
    echo "  Ask him for: (1) a population reference FASTA (PA3 or MM26 HiFi assembly),"
    echo "    (2) an aligned BAM (-> frag/gc/bam-models), (3) the VCF of pangenome variants"
    echo "    (9,495 SNVs / 4,567 indels / 3,344 large PAVs) AND which reference it was called against."
    echo "  Raw reads (if pulling yourself): the paper's Supplementary Table S1 lists the BioProject IDs."
    ;;

*)
    note "usage: DATA=hg002|soy|scn  bash scripts/delta/fetch_validation_data.sh"
    echo "  hg002 — real human, all 5 model-builders + round-trip (open; auto-fetches truth VCF)"
    echo "  soy   — messy polyploid, scale/sharding (open; auto-fetches reference)"
    echo "  scn   — soybean cyst nematode, weird genetics (mostly via colleague João Gomes Viana)"
    echo "  NEAT_DATA=$NEAT_DATA"
    exit 0 ;;
esac

# Index any FASTA we landed so it's ready for gen-reads / the builders.
for fa in "$NEAT_DATA"/*/*.fa; do [[ -f "$fa" && ! -f "$fa.fai" ]] && { echo "  indexing $(basename "$fa")"; samtools faidx "$fa" 2>/dev/null || true; }; done
note "done — staged under $NEAT_DATA/$DATA"

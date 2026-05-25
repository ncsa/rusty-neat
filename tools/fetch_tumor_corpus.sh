#!/usr/bin/env bash
# Convert the TCGA MC3 PUBLIC MAF into a VCF that rneat gen-mut-model
# can consume directly, then (optionally) train a pan-cancer tumor model.
#
# What MC3 is:
#   The Multi-Center Mutation Calling (MC3) project published consensus
#   somatic SNV/indel calls across every TCGA cancer project as a single
#   open-access MAF — ~3.6M mutations from ~10k tumors. It's the
#   community-standard public-tier source for pan-cancer mutation
#   training corpora.
#
# Where to get it:
#   gdc.cancer.gov/about-data/publications/pancanatlas → MC3 PUBLIC MAF
#   (typical filename: mc3.v0.2.8.PUBLIC.maf.gz, ~750 MB)
#
# What this script does:
#   1. Filter to PASS records (MC3 uses "." as the PASS equivalent; 80% of rows).
#   2. Keep SNPs only (drops the 5% indel records — see "Known limitations" below).
#   3. Re-shape MAF columns into VCF columns.
#   4. Add `chr` prefix to chromosomes by default (rneat's bundled defaults
#      are gnomAD-SV-derived which uses chr-prefixed contigs).
#   5. bgzip-compress the output (rneat handles BGZF via MultiGzDecoder).
#   6. Optional: drive `rneat gen-mut-model` against the converted VCF.
#
# Known limitations (intentional for v1):
#   - SNPs only. Indels need reference-FASTA anchor lookup (MAF uses "-"
#     for the absent allele; VCF requires both REF and ALT to be present
#     bases with a shared anchor). Adding indel support is a v2 follow-up;
#     see issue #186 thread.
#   - Pan-cancer model. No per-tissue split. The Tumor_Sample_Barcode TSS
#     code maps to project ID via a separate GDC code table; bundling
#     that mapping cleanly is also v2.
#
# Usage:
#   tools/fetch_tumor_corpus.sh --mc3-maf ~/code/data/mc3.v0.2.8.PUBLIC.maf.gz
#
# Run with --help for the full flag list.

set -euo pipefail

# ── Defaults ────────────────────────────────────────────────────────────
MC3_MAF=""
OUT_PREFIX="tcga_mc3_pancancer"
OUT_DIR="."
CHR_PREFIX="true"
REFERENCE=""
TRAIN_AFTER="false"
RNEAT_BIN="rneat"

usage() {
    cat <<'EOF'
Usage: fetch_tumor_corpus.sh --mc3-maf <path> [options]

Required:
  --mc3-maf      Path to mc3.v0.2.8.PUBLIC.maf.gz (downloaded from
                 gdc.cancer.gov/about-data/publications/pancanatlas)

Output:
  --out-dir      Directory for outputs (default: current directory)
  --out-prefix   Stem for output filenames (default: tcga_mc3_pancancer)
                 Produces <out-prefix>.vcf.gz

Conversion knobs:
  --no-chr-prefix   Emit bare chromosome names (1/2/X/Y) instead of
                    chr-prefixed. Default is to add `chr` so the output
                    matches rneat's gnomAD-SV-derived bundled defaults.

Optional training step:
  --train           After conversion, run `rneat gen-mut-model` to produce
                    a tumor mutation model JSON. Requires --reference.
  --reference       Path to the GRCh38 reference FASTA (only used when
                    --train is set; needs to use the same chr-prefix
                    convention as the converted VCF).
  --rneat-bin       Path to the rneat binary (default: "rneat" on PATH)

  -h, --help        Print this message
EOF
}

# ── Arg parsing ─────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --mc3-maf)         MC3_MAF="$2"; shift 2 ;;
        --out-dir)         OUT_DIR="$2"; shift 2 ;;
        --out-prefix)      OUT_PREFIX="$2"; shift 2 ;;
        --no-chr-prefix)   CHR_PREFIX="false"; shift ;;
        --train)           TRAIN_AFTER="true"; shift ;;
        --reference)       REFERENCE="$2"; shift 2 ;;
        --rneat-bin)       RNEAT_BIN="$2"; shift 2 ;;
        -h|--help)         usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

# ── Validation ──────────────────────────────────────────────────────────
[[ -z "$MC3_MAF"  ]] && { echo "Missing --mc3-maf" >&2; exit 2; }
[[ ! -f "$MC3_MAF" ]] && { echo "MC3 MAF not found: $MC3_MAF" >&2; exit 2; }

if [[ "$TRAIN_AFTER" == "true" ]]; then
    [[ -z "$REFERENCE"  ]] && { echo "--train requires --reference" >&2; exit 2; }
    [[ ! -f "$REFERENCE" ]] && { echo "Reference not found: $REFERENCE" >&2; exit 2; }
    if ! command -v "$RNEAT_BIN" >/dev/null 2>&1; then
        echo "rneat binary not found on PATH (override with --rneat-bin)" >&2
        exit 2
    fi
fi

if ! command -v bgzip >/dev/null 2>&1; then
    echo "bgzip not found on PATH (install htslib / samtools)" >&2
    exit 2
fi

mkdir -p "$OUT_DIR"
OUT_VCF="${OUT_DIR}/${OUT_PREFIX}.vcf.gz"

# ── Convert MAF → VCF ───────────────────────────────────────────────────
# MC3 MAF columns we care about (1-based):
#   5  Chromosome           e.g. "10", "X"
#   6  Start_Position
#   10 Variant_Type         SNP / DEL / INS / DNP / TNP / ONP
#   11 Reference_Allele     A/C/G/T (or "-" for INS, ignored in this recipe)
#   13 Tumor_Seq_Allele2    A/C/G/T for SNP
#   113 FILTER              "." = pass; any other value = some PanCanAtlas flag
#
# Skipped variant types in v1: DEL, INS (would need reference anchor lookup),
# DNP/TNP/ONP (multi-base substitutions, rare; not yet supported by rneat).

echo ">> Converting MC3 MAF → VCF (SNP-only, FILTER=. only)..."

# Pass the chr-prefix as a plain string to awk (no shell quoting around the
# value — `awk -v PFX=chr` is correct, `awk -v PFX='"chr"'` makes PFX equal
# to the literal six-character string "chr").
if [[ "$CHR_PREFIX" == "true" ]]; then
    CHR_PFX_VAL="chr"
else
    CHR_PFX_VAL=""
fi

# Two-step: emit header + body to a temp file, then sort the body and
# bgzip. We sort because MC3 is grouped by donor (Tumor_Sample_Barcode),
# not by genomic position — without sorting the output VCF can't be
# tabix-indexed and chunked tools choke. The sort is the slowest step
# (~30s for ~2.8M records) but the resulting file is clean.
TMP_BODY="$(mktemp)"
trap 'rm -f "$TMP_BODY"' EXIT

zcat "$MC3_MAF" | awk -F'\t' -v PFX="$CHR_PFX_VAL" '
    NR == 1 { next }                                # header line
    $113 != "."  { next }                           # not PASS
    $10  != "SNP" { next }                          # drop indels & multi-base subs
    $5 == "" || $6 == "" || $11 == "" || $13 == "" { next }  # malformed
    $11 == "-" || $13 == "-" { next }               # defense in depth — should not happen for SNP
    {
        chrom = PFX $5
        pos   = $6
        ref   = $11
        alt   = $13
        # Use "." for ID (TCGA barcodes are not VCF IDs and gen-mut-model
        # ignores them anyway). QUAL is also "." — MC3 does not preserve
        # per-call quality in a portable way.
        printf "%s\t%s\t.\t%s\t%s\t.\tPASS\t.\tGT\t0/1\n", chrom, pos, ref, alt
    }
' > "$TMP_BODY"

echo "    Sorting by (chrom, pos) for tabix-compatibility..."
{
    cat <<'HDR'
##fileformat=VCFv4.2
##source=tools/fetch_tumor_corpus.sh
##reference=TCGA-MC3-PUBLIC-v0.2.8
##INFO=<ID=MC3_FILTER,Number=1,Type=String,Description="Original MC3 FILTER value (always . in this VCF since we kept only those)">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype (synthetic 0/1 — somatic mutations are conventionally heterozygous in the tumor cell)">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NEAT_tumor_corpus
HDR
    LC_ALL=C sort -t$'\t' -k1,1 -k2,2n "$TMP_BODY"
} | bgzip -c > "$OUT_VCF"

# ── Tabix index (cheap; lets users bcftools-view -r quickly) ───────────
if command -v tabix >/dev/null 2>&1; then
    tabix -p vcf "$OUT_VCF"
    echo "    Indexed: ${OUT_VCF}.tbi"
fi

# ── Summary ─────────────────────────────────────────────────────────────
N_RECORDS=$(zcat "$OUT_VCF" | grep -cv '^#' || true)
OUT_SIZE=$(du -h "$OUT_VCF" | cut -f1)
cat <<EOF

────────────────────────────────────────────────────────────────
TCGA MC3 → VCF conversion complete

Input:     $MC3_MAF
Output:    $OUT_VCF
           ${OUT_SIZE} compressed, ${N_RECORDS} SNP records (PASS-filter)

Chromosome convention: $([ "$CHR_PREFIX" == "true" ] && echo "chr-prefixed (chr1, chr2, ..., chrX)" || echo "bare (1, 2, ..., X)")

Excluded from this VCF (intentional v1 simplifications):
  - Indels (DEL/INS):    ~173k records       — need reference anchor lookup
  - Multi-base subs:     ~174 records        — DNP/TNP/ONP, rarely modeled
  - Non-PASS records:    ~700k records       — population-AF flags etc.
EOF

# ── Optional training step ──────────────────────────────────────────────
if [[ "$TRAIN_AFTER" == "true" ]]; then
    echo ""
    echo ">> Training pan-cancer tumor mutation model..."
    MODEL_OUT="${OUT_DIR}/${OUT_PREFIX}_model.json.gz"
    TRAIN_CFG="${OUT_DIR}/${OUT_PREFIX}_train.yml"
    cat > "$TRAIN_CFG" <<EOF
reference: $REFERENCE
vcf_file: $OUT_VCF
output_file: $MODEL_OUT
overwrite_output: true
EOF
    "$RNEAT_BIN" gen-mut-model -c "$TRAIN_CFG"
    cat <<EOF

  Trained model:  $MODEL_OUT
  Training config: $TRAIN_CFG

  Plug into the cancer simulator with:
    tools/cancer_simulate.sh \\
        --reference $REFERENCE \\
        --output-prefix /path/to/out/sample \\
        --total-coverage 30 \\
        --purity 0.7 \\
        --tumor-model $MODEL_OUT
EOF
fi

echo "────────────────────────────────────────────────────────────────"
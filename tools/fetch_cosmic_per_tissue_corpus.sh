#!/usr/bin/env bash
# Build a per-tissue COSMIC SNV/indel corpus VCF, then (optionally) train a
# per-tissue tumor mutation model. The per-tissue counterpart to
# tools/fetch_cosmic_corpus.sh (which is pan-cancer).
#
# Why three input files (and how they join):
#   The Cosmic_GenomeScreensMutant *VCF* is tissue-aggregated — its ID column is
#   the COSV mutation ID, but it carries no sample/tissue/phenotype field, so it
#   alone can't be split by tissue. The per-sample TSV export does carry tissue,
#   indirectly:
#
#     Cosmic_GenomeScreensMutant_<v>_GRCh38.tsv.gz   COSMIC_PHENOTYPE_ID (col 6)
#                                                    + GENOMIC_MUTATION_ID = COSV (col 7)
#         joined on COSMIC_PHENOTYPE_ID to
#     Cosmic_Classification_<v>_GRCh38.tsv.gz        PRIMARY_SITE (col 2) = tissue
#
#   So we read the classification file to get the set of phenotype IDs for the
#   requested tissue, stream the (large) mutant TSV to collect every COSV ID
#   observed in those phenotypes, and finally filter the *VCF* down to those COSV
#   IDs. We take the alleles from the VCF — not the TSV — because the VCF already
#   carries proper VCF-anchored alleles, whereas the TSV uses COSMIC's raw allele
#   convention (indels are not left-anchored the way gen-mut-model / VCF expect).
#
# Downstream filtering matches fetch_cosmic_corpus.sh exactly:
#   - drop complex variants (REF>1 AND ALT>1; gen-mut-model has no Complex arm)
#   - drop ambiguous (non-ACGT) bases
#   - add `chr` prefix (MT -> chrM) by default
#   - sort + dedupe by (chrom,pos,ref,alt) — collapses multi-transcript rows
#
# Usage:
#   tools/fetch_cosmic_per_tissue_corpus.sh \
#       --tissue breast \
#       --cosmic-vcf     ~/code/data/cosmic_v104/Cosmic_GenomeScreensMutant_v104_GRCh38.vcf.gz \
#       --mutant-tsv     ~/code/data/cosmic_v104/Cosmic_GenomeScreensMutant_Tsv_v104_GRCh38.tar \
#       --classification ~/code/data/cosmic_v104/Cosmic_Classification_v104_GRCh38.tsv.gz \
#       --out-dir /tmp/cosmic_tissue --out-prefix cosmic_breast
#
# Run with --help for the full flag list.

set -euo pipefail

# ── Defaults ────────────────────────────────────────────────────────────
TISSUE=""
COSMIC_VCF=""
MUTANT_TSV=""
CLASSIFICATION=""
OUT_PREFIX=""
OUT_DIR="."
CHR_PREFIX="true"
REFERENCE=""
TRAIN_AFTER="false"
RNEAT_BIN="rneat"

usage() {
    cat <<'EOF'
Usage: fetch_cosmic_per_tissue_corpus.sh --tissue <site> --cosmic-vcf <vcf> \
           --mutant-tsv <tsv|tar> --classification <tsv.gz> [options]

Required:
  --tissue          COSMIC PRIMARY_SITE to keep (e.g. breast, skin, lung).
                    Must match a value in the Classification file's PRIMARY_SITE
                    column exactly (lowercase, underscores).
  --cosmic-vcf      Cosmic_GenomeScreensMutant_<v>_GRCh38.vcf.gz (alleles).
  --mutant-tsv      Cosmic_GenomeScreensMutant_<v>_GRCh38.tsv.gz, OR the .tar that
                    wraps it (the COSV<->phenotype map). Streamed, never extracted.
  --classification  Cosmic_Classification_<v>_GRCh38.tsv.gz (phenotype -> site).

Output:
  --out-dir         Directory for outputs (default: current directory)
  --out-prefix      Stem for output filenames (default: cosmic_<tissue>)
                    Produces <out-prefix>.vcf.gz

Conversion knobs:
  --no-chr-prefix   Emit bare chromosome names (1/2/X/Y/MT) instead of chr-prefixed.
                    Default adds `chr` (MT -> chrM), matching rneat's bundled models.

Optional training step:
  --train           After conversion, run `rneat gen-mut-model`. Requires --reference.
  --reference       GRCh38 reference FASTA (same chr-prefix convention as the VCF).
  --rneat-bin       Path to the rneat binary (default: "rneat" on PATH)

  -h, --help        Print this message
EOF
}

# ── Arg parsing ─────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --tissue)          TISSUE="$2"; shift 2 ;;
        --cosmic-vcf)      COSMIC_VCF="$2"; shift 2 ;;
        --mutant-tsv)      MUTANT_TSV="$2"; shift 2 ;;
        --classification)  CLASSIFICATION="$2"; shift 2 ;;
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
[[ -z "$TISSUE"         ]] && { echo "Missing --tissue" >&2; exit 2; }
[[ -z "$COSMIC_VCF"     ]] && { echo "Missing --cosmic-vcf" >&2; exit 2; }
[[ ! -f "$COSMIC_VCF"   ]] && { echo "COSMIC VCF not found: $COSMIC_VCF" >&2; exit 2; }
[[ -z "$MUTANT_TSV"     ]] && { echo "Missing --mutant-tsv" >&2; exit 2; }
[[ ! -f "$MUTANT_TSV"   ]] && { echo "Mutant TSV not found: $MUTANT_TSV" >&2; exit 2; }
[[ -z "$CLASSIFICATION" ]] && { echo "Missing --classification" >&2; exit 2; }
[[ ! -f "$CLASSIFICATION" ]] && { echo "Classification not found: $CLASSIFICATION" >&2; exit 2; }
[[ -z "$OUT_PREFIX"     ]] && OUT_PREFIX="cosmic_${TISSUE}"

if [[ "$TRAIN_AFTER" == "true" ]]; then
    [[ -z "$REFERENCE"  ]] && { echo "--train requires --reference" >&2; exit 2; }
    [[ ! -f "$REFERENCE" ]] && { echo "Reference not found: $REFERENCE" >&2; exit 2; }
    command -v "$RNEAT_BIN" >/dev/null 2>&1 || { echo "rneat binary not found (override with --rneat-bin)" >&2; exit 2; }
fi
command -v bgzip >/dev/null 2>&1 || { echo "bgzip not found on PATH (install htslib / samtools)" >&2; exit 2; }

mkdir -p "$OUT_DIR"
OUT_VCF="${OUT_DIR}/${OUT_PREFIX}.vcf.gz"

if [[ "$CHR_PREFIX" == "true" ]]; then CHR_PFX_VAL="chr"; else CHR_PFX_VAL=""; fi

# Stream the mutant TSV whether it's a bare .tsv.gz or the COSMIC .tar wrapper.
mutant_stream() {
    if [[ "$MUTANT_TSV" == *.tar ]]; then
        tar xOf "$MUTANT_TSV" --wildcards '*GenomeScreensMutant*.tsv.gz' | zcat
    else
        zcat "$MUTANT_TSV"
    fi
}

PHENO_FILE="$(mktemp)"
COSV_FILE="$(mktemp)"
TMP_BODY="$(mktemp)"
trap 'rm -f "$PHENO_FILE" "$COSV_FILE" "$TMP_BODY"' EXIT

# ── 1. phenotype IDs for the requested tissue ────────────────────────────
echo ">> [1/3] Selecting COSMIC_PHENOTYPE_IDs with PRIMARY_SITE='${TISSUE}'..."
zcat "$CLASSIFICATION" \
  | awk -F'\t' -v site="$TISSUE" 'NR>1 && $2==site {print $1}' \
  | sort -u > "$PHENO_FILE"
N_PHENO=$(wc -l < "$PHENO_FILE")
if [[ "$N_PHENO" -eq 0 ]]; then
    echo "No phenotypes with PRIMARY_SITE='${TISSUE}'. Valid sites:" >&2
    zcat "$CLASSIFICATION" | awk -F'\t' 'NR>1{print $2}' | sort -u | paste -sd' ' >&2
    exit 3
fi
echo "    ${N_PHENO} phenotype(s)."

# ── 2. COSV mutation IDs observed in those phenotypes ────────────────────
echo ">> [2/3] Scanning mutant TSV for COSV IDs in those phenotypes (one pass)..."
# col 6 = COSMIC_PHENOTYPE_ID, col 7 = GENOMIC_MUTATION_ID (COSV, = VCF ID col).
mutant_stream \
  | awk -F'\t' 'FNR==NR{ph[$1]=1; next} FNR>1 && ($6 in ph){print $7}' "$PHENO_FILE" - \
  | sort -u > "$COSV_FILE"
N_COSV=$(wc -l < "$COSV_FILE")
echo "    ${N_COSV} unique COSV mutation ID(s) for ${TISSUE}."
[[ "$N_COSV" -eq 0 ]] && { echo "No mutations for tissue '${TISSUE}'." >&2; exit 3; }

# ── 3. filter the VCF to those COSV IDs, then convert as fetch_cosmic_corpus.sh ──
echo ">> [3/3] Filtering VCF to ${TISSUE} mutations + converting..."
zcat "$COSMIC_VCF" \
  | awk -F'\t' -v PFX="$CHR_PFX_VAL" '
        FNR==NR { keep[$1]=1; next }            # COSV allowlist (first file)
        /^#/    { next }
        !($3 in keep) { next }                  # tissue filter on ID (COSV)
        {
            rl = length($4); al = length($5)
            if (rl > 1 && al > 1) next                          # drop complex
            if ($4 ~ /[^ACGTacgt]/ || $5 ~ /[^ACGTacgt]/) next  # ambiguous bases

            chrom = $1
            if (PFX == "chr") {
                if (chrom == "MT") chrom = "chrM"
                else               chrom = "chr" chrom
            }
            printf "%s\t%s\t.\t%s\t%s\t.\tPASS\t.\tGT\t0/1\n", chrom, $2, $4, $5
        }
    ' "$COSV_FILE" - \
  | LC_ALL=C sort -t$'\t' -k1,1 -k2,2n -k4,4 -k5,5 -S 2G \
  | awk -F'\t' '{ key=$1 FS $2 FS $4 FS $5; if (key != prev) { print; prev=key } }' \
  > "$TMP_BODY"

{
    cat <<HDR
##fileformat=VCFv4.2
##source=tools/fetch_cosmic_per_tissue_corpus.sh
##reference=COSMIC-GenomeScreensMutant-GRCh38
##cosmic_primary_site=${TISSUE}
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype (synthetic 0/1 — COSMIC mutations are conventionally heterozygous in tumor cell)">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NEAT_tumor_corpus
HDR
    cat "$TMP_BODY"
} | bgzip -c > "$OUT_VCF"

if command -v tabix >/dev/null 2>&1; then
    tabix -p vcf "$OUT_VCF" && echo "    Indexed: ${OUT_VCF}.tbi"
fi

N_RECORDS=$(zcat "$OUT_VCF" | grep -cv '^#' || true)
N_SNV=$(zcat "$OUT_VCF" | awk -F'\t' '!/^#/ && length($4)==1 && length($5)==1' | wc -l)
N_INDEL=$(( N_RECORDS - N_SNV ))
cat <<EOF

────────────────────────────────────────────────────────────────
COSMIC → per-tissue VCF complete  (tissue: ${TISSUE})

Output:  $OUT_VCF
         ${N_RECORDS} unique records (${N_SNV} SNV, ${N_INDEL} indel)
Chromosome convention: $([ "$CHR_PREFIX" == "true" ] && echo "chr-prefixed" || echo "bare")
────────────────────────────────────────────────────────────────
EOF

if [[ "$TRAIN_AFTER" == "true" ]]; then
    echo ""
    echo ">> Training ${TISSUE} tumor SNV/indel model..."
    MODEL_OUT="${OUT_DIR}/${OUT_PREFIX}_model.json.gz"
    TRAIN_CFG="${OUT_DIR}/${OUT_PREFIX}_train.yml"
    cat > "$TRAIN_CFG" <<EOF
reference: $REFERENCE
vcf_file: $OUT_VCF
output_file: $MODEL_OUT
overwrite_output: true
EOF
    "$RNEAT_BIN" gen-mut-model -c "$TRAIN_CFG"
    echo ""
    echo "  Trained SNV/indel model: $MODEL_OUT"
    echo "  (Graft a per-tissue sv_model on top to get a full per-tissue model.)"
fi

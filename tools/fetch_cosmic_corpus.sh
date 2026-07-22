#!/usr/bin/env bash
# Convert the COSMIC GenomeScreensMutant VCF into a deduplicated, chr-prefixed
# VCF that eidolon gen-mut-model can consume directly, then (optionally) train
# a pan-cancer tumor model.
#
# What COSMIC GenomeScreensMutant is:
#   COSMIC's "coding mutations from genome-wide screens" tier — a single
#   gzipped VCF covering all curated tumor mutations from WES/WGS studies.
#   v104 is ~50.6M records, ~1.1 GB compressed, GRCh38. Source:
#     https://cancer.sanger.ac.uk/cosmic/download/cosmic
#   Requires a free academic-license registration to download.
#
# A pre-trained model fit by this exact script against v104 is bundled at
# `tools/cosmic_v104_pancancer_model.json.gz` (9.2 KB). Users who only need
# a tumor model for `tools/cancer_simulate.sh` can skip the COSMIC download
# entirely and pass that file as `--tumor-model`. The bundled artifact's
# long-term placement is tracked at #186; expect it to move once the
# release/bundle/host decision is finalized.
#
# What this script does:
#   1. Drop complex variants (REF and ALT both >1 bp: multi-base substitutions
#      and delins). gen-mut-model classifies these as VariantType::Complex and
#      its training loop has no Complex arm — they would just log "Unknown
#      variant type, skipping". Pre-filtering saves ~150k log lines.
#   2. Drop records with ambiguous (non-ACGT) bases — defense in depth;
#      COSMIC v104 has zero such records but future releases may differ.
#   3. Add `chr` prefix to chromosomes by default. Maps `MT` → `chrM`
#      (UCSC convention used by GRCh38 references).
#   4. Sort by (chrom, pos, ref, alt) and **deduplicate** — COSMIC reports
#      the same locus once per transcript annotation, so v104 has ~11.8M
#      duplicate (chrom,pos,ref,alt) rows that would otherwise be counted
#      multiple times (gen-mut-model does no dedup of its own). Sample-level
#      recurrence lives in INFO/GENOME_SCREEN_SAMPLE_COUNT and is not used
#      by gen-mut-model anyway.
#   5. bgzip-compress the output (eidolon handles BGZF via MultiGzDecoder).
#   6. Optional: drive `eidolon gen-mut-model` against the converted VCF.
#
# Known limitations (intentional for v1):
#   - Pan-cancer model. No per-tissue split. (COSMIC encodes tissue via the
#     separate ClassificationsAndTissues tier, not the GenomeScreensMutant
#     VCF.)
#   - No recurrence weighting. Each unique (chrom,pos,ref,alt) is counted
#     once regardless of GENOME_SCREEN_SAMPLE_COUNT. v2 could expand by count.
#   - Multi-base substitutions and complex delins (~158k of 50.6M) are
#     dropped — see step 1.
#
# Usage:
#   tools/fetch_cosmic_corpus.sh \
#       --cosmic-vcf ~/code/data/cosmic_v104/Cosmic_GenomeScreensMutant_v104_GRCh38.vcf.gz
#
# Run with --help for the full flag list.

set -euo pipefail

# ── Defaults ────────────────────────────────────────────────────────────
COSMIC_VCF=""
OUT_PREFIX="cosmic_pancancer"
OUT_DIR="."
CHR_PREFIX="true"
REFERENCE=""
TRAIN_AFTER="false"
EIDOLON_BIN="eidolon"

usage() {
    cat <<'EOF'
Usage: fetch_cosmic_corpus.sh --cosmic-vcf <path> [options]

Required:
  --cosmic-vcf   Path to Cosmic_GenomeScreensMutant_<vN>_GRCh38.vcf.gz
                 (download from cancer.sanger.ac.uk/cosmic/download/cosmic
                 — requires free academic-license registration)

Output:
  --out-dir      Directory for outputs (default: current directory)
  --out-prefix   Stem for output filenames (default: cosmic_pancancer)
                 Produces <out-prefix>.vcf.gz

Conversion knobs:
  --no-chr-prefix   Emit bare chromosome names (1/2/X/Y/MT) instead of
                    chr-prefixed. Default is to add `chr` so the output
                    matches eidolon's gnomAD-SV-derived bundled defaults.
                    With chr prefix, MT is renamed chrM (UCSC convention).

Optional training step:
  --train           After conversion, run `eidolon gen-mut-model` to produce
                    a tumor mutation model JSON. Requires --reference.
  --reference       Path to the GRCh38 reference FASTA (only used when
                    --train is set; needs to use the same chr-prefix
                    convention as the converted VCF).
  --eidolon-bin       Path to the eidolon binary (default: "eidolon" on PATH)

  -h, --help        Print this message
EOF
}

# ── Arg parsing ─────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --cosmic-vcf)      COSMIC_VCF="$2"; shift 2 ;;
        --out-dir)         OUT_DIR="$2"; shift 2 ;;
        --out-prefix)      OUT_PREFIX="$2"; shift 2 ;;
        --no-chr-prefix)   CHR_PREFIX="false"; shift ;;
        --train)           TRAIN_AFTER="true"; shift ;;
        --reference)       REFERENCE="$2"; shift 2 ;;
        --eidolon-bin)       EIDOLON_BIN="$2"; shift 2 ;;
        -h|--help)         usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

# ── Validation ──────────────────────────────────────────────────────────
[[ -z "$COSMIC_VCF"  ]] && { echo "Missing --cosmic-vcf" >&2; exit 2; }
[[ ! -f "$COSMIC_VCF" ]] && { echo "COSMIC VCF not found: $COSMIC_VCF" >&2; exit 2; }

if [[ "$TRAIN_AFTER" == "true" ]]; then
    [[ -z "$REFERENCE"  ]] && { echo "--train requires --reference" >&2; exit 2; }
    [[ ! -f "$REFERENCE" ]] && { echo "Reference not found: $REFERENCE" >&2; exit 2; }
    if ! command -v "$EIDOLON_BIN" >/dev/null 2>&1; then
        echo "eidolon binary not found on PATH (override with --eidolon-bin)" >&2
        exit 2
    fi
fi

if ! command -v bgzip >/dev/null 2>&1; then
    echo "bgzip not found on PATH (install htslib / samtools)" >&2
    exit 2
fi

mkdir -p "$OUT_DIR"
OUT_VCF="${OUT_DIR}/${OUT_PREFIX}.vcf.gz"

# ── Filter + chr-prefix + sort + dedupe ─────────────────────────────────
# gen-mut-model classifies by REF/ALT length:
#   ref=1, alt=1  → SNP        (kept; trinucleotide transition recorded)
#   ref=1, alt>1  → Insertion  (kept; length-binned)
#   ref>1, alt=1  → Deletion   (kept; length-binned)
#   ref>1, alt>1  → Complex    (silently skipped by runner — drop here)

if [[ "$CHR_PREFIX" == "true" ]]; then
    CHR_PFX_VAL="chr"
else
    CHR_PFX_VAL=""
fi

TMP_BODY="$(mktemp)"
trap 'rm -f "$TMP_BODY"' EXIT

echo ">> Filtering COSMIC VCF (SNV + 1-bp-anchored indel only, ACGT only)..."

# Single pipeline: awk-filter | sort by (chrom,pos,ref,alt) | stream dedup.
# Stream dedup needs constant memory because sort already groups identicals.
zcat "$COSMIC_VCF" \
  | awk -F'\t' -v PFX="$CHR_PFX_VAL" '
        /^#/ { next }
        {
            rl = length($4); al = length($5)
            if (rl > 1 && al > 1) next                       # drop complex
            if ($4 ~ /[^ACGTacgt]/ || $5 ~ /[^ACGTacgt]/) next  # ambiguous bases

            chrom = $1
            if (PFX == "chr") {
                if (chrom == "MT") chrom = "chrM"
                else               chrom = "chr" chrom
            }
            # GT 0/1 makes gen-mut-model treat this as a het variant in a
            # synthetic single-sample column. COSMIC does not carry per-sample
            # genotype state; this is a conventional placeholder.
            printf "%s\t%s\t.\t%s\t%s\t.\tPASS\t.\tGT\t0/1\n", chrom, $2, $4, $5
        }
    ' \
  | LC_ALL=C sort -t$'\t' -k1,1 -k2,2n -k4,4 -k5,5 -S 2G \
  | awk -F'\t' '{
        key = $1 FS $2 FS $4 FS $5
        if (key != prev) { print; prev = key }
    }' \
  > "$TMP_BODY"

echo "    Assembling final VCF..."
{
    cat <<'HDR'
##fileformat=VCFv4.2
##source=tools/fetch_cosmic_corpus.sh
##reference=COSMIC-GenomeScreensMutant-GRCh38
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype (synthetic 0/1 — COSMIC mutations are conventionally heterozygous in tumor cell)">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NEAT_tumor_corpus
HDR
    cat "$TMP_BODY"
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
COSMIC → VCF conversion complete

Input:     $COSMIC_VCF
Output:    $OUT_VCF
           ${OUT_SIZE} compressed, ${N_RECORDS} unique records (post-dedup)

Chromosome convention: $([ "$CHR_PREFIX" == "true" ] && echo "chr-prefixed (chr1, ..., chrX, chrM)" || echo "bare (1, ..., X, MT)")

Excluded from this VCF (intentional v1 simplifications):
  - Complex substitutions/delins:  ~158k records — gen-mut-model skips them
  - Duplicate (chrom,pos,ref,alt):  ~11.8M rows  — multi-transcript noise
  - Recurrence weighting:           not applied   — each unique locus counted once
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
    "$EIDOLON_BIN" gen-mut-model -c "$TRAIN_CFG"
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

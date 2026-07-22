#!/usr/bin/env bash
# Build the complete per-tissue tumor models (#202): per-tissue SNP/indel spectrum
# (from COSMIC, this script) grafted with the per-tissue sv_model (from PCAWG, via
# tools/build_per_tissue_sv_models.sh — run that first).
#
# Pipeline per tissue:
#   fetch_cosmic_per_tissue_corpus.sh --train   -> per-tissue SNP/indel model
#   graft_sv_model.py                           -> graft the per-tissue sv_model on
#   => tools/cosmic_per_tissue_<label>.json.gz  (the bundled, user-facing model)
#
# The sv_model source files tools/cosmic_pancancer_sv_<label>.json.gz are produced
# by tools/build_per_tissue_sv_models.sh (a pan-cancer SNP/indel base + per-tissue
# sv_model); here they serve only as the sv_model donor.
#
# Tissue label mapping (COSMIC PRIMARY_SITE -> bundled model label):
#   breast -> BRCA,  skin -> skin,  lung -> lung
#
# Usage:
#   tools/build_per_tissue_models.sh \
#       --cosmic-dir ~/code/data/cosmic_v104 \
#       --reference  ~/code/data/hg38.fa.gz \
#       --out-dir    tools \
#       [--eidolon-bin ./target/release/eidolon]

set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

COSMIC_DIR=""
REFERENCE=""
OUT_DIR="${HERE}"
EIDOLON_BIN="eidolon"

usage() { grep '^#' "$0" | sed 's/^# \{0,1\}//'; }

while [[ $# -gt 0 ]]; do
    case "$1" in
        --cosmic-dir) COSMIC_DIR="$2"; shift 2 ;;
        --reference)  REFERENCE="$2"; shift 2 ;;
        --out-dir)    OUT_DIR="$2"; shift 2 ;;
        --eidolon-bin)  EIDOLON_BIN="$2"; shift 2 ;;
        -h|--help)    usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; exit 2 ;;
    esac
done

[[ -z "$COSMIC_DIR" ]] && { echo "Missing --cosmic-dir" >&2; exit 2; }
[[ -z "$REFERENCE"  ]] && { echo "Missing --reference" >&2; exit 2; }

VCF="${COSMIC_DIR}/Cosmic_GenomeScreensMutant_v104_GRCh38.vcf.gz"
TSV="${COSMIC_DIR}/Cosmic_GenomeScreensMutant_Tsv_v104_GRCh38.tar"
CLS="${COSMIC_DIR}/Cosmic_Classification_v104_GRCh38.tsv.gz"
for f in "$VCF" "$TSV" "$CLS" "$REFERENCE"; do
    [[ -f "$f" ]] || { echo "Required input missing: $f" >&2; exit 2; }
done

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

# COSMIC PRIMARY_SITE : bundled label
TISSUES=("breast:BRCA" "skin:skin" "lung:lung")

for pair in "${TISSUES[@]}"; do
    site="${pair%%:*}"
    label="${pair##*:}"
    sv_src="${OUT_DIR}/cosmic_pancancer_sv_${label}.json.gz"
    [[ -f "$sv_src" ]] || { echo "sv_model donor missing: $sv_src (run build_per_tissue_sv_models.sh first)" >&2; exit 3; }

    echo ""
    echo "==== ${label} (COSMIC site: ${site}) ===="
    "${HERE}/fetch_cosmic_per_tissue_corpus.sh" \
        --tissue "$site" \
        --cosmic-vcf "$VCF" --mutant-tsv "$TSV" --classification "$CLS" \
        --out-dir "$WORK" --out-prefix "cosmic_${label}" \
        --train --reference "$REFERENCE" --eidolon-bin "$EIDOLON_BIN"

    python3 "${HERE}/graft_sv_model.py" \
        --snv-model "${WORK}/cosmic_${label}_model.json.gz" \
        --sv-source "$sv_src" \
        --out "${OUT_DIR}/cosmic_per_tissue_${label}.json.gz"
done

echo ""
echo "Per-tissue models written to ${OUT_DIR}/cosmic_per_tissue_{BRCA,skin,lung}.json.gz"

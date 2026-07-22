#!/usr/bin/env bash
# Build per-tissue cancer SvModels from the PCAWG corpus (#202).
#
# For each tissue, this drives the same pipeline as the pan-cancer v1.14.0 refit
# (#218) but restricted to that tissue's donors:
#
#   build_pcawg_sv_vcf.py --projects <codes>   # tissue-filtered symbolic-SV VCF
#     -> eidolon gen-mut-model                    # fit length / CN distributions
#       -> normalize_pcawg_sv_model.py          # per-tumor rate/type-mix + splice
#
# INS length is germline (tissue-agnostic), so gnomAD is intentionally NOT
# re-streamed per tissue — normalize_pcawg_sv_model.py falls back to the
# --base-model's (gnomAD-derived) INS length, and the INS rate stays literature.
# Likewise the SNP/indel side of each model stays the pan-cancer COSMIC base;
# per-tissue SNP/indel is the separate MC3/COSMIC-classification work in #202.
#
# Tissue → PCAWG project-code mapping (donor counts among SV donors, v1.6):
#   BRCA  -> BRCA            (211)  breast; tandem-duplicator-enriched
#   skin  -> SKCM-US,MELA-AU (106)  melanoma; high UV mutation burden
#   lung  -> LUAD-US,LUSC-US ( 85)  lung adeno+squamous
# PCAWG's project codes diverge from TCGA's, and SKCM/LUAD alone are underpowered
# (~37 donors each), so skin/lung group the related project codes for stable fits.
#
# Usage:
#   tools/build_per_tissue_sv_models.sh \
#       --pcawg-dir   ~/code/data/pcawg \
#       --reference   ~/code/data/hg38.fa.gz \
#       --base-model  tools/cosmic_v104_pancancer_model.json.gz \
#       --out-dir     tools \
#       --eidolon-bin   ./target/release/eidolon

set -euo pipefail

PCAWG_DIR=""
REFERENCE=""
BASE_MODEL="tools/cosmic_v104_pancancer_model.json.gz"
OUT_DIR="tools"
EIDOLON_BIN="eidolon"
# tissue:label=project-codes (comma-separated). Override with --tissues.
TISSUES=("BRCA=BRCA" "skin=SKCM-US,MELA-AU" "lung=LUAD-US,LUSC-US")

usage() { sed -n '2,40p' "$0"; }

while [[ $# -gt 0 ]]; do
    case "$1" in
        --pcawg-dir)  PCAWG_DIR="$2"; shift 2 ;;
        --reference)  REFERENCE="$2"; shift 2 ;;
        --base-model) BASE_MODEL="$2"; shift 2 ;;
        --out-dir)    OUT_DIR="$2"; shift 2 ;;
        --eidolon-bin)  EIDOLON_BIN="$2"; shift 2 ;;
        --tissues)    IFS=' ' read -r -a TISSUES <<< "$2"; shift 2 ;;
        -h|--help)    usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

[[ -z "$PCAWG_DIR" ]] && { echo "Missing --pcawg-dir" >&2; exit 2; }
[[ -z "$REFERENCE" ]] && { echo "Missing --reference" >&2; exit 2; }
[[ -f "$REFERENCE" ]] || { echo "Reference not found: $REFERENCE" >&2; exit 2; }
[[ -f "$BASE_MODEL" ]] || { echo "Base model not found: $BASE_MODEL" >&2; exit 2; }

HERE="$(cd "$(dirname "$0")" && pwd)"
SAMPLE_SHEET="${PCAWG_DIR}/metadata/pcawg_sample_sheet.tsv"
PLOIDY="${PCAWG_DIR}/consensus_cnv/consensus.20170217.purity.ploidy.txt.gz"
[[ -f "$SAMPLE_SHEET" ]] || { echo "Sample sheet not found: $SAMPLE_SHEET (run fetch_pcawg_sv_corpus.sh)" >&2; exit 2; }

mkdir -p "$OUT_DIR"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

for entry in "${TISSUES[@]}"; do
    label="${entry%%=*}"
    projects="${entry#*=}"
    echo ""
    echo "════════════ ${label}  (projects: ${projects}) ════════════"

    corpus="${WORK}/${label}.vcf.gz"
    fit="${WORK}/${label}_fit.json.gz"
    model_out="${OUT_DIR}/cosmic_pancancer_sv_${label}.json.gz"

    echo ">> [1/3] build tissue-filtered symbolic-SV VCF"
    python3 "${HERE}/build_pcawg_sv_vcf.py" \
        --bedpe-dir "${PCAWG_DIR}/consensus_sv" \
        --cna-dir "${PCAWG_DIR}/consensus_cnv" \
        --purity-ploidy "$PLOIDY" \
        --sample-sheet "$SAMPLE_SHEET" \
        --projects "$projects" --tissue-label "$label" \
        --out "$corpus"

    echo ">> [2/3] fit sv_model with gen-mut-model"
    fit_cfg="${WORK}/${label}_fit.yml"
    cat > "$fit_cfg" <<EOF
reference: $REFERENCE
vcf_file: $corpus
output_file: $fit
overwrite_output: true
EOF
    "$EIDOLON_BIN" gen-mut-model -c "$fit_cfg"

    echo ">> [3/3] normalize + splice into ${model_out}"
    python3 "${HERE}/normalize_pcawg_sv_model.py" \
        --fitted-model "$fit" \
        --sidecar "${corpus%.vcf.gz}.counts.json" \
        --base-model "$BASE_MODEL" \
        --out "$model_out"
    echo ">> wrote ${model_out}"
done

echo ""
echo "Per-tissue SV models written to ${OUT_DIR}/cosmic_pancancer_sv_<tissue>.json.gz"
echo "Use with: tools/cancer_simulate.sh --tumor-model <model> ..."

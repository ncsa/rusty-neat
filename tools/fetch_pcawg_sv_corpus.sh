#!/usr/bin/env bash
# Download the open-tier PCAWG consensus structural-variant and copy-number
# callsets, verify checksums, and extract them so tools/build_pcawg_sv_vcf.py
# can convert them into a symbolic-SV VCF that `eidolon gen-mut-model` fits.
#
# This is the data-acquisition half of the #218 cancer-SvModel refit: it
# replaces the heuristic literature rates currently injected by
# tools/inject_cancer_sv_model.py with values counted from real PCAWG donors.
#
# What gets fetched (all open-access, no credentials):
#
#   consensus_sv/  — per-donor PASS-only breakpoint calls in BEDPE
#       final_consensus_sv_bedpe_passonly.icgc.public.tgz  (1926 donors)
#       final_consensus_sv_bedpe_passonly.tcga.public.tgz  ( 822 donors)
#     svclass ∈ {DEL, DUP, TRA, h2hINV, t2tINV}; no INS, no CNV.
#     → maps to SvType Del / Dup / Bnd(TRA) / Inv(h2hINV+t2tINV).
#
#   consensus_cnv/ — per-donor copy-number segments + per-sample purity/ploidy
#       consensus.20170119.somatic.cna.icgc.public.tar.gz (1950 donors)
#       consensus.20170119.somatic.cna.tcga.public.tar.gz ( 828 donors)
#       consensus.20170217.purity.ploidy.txt.gz
#     segment cols: chromosome start end total_cn major_cn minor_cn star
#     → CNV CN distribution + focal CNV length distribution (non-neutral
#       segments, where total_cn != round(per-sample ploidy)).
#
# NOT fetched here:
#   - INS: PCAWG retrotransposition (MEI) calls are CONTROLLED-access only.
#     The v1.14.0 refit takes INS *length* from the on-disk gnomAD-SV germline
#     callset instead (its one strongly-populated type), with INS *rate* left
#     as the somatic literature estimate. See docs/cancer_simulator.md / #218.
#   - The consensus SV/MEI VCFs (controlled) — the open BEDPE is sufficient.
#
# The mirror lives behind the ICGC-ARGO object store, NOT the default AWS S3
# endpoint (which 404s for this bucket). Hence --endpoint-url below.
# Ref: https://docs.icgc-argo.org/docs/data-access/icgc-25k-data
#
# Usage:
#   tools/fetch_pcawg_sv_corpus.sh --out-dir ~/code/data/pcawg
#
# Run with --help for the full flag list.

set -euo pipefail

# ── Defaults ────────────────────────────────────────────────────────────
OUT_DIR="."
ENDPOINT="https://object.genomeinformatics.org"
BUCKET="s3://icgc25k-open/PCAWG"
SKIP_EXISTING="true"

usage() {
    cat <<'EOF'
Usage: fetch_pcawg_sv_corpus.sh [options]

Output:
  --out-dir       Directory to download into (default: current directory).
                  Creates <out-dir>/consensus_sv/ and <out-dir>/consensus_cnv/,
                  each with the tarballs, extracted icgc/ + tcga/ trees, and
                  (for CNV) the purity/ploidy table.

Source overrides (rarely needed):
  --endpoint-url  Object-store endpoint (default: https://object.genomeinformatics.org)
  --bucket        Bucket/prefix    (default: s3://icgc25k-open/PCAWG)
  --force         Re-download even if the tarball already exists and matches md5.

  -h, --help      Print this message
EOF
}

# ── Arg parsing ─────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --out-dir)      OUT_DIR="$2"; shift 2 ;;
        --endpoint-url) ENDPOINT="$2"; shift 2 ;;
        --bucket)       BUCKET="$2"; shift 2 ;;
        --force)        SKIP_EXISTING="false"; shift ;;
        -h|--help)      usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

if ! command -v aws >/dev/null 2>&1; then
    echo "aws CLI not found on PATH (needed to reach the ICGC object store)" >&2
    exit 2
fi

# fname  md5  subdir
# (md5s from the per-folder README.md manifests on the mirror, 2024-07-18.)
MANIFEST=(
  "consensus_sv/final_consensus_sv_bedpe_passonly.icgc.public.tgz   99a20a2f1828e94275c92612930697c2 consensus_sv"
  "consensus_sv/final_consensus_sv_bedpe_passonly.tcga.public.tgz   650981d737522350503bba74e6d2c65d consensus_sv"
  "consensus_cnv/consensus.20170119.somatic.cna.icgc.public.tar.gz  b16a44e5b52302136b2dfa2d2f4b2000 consensus_cnv"
  "consensus_cnv/consensus.20170119.somatic.cna.tcga.public.tar.gz  0d05a8ec9f07aed1c70c6ebbe1db9aad consensus_cnv"
  "consensus_cnv/consensus.20170217.purity.ploidy.txt.gz            a3cf8c9d0521b23b31a57b0a35379d8e consensus_cnv"
)

aws_get() {  # remote-relative-path  local-path
    aws s3 cp "${BUCKET}/$1" "$2" \
        --endpoint-url "$ENDPOINT" --no-sign-request
}

mkdir -p "$OUT_DIR"

echo ">> Fetching PCAWG open-tier SV + CNV callsets into $OUT_DIR"
for entry in "${MANIFEST[@]}"; do
    read -r remote md5 subdir <<<"$entry"
    fname="$(basename "$remote")"
    dest_dir="${OUT_DIR}/${subdir}"
    dest="${dest_dir}/${fname}"
    mkdir -p "$dest_dir"

    if [[ "$SKIP_EXISTING" == "true" && -f "$dest" ]] \
       && echo "${md5}  ${dest}" | md5sum -c --status 2>/dev/null; then
        echo "   [have] $fname (md5 ok)"
    else
        echo "   [get ] $fname"
        aws_get "$remote" "$dest"
        echo "${md5}  ${dest}" | md5sum -c \
            || { echo "md5 mismatch for $dest" >&2; exit 3; }
    fi
done

# ── Extract the tarballs (purity/ploidy stays gzipped; the adapter zcats it) ──
echo ">> Extracting tarballs..."
extract_into() {  # tarball  target_subdir
    local tarball="$1" target="$2"
    mkdir -p "$target"
    # Each tarball nests its own dir layout; extract flat-ish into target/.
    tar xzf "$tarball" -C "$target"
}

SV_DIR="${OUT_DIR}/consensus_sv"
CNV_DIR="${OUT_DIR}/consensus_cnv"
extract_into "${SV_DIR}/final_consensus_sv_bedpe_passonly.icgc.public.tgz" "${SV_DIR}/icgc"
extract_into "${SV_DIR}/final_consensus_sv_bedpe_passonly.tcga.public.tgz" "${SV_DIR}/tcga"
extract_into "${CNV_DIR}/consensus.20170119.somatic.cna.icgc.public.tar.gz" "${CNV_DIR}/icgc"
extract_into "${CNV_DIR}/consensus.20170119.somatic.cna.tcga.public.tar.gz" "${CNV_DIR}/tcga"

# ── Sample sheet: aliquot_id → dcc_project_code, for per-tissue filtering (#202) ──
# (No published md5 in the consensus READMEs, so fetched without a checksum.)
META_DIR="${OUT_DIR}/metadata"
mkdir -p "$META_DIR"
SAMPLE_SHEET="${META_DIR}/pcawg_sample_sheet.tsv"
if [[ "$SKIP_EXISTING" == "true" && -s "$SAMPLE_SHEET" ]]; then
    echo "   [have] pcawg_sample_sheet.tsv"
else
    echo "   [get ] pcawg_sample_sheet.tsv (donor→tissue metadata)"
    aws_get "donors_and_biospecimens/pcawg_sample_sheet.tsv" "$SAMPLE_SHEET"
fi

N_BEDPE=$(find "${SV_DIR}/icgc" "${SV_DIR}/tcga" -name "*.bedpe.gz" 2>/dev/null | wc -l)
N_CNA=$(find "${CNV_DIR}/icgc" "${CNV_DIR}/tcga" -name "*.cna.txt" -o -name "*.somatic.cna.txt" 2>/dev/null | wc -l)

cat <<EOF

────────────────────────────────────────────────────────────────
PCAWG open-tier corpus ready in $OUT_DIR

  consensus_sv/  : ${N_BEDPE} per-donor BEDPE files (DEL/DUP/TRA/INV)
  consensus_cnv/ : CNA segments + purity.ploidy table
                   (consensus.20170217.purity.ploidy.txt.gz)

Next: convert to a symbolic-SV VCF gen-mut-model can fit:

  tools/build_pcawg_sv_vcf.py \\
      --bedpe-dir   ${SV_DIR} \\
      --cna-dir     ${CNV_DIR} \\
      --purity-ploidy ${CNV_DIR}/consensus.20170217.purity.ploidy.txt.gz \\
      --gnomad-ins-vcf ~/code/data/gnomad/gnomad_v4_sv.sites.vcf.gz \\
      --out pcawg_sv_corpus.vcf.gz

(INS length is sourced from gnomAD-SV; see #218 / docs/cancer_simulator.md.)
────────────────────────────────────────────────────────────────
EOF

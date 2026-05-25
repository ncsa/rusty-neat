#!/usr/bin/env bash
# Simulate a tumor / normal mixture sample by running rneat gen-reads twice
# over the same reference and concatenating the resulting FASTQs.
#
# Architecture (per docs/cancer_simulator.md):
#
#   Pass 1 (normal): coverage = (1 - purity) × total_coverage.
#     Produces a "germline" FASTQ + a golden VCF of randomly-generated
#     germline variants (the "patient's genome").
#
#   Pass 2 (tumor):  coverage = purity × total_coverage.
#     Re-uses pass 1's golden VCF as input_vcf — guaranteeing the tumor
#     cells carry the same germline variants as the normal cells. Adds
#     tumor-specific somatic mutations on top (via the tumor mutation
#     model and optional sv_rate_scale).
#
#   Merge: concatenate FASTQs (independent reads, just append). The two
#     per-pass golden VCFs survive for downstream truth-set work; merging
#     them with origin annotation is rneat issue #185, not this script.
#
# Ported from NEAT 2.1's `genReadsTumorTutorial/simulate.sh` idiom.
#
# Usage:
#   tools/cancer_simulate.sh \
#       --reference /path/to/ref.fa \
#       --output-prefix /path/to/out/sample \
#       --total-coverage 30 \
#       --purity 0.7
#
# Run with --help for the full flag list.

set -euo pipefail

# ── Defaults ────────────────────────────────────────────────────────────
REFERENCE=""
OUTPUT_PREFIX=""
TOTAL_COVERAGE="30"
PURITY="0.5"
NORMAL_MODEL=""
TUMOR_MODEL=""
GERMLINE_VCF=""
READ_LEN="151"
PAIRED_END="false"
FRAGMENT_MEAN=""
FRAGMENT_ST_DEV=""
SV_RATE_SCALE="0.0"
RNG_SEED_ROOT="cancer-simulate"
RNEAT_BIN="rneat"

usage() {
    cat <<'EOF'
Usage: cancer_simulate.sh [options]

Required:
  --reference        FASTA reference (.fa or .fa.gz)
  --output-prefix    Path stem for outputs (directory must be writable; created if absent).
                     Produces:
                       <prefix>_normal_r1.fastq.gz  (and _r2 if paired)
                       <prefix>_normal.vcf.gz       (germline truth)
                       <prefix>_tumor_r1.fastq.gz   (and _r2 if paired)
                       <prefix>_tumor.vcf.gz        (germline + somatic truth)
                       <prefix>_merged_r1.fastq.gz  (and _r2 if paired)

Coverage / purity:
  --total-coverage   Combined coverage of the merged output (default: 30)
  --purity           Tumor cell fraction in [0.0, 1.0]                (default: 0.5)
                     Tumor pass gets purity × total, normal pass gets the rest.

Mutation models (optional; both default to rneat's bundled MutationModel::default()
which is germline-derived. Supplying a cancer-trained tumor model is strongly
recommended for any non-toy run — see issue #186):
  --normal-model     Path to a .json.gz mutation model for the normal pass
  --tumor-model      Path to a .json.gz mutation model for the tumor pass

Germline VCF (optional):
  --germline-vcf     A VCF to use as the shared germline. If omitted, pass 1
                     generates one de novo and pass 2 consumes it. Supply this
                     when you want a specific germline genotype (or when running
                     multiple tumor scenarios over the same germline).

Read parameters:
  --read-len         Read length in bp (default: 151)
  --paired-ended     Generate paired-end reads. Requires --fragment-mean and
                     --fragment-st-dev. (default: single-ended)
  --fragment-mean    Mean insert length, paired-end only
  --fragment-st-dev  Insert-length stddev, paired-end only

Tumor-pass SV generation:
  --sv-rate-scale    Multiplier for the tumor model's per_base_rate (default: 0,
                     SVs from the model are disabled). Set to 1.0 to use the
                     model's nominal rate; higher values stress-test SV callers.

Reproducibility / wiring:
  --rng-seed         Seed root (default: "cancer-simulate"); per-pass seeds
                     append "-normal" and "-tumor" so the two passes don't
                     share RNG state.
  --rneat-bin        Path to the rneat binary (default: "rneat" on PATH)

  -h, --help         Print this message
EOF
}

# ── Arg parsing ─────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --reference)        REFERENCE="$2"; shift 2 ;;
        --output-prefix)    OUTPUT_PREFIX="$2"; shift 2 ;;
        --total-coverage)   TOTAL_COVERAGE="$2"; shift 2 ;;
        --purity)           PURITY="$2"; shift 2 ;;
        --normal-model)     NORMAL_MODEL="$2"; shift 2 ;;
        --tumor-model)      TUMOR_MODEL="$2"; shift 2 ;;
        --germline-vcf)     GERMLINE_VCF="$2"; shift 2 ;;
        --read-len)         READ_LEN="$2"; shift 2 ;;
        --paired-ended)     PAIRED_END="true"; shift ;;
        --fragment-mean)    FRAGMENT_MEAN="$2"; shift 2 ;;
        --fragment-st-dev)  FRAGMENT_ST_DEV="$2"; shift 2 ;;
        --sv-rate-scale)    SV_RATE_SCALE="$2"; shift 2 ;;
        --rng-seed)         RNG_SEED_ROOT="$2"; shift 2 ;;
        --rneat-bin)        RNEAT_BIN="$2"; shift 2 ;;
        -h|--help)          usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

# ── Validation ──────────────────────────────────────────────────────────
[[ -z "$REFERENCE" ]] && { echo "Missing --reference" >&2; exit 2; }
[[ ! -f "$REFERENCE" ]] && { echo "Reference file not found: $REFERENCE" >&2; exit 2; }
[[ -z "$OUTPUT_PREFIX" ]] && { echo "Missing --output-prefix" >&2; exit 2; }

# Purity in (0, 1) — endpoints would zero one of the passes; refuse rather
# than silently skip a pass (use plain rneat gen-reads for purity=0 / 1).
if ! awk -v p="$PURITY" 'BEGIN { exit !(p > 0 && p < 1) }'; then
    echo "Purity must be strictly between 0 and 1; got: $PURITY" >&2
    echo "  (For purity=0 or 1 use rneat gen-reads directly — no mixing needed.)" >&2
    exit 2
fi

if [[ "$PAIRED_END" == "true" ]]; then
    [[ -z "$FRAGMENT_MEAN"   ]] && { echo "--paired-ended requires --fragment-mean"   >&2; exit 2; }
    [[ -z "$FRAGMENT_ST_DEV" ]] && { echo "--paired-ended requires --fragment-st-dev" >&2; exit 2; }
fi

if ! command -v "$RNEAT_BIN" >/dev/null 2>&1; then
    echo "rneat binary not found on PATH (override with --rneat-bin): $RNEAT_BIN" >&2
    exit 2
fi

# ── Per-pass coverage ────────────────────────────────────────────────────
# Round to nearest integer; reject if either pass would round to zero.
NORMAL_COVERAGE=$(awk -v t="$TOTAL_COVERAGE" -v p="$PURITY" 'BEGIN { printf "%d", t * (1 - p) + 0.5 }')
TUMOR_COVERAGE=$(awk -v t="$TOTAL_COVERAGE" -v p="$PURITY" 'BEGIN { printf "%d", t * p + 0.5 }')

if [[ "$NORMAL_COVERAGE" -lt 1 || "$TUMOR_COVERAGE" -lt 1 ]]; then
    echo "Per-pass coverage rounds to zero (normal=$NORMAL_COVERAGE, tumor=$TUMOR_COVERAGE)." >&2
    echo "  Bump --total-coverage or pick a less extreme --purity." >&2
    exit 2
fi

# ── Output layout ────────────────────────────────────────────────────────
OUTPUT_DIR="$(dirname "$OUTPUT_PREFIX")"
PREFIX_NAME="$(basename "$OUTPUT_PREFIX")"
mkdir -p "$OUTPUT_DIR"

NORMAL_PREFIX="${PREFIX_NAME}_normal"
TUMOR_PREFIX="${PREFIX_NAME}_tumor"
MERGED_PREFIX="${PREFIX_NAME}_merged"

# ── Pass-1 config writer ─────────────────────────────────────────────────
write_config() {
    local config_path="$1"
    local prefix_name="$2"
    local coverage="$3"
    local model_path="$4"
    local input_vcf="$5"
    local seed_suffix="$6"
    local sv_scale="$7"

    cat > "$config_path" <<EOF
reference: $REFERENCE
read_len: $READ_LEN
coverage: $coverage
ploidy: 2
paired_ended: $PAIRED_END
EOF
    if [[ "$PAIRED_END" == "true" ]]; then
        cat >> "$config_path" <<EOF
fragment_mean: $FRAGMENT_MEAN
fragment_st_dev: $FRAGMENT_ST_DEV
EOF
    fi
    cat >> "$config_path" <<EOF
produce_fastq: true
produce_vcf: true
overwrite_output: true
output_dir: $OUTPUT_DIR
output_filename: $prefix_name
rng_seed: ${RNG_SEED_ROOT}-${seed_suffix}
sv_rate_scale: $sv_scale
EOF
    # NOTE: don't use `[[ -n "$x" ]] && echo ...` here. A false `[[ ]]` test
    # returns exit 1, and with `set -e` a function whose last command exits
    # non-zero aborts the whole script — even if the && short-circuit was the
    # intended behavior. Spell it out with `if`.
    if [[ -n "$model_path" ]]; then
        echo "mutation_model: $model_path" >> "$config_path"
    fi
    if [[ -n "$input_vcf" ]]; then
        echo "input_vcf: $input_vcf" >> "$config_path"
    fi
    return 0
}

# ── Pass 1: normal ───────────────────────────────────────────────────────
echo ">> Pass 1: normal at ${NORMAL_COVERAGE}× coverage"
NORMAL_CONFIG="${OUTPUT_DIR}/${NORMAL_PREFIX}.yml"
write_config "$NORMAL_CONFIG" "$NORMAL_PREFIX" "$NORMAL_COVERAGE" \
    "$NORMAL_MODEL" "$GERMLINE_VCF" "normal" "0.0"
"$RNEAT_BIN" gen-reads -c "$NORMAL_CONFIG"

# If the user didn't supply a germline VCF, use pass 1's golden as the
# input for pass 2 — that's the whole point of running two passes.
if [[ -z "$GERMLINE_VCF" ]]; then
    GERMLINE_VCF="${OUTPUT_DIR}/${NORMAL_PREFIX}.vcf.gz"
    if [[ ! -f "$GERMLINE_VCF" ]]; then
        echo "Expected pass-1 golden VCF at $GERMLINE_VCF but file not found." >&2
        echo "  Check that produce_vcf is honored in your gen-reads build." >&2
        exit 3
    fi
fi

# ── Pass 2: tumor ────────────────────────────────────────────────────────
echo ">> Pass 2: tumor at ${TUMOR_COVERAGE}× coverage (germline from ${GERMLINE_VCF})"
TUMOR_CONFIG="${OUTPUT_DIR}/${TUMOR_PREFIX}.yml"
write_config "$TUMOR_CONFIG" "$TUMOR_PREFIX" "$TUMOR_COVERAGE" \
    "$TUMOR_MODEL" "$GERMLINE_VCF" "tumor" "$SV_RATE_SCALE"
"$RNEAT_BIN" gen-reads -c "$TUMOR_CONFIG"

# ── Merge FASTQs ─────────────────────────────────────────────────────────
# Gzipped FASTQ files are concatenable as-is — multi-stream gzip is part
# of the spec and every standard reader (fastp, samtools, MultiGzDecoder,
# bgzf-aware tools) handles it. cat is the right tool here, not a re-gzip.
NORMAL_R1="${OUTPUT_DIR}/${NORMAL_PREFIX}_r1.fastq.gz"
TUMOR_R1="${OUTPUT_DIR}/${TUMOR_PREFIX}_r1.fastq.gz"
MERGED_R1="${OUTPUT_DIR}/${MERGED_PREFIX}_r1.fastq.gz"
cat "$NORMAL_R1" "$TUMOR_R1" > "$MERGED_R1"

if [[ "$PAIRED_END" == "true" ]]; then
    NORMAL_R2="${OUTPUT_DIR}/${NORMAL_PREFIX}_r2.fastq.gz"
    TUMOR_R2="${OUTPUT_DIR}/${TUMOR_PREFIX}_r2.fastq.gz"
    MERGED_R2="${OUTPUT_DIR}/${MERGED_PREFIX}_r2.fastq.gz"
    cat "$NORMAL_R2" "$TUMOR_R2" > "$MERGED_R2"
fi

# ── Summary ──────────────────────────────────────────────────────────────
cat <<EOF

────────────────────────────────────────────────────────────────
Cancer simulation complete

Configuration:
  Reference:        $REFERENCE
  Total coverage:   ${TOTAL_COVERAGE}×
  Purity:           $PURITY  (tumor cell fraction)
    → normal pass:  ${NORMAL_COVERAGE}× coverage
    → tumor pass:   ${TUMOR_COVERAGE}× coverage
  Paired-ended:     $PAIRED_END
  Tumor SV scale:   $SV_RATE_SCALE

Outputs (in ${OUTPUT_DIR}):
  Normal FASTQ:     ${NORMAL_PREFIX}_r1.fastq.gz $( [[ "$PAIRED_END" == "true" ]] && echo ", ${NORMAL_PREFIX}_r2.fastq.gz" )
  Normal truth:     ${NORMAL_PREFIX}.vcf.gz             (germline only)
  Tumor FASTQ:      ${TUMOR_PREFIX}_r1.fastq.gz $( [[ "$PAIRED_END" == "true" ]] && echo ", ${TUMOR_PREFIX}_r2.fastq.gz" )
  Tumor truth:      ${TUMOR_PREFIX}.vcf.gz              (germline + somatic)
  Merged FASTQ:     ${MERGED_PREFIX}_r1.fastq.gz$( [[ "$PAIRED_END" == "true" ]] && echo ", ${MERGED_PREFIX}_r2.fastq.gz" )

The merged FASTQ is what a downstream somatic-variant pipeline should
consume. Both per-pass golden VCFs are preserved for benchmarking; a
proper unified truth set with germline/somatic origin tags is tracked
at rneat issue #185.
────────────────────────────────────────────────────────────────
EOF

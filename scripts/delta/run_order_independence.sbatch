#!/bin/bash
# SLURM job: order-independence + determinism validation (report §4 item "#3 SHUFFLE").
#
# Confirms that rneat's output is reproducible and invariant to things that MUST
# NOT change a correct simulation, while being explicit about what legitimately
# does change. All runs use ONE genome (default yeast: multi-contig + fast) at a
# FIXED seed; only the variable under test changes between runs. Each run emits a
# truth VCF (reads off — variants are the determinism signal); we compare the
# header-less, fully-sorted variant bodies by md5 (set-equality) and use position
# keys to check shard disjointness.
#
# Checks (each PASS/FAIL printed in a parseable RESULTS block at the end):
#   A determinism_1thread       same config+seed, run twice, 1 thread  -> identical
#   B determinism_Nthread       same config+seed, run twice, N threads -> identical
#   C contig_order_independent  shuffle contig order in the FASTA      -> identical set
#   D shard_order_independent   same per-shard seeds, merged in two
#                               different shard orders                 -> identical set
#   E shard_disjoint            no variant emitted by >1 shard window  -> 0 duplicates
# Informational (NOT pass/fail — expected to differ by design):
#   - multithread vs single-thread variant set (RNG stream may differ)
#   - sharded union vs monolithic (shards use distinct per-shard seeds by design)
#
# Why this matters: it validates that the HPC region-sharding (genome_array.sbatch
# + merge_shards.sh) is correct and reproducible, and pins down whether monolithic
# runs are themselves contig-order-independent. Fast (~few min on yeast).
#
# Prereqs: setup.sh (rneat build + bioinf env w/ bcftools, samtools). Submit from repo root.
# Usage:
#   sbatch scripts/delta/run_order_independence.sbatch
#   REFERENCE=$SCRATCH/neat_data/yeast.fa COVERAGE=30 MT_THREADS=8 NSHARDS=4 \
#     sbatch scripts/delta/run_order_independence.sbatch

#SBATCH --job-name=rneat-orderindep
#SBATCH --partition=cpu
#SBATCH --account=bhrd-delta-cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -uo pipefail

REPO_ROOT="${RNEAT_REPO:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/../.." && pwd)}}"
source "$REPO_ROOT/scripts/delta/lib_report.sh"   # setup_conda/conda_activate + $SCRATCH

REFERENCE="${REFERENCE:-$SCRATCH/neat_data/yeast.fa}"
OUTDIR="${OUTDIR:-$SCRATCH/orderindep_$SLURM_JOB_ID}"
SEED="${SEED:-rneat order independence}"   # one fixed seed for every monolithic run
COVERAGE="${COVERAGE:-30}"
READ_LEN="${READ_LEN:-151}"
MT_THREADS="${MT_THREADS:-8}"              # threads for the multi-thread determinism check
NSHARDS="${NSHARDS:-4}"                     # target shard count for the sharding checks
CARGO_TARGET_DIR="${CARGO_TARGET_DIR:-$SCRATCH/cargo-target/rusty-neat}"
RNEAT_BIN="$CARGO_TARGET_DIR/release/rneat"

source "$HOME/.cargo/env" 2>/dev/null || true
module load samtools/1.22-cce19.0.0 2>/dev/null || true
setup_conda
conda_activate bioinf   # bcftools, samtools, seqkit

[[ -x "$RNEAT_BIN" ]] || { echo "rneat binary not found: $RNEAT_BIN (run setup.sh)" >&2; exit 1; }
[[ -f "$REFERENCE" ]] || { echo "reference not found: $REFERENCE" >&2; exit 1; }
mkdir -p "$OUTDIR"
[[ -f "${REFERENCE}.fai" ]] || samtools faidx "$REFERENCE"

echo "=== order-independence / determinism check: $SLURM_JOB_ID ==="
echo "ref=$REFERENCE  seed='$SEED'  cov=${COVERAGE}x  mt_threads=$MT_THREADS  nshards=$NSHARDS"
echo "out=$OUTDIR"

# ── one rneat run -> $1/sample.vcf.gz ────────────────────────────────────
# run_rneat <outdir> <reference> <seed> <threads> [target_bed]
run_rneat() {
    local wd="$1" ref="$2" seed="$3" threads="$4" bed="${5:-}"
    mkdir -p "$wd"
    {
        echo "reference: $ref"
        echo "read_len: $READ_LEN"
        echo "coverage: $COVERAGE"
        echo "ploidy: 2"
        echo "paired_ended: true"
        echo "fragment_mean: 350"
        echo "fragment_st_dev: 50"
        echo "produce_fastq: false"   # variants only — the determinism signal
        echo "produce_vcf: true"
        echo "produce_bam: false"
        echo "overwrite_output: true"
        echo "output_dir: $wd"
        echo "output_filename: sample"
        echo "rng_seed: $seed"
        echo "sv_rate_scale: 0"
        [[ -n "$bed" ]] && echo "target_bed: $bed"
        echo "num_threads: $threads"
    } > "$wd/sim.yml"
    "$RNEAT_BIN" gen-reads -c "$wd/sim.yml" > "$wd/run.log" 2>&1
    [[ -f "$wd/sample.vcf.gz" ]] || { echo "  no VCF from $wd (see run.log)" >&2; return 1; }
}

# header-less, fully-sorted body -> stable md5 (set-equality across runs)
md5of() { bcftools view -H "$1" 2>/dev/null | LC_ALL=C sort | md5sum | cut -d' ' -f1; }
nvar()  { bcftools view -H "$1" 2>/dev/null | wc -l; }

# ── A/B: determinism, 1 thread and N threads ─────────────────────────────
echo; echo "[A] determinism @ 1 thread (run twice)..."
run_rneat "$OUTDIR/a1" "$REFERENCE" "$SEED" 1
run_rneat "$OUTDIR/a2" "$REFERENCE" "$SEED" 1
echo "[B] determinism @ $MT_THREADS threads (run twice)..."
run_rneat "$OUTDIR/b1" "$REFERENCE" "$SEED" "$MT_THREADS"
run_rneat "$OUTDIR/b2" "$REFERENCE" "$SEED" "$MT_THREADS"

# ── C: contig-order independence ─────────────────────────────────────────
echo "[C] contig-order independence (reverse contig order in FASTA)..."
SHUF_REF="$OUTDIR/ref_shuffled.fa"
# reverse the contig order, then re-emit the FASTA in that order
mapfile -t CONTIGS < <(cut -f1 "${REFERENCE}.fai")
REV=(); for (( i=${#CONTIGS[@]}-1; i>=0; i-- )); do REV+=("${CONTIGS[$i]}"); done
samtools faidx "$REFERENCE" "${REV[@]}" > "$SHUF_REF"
samtools faidx "$SHUF_REF"
run_rneat "$OUTDIR/c1" "$SHUF_REF" "$SEED" 1

# ── D/E: shard-order independence + disjointness ─────────────────────────
echo "[D/E] sharding (~$NSHARDS windows, SAME per-shard seeds, two merge orders)..."
GENOME_BP=$(awk '{s+=$2} END{print s}' "${REFERENCE}.fai")
SHARD_BP=$(( (GENOME_BP + NSHARDS - 1) / NSHARDS ))
SHARD_DIR="$OUTDIR/shards"
N=$(bash "$REPO_ROOT/scripts/delta/make_shard_beds.sh" "${REFERENCE}.fai" "$SHARD_BP" "$SHARD_DIR")
echo "  made $N shard bed(s)"
SHVCF=()
for (( id=0; id<N; id++ )); do
    bed="$(printf '%s/shard_%04d.bed' "$SHARD_DIR" "$id")"
    sd="$OUTDIR/shard_$(printf '%04d' "$id")"
    run_rneat "$sd" "$REFERENCE" "$SEED shard $id" 1 "$bed"   # fixed per-shard seed
    v="$sd/sample.vcf.gz"; bcftools index -f -t "$v" 2>/dev/null || true
    SHVCF+=("$v")
done
# merge in natural order and in reversed order; both must give the same set
REVVCF=(); for (( i=${#SHVCF[@]}-1; i>=0; i-- )); do REVVCF+=("${SHVCF[$i]}"); done
bcftools concat -a -Ou "${SHVCF[@]}"  2>/dev/null | bcftools sort -Oz -o "$OUTDIR/merge_fwd.vcf.gz" 2>/dev/null
bcftools concat -a -Ou "${REVVCF[@]}" 2>/dev/null | bcftools sort -Oz -o "$OUTDIR/merge_rev.vcf.gz" 2>/dev/null
# disjointness: a CHROM:POS:REF:ALT key emitted by >1 shard window = double-count
DUP=$(for v in "${SHVCF[@]}"; do bcftools view -H "$v" | awk '{print $1":"$2":"$4":"$5}'; done \
        | LC_ALL=C sort | uniq -d | wc -l)

# ── verdicts ─────────────────────────────────────────────────────────────
A1=$(md5of "$OUTDIR/a1/sample.vcf.gz"); A2=$(md5of "$OUTDIR/a2/sample.vcf.gz")
B1=$(md5of "$OUTDIR/b1/sample.vcf.gz"); B2=$(md5of "$OUTDIR/b2/sample.vcf.gz")
C1=$(md5of "$OUTDIR/c1/sample.vcf.gz")
MF=$(md5of "$OUTDIR/merge_fwd.vcf.gz"); MR=$(md5of "$OUTDIR/merge_rev.vcf.gz")
verdict() { [[ "$2" == "$3" ]] && echo "PASS" || echo "FAIL"; }

echo
echo "════════════════════════════════════════════════════════════════"
echo "ORDER-INDEPENDENCE RESULTS — $SLURM_JOB_ID ($(basename "$REFERENCE"))"
echo "════════════════════════════════════════════════════════════════"
printf '  %-28s %s\n' "determinism_1thread:"      "$(verdict x "$A1" "$A2")  ($A1 vs $A2)"
printf '  %-28s %s\n' "determinism_${MT_THREADS}thread:" "$(verdict x "$B1" "$B2")  ($B1 vs $B2)"
printf '  %-28s %s\n' "contig_order_independent:"  "$(verdict x "$A1" "$C1")  (base $A1 vs shuffled $C1)"
printf '  %-28s %s\n' "shard_order_independent:"   "$(verdict x "$MF" "$MR")  (fwd $MF vs rev $MR)"
printf '  %-28s %s\n' "shard_disjoint:"            "$([[ "$DUP" -eq 0 ]] && echo PASS || echo FAIL)  ($DUP duplicate keys)"
echo "  ── informational (differences here are expected, not failures) ──"
printf '  %-28s %s\n' "multithread_vs_1thread:"    "$([[ "$A1" == "$B1" ]] && echo same || echo differ)"
printf '  %-28s %s\n' "variant_counts:"            "base=$(nvar "$OUTDIR/a1/sample.vcf.gz") shuffled=$(nvar "$OUTDIR/c1/sample.vcf.gz") shard_union=$(nvar "$OUTDIR/merge_fwd.vcf.gz")"
echo "════════════════════════════════════════════════════════════════"

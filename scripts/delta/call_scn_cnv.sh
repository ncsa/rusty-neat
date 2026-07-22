#!/bin/bash
# SLURM job: call COPY-NUMBER VARIANTS on the staged SCN data with `delly cnv`, to populate the
# biologically-central CNV signal that `delly sr` (call_scn_sv.sh) does NOT produce — its output
# is DEL/DUP/INV/BND, so the built sv_model's `cnv_copy_number_distribution` comes out empty.
# CNV / gene copy-number is at the heart of SCN virulence (host Rhg1/Rhg4 resolve via copy-number;
# the SCN genome has extensive parasitism-gene duplications — see call_scn_sv.sh refs), so a
# read-depth CNV pass is the highest-value biology add for the SV model.
#
# ⚠️ `delly cnv` REQUIRES a mappability map (`-m`). There's no prebuilt map for a non-model genome
# like SCN, so we generate one with `dicey` (dellytools), following Delly's documented procedure:
# chop the genome into reads → realign → measure per-base mappability. This is the heavy/fiddly
# part and dicey's CLI is VERSION-SENSITIVE (we learned the hard way that Delly 2.x renamed
# `call`→`sr`) — if the chop/mappability2 calls below fail, check `dicey --help` for your version,
# or pass a prebuilt map with MAP=<path> to skip generation entirely.
#
# Prereqs: strain staged by stage_scn.sh (REF + BAM); delly (static v2.3.1, DELLY=<path>) + dicey
# (conda install -n bioinf -c bioconda dicey) + bwa-mem2 + bcftools; samtools/htslib modules.
#
# Usage:
#   DELLY=$SCRATCH/bin/delly sbatch scripts/delta/call_scn_cnv.sh                 # MM26/MM26A
#   DELLY=… ACC=GCA_040805705.1 SRR=SRR27329600 sbatch scripts/delta/call_scn_cnv.sh   # PA3/BD3
#   DELLY=… MAP=/path/to/prebuilt.map.fa.gz sbatch scripts/delta/call_scn_cnv.sh   # skip dicey

#SBATCH --job-name=eidolon-scncnv
#SBATCH --partition=cpu
#SBATCH --account=bhrd-delta-cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -euo pipefail

REPO_ROOT="${EIDOLON_REPO:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/../.." && pwd)}}"
source "$REPO_ROOT/scripts/delta/lib_report.sh"

D="${DATA_DIR:-$SCRATCH/neat_data/scn}"
ACC="${ACC:-GCA_040805935.1}"
SRR="${SRR:-SRR27329602}"
REF="${REFERENCE:-$D/${ACC}.fa}"
BAM="${BAM:-$D/${SRR}.bam}"
MAP="${MAP:-$D/${ACC}.map.fa.gz}"           # mappability map; auto-generated with dicey if absent
DELLY="${DELLY:-delly}"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

module load samtools/1.22-cce19.0.0
module load htslib/1.22-gcc13.3.1
setup_conda
conda_activate bioinf

"$DELLY" --version >/dev/null 2>&1 || { echo "delly '$DELLY' not runnable — pass DELLY=<static binary>." >&2; exit 1; }
"$DELLY" --help 2>&1 | grep -qE '^[[:space:]]*cnv[[:space:]]' || { echo "this delly has no 'cnv' subcommand." >&2; exit 1; }
[[ -s "$REF" ]] || { echo "reference not staged: $REF (run stage_scn.sh)" >&2; exit 1; }
[[ -s "$BAM" ]] || { echo "BAM not staged: $BAM (run stage_scn.sh)" >&2; exit 1; }
[[ -f "$REF.fai" ]] || samtools faidx "$REF"
[[ -f "$BAM.bai" ]] || samtools index "$BAM"

echo "=== SCN CNV calling: ref=$ACC  bam=$(basename "$BAM") ==="

# ── 1. mappability map (dicey) — VERSION-SENSITIVE; skip if MAP already exists ────
if [[ ! -s "$MAP" ]]; then
    command -v dicey >/dev/null 2>&1 || {
        echo "dicey not found (conda install -n bioinf -c bioconda dicey), or pass MAP=<prebuilt>." >&2
        exit 1; }
    echo "generating mappability map with dicey (Delly's documented procedure)..."
    MW="$D/cnv_map_${ACC}"; mkdir -p "$MW"
    ( cd "$MW"
      # chop the genome into simulated reads, realign to itself, measure mappability.
      dicey chop "$REF"
      [[ -s "$REF.bwt.2bit.64" ]] || bwa-mem2 index "$REF"
      bwa-mem2 mem -t "$THREADS" "$REF" read1.fq.gz read2.fq.gz 2>bwa.log \
          | samtools sort -@ "$THREADS" -o srt.bam -
      samtools index srt.bam
      dicey mappability2 srt.bam            # → map.fa.gz (+ any .fai/.gzi)
      # dicey's CLI is version-sensitive; fail loudly if the expected output is absent
      # rather than letting the next step operate on a missing/renamed file.
      [[ -s map.fa.gz ]] || { echo "ERROR: dicey mappability2 did not produce map.fa.gz (CLI version mismatch?)" >&2; exit 1; }
      # normalize to a bgzipped, faidx'd map at $MAP
      gunzip -f map.fa.gz
      bgzip -f map.fa
      mv -f map.fa.gz "$MAP"
      samtools faidx "$MAP"
    )
    rm -rf "$MW"
fi
[[ -s "$MAP" ]] || { echo "mappability map missing after generation: $MAP" >&2; exit 1; }

# ── 2. delly cnv (read-depth copy-number) ────────────────────────────────────────
CNV_BCF="$D/${SRR}.cnv.bcf"
CNV_VCF="$D/${SRR}.cnv.vcf.gz"
if [[ ! -s "$CNV_VCF" ]]; then
    echo "delly cnv..."
    OMP_NUM_THREADS="$THREADS" "$DELLY" cnv -g "$REF" -m "$MAP" -o "$CNV_BCF" "$BAM"
    bcftools view -f PASS -O z -o "$CNV_VCF" "$CNV_BCF"
    bcftools index -t "$CNV_VCF"
fi

ncnv=$(bcftools view -H "$CNV_VCF" | wc -l)
echo "── copy-number spread (PASS; CN from FORMAT/RDCN where available) ──"
bcftools view -H "$CNV_VCF" | grep -oE 'SVTYPE=[A-Z]+' | sort | uniq -c | sort -rn || true
echo "  total PASS CNV segments: $ncnv"

echo
echo "════════════════════════════════════════════════════════════════"
echo "SCN CNVs: $ncnv PASS segments → $CNV_VCF"
echo "To fold CNV into the model, concat it into the all-variants VCF (sample-normalized, like"
echo "call_scn_sv.sh does) and rebuild — the sv_model's cnv_copy_number_distribution then populates:"
echo "  # reheader $CNV_VCF + $D/${SRR}.all.vcf.gz to a common sample, bcftools concat -a | sort,"
echo "  # then: REFERENCE=$REF INPUT_BAM=$BAM INPUT_FASTQ=$D/${SRR}.fastq.gz INPUT_VCF=<merged> \\"
echo "  #       sbatch scripts/delta/model_builders.sbatch"
echo "════════════════════════════════════════════════════════════════"

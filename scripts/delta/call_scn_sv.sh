#!/bin/bash
# SLURM job: call STRUCTURAL VARIANTS on the staged SCN data with Delly, so gen-mut-model's
# SV model (CNV / BND / INV / DEL / DUP) can be exercised on a real complex organism — the
# short-variant staging (stage_scn.sh, bcftools) only yields SNPs + small indels, leaving the
# SV model unpopulated (sv_model: None).
#
# WHY SVs are worth chasing in SCN (Heterodera glycines): the SCN system is structurally rich.
# CNV / tandem duplication is central to it — e.g. the host Rhg1/Rhg4 resistance loci resolve
# largely through copy-number (Patil et al. 2019, Plant Biotech J, doi:10.1111/pbi.13086;
# Grunwald et al. 2021, Plant Genome, doi:10.1002/tpg2.20152), and effector genes differ
# between virulent vs avirulent SCN inbred populations like ours (Kwon et al. 2024,
# Phytopathology, doi:10.1094/PHYTO-03-24-0095-R). The SCN genome itself is repeat/TE-rich with
# extensive parasitism-gene duplications (Masonbrink et al. 2019, BMC Genomics). So DELs/DUPs
# and repeat/TE-driven rearrangements are expected.
#
# ⚠️ EXPECTATIONS / CAVEATS:
#   - Short-read SV calling on a repeat/TE-rich genome is noisy — expect many calls in
#     repetitive regions; keep only Delly PASS and eyeball the SVTYPE spread.
#   - TE *insertions* specifically are hard from short reads (they live in repeats) — DEL/DUP/
#     INV are more reliably recovered than novel INS. Long reads would do better; we work with
#     what's on NCBI.
#   - Pool-seq: the sample is 150 pooled individuals, so SV genotypes/AF are pooled pseudo-
#     values (as with the short VCF). Fine for populating + exercising the SV model.
#
# Pipeline: delly short-read SV call (sr on 2.x, call on 1.x) on the staged BAM -> keep PASS -> report SVTYPE counts -> merge with the
# short VCF (sample-name normalized) into an all-variants VCF -> print the gen-mut-model command
# so the built model carries SVs. Delly emits symbolic ALTs (<DEL>/<DUP>/<INV>/<INS>/BND) that
# gen-mut-model's SvType::from_alt parses directly.
#
# Prereqs: the strain already staged by stage_scn.sh (BAM + REF present); delly + bcftools
# (bioinf conda env), samtools/htslib modules.
#
# Usage (defaults = MM26 strain, matching stage_scn.sh defaults):
#   sbatch scripts/delta/call_scn_sv.sh
#   ACC=GCA_040805705.1 SRR=SRR27329600 sbatch scripts/delta/call_scn_sv.sh   # PA3 / MM-BD3A

#SBATCH --job-name=rneat-scnsv
#SBATCH --partition=cpu
#SBATCH --account=bhrd-delta-cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=8:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -euo pipefail

REPO_ROOT="${RNEAT_REPO:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/../.." && pwd)}}"
source "$REPO_ROOT/scripts/delta/lib_report.sh"

D="${DATA_DIR:-$SCRATCH/neat_data/scn}"
ACC="${ACC:-GCA_040805935.1}"            # MM26 (default); PA3 = GCA_040805705.1
SRR="${SRR:-SRR27329602}"                # MM26A (default); MM-BD3A = SRR27329600
REF="${REFERENCE:-$D/${ACC}.fa}"
BAM="${BAM:-$D/${SRR}.bam}"
SHORT_VCF="${SHORT_VCF:-$D/${SRR}.vcf.gz}"   # the bcftools SNP/indel VCF from stage_scn.sh

module load samtools/1.22-cce19.0.0
module load htslib/1.22-gcc13.3.1
setup_conda
conda_activate bioinf

# Delly: default to `delly` on PATH, override with DELLY=<path>. bioconda's delly often can't
# load its Boost at runtime (libboost_iostreams.so.1.85.0 — and Delta's `module load boost` is
# 1.88, the WRONG soname, so it won't satisfy it). Robust fix = Delly's statically-linked
# release binary (no shared-lib deps):
#   wget -O $SCRATCH/bin/delly https://github.com/dellytools/delly/releases/download/v2.3.1/delly-v2.3.1-linux-amd64
#   chmod +x $SCRATCH/bin/delly     # then run this job with  DELLY=$SCRATCH/bin/delly
DELLY="${DELLY:-delly}"
# Actually RUN it (not just `command -v`) so a binary that resolves but can't load its libs
# fails here with a clear message instead of mid-call.
"$DELLY" --version >/dev/null 2>&1 || {
    echo "delly '$DELLY' not runnable (missing on PATH, or a shared-lib load error like the Boost one)." >&2
    echo "Use the static binary and pass DELLY=<path> — see the comment above this line." >&2
    exit 1; }
[[ -s "$REF" ]] || { echo "reference not staged: $REF (run stage_scn.sh first)" >&2; exit 1; }
[[ -s "$BAM" ]] || { echo "BAM not staged: $BAM (run stage_scn.sh first)" >&2; exit 1; }
[[ -f "$REF.fai" ]] || samtools faidx "$REF"
[[ -f "$BAM.bai" ]] || samtools index "$BAM"

echo "=== SCN SV calling: ref=$ACC  bam=$(basename "$BAM") ==="

# ── 1. Delly call → PASS-only SV VCF ─────────────────────────────────────────────
BCF="$D/${SRR}.delly.bcf"
SV_VCF="$D/${SRR}.sv.vcf.gz"
if [[ ! -s "$SV_VCF" ]]; then
    # Delly 2.x renamed short-read SV calling from `call` to `sr` (and added `lr` for
    # long reads — useful when we get the PacBio/ONT data). Auto-detect so this works on
    # both delly 1.x (`call`) and 2.x (`sr`); same -g/-o/BAM interface either way.
    if "$DELLY" --help 2>&1 | grep -qE '^[[:space:]]*sr[[:space:]]'; then
        SV_CMD=sr
    else
        SV_CMD=call
    fi
    echo "delly $SV_CMD (single sample; DEL/DUP/INV/INS/BND)..."
    # OMP threads via env; delly parallelizes across SV types.
    OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}" "$DELLY" "$SV_CMD" -g "$REF" -o "$BCF" "$BAM"
    bcftools view -f PASS -O z -o "$SV_VCF" "$BCF"
    bcftools index -t "$SV_VCF"
fi

echo "── SVTYPE spread (PASS) ──"
bcftools view -H "$SV_VCF" | grep -oE 'SVTYPE=[A-Z]+' | sort | uniq -c | sort -rn || true
nsv=$(bcftools view -H "$SV_VCF" | wc -l)
echo "  total PASS SVs: $nsv"

# ── 2. merge short (SNP/indel) + SV into one all-variants VCF ─────────────────────
# gen-mut-model takes ONE VCF; concat needs identical sample names, so reheader both to a
# common sample (real file, NOT <(...) process-sub — that's a known Delta srun/sbatch trap).
ALL_VCF="$D/${SRR}.all.vcf.gz"
if [[ -s "$SHORT_VCF" && "$nsv" -gt 0 ]]; then
    printf 'scn\n' > "$D/${SRR}.sample.txt"
    bcftools reheader -s "$D/${SRR}.sample.txt" "$SHORT_VCF" -o "$D/${SRR}.short.rh.vcf.gz"
    bcftools reheader -s "$D/${SRR}.sample.txt" "$SV_VCF"    -o "$D/${SRR}.sv.rh.vcf.gz"
    bcftools index -t "$D/${SRR}.short.rh.vcf.gz"
    bcftools index -t "$D/${SRR}.sv.rh.vcf.gz"
    bcftools concat -a "$D/${SRR}.short.rh.vcf.gz" "$D/${SRR}.sv.rh.vcf.gz" \
        | bcftools sort -O z -o "$ALL_VCF"
    bcftools index -t "$ALL_VCF"
    rm -f "$D/${SRR}.short.rh.vcf.gz"* "$D/${SRR}.sv.rh.vcf.gz"* "$D/${SRR}.sample.txt"
    nall=$(bcftools view -H "$ALL_VCF" | wc -l)
else
    ALL_VCF="$SV_VCF"; nall="$nsv"
fi

echo
echo "════════════════════════════════════════════════════════════════"
echo "SCN SVs called: $nsv PASS SVs; all-variants VCF: $nall records → $ALL_VCF"
echo "Rebuild the models with the SV-inclusive VCF (harness handles the gen-mut-model config):"
echo "  REFERENCE=$REF INPUT_BAM=$D/${SRR}.bam \\"
echo "  INPUT_FASTQ=$D/${SRR}.fastq.gz INPUT_VCF=$ALL_VCF \\"
echo "    sbatch scripts/delta/model_builders.sbatch"
echo "Then confirm the SV model populated (expect True):"
echo "  zcat <modelbuild_JOB>/models/mut_model.json.gz | \\"
echo "    python3 -c 'import json,sys;print(\"sv_model populated:\", json.load(sys.stdin)[\"sv_model\"] is not None)'"
echo "════════════════════════════════════════════════════════════════"

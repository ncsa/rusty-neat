#!/usr/bin/env bash
# One-time environment setup for rusty-neat on Delta (NCSA ACCESS).
# Run this interactively before submitting any SLURM jobs.
#
# What this does:
#   1. Installs rustup + a stable toolchain into $HOME/.cargo (persists across jobs)
#   2. Creates a conda env with NEAT 4 (the Python predecessor) for benchmarking
#   3. Creates a conda env (bioinf) with tools not available as Delta modules:
#        bcftools, bwa-mem2, gatk4
#      Note: samtools and htslib (tabix/bgzip) ARE available as Delta modules
#        (samtools/1.22-cce19.0.0 and htslib/1.22-gcc13.3.1) and are loaded
#        by the SLURM jobs directly. bwa-mem2 and gatk are not available.
#   4. Pulls the hap.py Apptainer image used by cancer_pipeline.sbatch
#
# Usage (from the repo root):
#   bash scripts/delta/setup.sh
#
# Re-running is safe — each step skips if already done.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
# Resolve Delta paths ($SCRATCH isn't exported on Delta). setup.sh runs
# interactively, so $0 is reliable here (unlike the spooled sbatch jobs).
source "$(dirname "$0")/lib_report.sh"
SIF_DIR="${SIF_DIR:-$SCRATCH/sif}"          # where to cache Apptainer .sif files
CARGO_TARGET_DIR="${CARGO_TARGET_DIR:-$SCRATCH/cargo-target/rusty-neat}"
CONDA_ENV_NAME="neat4"
BIOINF_ENV_NAME="bioinf"

echo "=== rusty-neat Delta setup ==="
echo "Repo:       $REPO_ROOT"
echo "SIF dir:    $SIF_DIR"
echo "Cargo target: $CARGO_TARGET_DIR"

# ── 1. Rust toolchain + rneat binary ────────────────────────────────────
# Check for a module first; fall back to rustup if absent.
if module load rust 2>/dev/null; then
    echo "[1/4] Rust module loaded: $(rustc --version)"
elif command -v cargo >/dev/null 2>&1; then
    echo "[1/4] Rust already installed via rustup: $(rustc --version)"
else
    echo "[1/4] Installing rustup (no rust module found)..."
    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y --default-toolchain stable
    source "$HOME/.cargo/env"
    echo "      $(rustc --version)"
fi
[[ -f "$HOME/.cargo/env" ]] && source "$HOME/.cargo/env"

# Redirect build artifacts to $SCRATCH — a full release build generates
# 10-15 GB of crate artifacts that would fill a typical HPC home quota.
export CARGO_TARGET_DIR
mkdir -p "$CARGO_TARGET_DIR"
echo "[1/4] Building rneat release binary (artifacts → \$SCRATCH)..."
cd "$REPO_ROOT"
cargo build --release 2>&1 | tail -5
echo "      Binary: $CARGO_TARGET_DIR/release/rneat"

# ── 2. NEAT 4 conda env ─────────────────────────────────────────────────
setup_conda   # load the conda module + bootstrap `conda` (Delta: miniforge3-python)

# NEAT 4 + bcftools. NEAT shells out to `bcftools sort` at runtime (its
# environment.yml pins bcftools==1.22*), but the bioconda `neat` package does
# NOT pull it — without bcftools NEAT dies with FileNotFoundError: 'bcftools'.
# Used by benchmark.sbatch (throughput) and germline_e2e.sbatch (NEAT fidelity
# arm). If the env already exists, RETROFIT bcftools (a pre-existing neat4 from
# before this was added wouldn't otherwise get it).
if conda env list | grep -q "^${CONDA_ENV_NAME} "; then
    echo "[2/4] Conda env '$CONDA_ENV_NAME' exists — ensuring bcftools is present..."
    conda run -n "$CONDA_ENV_NAME" bcftools --version >/dev/null 2>&1 \
        || conda install -y -n "$CONDA_ENV_NAME" -c conda-forge -c bioconda bcftools
else
    echo "[2/4] Creating conda env for NEAT 4 (bioconda)..."
    conda create -y -n "$CONDA_ENV_NAME" -c conda-forge -c bioconda neat bcftools
fi
# `neat` has no --version flag; verify the tools resolve instead.
echo "      neat4: $(conda run -n "$CONDA_ENV_NAME" which neat 2>/dev/null || echo 'neat MISSING')"
echo "             $(conda run -n "$CONDA_ENV_NAME" bcftools --version 2>/dev/null | head -1 || echo 'bcftools MISSING')"

# ── 3. bioinf conda env (bwa-mem2, gatk4, bcftools) ────────────────────
# samtools and htslib (tabix/bgzip) are available as Delta modules.
# bcftools, bwa-mem2, and gatk4 are not — install via bioconda.
if conda env list | grep -q "^${BIOINF_ENV_NAME} "; then
    echo "[3/4] Conda env '$BIOINF_ENV_NAME' already exists — skipping."
else
    echo "[3/4] Creating bioinf conda env (bcftools, bwa-mem2, gatk4)..."
    conda create -y -n "$BIOINF_ENV_NAME" \
        -c bioconda -c conda-forge \
        bcftools bwa-mem2 gatk4
    echo "      Tools installed:"
    # `|| true`: head -1 closes the pipe after one line, so the tool gets
    # SIGPIPE (exit 141) and `set -o pipefail` would otherwise abort setup.
    conda run -n "$BIOINF_ENV_NAME" bcftools  --version 2>&1 | head -1 || true
    conda run -n "$BIOINF_ENV_NAME" bwa-mem2  version   2>&1 | head -1 || true
    conda run -n "$BIOINF_ENV_NAME" gatk      --version 2>&1 | head -1 || true
fi

# ── 4. hap.py Apptainer image ───────────────────────────────────────────
mkdir -p "$SIF_DIR"
HAPPY_SIF="$SIF_DIR/happy_v0.3.12.sif"
if [[ -f "$HAPPY_SIF" ]]; then
    echo "[4/4] hap.py SIF already present: $HAPPY_SIF"
else
    echo "[4/4] Pulling hap.py image (this takes a few minutes)..."
    apptainer pull "$HAPPY_SIF" docker://jmcdani20/hap.py:v0.3.12
    echo "      Saved to: $HAPPY_SIF"
fi

cat <<EOF

Setup complete.

Next steps:
  1. Place reference genome(s) and input data in \$SCRATCH
  2. ACCESS account is preset to bhrd-delta-cpu in the *.sbatch files.
     Override per-submit if needed:  sbatch --account=<other> scripts/delta/<job>.sbatch
  3. (optional) Point results at durable storage for the ACCESS final report —
     defaults to \${WORK:-\$HOME}/rneat-access-results (NOT scratch, which is purged):
       export RESULTS_DIR=/projects/<account>/rneat-access-results
  4. Submit jobs:
       sbatch scripts/delta/baseline_capture.sbatch   # simulator baseline
       sbatch scripts/delta/germline_e2e.sbatch       # SNP/indel recall
       sbatch scripts/delta/benchmark.sbatch          # NEAT-vs-rneat throughput
       sbatch scripts/delta/cancer_pipeline.sbatch    # somatic SNV/indel recall
  5. Build the ACCESS final-report document once runs finish:
       bash scripts/delta/collect_report.sh           # -> \$RESULTS_DIR/REPORT.md

Each job persists its result artifacts + a provenance/resource manifest under
\$RESULTS_DIR (off scratch); collect_report.sh aggregates them (incl. total
core-hours ≈ SUs) into REPORT.md.
EOF
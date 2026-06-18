#!/usr/bin/env bash
# One-time environment setup for rusty-neat on Delta (NCSA ACCESS).
# Run this interactively before submitting any SLURM jobs.
#
# What this does:
#   1. Installs rustup + a stable toolchain into $HOME/.cargo (persists across jobs)
#   2. Creates a conda env with NEAT 4 (the Python predecessor) for benchmarking
#   3. Pulls the hap.py Apptainer image used by cancer_pipeline.sbatch
#
# Usage (from the repo root):
#   bash scripts/delta/setup.sh
#
# Re-running is safe — each step skips if already done.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
SIF_DIR="${SIF_DIR:-$SCRATCH/sif}"          # where to cache Apptainer .sif files
CONDA_ENV_NAME="neat4"

echo "=== rusty-neat Delta setup ==="
echo "Repo:    $REPO_ROOT"
echo "SIF dir: $SIF_DIR"

# ── 1. Rust toolchain ───────────────────────────────────────────────────
if command -v cargo >/dev/null 2>&1; then
    echo "[1/3] Rust already installed: $(rustc --version)"
else
    echo "[1/3] Installing rustup..."
    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y --default-toolchain stable
    source "$HOME/.cargo/env"
    echo "      $(rustc --version)"
fi

echo "[1/3] Building rneat release binary..."
cd "$REPO_ROOT"
source "$HOME/.cargo/env"
cargo build --release 2>&1 | tail -5
echo "      Binary: $REPO_ROOT/target/release/rneat"

# ── 2. NEAT 4 conda env ─────────────────────────────────────────────────
module load anaconda3_cpu 2>/dev/null || module load miniforge 2>/dev/null || true

if conda env list | grep -q "^${CONDA_ENV_NAME} "; then
    echo "[2/3] Conda env '$CONDA_ENV_NAME' already exists — skipping."
else
    echo "[2/3] Creating conda env for NEAT 4..."
    conda create -y -n "$CONDA_ENV_NAME" python=3.11
    conda run -n "$CONDA_ENV_NAME" pip install neat-genreads
    echo "      NEAT 4 installed: $(conda run -n $CONDA_ENV_NAME neat --version 2>&1 || echo 'check with: conda run -n $CONDA_ENV_NAME neat --version')"
fi

# ── 3. hap.py Apptainer image ───────────────────────────────────────────
mkdir -p "$SIF_DIR"
HAPPY_SIF="$SIF_DIR/happy_v0.3.12.sif"
if [[ -f "$HAPPY_SIF" ]]; then
    echo "[3/3] hap.py SIF already present: $HAPPY_SIF"
else
    echo "[3/3] Pulling hap.py image (this takes a few minutes)..."
    apptainer pull "$HAPPY_SIF" docker://jmcdani20/hap.py:v0.3.12
    echo "      Saved to: $HAPPY_SIF"
fi

cat <<EOF

Setup complete.

Next steps:
  1. Place reference genome(s) and input data in \$SCRATCH
  2. Edit scripts/delta/benchmark.sbatch  — set ACCOUNT and DATA_DIR
  3. Edit scripts/delta/cancer_pipeline.sbatch — set ACCOUNT, REFERENCE, etc.
  4. Submit:
       sbatch scripts/delta/benchmark.sbatch
       sbatch scripts/delta/cancer_pipeline.sbatch

EOF
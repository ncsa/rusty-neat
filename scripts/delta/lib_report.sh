#!/usr/bin/env bash
# Shared reporting/persistence helpers for the Delta SLURM jobs.
#
# Source this from each *.sbatch right after REPO_ROOT is set:
#     source "$REPO_ROOT/scripts/delta/lib_report.sh"
# (NOT via $(dirname "$0") — under sbatch $0 is a spooled copy). Sourcing also
# resolves the Delta filesystem paths below, so do it before any $SCRATCH use.
# Call archive_run() once near the end of the job.
#
# WHY THIS EXISTS: the jobs do their heavy work in $SCRATCH, which Delta PURGES
# (files untouched ~30 days are deleted). For the ACCESS final report we need
# outcomes preserved on DURABLE storage, consolidated, with per-run provenance
# and resource (core-hour / SU) accounting. archive_run() copies the small
# report artifacts off scratch and writes a run_manifest.tsv; collect_report.sh
# later aggregates all runs into a single REPORT.md.

# ── Delta filesystem paths ───────────────────────────────────────────────
# Delta does NOT export $SCRATCH/$WORK to jobs, and these scripts run under
# `set -u`, so resolve them here (this file is sourced before any $SCRATCH use).
# Layout (confirmed): scratch = /scratch/<project>/$USER (purged ~30 days);
# durable project space = /projects/<project>/$USER.
# Override by exporting SCRATCH / RESULTS_DIR, or by setting ACCESS_PROJECT.
: "${USER:=$(id -un)}"
: "${ACCESS_PROJECT:=bhrd}"
: "${SCRATCH:=/scratch/$ACCESS_PROJECT/$USER}"
if [[ -n "${SLURM_JOB_ID:-}" && ! -d "$SCRATCH" ]]; then
    echo "ERROR: scratch dir '$SCRATCH' not found — export SCRATCH=... or set ACCESS_PROJECT=..." >&2
    exit 1
fi

# Durable results root for the ACCESS final report — the project filesystem
# (persists), NEVER scratch. Override with RESULTS_DIR=...
RESULTS_DIR="${RESULTS_DIR:-/projects/$ACCESS_PROJECT/$USER/rneat-access-results}"

# ── conda activation (Delta) ─────────────────────────────────────────────
# On Delta the miniforge module puts conda's legacy `activate` on PATH but NOT
# `conda` itself, so `source activate` bootstraps the base env + shell
# functions. Override the module with CONDA_MODULE=...; no-op if conda is
# already available. The `set +u` is required: conda's activate scripts read
# unbound vars (e.g. $PS1) and would abort under our `set -u`.
setup_conda() {
    command -v conda >/dev/null 2>&1 && return 0
    module load "${CONDA_MODULE:-miniforge3-python}" 2>/dev/null || true
    set +u
    command -v conda >/dev/null 2>&1 || source activate 2>/dev/null || true
    if ! command -v conda >/dev/null 2>&1 && [[ -n "${CONDA_PREFIX:-}" ]]; then
        source "$CONDA_PREFIX/etc/profile.d/conda.sh" 2>/dev/null || true
    fi
    set -u
    command -v conda >/dev/null 2>&1 || {
        echo "ERROR: conda unavailable after 'module load ${CONDA_MODULE:-miniforge3-python}' + 'source activate'." >&2
        echo "       Set CONDA_MODULE=<your module> or initialize conda before running." >&2
        return 1
    }
}

# set-u-safe environment activation (same unbound-var reason as above).
conda_activate() {
    set +u
    conda activate "$1" || source activate "$1"
    local rc=$?
    set -u
    return "$rc"
}

# Emit five space-separated resource values for the current SLURM job:
#   elapsed_s alloc_cpus alloc_nodes maxrss_kb core_hours
# Best-effort: prints zeros when not under SLURM or sacct is unavailable, so the
# caller can always `read` exactly five fields.
_resource_values() {
    local jid="${SLURM_JOB_ID:-}"
    if [[ -z "$jid" ]] || ! command -v sacct >/dev/null 2>&1; then
        echo "0 0 0 0 0.0"; return
    fi
    # First (allocation) line; ElapsedRaw is seconds, MaxRSS like 1234K/M/G.
    local line elapsed cpus nodes maxrss
    line=$(sacct -j "$jid" --noheader --parsable2 \
        --format=ElapsedRaw,AllocCPUS,AllocNodes,MaxRSS 2>/dev/null | head -1)
    elapsed=$(echo "$line" | cut -d'|' -f1); cpus=$(echo "$line" | cut -d'|' -f2)
    nodes=$(echo "$line" | cut -d'|' -f3); maxrss=$(echo "$line" | cut -d'|' -f4)
    elapsed=${elapsed:-0}; cpus=${cpus:-0}; nodes=${nodes:-0}
    local rss_kb=0
    case "$maxrss" in
        *K) rss_kb=${maxrss%K} ;;
        *M) rss_kb=$(awk -v v="${maxrss%M}" 'BEGIN{printf "%.0f", v*1024}') ;;
        *G) rss_kb=$(awk -v v="${maxrss%G}" 'BEGIN{printf "%.0f", v*1048576}') ;;
    esac
    # Core-hours = elapsed_hours * allocated CPUs (≈ Delta CPU SUs).
    local core_hours
    core_hours=$(awk -v e="$elapsed" -v c="$cpus" 'BEGIN{printf "%.3f", (e/3600.0)*c}')
    echo "${elapsed:-0} ${cpus:-0} ${nodes:-0} ${rss_kb:-0} ${core_hours:-0.0}"
}

# archive_run <kind> <source_outdir> [artifact-file ...]
# Copies the named (small) report artifacts off scratch into the durable
# results dir and writes run_manifest.tsv with provenance + resource usage.
# Large data (FASTQ/BAM) is intentionally left in scratch — only report-
# relevant files (csv/tsv/json/logs) are persisted.
archive_run() {
    local kind="$1" outdir="$2"; shift 2
    local jid="${SLURM_JOB_ID:-local}"
    local dest="$RESULTS_DIR/$kind/job_${jid}"
    mkdir -p "$dest"

    local f
    for f in "$@"; do
        [[ -e "$outdir/$f" ]] && cp -f "$outdir/$f" "$dest/" 2>/dev/null || true
    done
    # SLURM stdout/err (written to the submit dir as <jobname>_<jobid>.out/.err).
    local base="${SLURM_JOB_NAME:-job}_${jid}"
    [[ -f "${base}.out" ]] && cp -f "${base}.out" "$dest/" 2>/dev/null || true
    [[ -f "${base}.err" ]] && cp -f "${base}.err" "$dest/" 2>/dev/null || true

    local ver git_desc
    ver="$("${RNEAT_BIN:-rneat}" --version 2>/dev/null || echo unknown)"
    git_desc="$(git -C "${REPO_ROOT:-.}" describe --tags --always --dirty 2>/dev/null || echo unknown)"
    local elapsed_s alloc_cpus alloc_nodes maxrss_kb core_hours
    read -r elapsed_s alloc_cpus alloc_nodes maxrss_kb core_hours < <(_resource_values)

    # Tab-separated key/value manifest — robust to parse with awk, no jq needed.
    printf '%s\t%s\n' \
        kind          "$kind" \
        date_utc      "$(date -u +%Y-%m-%dT%H:%M:%SZ)" \
        slurm_job_id  "$jid" \
        rneat_version "$ver" \
        git           "$git_desc" \
        reference     "$(basename "${REFERENCE:-NA}")" \
        elapsed_s     "${elapsed_s:-0}" \
        alloc_cpus    "${alloc_cpus:-0}" \
        alloc_nodes   "${alloc_nodes:-0}" \
        maxrss_kb     "${maxrss_kb:-0}" \
        core_hours    "${core_hours:-0.0}" \
        artifacts     "$*" \
        > "$dest/run_manifest.tsv"

    echo "[archive] $kind -> $dest  (core_hours=${core_hours:-0.0}, files: $*)"
}

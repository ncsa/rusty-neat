#!/usr/bin/env bash
# Re-freeze the PERFORMANCE rows of baseline_metrics.tsv from a fresh benchmark
# run on the CURRENT build — without disturbing the fidelity/determinism rows.
#
# WHY THIS EXISTS
#   The v1.19.1 logging fix (#340) is the motivating case. Every Phase-1/2 Delta
#   run used the old default log level `trace`, whose per-base debug events fire
#   ~coverage*read_len*ref_bp times and burned most of the wall-clock on log I/O.
#   So the frozen `*_wall_s` / `*_rss_mb` baselines were captured UNDER that
#   handicap and are now stale (too slow). The fix does NOT change output bytes,
#   so the fidelity/determinism/recall baselines stay valid — only the perf rows
#   must be re-measured on the fixed build and re-frozen. This script does exactly
#   that: it rewrites only the `*_wall_s` and `*_rss_mb` rows and leaves every
#   other row untouched.
#
# USAGE
#   1) Build the fixed binary on Delta and run the perf harness (n=3 median):
#        bash scripts/delta/setup.sh
#        MODE=submit TIER=0 LABEL=logfix bash scripts/delta/regression_suite.sh
#        # (TIER=1 if you also want chr22 perf re-frozen; TIER=0 = ecoli only)
#   2) When it finishes, collect the candidate metrics:
#        MODE=collect MANIFEST=$HOME/regression_logfix.manifest \
#          bash scripts/delta/regression_suite.sh > $HOME/logfix.candidate.tsv
#   3) Splice the fresh perf rows into the baseline (writes a .bak first):
#        bash scripts/delta/rebaseline_perf.sh $HOME/logfix.candidate.tsv
#   4) Review the printed diff, then commit baseline_metrics.tsv with a note that
#      the perf rows were re-frozen on the fixed build (git SHA / eidolon --version).
#
# SAFETY CHECK (recommended, separate): confirm the logging fix left OUTPUT
#   byte-identical. Run baseline_capture.sbatch on the fixed build and compare its
#   baseline.json `fastq_md5` / `vcf_records_md5` to a pre-fix run's. This script
#   can diff two such files for you:
#        bash scripts/delta/rebaseline_perf.sh --check-identity OLD.json NEW.json
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
BASELINE="${BASELINE:-$REPO_ROOT/scripts/delta/baseline_metrics.tsv}"

# ── --check-identity mode: assert two baseline.json files agree on output md5 ──
if [[ "${1:-}" == "--check-identity" ]]; then
    old="${2:?usage: --check-identity OLD.json NEW.json}"
    new="${3:?usage: --check-identity OLD.json NEW.json}"
    get() { grep -oE "\"$1\"[[:space:]]*:[[:space:]]*\"[0-9a-f]+\"" "$2" | grep -oE '[0-9a-f]{32}' | head -1; }
    rc=0
    for key in fastq_md5 vcf_records_md5; do
        o=$(get "$key" "$old"); n=$(get "$key" "$new")
        if [[ -n "$o" && "$o" == "$n" ]]; then
            echo "  $key: IDENTICAL ($o)"
        else
            echo "  $key: MISMATCH  old=$o  new=$n" >&2; rc=1
        fi
    done
    if [[ "$rc" -eq 0 ]]; then
        echo "OUTPUT-IDENTITY: PASS — the change did not alter simulator output."
    else
        echo "OUTPUT-IDENTITY: FAIL — the change altered output; this is NOT a pure perf/logging change." >&2
    fi
    exit "$rc"
fi

# ── perf-row splice mode ───────────────────────────────────────────────────────
CAND="${1:?usage: rebaseline_perf.sh CANDIDATE.tsv   (from regression_suite.sh MODE=collect)}"
[[ -f "$BASELINE" ]] || { echo "baseline not found: $BASELINE" >&2; exit 1; }
[[ -f "$CAND" ]]     || { echo "candidate not found: $CAND" >&2; exit 1; }

# A perf metric is any row whose name ends in _wall_s or _rss_mb.
is_perf() { [[ "$1" == *_wall_s || "$1" == *_rss_mb ]]; }

# Pull the fresh value for a perf metric out of the candidate TSV (col 2).
cand_val() { awk -F'\t' -v m="$1" '$1==m {print $2; exit}' "$CAND"; }

TMP="$(mktemp)"
CHANGED=0; MISSING=""
while IFS= read -r line || [[ -n "$line" ]]; do
    # pass through comments / blanks unchanged
    if [[ "$line" == \#* || -z "${line// }" ]]; then echo "$line" >> "$TMP"; continue; fi
    IFS=$'\t' read -r metric value sd gate tier <<< "$line"
    if is_perf "$metric"; then
        new_val="$(cand_val "$metric")"
        if [[ -n "$new_val" && "$new_val" != "NA" ]]; then
            printf '%s\t%s\t%s\t%s\t%s\n' "$metric" "$new_val" "$sd" "$gate" "$tier" >> "$TMP"
            [[ "$new_val" != "$value" ]] && { echo "  $metric: $value -> $new_val"; CHANGED=$((CHANGED+1)); }
        else
            echo "$line" >> "$TMP"                 # no fresh value → keep old, warn
            MISSING="$MISSING $metric"
        fi
    else
        echo "$line" >> "$TMP"                     # non-perf row: never touch
    fi
done < "$BASELINE"

echo
if [[ -n "$MISSING" ]]; then
    echo "WARNING: no fresh candidate value for:$MISSING (kept old value)." >&2
    echo "         (run a TIER that produces these — TIER=1 adds chr22_*.)" >&2
fi
if [[ "$CHANGED" -eq 0 ]]; then
    echo "No perf rows changed — baseline left as-is. (Nothing written.)"
    rm -f "$TMP"; exit 0
fi
cp -p "$BASELINE" "$BASELINE.bak"
mv "$TMP" "$BASELINE"
echo
echo "Re-froze $CHANGED perf row(s) in $BASELINE  (backup: $BASELINE.bak)"
echo "Fidelity/determinism/recall rows were NOT touched."
echo "Now: git diff scripts/delta/baseline_metrics.tsv ; commit with the build's git SHA."

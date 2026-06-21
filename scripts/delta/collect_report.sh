#!/usr/bin/env bash
# Aggregate all archived Delta runs into a single report-ready REPORT.md for the
# ACCESS final report. Reads the run_manifest.tsv files (and the key result
# artifacts beside them) that archive_run() persisted under RESULTS_DIR.
#
# Usage:
#   bash scripts/delta/collect_report.sh                 # writes $RESULTS_DIR/REPORT.md
#   RESULTS_DIR=/path bash scripts/delta/collect_report.sh [output.md]
#
# No external deps beyond coreutils/awk — safe to run on a login node.

set -euo pipefail

RESULTS_DIR="${RESULTS_DIR:-${WORK:-$HOME}/rneat-access-results}"
OUT="${1:-$RESULTS_DIR/REPORT.md}"

[[ -d "$RESULTS_DIR" ]] || { echo "No results dir: $RESULTS_DIR (run some jobs first)" >&2; exit 1; }
mapfile -t MANIFESTS < <(find "$RESULTS_DIR" -name run_manifest.tsv | sort)
[[ "${#MANIFESTS[@]}" -gt 0 ]] || { echo "No run_manifest.tsv under $RESULTS_DIR" >&2; exit 1; }

# Read one value from a TSV manifest by key.
mval() { awk -F'\t' -v k="$2" '$1==k{print $2; exit}' "$1"; }

{
    echo "# rneat on Delta (NCSA ACCESS) — results"
    echo
    echo "_Generated $(date -u +%Y-%m-%dT%H:%M:%SZ) from \`$RESULTS_DIR\` — ${#MANIFESTS[@]} run(s)._"
    echo
    echo "## Run summary"
    echo
    echo "| Kind | Date (UTC) | rneat | git | Reference | Core-hours | Max RSS (GB) | Job |"
    echo "|------|------------|-------|-----|-----------|-----------:|-------------:|-----|"
    total_ch="0.0"
    for m in "${MANIFESTS[@]}"; do
        kind=$(mval "$m" kind);        date=$(mval "$m" date_utc)
        ver=$(mval "$m" rneat_version); gitd=$(mval "$m" git)
        ref=$(mval "$m" reference);     ch=$(mval "$m" core_hours)
        rss=$(mval "$m" maxrss_kb);     jid=$(mval "$m" slurm_job_id)
        # sacct often fails to record MaxRSS on Delta (empty/0); show "—" rather
        # than a misleading "0.00 GB". Per-tool memory is in the benchmark's
        # median.tsv (peak_rss_mb, from /usr/bin/time).
        if [[ "${rss:-0}" =~ ^0*$ || -z "${rss:-}" ]]; then
            rss_gb="—"
        else
            rss_gb=$(awk -v r="$rss" 'BEGIN{printf "%.2f", r/1048576}')
        fi
        echo "| $kind | $date | $ver | $gitd | $ref | ${ch:-0} | $rss_gb | $jid |"
        total_ch=$(awk -v t="$total_ch" -v c="${ch:-0}" 'BEGIN{printf "%.3f", t+c}')
    done
    echo
    echo "**Total core-hours (≈ Delta CPU SUs): $total_ch**"
    echo

    # Inline the key result artifact for each run, grouped by kind.
    echo "## Outcomes by run"
    for m in "${MANIFESTS[@]}"; do
        d=$(dirname "$m"); kind=$(mval "$m" kind); jid=$(mval "$m" slurm_job_id)
        echo
        echo "### $kind — job $jid"
        local_emitted=0
        for art in baseline.json median.tsv tune.tsv \
                   rneat_scored.summary.csv neat_scored.summary.csv \
                   scored.summary.csv scored.stats.csv het_af_histogram.tsv; do
            if [[ -f "$d/$art" ]]; then
                echo
                echo "<details><summary>\`$art\`</summary>"
                echo
                echo '```'
                # Cap very long files so the report stays readable.
                head -200 "$d/$art"
                [[ $(wc -l < "$d/$art") -gt 200 ]] && echo "... (truncated)"
                echo '```'
                echo
                echo "</details>"
                local_emitted=1
            fi
        done
        [[ "$local_emitted" -eq 0 ]] && echo "_(no inlined artifacts; see $d)_"
    done

    echo
    echo "---"
    echo "_Max RSS via sstat/sacct (\"—\" = not recorded by Slurm accounting); the"
    echo "benchmark's median.tsv carries per-tool peak RSS from /usr/bin/time._"
    echo "_Artifacts and SLURM logs for each run live under \`$RESULTS_DIR/<kind>/job_<id>/\`._"
} > "$OUT"

echo "Wrote $OUT (${#MANIFESTS[@]} runs, total core-hours above)."

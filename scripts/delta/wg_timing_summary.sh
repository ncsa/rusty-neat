#!/usr/bin/env bash
# Summarize per-window wall-clock for a sharded whole-genome run — the
# "simulation-only timing" number (no downstream align/call). Reads the
# /usr/bin/time -v `time.txt` each shard writes under SHARD_OUTROOT.
#
#   bash scripts/delta/wg_timing_summary.sh SHARD_OUTROOT [EXPECT_SHARDS]
#
# Prints: completeness, per-window min/median/mean/max, total compute (sum of
# window wall-clocks = the single-core work), and the heaviest few windows with
# their region (the tail that sets the per-task walltime).
set -uo pipefail

OUTROOT="${1:?usage: wg_timing_summary.sh SHARD_OUTROOT [EXPECT_SHARDS]}"
EXPECT="${2:-}"
[[ -d "$OUTROOT" ]] || { echo "not a directory: $OUTROOT" >&2; exit 1; }

mapfile -t DIRS < <(find "$OUTROOT" -maxdepth 1 -type d -name 'shard_*' | sort)
n_dirs=${#DIRS[@]}
complete=0
TMP="$(mktemp)"                         # rows: seconds<TAB>shard<TAB>region
for d in "${DIRS[@]}"; do
    sh=$(basename "$d")
    [[ -f "$d/s_r1.fastq.gz" && -f "$d/s.vcf.gz" ]] && complete=$(( complete + 1 ))
    t="$d/time.txt"; [[ -f "$t" ]] || continue
    el=$(awk -F': ' '/Elapsed \(wall clock\)/{print $NF}' "$t")
    [[ -n "$el" ]] || continue
    region=$(cat "$d"/*.bed 2>/dev/null | head -1 | tr '\t' ':')
    awk -v e="$el" -v s="$sh" -v r="$region" '
        BEGIN{ n=split(e,a,":"); sec=(n==3)?a[1]*3600+a[2]*60+a[3]:(n==2)?a[1]*60+a[2]:a[1];
               printf "%.1f\t%s\t%s\n", sec, s, r }' >> "$TMP"
done

echo "════════════════════════════════════════════════════════════════"
echo "Whole-genome simulation timing — $OUTROOT"
echo "  shards found: $n_dirs   complete (FASTQ+VCF): $complete${EXPECT:+ / $EXPECT expected}"
if [[ -n "$EXPECT" && "$complete" -ne "$EXPECT" ]]; then
    echo "  WARNING: $(( EXPECT - complete )) shard(s) incomplete — rerun the array before merging." >&2
fi

if [[ -s "$TMP" ]]; then
    sort -n "$TMP" | awk -F'\t' '
        { v[NR]=$1; sum+=$1; sh[NR]=$2; rg[NR]=$3 }
        END{
            n=NR; if(n==0) exit
            med=(n%2)?v[(n+1)/2]:(v[n/2]+v[n/2+1])/2
            printf "  windows timed : %d\n", n
            printf "  per-window sec: min=%.1f  median=%.1f  mean=%.1f  max=%.1f\n", v[1], med, sum/n, v[n]
            printf "  total compute : %.1f core-sec  (%.2f core-hours, single-core sum)\n", sum, sum/3600
            print  "  heaviest windows (set the per-task walltime):"
            for(i=n; i>n-5 && i>=1; i--) printf "     %8.1fs  %s  %s\n", v[i], sh[i], rg[i]
        }'
else
    echo "  (no time.txt files found — was the array run with /usr/bin/time?)"
fi
echo "════════════════════════════════════════════════════════════════"
rm -f "$TMP"

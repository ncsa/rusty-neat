#!/usr/bin/env bash
# Turn a run_genome_reps.sh manifest into the wall-clock table for report 3.6:
# per (K, rep) the TRUE array span (first task Start -> last task End) and the
# per-shard rneat time distribution, then mean +/- sd of wall-clock per K. Pass
# a baseline shared-partition run to print the realistic-vs-optimal contrast.
#
# Usage:
#   bash scripts/delta/aggregate_reps.sh $SCRATCH/wg_sweep_hg38.manifest
#   # include the realistic (shared-partition) baseline run for contrast:
#   BASELINE_JOB=19539584 BASELINE_OUTROOT=$SCRATCH/wg_hg38b \
#     bash scripts/delta/aggregate_reps.sh $SCRATCH/wg_sweep_hg38.manifest
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
source "$REPO_ROOT/scripts/delta/lib_report.sh"

MANIFEST="${1:?usage: aggregate_reps.sh <manifest>}"
[[ -f "$MANIFEST" ]] || { echo "manifest not found: $MANIFEST" >&2; exit 1; }

# wall-clock (seconds) of an array job = last finite End - first finite Start.
job_wall() {
    local jid="$1" start end
    start=$(sacct -j "$jid" -X -n -P -o Start | grep -vE 'Unknown|None' | sort | head -1)
    end=$(sacct  -j "$jid" -X -n -P -o End   | grep -vE 'Unknown|None' | sort | tail -1)
    [[ -n "$start" && -n "$end" ]] || { echo ""; return; }
    echo $(( $(date -d "$end" +%s) - $(date -d "$start" +%s) ))
}

# per-shard rneat time distribution (from an OUTROOT's time.txt files).
shard_dist() {
    local outroot="$1"
    for f in "$outroot"/shard_*/time.txt; do
        [[ -f "$f" ]] || continue
        local e; e=$(awk -F'): ' '/Elapsed/{print $2}' "$f")
        [[ -n "$e" ]] && awk -F: -v x="$e" 'BEGIN{n=split(x,a,":"); print (n==3)?a[1]*3600+a[2]*60+a[3]:a[1]*60+a[2]}'
    done | sort -n | awk '{v[NR]=$1; s+=$1} END{ if(NR==0){printf "n/a"; next}
        printf "n=%d median=%.1fm p90=%.1fm max=%.1fm", NR, v[int(NR/2)]/60, v[int(NR*0.9)]/60, v[NR]/60 }'
}

fmt() { awk -v s="$1" 'BEGIN{ if(s==""){print "?"} else printf "%dm%02ds", int(s/60), s%60 }'; }

printf '%-4s %-4s %-10s %-6s %-10s  %s\n' K rep job nodes wall "per-shard rneat time"
printf '%s\n' "--------------------------------------------------------------------------------"

declare -A WALLS
while read -r K rep jid nodes outroot; do
    [[ "$K" == \#* || -z "${K:-}" ]] && continue
    w=$(job_wall "$jid")
    WALLS[$K]+="$w "
    printf '%-4s %-4s %-10s %-6s %-10s  %s\n' "$K" "$rep" "$jid" "$nodes" "$(fmt "$w")" "$(shard_dist "$outroot")"
done < "$MANIFEST"

echo
echo "Wall-clock per K (the optimal-config sweep, exclusive nodes):"
for K in $(echo "${!WALLS[@]}" | tr ' ' '\n' | sort -n); do
    echo "${WALLS[$K]}" | awk -v k="$K" '{ n=0; mn=1e18;
        for(i=1;i<=NF;i++){ if($i!=""){x[++n]=$i; s+=$i; if($i<mn)mn=$i} }
        if(n==0){ printf "  K=%-3s no data\n", k; next }
        m=s/n; for(i=1;i<=n;i++) v+=(x[i]-m)^2; sd=(n>1)?sqrt(v/(n-1)):0;
        printf "  K=%-3s n=%d  mean=%dm%02ds  sd=%dm%02ds  best=%dm%02ds\n",
            k, n, int(m/60), int(m)%60, int(sd/60), int(sd)%60, int(mn/60), int(mn)%60 }'
done

if [[ -n "${BASELINE_JOB:-}" ]]; then
    echo
    echo "Realistic shared-partition baseline (for contrast):"
    bw=$(job_wall "$BASELINE_JOB")
    printf '  job=%s  wall=%s  %s\n' "$BASELINE_JOB" "$(fmt "$bw")" \
        "$([[ -n "${BASELINE_OUTROOT:-}" ]] && shard_dist "$BASELINE_OUTROOT")"
fi

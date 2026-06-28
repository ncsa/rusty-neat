#!/usr/bin/env bash
# Aggregate input-variety germline_e2e results into a TSV for the fidelity chart:
# one row per genome with size + SNP/indel recall/precision + query Ts/Tv (PASS).
# Usage: bash scripts/delta/collect_variety.sh <manifest>   (from run_input_variety.sh)
set -euo pipefail

MANIFEST="${1:?usage: collect_variety.sh <manifest>}"
[[ -f "$MANIFEST" ]] || { echo "manifest not found: $MANIFEST" >&2; exit 1; }

printf 'genome\tsize_Mb\tSNP_recall\tSNP_prec\tINDEL_recall\tINDEL_prec\tTsTv\n'
while read -r name jid outdir ref _; do
    [[ "$name" == \#* || -z "${name:-}" ]] && continue
    f="$outdir/rneat_scored.summary.csv"
    [[ -f "$f" ]] || { printf '%s\tNA\tNA\tNA\tNA\tNA\tNA\n' "$name"; continue; }
    size=$(awk '{s+=$2} END{printf "%.1f", s/1e6}' "${ref}.fai" 2>/dev/null || echo NA)
    python3 - "$f" "$name" "$size" <<'PY'
import csv, sys
f, name, size = sys.argv[1], sys.argv[2], sys.argv[3]
rows = {r["Type"]: r for r in csv.DictReader(open(f)) if r.get("Filter") == "PASS"}
def num(t, k):
    try: return f"{float(rows[t][k]):.4f}"
    except Exception: return "NA"
tstv = rows.get("SNP", {}).get("QUERY.TOTAL.TiTv_ratio", "")
try: tstv = f"{float(tstv):.3f}"
except Exception: tstv = "NA"
print("\t".join([name, size, num("SNP","METRIC.Recall"), num("SNP","METRIC.Precision"),
                 num("INDEL","METRIC.Recall"), num("INDEL","METRIC.Precision"), tstv]))
PY
done < "$MANIFEST"

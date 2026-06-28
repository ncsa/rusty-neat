#!/usr/bin/env bash
# Tabulate a cancer-sweep manifest (from run_cancer_sweeps.sh) into a TSV for charting.
#
#   purity/coverage/tmr : per point, som.py somatic SNV + indel recall/precision
#                         (Mutect2, $outdir/scored.stats.csv) + derived normal cov.
#   tissue              : per tissue, per-type SV recall from truvari summaries
#                         ($outdir/truvari_<caller>_<TYPE>/summary.json), Manta + Delly.
#
# Best-effort + defensive: a missing/locked output yields NA rather than aborting.
# Usage: bash scripts/delta/collect_cancer_sweeps.sh <axis manifest>
set -euo pipefail

MANIFEST="${1:?usage: collect_cancer_sweeps.sh <manifest>}"
[[ -f "$MANIFEST" ]] || { echo "manifest not found: $MANIFEST" >&2; exit 1; }
axis=$(awk '!/^#/ && NF {print $1; exit}' "$MANIFEST")

if [[ "$axis" == "tissue" ]]; then
    printf 'tissue\tcaller\tDEL_recall\tDUP_recall\tINV_recall\tBND_recall\n'
    while read -r ax point jid outdir total purity tmr model; do
        [[ "$ax" == \#* || -z "${ax:-}" ]] && continue
        for caller in manta delly; do
            python3 - "$outdir" "$point" "$caller" <<'PY'
import json, sys, os
outdir, tissue, caller = sys.argv[1:4]
def rec(t):
    p = os.path.join(outdir, f"truvari_{caller}_{t}", "summary.json")
    try:
        d = json.load(open(p))
        v = d.get("recall")
        return f"{float(v):.3f}" if v is not None else "NA"
    except Exception:
        return "NA"
print("\t".join([tissue, caller, rec("DEL"), rec("DUP"), rec("INV"), rec("BND")]))
PY
        done
    done < "$MANIFEST"
else
    printf 'point\ttotal_cov\tpurity\tnormal_cov\tSNV_recall\tSNV_prec\tINDEL_recall\tINDEL_prec\n'
    while read -r ax point jid outdir total purity tmr model; do
        [[ "$ax" == \#* || -z "${ax:-}" ]] && continue
        f="$outdir/scored.stats.csv"
        if [[ ! -f "$f" ]]; then
            printf '%s\t%s\t%s\tNA\tNA\tNA\tNA\tNA\n' "$point" "$total" "$purity"; continue
        fi
        python3 - "$f" "$point" "$total" "$purity" <<'PY'
import csv, sys
f, point, total, purity = sys.argv[1:5]
try: nc = f"{float(total)*(1-float(purity)):.0f}"
except Exception: nc = "NA"
rows = {}
for r in csv.DictReader(open(f)):
    rows[(r.get("type") or "").strip().lower()] = r
def cell(type_aliases, metric):
    for t in type_aliases:
        r = rows.get(t)
        if not r: continue
        for c in r:
            if c.strip().lower() == metric:
                try: return f"{float(r[c]):.4f}"
                except Exception: return "NA"
    return "NA"
print("\t".join([point, total, purity, nc,
                 cell(("snvs","snp","snps"), "recall"),
                 cell(("snvs","snp","snps"), "precision"),
                 cell(("indels","indel"),    "recall"),
                 cell(("indels","indel"),    "precision")]))
PY
    done < "$MANIFEST"
fi

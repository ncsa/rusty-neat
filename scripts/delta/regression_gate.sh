#!/usr/bin/env bash
# Compare candidate metrics against a frozen baseline and apply regression gates.
# Exits non-zero if any gate FAILs (CI / --dependency friendly).
#
# Both files are TSV: metric <tab> value [<tab> sd <tab> gate <tab> tier].
# The BASELINE carries sd/gate/tier; the CANDIDATE needs only metric <tab> value.
# Gates (from baseline): exact | ge2sd (>= value-2*sd) | band:LO:HI | lepct:P.
#
# Usage:
#   regression_gate.sh baseline_metrics.tsv candidate_metrics.tsv
#   regression_gate.sh baseline_metrics.tsv candidate_metrics.tsv 1   # only tier<=1
set -euo pipefail

BASE="${1:?usage: regression_gate.sh <baseline.tsv> <candidate.tsv> [max_tier]}"
CAND="${2:?need candidate.tsv}"
MAXTIER="${3:-9}"

python3 - "$BASE" "$CAND" "$MAXTIER" <<'PY'
import sys
base, cand, maxtier = sys.argv[1], sys.argv[2], int(sys.argv[3])

def load(f):
    d = {}
    for ln in open(f):
        ln = ln.rstrip("\n")
        if not ln.strip() or ln.lstrip().startswith("#"):
            continue
        p = [x.strip() for x in ln.split("\t")]
        d[p[0]] = p[1:]
    return d

B, C = load(base), load(cand)
def num(x):
    try: return float(x)
    except Exception: return None

rows, fails, missing = [], 0, 0
for m, v in B.items():
    if len(v) < 4:
        continue
    val, sd, gate, tier = v[0], v[1], v[2], v[3]
    if int(tier) > maxtier:
        continue
    if m not in C:
        rows.append((m, tier, "SKIP", "absent in candidate")); missing += 1; continue
    cv = C[m][0]
    if gate == "exact":
        ok = (cv == val); detail = f"{cv}  (baseline {val})"
    elif gate == "ge2sd":
        thr = num(val) - 2*num(sd); ok = num(cv) >= thr
        detail = f"{cv} ≥ {thr:.3f}   (baseline {val} ± {sd})"
    elif gate.startswith("band:"):
        _, lo, hi = gate.split(":"); ok = num(lo) <= num(cv) <= num(hi)
        detail = f"{cv} ∈ [{lo}, {hi}]"
    elif gate.startswith("lepct:"):
        p = num(gate.split(":")[1]); thr = num(val)*(1+p/100); ok = num(cv) <= thr
        detail = f"{cv} ≤ {thr:.2f}   (baseline {val} +{p:.0f}%)"
    else:
        ok = False; detail = f"unknown gate '{gate}'"
    rows.append((m, tier, "PASS" if ok else "FAIL", detail))
    if not ok: fails += 1

w = max((len(r[0]) for r in rows), default=6)
print(f"{'METRIC':<{w}}  T  VERDICT  DETAIL")
print("-" * (w + 28))
for m, t, verdict, d in rows:
    mark = {"PASS": "✓", "FAIL": "✗", "SKIP": "·"}[verdict]
    print(f"{m:<{w}}  {t}  {mark} {verdict:<5} {d}")
npass = len(rows) - fails - missing
print("-" * (w + 28))
print(f"{npass} pass · {fails} FAIL · {missing} skipped   (gate: regression = >2sd below baseline)")
sys.exit(1 if fails else 0)
PY

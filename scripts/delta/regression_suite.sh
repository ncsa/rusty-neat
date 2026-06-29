#!/usr/bin/env bash
# Run a regression tier on the CURRENT build and collect metrics for regression_gate.sh.
# Thin orchestrator over the existing harnesses — no new simulation logic.
#
#   MODE=submit  TIER=0|1   — submit the tier's jobs, write a manifest of (kind prefix jobid outdir)
#   MODE=collect MANIFEST=… — after the jobs finish, extract metrics → stdout (the candidate TSV)
#
# Typical flow (on a feature branch, from repo root):
#   MODE=submit TIER=1 LABEL=cand bash scripts/delta/regression_suite.sh      # submits, prints manifest path
#   # ...wait for the jobs (squeue --me)...
#   MODE=collect MANIFEST=<path> bash scripts/delta/regression_suite.sh > candidate_metrics.tsv
#   bash scripts/delta/regression_gate.sh scripts/delta/baseline_metrics.tsv candidate_metrics.tsv "$TIER"
#
# Tier 2 (variety, SV per-type reps, sweeps) is heavier and seam-specific — drive it
# with run_input_variety.sh / run_cancer_sweeps.sh (REPS=3) + their collect_* scripts,
# then append rows to the candidate TSV (sv_<type>_recall, etc.). See docs/regression_protocol.md.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
source "$REPO_ROOT/scripts/delta/lib_report.sh"   # $SCRATCH + RESULTS_DIR

MODE="${MODE:-submit}"
TIER="${TIER:-1}"
LABEL="${LABEL:-cand}"
DATA="${DATA:-$SCRATCH/neat_data}"
MANIFEST="${MANIFEST:-${RESULTS_DIR:-$HOME}/regression_${LABEL}.manifest}"
D="$REPO_ROOT/scripts/delta"

# ── submit ───────────────────────────────────────────────────────────────
if [[ "$MODE" == "submit" ]]; then
    mkdir -p "$(dirname "$MANIFEST")"
    echo "# kind prefix jobid outdir" > "$MANIFEST"
    sub() { # kind prefix outdir sbatch-script env...
        local kind="$1" prefix="$2" out="$3" script="$4"; shift 4
        local jid; jid=$(env "$@" OUTDIR="$out" sbatch --parsable "$D/$script")
        echo "$kind $prefix $jid $out" | tee -a "$MANIFEST"
    }
    # Tier 0 — determinism + ecoli fidelity + perf (perf job covers all its genomes)
    sub order   -     "$SCRATCH/reg_${LABEL}_order"      run_order_independence.sbatch REFERENCE="$DATA/yeast.fa"
    sub germline ecoli "$SCRATCH/reg_${LABEL}_germ_ecoli" germline_e2e.sbatch          REFERENCE="$DATA/ecoli.fa" TOOLS=rneat
    sub perf    -     "$SCRATCH/reg_${LABEL}_perf"       benchmark.sbatch
    if [[ "$TIER" -ge 1 ]]; then
        sub germline chr22 "$SCRATCH/reg_${LABEL}_germ_chr22" germline_e2e.sbatch   REFERENCE="$DATA/chr22.fa" TOOLS=rneat
        sub cancer  somatic "$SCRATCH/reg_${LABEL}_cancer"    cancer_pipeline.sbatch REFERENCE="$DATA/chr22.fa" TOTAL_COVERAGE=30 PURITY=0.6 PRUNE=1
    fi
    echo
    echo "Submitted TIER=$TIER. Manifest: $MANIFEST   (watch: squeue --me)"
    echo "When done: MODE=collect MANIFEST=$MANIFEST bash $0 > candidate_metrics.tsv"
    exit 0
fi

# ── collect ──────────────────────────────────────────────────────────────
[[ -f "$MANIFEST" ]] || { echo "manifest not found: $MANIFEST" >&2; exit 1; }

while read -r kind prefix jobid outdir; do
    [[ "$kind" == \#* || -z "${kind:-}" ]] && continue
    case "$kind" in
      order)
        # parse the RESULTS block from the job's .out (in the submit dir)
        f=$(ls -t "$D/../../rneat-orderindep_${jobid}.out" rneat-orderindep_${jobid}.out 2>/dev/null | head -1 || true)
        [[ -f "$f" ]] || f="rneat-orderindep_${jobid}.out"
        get(){ grep -m1 "$1" "$f" 2>/dev/null | grep -oE 'PASS|FAIL' | head -1; }
        # thread-invariance: PASS iff 1- and N-thread determinism hold and mt==1thread
        ti=FAIL; grep -q 'multithread_vs_1thread:.*same' "$f" 2>/dev/null && \
                 [[ "$(get determinism_1thread)" == PASS ]] && ti=PASS
        printf '%s\t%s\n' determinism_thread_invariant "$ti"
        printf '%s\t%s\n' shard_order_independent "$(get shard_order_independent)"
        printf '%s\t%s\n' shard_disjoint          "$(get shard_disjoint)"
        printf '%s\t%s\n' contig_order_independent "$(get contig_order_independent)"
        ;;
      germline)
        python3 - "$outdir/rneat_scored.summary.csv" "$prefix" <<'PY'
import csv, sys
f, pre = sys.argv[1], sys.argv[2]
try: rows = {r["Type"]: r for r in csv.DictReader(open(f)) if r.get("Filter") == "PASS"}
except Exception: rows = {}
def g(t, k):
    try: return f"{float(rows[t][k]):.4f}"
    except Exception: return "NA"
print(f"{pre}_snp_recall\t{g('SNP','METRIC.Recall')}")
print(f"{pre}_snp_precision\t{g('SNP','METRIC.Precision')}")
print(f"{pre}_indel_recall\t{g('INDEL','METRIC.Recall')}")
print(f"{pre}_indel_precision\t{g('INDEL','METRIC.Precision')}")
try: print(f"{pre}_tstv\t{float(rows['SNP']['QUERY.TOTAL.TiTv_ratio']):.3f}")
except Exception: print(f"{pre}_tstv\tNA")
PY
        ;;
      perf)
        # median.tsv: genome size_mb threads tool wall_s peak_rss_mb  (rneat, 1 thread)
        m="$outdir/median.tsv"
        for g in ecoli chr22; do
          awk -F'\t' -v g="$g" '$1==g && $4=="rneat" && $3==1 {print g"_wall_s\t"$5"\n"g"_rss_mb\t"$6}' "$m" 2>/dev/null
        done
        ;;
      cancer)
        python3 - "$outdir/scored.stats.csv" <<'PY'
import csv, sys
f = sys.argv[1]
try: rows = {(r.get("type") or "").strip().lower(): r for r in csv.DictReader(open(f))}
except Exception: rows = {}
def g(types, k):
    for t in types:
        r = rows.get(t)
        if not r: continue
        for c in r:
            if c.strip().lower() == k:
                try: return f"{float(r[c]):.4f}"
                except Exception: return "NA"
    return "NA"
print(f"somatic_snv_recall\t{g(('snvs','snp','snps'),'recall')}")
print(f"somatic_indel_recall\t{g(('indels','indel'),'recall')}")
PY
        ;;
    esac
done < "$MANIFEST"

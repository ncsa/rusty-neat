# Replication coverage audit

**Question:** before moving on from validation, which results rest on enough runs
to trust, and which are single-draw point estimates that warrant replication —
even where the single number looks good?

**Short answer:** there are **no qualitative gaps** (every planned validation ran),
but most *fidelity/recall* numbers are **single draws**. Performance and the
scaling/allocator work are already replicated. The one statistically *fragile*
headline is the **§3.7 SV per-type recall** (single run **and** small event
counts). Determinism/order-independence need no reps (they are exact).

## Classification

Each result falls into one of four buckets:

- **EXACT** — deterministic; repeating yields identical output, so n=1 is definitive.
- **REPLICATED** — already run ≥2× (median or mean±sd reported).
- **SINGLE / TIGHT** — one run, but over a large denominator (≥10⁴ variants), so the
  sampling CI is already tight; reps only add error bars.
- **SINGLE / FRAGILE** — one run over a small denominator or a noisy measurement;
  the point estimate could move on a re-draw → **reps needed**.

| § | Result | n | Denominator / nature | Bucket |
|---|---|---|---|---|
| 3.1 | Determinism (1 vs 8 thread, byte-identical) | 1 | exact (md5) | **EXACT** |
| 3.10 | Order-independence (determ., thread-inv., shard-order, disjoint) | 1 | exact (md5) | **EXACT** |
| 3.1 | Coverage / breadth / Ts/Tv / het-hom | 1 | whole-chr aggregates | SINGLE / TIGHT |
| 3.2 | Performance vs NEAT (wall, RSS) | **3** | median of 3 (`benchmark.sbatch REPS=3`) | **REPLICATED** |
| 3.6 | Thread × allocator × NUMA isolation | **2** | 2 reps | **REPLICATED** |
| 3.6 | Procs-per-node throughput, whole-genome wall | ≥3* | `run_genome_reps REPS=3` + mean±sd | **REPLICATED*** |
| 3.3 | Germline fidelity chr22 (SNP/indel recall) | 1 | ~10⁵ SNVs → tight CI | SINGLE / TIGHT |
| 4 | Input variety (5 genomes) | 1 ea | 10⁵–4·10⁵ variants → tight | SINGLE / TIGHT |
| 3.4 | Cancer end-to-end (somatic SNV/indel) | 1 | moderate somatic count | SINGLE / FRAGILE |
| 3.9 | Cross-caller (Delly, VarScan2) | 1 | moderate | SINGLE / FRAGILE |
| 4 | Cancer sweeps (purity / coverage / tmr) | 1/pt | moderate; some small | SINGLE / FRAGILE |
| 3.7 | **SV per-type recall (DEL/DUP/INV)** | 1 | **small event count** (e.g. INV recall 1.000 may be 3/3) | **SINGLE / FRAGILE ← top priority** |
| 4 | Per-tissue SV spectrum | 1 | ~200–280 SVs → adequate | SINGLE / TIGHT |

\* The reps *machinery* (`run_genome_reps.sh` + `aggregate_reps.sh`, mean±sd per K)
exists and was used for the K-sweep; confirm the headline whole-genome wall-clock
(3 h 41 m shared / 2 h 30 m exclusive) is a rep mean, not a single run — and note
that the **shared-partition** number is inherently high-variance (cross-tenant
contention) and should be reported as a range, not a point.

## What this means

- **The conclusions don't change** with reps — the qualitative findings (flat
  fidelity across genomes, the purity-collapse mechanism, coverage saturation,
  tissue spectra, thread-vs-process scaling) are mechanistically explained and
  robust. Reps buy **error bars and fluke-protection**, not new conclusions.
- **The one number not to publish bare is §3.7 SV per-type recall.** With only a
  handful of each SV type in the chr22/60×/0.8 truth, `INV = 1.000` and
  `DEL = 0.882` carry wide confidence intervals. This is the must-replicate item.

## Prioritized replication plan

1. **§3.7 SV per-type recall — MUST.** 3 reps at 60×/0.8 with distinct seeds *and*
   an elevated SV count (more events per run shrinks the CI directly, as the
   tissue-spectrum run showed). Report mean ± sd per type.
2. **Cancer recall — SHOULD.** 3 reps of the end-to-end point (§3.4) and the key
   sweep points (purity 0.5/0.9, coverage 60×, tmr 1e-5). Confirms the purity
   collapse and resolves the suspicious `0.7 ≡ 0.9` four-digit tie (saturation vs
   artifact).
3. **Germline fidelity + input variety — NICE.** 3 reps each for error bars; cheap,
   and turns point estimates into mean ± sd. Conclusions won't move.
4. **Whole-genome wall-clock — CONFIRM + CAVEAT.** Verify the headline used the
   3-rep machinery; report the shared-partition figure as a range.
5. **Determinism / order-independence — NONE.** Exact by construction; a re-run is
   identical. (Optional one-off: confirm the contig-order sensitivity on a second
   genome to show it generalizes — a correctness check, not a rep.)

## How to run the reps

The pipelines already derive a distinct seed per replicate; the sweep drivers now
take a `REPS` env that submits each condition `REPS` times into `…_rep<r>` output
dirs with distinct seeds (default `REPS=1` = unchanged behavior):

```bash
# input variety, 3 reps per genome
REPS=3 GENOMES="…" bash scripts/delta/run_input_variety.sh

# cancer sweeps, 3 reps per point (any axis)
REPS=3 AXIS=purity bash scripts/delta/run_cancer_sweeps.sh

# SV per-type (top priority): tissue axis is the SV path — 3 reps, elevated SV rate
REPS=3 AXIS=tissue SV_RATE_SCALE=20 bash scripts/delta/run_cancer_sweeps.sh
```

The collectors emit one row per rep (the point/genome label carries the `_rep<r>`
suffix); aggregate to **mean ± sd** by grouping on the base label. A future
`aggregate_*` helper can fold this automatically; for now the per-rep rows are
small enough to average on paste.

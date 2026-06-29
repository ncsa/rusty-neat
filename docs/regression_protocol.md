# rneat re-validation / regression protocol

When a feature changes rneat (e.g. adapter modeling, #125), this protocol confirms
the change does not:

- **(M) blow up memory** — peak RSS regression,
- **(S) slow down** — wall-clock regression,
- **(D) emit errant data** — fidelity or determinism regression.

It reuses the Delta validation harnesses already built, and gates a candidate
build against a frozen baseline with **statistically-grounded tolerances**.

## Why this is now possible: the replication study gives the noise floor

The Phase-2 replication runs measured the *run-to-run* standard deviation of every
metric (see [`replication_audit.md`](replication_audit.md)). That sd is exactly
what separates a **real regression** from normal scatter:

| metric | measured run-to-run sd | gate uses |
|---|---|---|
| SNV recall | ~0.01 | tight band |
| indel recall | ~0.05 | wider band |
| SV per-type recall | ~0.09–0.11 | widest band |
| Ts/Tv | ~0.14 | plausibility band 2.0–2.6 |
| determinism (md5) | 0 (exact) | zero tolerance |
| wall-clock | HPC-noisy → use n=3 median | % band |
| peak RSS | stable | % band + ceiling |

A gate fires when a candidate moves **more than k·sd below baseline** (default
k = 2 ≈ 95 % confidence). Without the reps you could not set these — you'd be
guessing whether a 0.97→0.94 dip is a bug or noise.

## The suite: failure mode → metric → harness → gate

| Mode | Metric | Harness | Gate |
|---|---|---|---|
| **M** | peak RSS | `benchmark.sbatch` (`/usr/bin/time -v`) | ≤ baseline × 1.15 **and** below an absolute ceiling |
| **S** | wall-clock (n=3 median) | `benchmark.sbatch` | median ≤ baseline × 1.20 (absorbs HPC noise) |
| **D** | determinism + thread-invariance | `run_order_independence.sbatch` | **byte-identical (exact)** — any drift fails |
| **D** | germline SNP/indel recall, precision, Ts/Tv | `germline_e2e.sbatch` | recall ≥ base − 2 sd; precision ≥ base − 2 sd; Ts/Tv ∈ [2.0, 2.6] |
| **D** | somatic SNV/indel recall, precision | `cancer_pipeline.sbatch` | ≥ base − 2 sd |
| **D** | SV per-type recall | `sv_pipeline.sbatch` (REPS) | ≥ base − 2 sd (only for SV-touching changes) |

All harnesses already emit the needed numbers via the `collect_*` scripts and
`archive_run` (which also records git SHA + core-hours for provenance).

## Tiers — run the cheapest tier that covers the change

| Tier | Scope | Cost | When |
|---|---|---|---|
| **0 — Smoke** | determinism (exact) + **E. coli** germline fidelity + wall/RSS vs baseline | minutes, 1 small node | **every commit / PR** (cheap enough to be near-CI) |
| **1 — Standard** | chr22 germline + cancer fidelity + perf (n=3) + order-independence | ~1 h | **every feature PR** — the default gate |
| **2 — Deep** | input variety (multi-genome) + SV per-type (REPS=3) + cancer sweeps | hours | **risky/perf/SV features, or pre-release** |

Tier 0 alone catches the big three failure modes grossly (crash, memory blowup,
non-determinism, fidelity collapse). Tiers 1–2 add resolution and breadth.

## Feature-touch matrix — which extra checks a change demands

| If the change touches… | Require |
|---|---|
| read generation / error / quality / **adapters** | Tier 0 + 1, plus read-level checks (below) |
| RNG / threading / sharding / seeding | **order-independence is mandatory** (Tier 1) |
| SV / CNV / cancer mutation model | **Tier 2 SV per-type (REPS=3)** |
| performance / memory / allocation | **benchmark Tier 1 + RSS ceiling mandatory** |
| VCF / truth-set output | fidelity + **truth-VCF determinism (md5)** |

## Worked example — adapter modeling (#125)

Adapters add read-level sequence (3′ adapter readthrough, etc.). Two question sets:

**Regression (must NOT change):**
- Determinism still byte-identical (exact).
- Wall-clock / peak RSS within tolerance — adapter logic must not inflate either.
- Truth VCF unchanged — adapters are a read artifact, not variants (md5 of truth
  VCF vs no-adapter run must match).
- **Post-trim fidelity == baseline:** simulate with adapters → trim (cutadapt /
  Trim Galore) → align → call; SNP/indel recall/precision must match the
  no-adapter baseline within 2 sd. Adapters that survive trimming and corrupt
  calls would show up here.

**New-feature checks (adapter-specific, added to Tier 0/1):**
- Adapters present in raw reads at the configured rate and 3′ position
  (cutadapt detection count ≈ expected).
- A standard trimmer recovers the pre-adapter read set (read count + base-level).
- No adapter bases leak into the aligned, properly soft-clipped fraction.

**Tier plan for adapters:** Tier 0 on every commit (incl. the adapter-presence
check), Tier 1 on the PR; Tier 2 only if the same change also touches SVs.

## Baseline management

- Freeze baseline metrics on `develop` into a committed `baseline_metrics.tsv`
  (one row per metric: value, sd, n, git SHA, date). Re-freeze on each release.
- A candidate run produces `candidate_metrics.tsv` in the same shape.
- `regression_gate.sh` diffs the two, applies the gates + tolerances above, prints
  a PASS/FAIL table, and **exits non-zero on any FAIL** (CI / `--dependency`-friendly).

## Automation (built)

Three files implement this, thin wrappers over the existing harnesses — no new
simulation logic:

- **`scripts/delta/baseline_metrics.tsv`** — the frozen baseline, *seeded from the
  Phase-1/2 validation + replication numbers* (value, sd, gate, tier). Re-freeze on
  each release by running the collect step on `develop`.
- **`scripts/delta/regression_suite.sh`** — `MODE=submit TIER=n` submits the tier's
  jobs (order-independence + germline + perf, plus chr22 germline/cancer at Tier 1)
  and writes a manifest; `MODE=collect MANIFEST=…` extracts metrics into a candidate
  TSV. (Tier 2 — variety, SV per-type reps, sweeps — is driven by
  `run_input_variety.sh` / `run_cancer_sweeps.sh REPS=3` and appended.)
- **`scripts/delta/regression_gate.sh baseline candidate [max_tier]`** — applies the
  gates, prints a PASS/FAIL table, and exits non-zero on any FAIL (CI-friendly).

Usage on a feature branch (from repo root):

```bash
MODE=submit TIER=1 LABEL=cand bash scripts/delta/regression_suite.sh   # submit; prints manifest path
# ...wait for jobs (squeue --me)...
MODE=collect MANIFEST=<path> bash scripts/delta/regression_suite.sh > candidate_metrics.tsv
bash scripts/delta/regression_gate.sh scripts/delta/baseline_metrics.tsv candidate_metrics.tsv 1
```

The gate logic and all collect parsers are unit-tested; the `sbatch` submit path
runs on Delta. Capturing a fresh baseline is just the same collect step pointed at
a `develop` run.

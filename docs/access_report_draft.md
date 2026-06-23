# Benchmarking rneat: Enhancing Performance of a Genome Sequencing Simulator for Cancer Genomics with Rust

**ACCESS allocation report — DRAFT**
Allocation: `bhrd-delta-cpu` (NCSA Delta). Core-hours used to date: ~2,500 of ~397,500.

---

## 1. Project summary

`rneat` is a Rust port and extension of NEAT, a next-generation-sequencing read
simulator. It generates FASTQ, a golden BAM, and a truth VCF whose statistical
properties match real data, and it adds a native tumor/normal cancer workflow
(`gen-cancer-reads`) with structural variants and an origin-tagged somatic truth
set. This project used Delta to (1) verify rneat produces high-quality data on a
small test set, (2) compare rneat to its predecessor NEAT 4 at scale on both
performance and fidelity, and (3) validate the cancer workflow end-to-end through
a real somatic-variant-calling pipeline.

A central outcome is methodological: **benchmarking rneat against NEAT 4 and
against GATK/Mutect2 surfaced — and we fixed — five real simulator defects** that
unit tests had not caught. The verification process is itself a result.

---

## 2. Methods

All runs on Delta CPU nodes (dual-socket AMD EPYC, 128 cores). Pipeline tools:
BWA-MEM2 (alignment), GATK HaplotypeCaller (germline calling), GATK Mutect2 +
FilterMutectCalls (somatic calling), hap.py / som.py (scoring against truth),
NEAT 4.6.1 (reference simulator). Every run archives its result artifacts and
SLURM resource usage (core-hours) to durable project storage; `collect_report.sh`
aggregates them.

Reference genomes: chr22 (GRCh38) for the verification tier; E. coli, yeast for
performance scaling. Default 30× coverage, 151 bp paired-end, ploidy 2 unless
noted.

---

## 3. Results (Phase 1 — complete)

### 3.1 Simulator data quality (chr22, 30×)

| Property | Result |
|---|---|
| Determinism (FASTQ + VCF, 1 vs 8 threads) | **byte-identical (PASS)** |
| Mean coverage (covered bases) | 29.96× (requested 30×) |
| Breadth | 0.77 (chr22 N-gaps correctly left uncovered) |
| Heterozygous allele fraction | 0.486 (clean diploid binomial ~0.5) |
| Ts/Tv | **2.35** (human-realistic) |

Determinism holding byte-for-byte across thread counts is a guarantee NEAT 4 does
not make.

### 3.2 Performance vs NEAT 4 (single thread)

| Genome | rneat wall | NEAT wall | **Speedup** | rneat RSS | NEAT RSS | **Memory** |
|---|---|---|---|---|---|---|
| E. coli (4.5 MB) | 33.0 s | 88.2 s | **2.7×** | 41 MB | 280 MB | **6.9×** |
| yeast (11 MB) | 77.0 s | 214.7 s | **2.8×** | 53 MB | 172 MB | **3.2×** |
| chr22 (49 MB) | 255.0 s | 749.0 s | **2.9×** | 227 MB | 1528 MB | **6.7×** |

rneat is consistently ~2.7–2.9× faster and uses 3–7× less, flatter memory.

**Multicore behavior (current default).** With its default per-contig
parallelism, rneat does not gain from added threads on Delta's dual-socket EPYC —
wall-clock is flat-to-slightly-higher with more threads (yeast: 1 thread 77 s →
4 threads 110 s, plateau ~115 s), while on a single-socket machine it scales
modestly (2 threads = 1.35×). rneat's single-thread efficiency, not thread
scaling, is the source of its speed advantage. The behavior is consistent with
memory-bandwidth saturation / allocator contention on NUMA hardware; Phase 2
investigates whether code changes can recover multicore scaling, in which case
these numbers will be revised.

### 3.3 Germline fidelity vs NEAT 4 (chr22, 30×, identical pipeline)

Caller recovery of each simulator's truth set through the same GATK pipeline:

| | rneat SNP | NEAT SNP | rneat INDEL | NEAT INDEL |
|---|---|---|---|---|
| Recall | 0.989 | 0.989 | 0.982 | 0.969 |
| Precision | 0.9996 | 0.991 | 0.989 | 0.972 |

**Statistically equivalent**, with rneat marginally cleaner on precision. rneat
faithfully reproduces NEAT 4's germline behavior.

### 3.4 Cancer end-to-end (chr22, 30×, purity 0.6)

Mutect2 tumor/normal recovery of rneat's somatic truth (after fixes + truth/query
normalization):

| | Recall | Precision |
|---|---|---|
| Somatic SNV | 0.925 | 0.872 |
| Somatic indel | 0.900 | 0.844 |

rneat's simulated tumor/normal data is recovered well by a standard somatic
pipeline. **Scope:** SNV/indel only (Mutect2 does not call SVs), single
purity/coverage point — both expanded in Phase 2.

### 3.5 Defects found and fixed through verification

Benchmarking surfaced five real bugs invisible to unit tests; each was diagnosed,
fixed, and re-verified end-to-end:

| # | Defect | Effect when fixed |
|---|---|---|
| #280 | Default mutation model near-random (uniform substitution matrix) | Ts/Tv 0.72 → 2.35 (NEAT parity) |
| #285 | Alt alleles written only to forward-strand reads | somatic SNV recall 0.33 → 0.94 |
| #287 | Indels garbled on reverse-strand reads | indel precision 0.07 → 0.49 (read-level) |
| #289 | Coverage regression introduced by #287 (R2 deletion truncation) | restored −14.5% lost read pairs |
| #290 | som.py indel normalization off by default (scoring) | cancer indel recall 0.60 → 0.90 |

All affect every paired-end run, not just cancer. The somatic-SNV journey
(0.33 → 0.94) and cancer-indel journey (0.20 → 0.90) trace the cumulative impact.

### 3.6 Multicore scaling and the whole-genome strategy

Benchmarking exposed a scaling characteristic that reshaped the whole-genome
plan. On Delta's dual-socket EPYC nodes, rneat's in-process (rayon) threading
initially *regressed* — yeast 30× took 218 s on 1 thread but 278 s at 2 and
337 s at 16, monotonically slower with more threads. A systematic isolation
sweep (thread count × {glibc, mimalloc, jemalloc} × {default, numactl
--interleave}, 2 reps) showed *swapping* the allocator or NUMA policy made no
difference — but that pointed at allocation *volume* rather than allocator
*choice*. Eliminating the read loop's per-base heap allocation removed the
at-scale degradation: 16-thread dropped from 337 s to **249 s (−26 %)**, now
within ~14 % of single-thread, and the curve recovers with thread count instead
of diverging.

A residual 1→2-thread step (+28 %, ~218→280 s) remains; it is allocator- and
NUMA-independent (identical under mimalloc/jemalloc), i.e. a memory-bandwidth /
cross-socket floor, not something further code removes. Net: in-process
threading is now **break-even — no longer harmful, but no net speedup** over
single-thread, which remains the fastest per-worker mode.

The decisive experiment for the whole-genome strategy was a **procs-per-node
sweep**: K independent single-threaded rneat *processes* concurrently on one
node. Processes **scale** — aggregate throughput rises ~3.5× from 1 to 64
(0.0118 → 0.0409 jobs/s) before the bandwidth ceiling plateaus it:

| processes/node | 1 | 2 | 4 | 8 | 16 | 32 | 64 | 128 |
|---|---|---|---|---|---|---|---|---|
| throughput (jobs/s) | 0.012 | 0.014 | 0.018 | 0.028 | 0.037 | 0.041 | **0.041** | 0.039 |
| per-process efficiency | 100% | 60% | 39% | 29% | 19% | 11% | 5% | 3% |

Separate processes share no thread pool, sequence-block, or counter state and no
allocator, yet still plateau at 32–64 — so the ceiling is genuine per-node
**memory bandwidth**, library-independent. But because each single-threaded
process runs at full speed, packing many per node and spreading across nodes
extracts that bandwidth far better than threading a single process does.

**Consequence — rneat's HPC parallelism model is multi-process region-sharding,
not threading.** A genome is partitioned into anchor-windows (rneat's
generation-time BED filter assigns each fragment to exactly one window in
absolute coordinates, so shards reconstruct a whole-genome run with no
double-counting or gaps), simulated as one single-threaded process per window
via a SLURM array spread across nodes, then concatenated.

**First full GRCh38 run — and the contention discovery.** The initial
whole-genome run (306 × 25 Mb shards, one single-threaded rneat per shard,
submitted as a SLURM array on Delta's *shared* `cpu` partition with a 16-core
slice per shard) completed in **3 h 41 m** — far above the ~20–26 min the
procs-per-node sweep projected. Per-shard timing was sharply **bimodal**: median
2.5 min, but p90 131 min and a worst shard of **170 min** for an ordinary 25 Mb
window that runs in ~8 min on an unloaded node (confirmed by re-running that exact
window on an isolated node: 7 m 53 s). Grouping shard time by compute node exposed
the cause: the slow shards clustered on nodes holding only 1–2 of our shards,
while nodes that ran 12–17 of our shards each were uniformly fast — the *inverse*
of what our own packing would produce. The slowdown therefore came from **other
tenants' jobs** co-scheduled on those nodes: because rneat is
memory-bandwidth-bound (above), a bandwidth-hungry neighbour starves it,
inflating an 8-min window 10–20×. On a shared partition this is uncontrollable
and makes the whole-genome wall-clock irreproducible.

**Fix — exclusive nodes, packed by us.** Requesting whole nodes (`--exclusive`)
and packing K shards per node ourselves removes every stranger: the only
remaining contention is our own K processes competing for the node's bandwidth.
This worked — a single isolated 25 Mb window dropped back to **7 m 53 s** — but
the path to a clean whole-genome number then exposed two further trade-offs that,
together with the shared-partition result, form the complete HPC story.

**Trade-off 1 — straggler ceiling vs. completeness.** The first exclusive sweep
paired packing with a *tight* per-task walltime (`10 + 6K` min) as
straggler-proofing. But packing slows each window (Trade-off 2), so the ceiling
silently *killed* the heavy windows — every K's slowest shard landed exactly on
its walltime (K=4: 34 min, K=8: 58 min, K=16: 103 min) and runs finished
incomplete (197–268 of 306 shards). Lesson: once `--exclusive` removes the
cross-tenant straggler, the tight ceiling is not only unnecessary, it is
*harmful*; the walltime must comfortably exceed the packed-window time, and
`--requeue` alone covers genuine node failure. (Fixed to `30 + 20K` min.)

**Trade-off 2 — packing density vs. per-shard speed.** Even with strangers gone,
rneat is hard memory-bandwidth-bound, so our *own* co-packed processes slow each
other roughly linearly: a window that runs ~8 min alone takes ≥34 min at K=4,
≥58 min at K=8, ≥103 min at K=16. The node's bandwidth saturates by ~K=4, so
**packing more shards per node buys no extra throughput** — it only inflates
per-shard latency. Wall-clock is therefore minimized by *low* K (each shard near
full speed) spread across *many* nodes, not by dense packing.

**Trade-off 3 — compute speed vs. schedulability.** But low-K-many-nodes is
exactly the hardest thing to schedule: `--exclusive` claims a full 128-core node
even for a 1-core (K=1) task, so minimizing compute wall-clock means requesting
dozens of whole nodes, which on a busy shared cluster queues behind everyone
else (observed: 77- and 153-task low-K arrays stuck `PD (Resources/Priority)`
with no estimable start). The bandwidth-optimal configuration is the
scheduling-pessimal one.

These three trade-offs resolve into one rule and three regimes:

| Approach | Compute speed | Schedulability | Reproducible? | Best when |
|---|---|---|---|---|
| Shared partition, any packing | unpredictable (3 h 41 m; 8-min windows → 170 min under co-tenants) | instant | no | never the right choice |
| Exclusive, **low** K, many nodes | fastest (~8 min/window) | poor on a busy cluster | yes | cluster is idle / reserved nodes |
| Exclusive, **moderate** K, fewer nodes | good | schedules quickly | yes | busy shared cluster (pragmatic) |

**Recommendation for HPC users.** Always run on **exclusive nodes** — a shared
per-core slice exposes a bandwidth-bound simulator to unbounded cross-tenant
contention and can inflate wall-clock 10–20× with no warning. Then **match
packing density to node availability**: pack *light* (1–2 shards/node) across as
many nodes as you can get when the cluster is idle or you hold a reservation
(minimizes compute wall-clock); pack *heavier* (8–16/node) when nodes are scarce
(fewer whole-node grabs → schedules sooner, at the cost of slower per-shard time).
Use a generous per-task walltime plus `--requeue`; never a tight ceiling. Tooling
(`make_shard_beds.sh`, `genome_array.sbatch` [exclusive packing],
`run_genome_reps.sh` [K-sweep × reps], `aggregate_reps.sh`, `merge_shards.sh`,
`tune_procs_per_node.sbatch`) is implemented and validated.

> **[PENDING — clean exclusive numbers]** A completeness-verified K=2/K=4 run
> (generous walltime, on whatever exclusive nodes schedule) will fill the exact
> optimal whole-genome wall-clock ± sd. The trade-off structure above is final
> regardless of the precise figure.

The honest framing for the report: rneat's advantage is **single-thread
efficiency + low, flat memory**, and HPC throughput is reclaimed by **process-
level sharding** — a defensible, mechanistically-grounded story rather than a
(false) claim of linear thread scaling. Since each shard runs single-threaded on
a (relatively slow) Delta core, ongoing single-thread optimization (reducing
per-read allocations) compounds directly across the whole genome.

---

## 4. Phase 2 — robustness at scale (planned, key deliverables)

Phase 1 established correctness and performance on a focused test set (chr22 and
small genomes). Phase 2 closes the four gaps that Phase 1 deliberately did not
cover, with an explicit goal of **avoiding over-tuning to chr22**: exercise rneat
on substantial and *varied* inputs to confirm robustness and surface any
remaining defects of the kind already found.

1. **Whole-genome scale.** Run the full germline and cancer pipelines on GRCh38
   (whole genome and large chromosome subsets), not just chr22. Confirms the
   Phase 1 metrics hold at genome scale and that nothing was overfit to one
   chromosome.

2. **Multicore scaling / HPC tuning — substantially complete (§3.6).** The
   thread regression was characterized (memory-bandwidth-bound, allocator- and
   NUMA-independent) and the strategy resolved: multi-process region-sharding.
   The first full GRCh38 sharded run surfaced a further finding — on a *shared*
   partition, cross-tenant memory-bandwidth contention inflated the wall-clock to
   3 h 41 m (vs an ~8 min/window unloaded cost). The remaining item is the
   **exclusive-node packing sweep** (K ∈ {4, 8, 16} shards/node × reps) to chart
   the reproducible optimal whole-genome wall-clock and quantify the
   shared-vs-exclusive gap that grounds the HPC-usage recommendation.

3. **Exercise the new cancer features deeply.** Validate structural-variant
   realism downstream with SV-aware callers (Manta / Delly / GRIDSS), which the
   SNV-focused Mutect2 path does not test. Exercise CNVs, BND/INV/INS, per-tissue
   cancer models (BRCA / skin / lung), and the purity-mixing machinery — the
   areas most likely to harbor bugs of the type found in Phase 1.

4. **Input variety and parameter sweeps.** Run across a range of genome sizes and
   compositions (bacterial, fungal, invertebrate, mammalian) and across cancer
   parameters (purity, coverage, tumor mutation rate) to ensure rneat handles
   diverse inputs rather than a single tuned case.

5. **Simulation-only timing for the community.** Benchmark *just the rneat
   simulation stage* (FASTQ + truth VCF generation, no downstream alignment or
   variant calling) across genome sizes, coverages, and the sharded whole-genome
   configuration, to publish concrete real-world wall-clock / throughput numbers
   users can plan around (e.g. "whole human genome at 30× in N minutes on K
   nodes"). This isolates the simulator's own performance — rneat's core
   contribution — from pipeline overhead, and is the headline deliverable for
   communicating rneat's HPC performance to the community. (The sharded GRCh38
   run in item 2 is the first such measurement at whole-genome scale.)

---

## 5. Phase 3 — stretch goals (time permitting)

- **Long-read simulation.** Extend testing to long-read-style data
  (longer fragments / single-molecule profiles) and validate against long-read
  aligners/callers.
- **Plant genomes and polyploidy.** Exercise large, repetitive plant genomes and
  scope tuning for polyploid simulation (rneat's `ploidy` parameter and the
  allele-dosage work tracked separately), where structural variation and high
  copy number stress the simulator differently than human data.

---

## 6. Resource usage

| | Core-hours |
|---|---|
| Documented production runs (9, archived) | **491** |
| Total charged incl. exploratory / failed / timed-out runs | ~2,500 |
| Remaining allocation | ~397,500 |

The 9 archived production runs (3 baseline, 1 benchmark, 4 cancer, 1 germline)
total 491 core-hours. The gap to the ~2,500 charged is dominated by a single
timed-out benchmark that ran `--exclusive` (charging all 128 cores for 8 h ≈
1,000 core-hours before being redesigned) plus failed/iterated runs during
bug-fixing. Phase 1 was deliberately economical (chr22-scale); Phase 2's
whole-genome runs and tuning sweeps are the primary consumers of the remaining
allocation, where the binding constraint is per-run wall-clock, not core-hours.

---

_Run-level metrics and provenance for every job (version, git commit, reference,
core-hours) are archived under the project results directory and assembled by
`scripts/delta/collect_report.sh`._

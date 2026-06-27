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
other **superlinearly**. From the clean, complete (uncapped) runs, the heaviest
25 Mb window took ~8 min alone (K=1), **29 min at K=2, and 101 min at K=4** — the
node's bandwidth saturates almost immediately (by K=2), so packing more shards per
node buys no extra throughput, only steeply rising per-shard latency. (The
human-genome penalty is harsher than the earlier yeast procs-per-node sweep
implied: a 25 Mb window's working set drives far more memory traffic than yeast.)
Wall-clock compute is therefore minimized by *low* K (each shard near full speed)
spread across *many* nodes, not by dense packing.

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

**Measured results (each run complete, 306/306 shards).** Two
completeness-verified exclusive runs against the shared baseline:

| Run | Whole-genome wall | Heaviest window | Median window | Effective nodes |
|---|---|---|---|---|
| Shared partition (16-core slice) | 3 h 41 m | 170 min | 2.4 min | ~varies (contended) |
| Exclusive, **K=2** | **2 h 30 m** | **29 min** | 0.8 min | ~30 (in waves) |
| Exclusive, **K=4** | 12 h 32 m | 101 min | 0.4 min | ~10 (starved) |

Two results stand out. First, **exclusive K=2 beat the shared run on every axis** —
faster wall-clock (2 h 30 m vs 3 h 41 m) and, more durably, a **6× tighter
per-window ceiling (29 min vs 170 min)** with no contention outliers; that ceiling
improvement is reproducible and scheduling-independent, where the wall-clock is
not. Second, **K=4's 12.5 h wall is almost entirely scheduling latency, not
compute**: its heaviest window was 101 min, but its 77 whole-node tasks trickled
through only ~10 concurrently-available exclusive nodes (vs ~30 for K=2, which was
submitted first and got nodes first). The K2/K4 wall gap is therefore
node-availability luck, not packing — Trade-off 3 made literal. With full node
availability, K=2's compute floor is a single ~29-min wave (153 nodes), i.e. a
whole human genome at 30× in **~30 min** when the nodes are actually there; the
2 h 30 m reflects getting only ~30 of them, in waves, on a busy day.

The honest framing for the report: rneat's advantage is **single-thread
efficiency + low, flat memory**, and HPC throughput is reclaimed by **process-
level sharding** — a defensible, mechanistically-grounded story rather than a
(false) claim of linear thread scaling. Since each shard runs single-threaded on
a (relatively slow) Delta core, ongoing single-thread optimization (reducing
per-read allocations) compounds directly across the whole genome.

### 3.7 Structural-variant validation (chr22, Manta + truvari)

The cancer workflow's structural variants (DEL / DUP / INV / BND / CNV) had only
ever been validated through Mutect2, which scores SNVs and indels — never through
an SV-aware caller. We closed that gap: simulate an SV-rich tumor/normal pair
(60×, purity 0.8), call somatic SVs with **Manta** (tumor/normal mode), and score
against rneat's SV truth with **truvari**.

On the classes Manta and truvari can score (DEL/DUP/INV), rneat's variants —
**including large events up to 1 Mb** — are recovered at **0.88–1.00 recall
(32/35 = 91%)**. The strongest evidence, though, looks *below* the caller, at the
reads themselves: rneat does not edit the reference; it emits breakpoint signal
through a dedicated chimeric-read pass that stitches flanks across each junction,
and that signal is present and correctly placed for **every** class — including
where the caller fails to recover it.

- **CNV — two ways.** Manta cannot call copy number. (i) *By depth:* a `CN=4`
  region showed 110.5× vs a 60.4× flank = **1.83×**, matching the expected
  `(0.8·4 + 0.2·2)/2 = 1.8×` at purity 0.8. (ii) *By a copy-number caller:* GATK
  `DenoiseReadCounts`→`ModelSegments`→`CallCopyRatioSegments` recovered **4/7 truth
  CNVs with the correct gain/loss direction** (the misses are the small/low-amplitude
  het CNVs — kb-scale events near the read-depth segmentation limit). Both confirm
  rneat's CNVs carry correct, caller-detectable copy-ratio signal.
- **Somatic specificity (no leak).** A homozygous somatic deletion was depleted
  5× in the tumor (60×→12×) while the normal stayed at full depth — the SV is real
  in the reads and correctly tumor-restricted.
- **BND signal is in the reads even when uncalled.** Every truth-BND locus
  inspected carried **8–73 discordant/split reads** — the same range as the
  *detected* BNDs. Manta calls only ~half of them (re-representing tandem-dup
  pairs as `DUP:TANDEM` at the correct coordinates); the misses are a *caller*
  limitation — translocations are the hardest short-read class — not absent signal.

Per-type recovery (Manta + truvari, 60×/0.8):

| SV type | Truth | Recovery | Bounded by |
|---|---|---|---|
| DEL | 17 | **0.882** (15/17) | — (incl. 376 / 645 kb) |
| DUP | 14 | **0.929** (13/14) | — (incl. 1 Mb) |
| INV | 4 | **1.000** (4/4) | — |
| BND | 22 | signal present at read level; ~50% Manta-called | truvari can't benchmark breakends; Manta translocation recall |
| CNV | 7 | depth 1.83× + GATK caller 4/7 by direction | Manta does not call CNV |

Overall figures across all 64 truth SVs (recall 0.500, precision 0.561) are bounded
by BND (unscoreable by truvari) and CNV (uncallable by Manta) — tool constraints,
not simulator defects.

**Confirmed by an independent second caller.** Delly (read-pair/split/depth-based,
no assembly) recovers the same SVs at the *same* rates — **DEL 0.882, DUP 0.929,
INV 1.000, identical to Manta** — with slightly higher precision (0.667 vs 0.561),
and likewise scores 0/22 BND and 0/7 CNV. Two independent callers agreeing this
closely is strong evidence rneat's SV data is not tuned to one tool, and that the
BND/CNV gaps are scorer/caller limitations (truvari cannot benchmark breakends;
neither caller emits truvari-matchable CNV) rather than simulator defects.

**Per-tissue models shift the spectrum correctly.** Swapping the pan-cancer model
for tissue-specific COSMIC SV models reproduces each tissue's dominant SV type in
the truth set, tracking the model's type probabilities, while per-type *recovery*
stays tissue-independent:

| Tissue | Dominant truth type (✓ model) | DEL / DUP / INV recall | BND frac | Overall |
|---|---|---|---|---|
| BRCA | **DUP** (37, > DEL 27) — model Dup 0.32 | 0.89 / 0.94 / 1.00 (0.94) | 32% | 0.585 |
| skin | **BND**-enriched (45) — model Bnd 0.31 | 0.87 / 0.94 / 0.80 (0.89) | 47% | 0.432 |
| lung | **DEL** (35) — model Del 0.35 | 0.86 / 0.93 / 1.00 (0.90) | 27% | 0.614 |

Each tissue's largest bucket matches its model's largest type probability.
Scoreable-class recall holds at **~0.89–0.94 regardless of tissue**; a tissue's
*overall* recall simply tracks how much of its spectrum is BND (skin 47% → 0.43,
lung 27% → 0.61) — a scorer limitation, not detection quality.

**Recall vs. tumor purity — robust, with a coverage-model caveat at the extreme.**
Holding model and total coverage fixed (60×) and sweeping purity, DEL+DUP+INV recall
is **0.91 / 0.94 / 0.94** at purity 0.3 / 0.5 / 0.7, then drops to **0.17 at 0.9**.
The collapse is *not* a detection failure: `cancer_simulate` splits a fixed coverage
budget by purity (`normal = (1−purity)·total`), so at purity 0.9 the matched normal
falls to **6×**, and Manta — unable to score somatic confidence against so thin a
normal — filters 37 of 68 calls as `MinSomaticScore` (which `--passonly` then drops).
Recall is stable as long as the matched normal stays ≳18× (purity ≤0.7 here); the
extreme-purity drop is a simulation coverage-model artifact (#315: decouple
matched-normal depth from purity), not SV generation. Real tumor/normal pairs
sequence the normal at an independent depth, avoiding this.

This is the first end-to-end validation of rneat's structural-variant output, and
the read-level signal is correct across all classes and sizes. (Detection needs
adequate depth and SV count: a 30×/0.6 chr22 run yields too few somatic SVs to
measure; 60×/0.8 gives a scoreable set.)

**Honest caveats — tracked as realism enhancements (epic #311).** The simulated SVs
are *idealized*: clean cuts in unique sequence, lacking the microhomology, imprecise
breakpoints, and repeat / segmental-duplication context that make real somatic SVs
hard to call. The 0.88–1.00 DEL/DUP/INV recall is therefore an **upper bound**
relative to real, repeat-embedded variants (#312). Two further items: rneat encodes
tandem-dup junctions as generic BND pairs rather than the canonical `DUP` that
callers and truth sets expect (#313), and SV-scale insertions did not appear in any
run (`INS: 0`) despite a non-zero model probability (#314). None is a defect in the
validated read generation — all are realism / interoperability improvements.

### 3.8 At-scale validation (chr1–3 subset, ~690 Mb)

To confirm the chr22 metrics are not overfit to one chromosome, both pipelines were
re-run on a **chr1–chr3 subset of GRCh38** — ~690 Mb, ~13× chr22, entirely different
sequence. Every metric held or slightly improved:

| Pipeline | Metric | chr1–3 | chr22 |
|---|---|---|---|
| Germline (rneat) | SNP recall / precision | **0.990 / 0.9997** | 0.989 / 0.9996 |
| Germline (rneat) | INDEL recall / precision | **0.988 / 0.990** | 0.982 / 0.989 |
| Cancer (Mutect2 T/N) | somatic SNV recall / precision | **0.944 / 0.910** | 0.925 / 0.872 |
| Cancer (Mutect2 T/N) | somatic INDEL recall / precision | **0.908 / 0.885** | 0.900 / 0.844 |

Germline Ts/Tv held at 2.34 (truth = query). rneat's fidelity is therefore not an
artifact of chr22 — it holds at 13× the scale on different sequence. Whole-genome
runs were not pursued (no divergence to investigate); the Mutect2 tumor/normal step
is the practical scale ceiling (~3.2 h on chr1–3, and ≫48 h projected genome-wide
without interval scatter).

*Tooling note:* hap.py/som.py hardcodes UCSC chr-prefixing of input VCFs, so an
Ensembl-named reference (GRCh38: `1/2/3`) needs a chr-prefixed copy for scoring —
now handled automatically by the cancer pipeline.

### 3.9 Cross-caller coverage and mutational-signature fidelity

To avoid judging rneat's cancer features by a single caller, second callers were
added across classes on the existing Delta harnesses — **Delly** (somatic SV, vs
Manta), **Strelka2** (somatic SNV/indel, vs Mutect2), and **GATK somatic-CNV** (vs
the depth check of §3.7) — plus a **mutational-signature** check
(SigProfilerAssignment). Cross-caller *agreement* is the strongest anti-overfit
signal; the signature check probes something variant-recall cannot: does rneat
reproduce the COSMIC mutational *signature* its tumor model encodes?

**At the signature level, it does not — a real fidelity limitation.** Fitting
COSMIC signatures to 6,225 simulated somatic SNVs (chr1–3) assigned them entirely
to flat, clock-like signatures (SBS5 45 %, SBS3 31 %, SBS41 17 %, SBS12 7 %), with
**no SBS1** — the near-universal C>T-at-CpG signature present in essentially every
real tumor. The 96-context spectrum confirms it: the four CpG `[C>T]G` contexts are
the *lowest* C>T contexts (16–23 counts vs ~65 average).

**Root cause (verified both ends).** The bundled COSMIC model is faithful — its
per-context substitution model encodes strong CpG C>T enrichment (conditional
0.78–0.88 at CpG vs 0.39–0.62 elsewhere). But rneat places mutations at a
**context-independent rate** and conditions only the *alt allele* on trinucleotide
context. The realized spectrum is therefore *(genome trinucleotide frequency) ×
(conditional alt)*, not the COSMIC signature; because CpG dinucleotides are ~10×
genome-depleted, the CpG C>T peak never forms. COSMIC signatures are **rate
patterns** (mutations-per-context), which a uniform-placement model cannot
represent.

**Interpretation.** rneat faithfully reproduces the substitution-*type*
distribution and transition bias (Ts/Tv 2.34, §3.1/3.3) and recovers variants well
across independent callers and at scale (§3.3/3.4/3.7/3.8) — but simulated somatic
SNVs do not reconstruct the input COSMIC *signature* under signature extraction.
This is a limitation of the underlying NEAT trinucleotide approach (the germline
default model shares it), and the fix is to weight mutation *placement* by the
signature's per-context probability rather than placing uniformly (tracked: #320).
It is exactly the class of gap recall-based metrics cannot surface — and the reason
the signature check was added.

**Cross-caller results (#317).** The broadened coverage reinforces the anti-overfit
picture across independent callers. **SVs:** Delly reproduces Manta's per-type recall
*exactly* (DEL 0.882 / DUP 0.929 / INV 1.000; §3.7). **SNV/indel:** VarScan2 — a
different statistical model from Mutect2 — calls rneat's somatic variants at high
**precision** (0.93 SNV / 0.94 indel, vs Mutect2's 0.87 / 0.84): the calls it makes
*match the truth*, so the data is not a Mutect2-specific artifact. Its lower recall
(0.63 / 0.53 vs Mutect2's 0.93 / 0.90) reflects VarScan2's more conservative model
and the high-confidence (`Somatic.hc`) filter, not a simulator property — Mutect2
(sensitive) and VarScan2 (precise) bracket the truth as expected for two callers.
(Strelka2 was dropped — its 2018 build SIGSEGVs on Delta's stack. GATK somatic-CNV
recovered **4/7 CNVs by direction** — no-PoN, since its panel step is a Spark tool
incompatible with the env's Java 25 — corroborating the depth check, §3.7.)

---

## 4. Phase 2 — robustness at scale (planned, key deliverables)

Phase 1 established correctness and performance on a focused test set (chr22 and
small genomes). Phase 2 closes the four gaps that Phase 1 deliberately did not
cover, with an explicit goal of **avoiding over-tuning to chr22**: exercise rneat
on substantial and *varied* inputs to confirm robustness and surface any
remaining defects of the kind already found.

1. **Whole-genome scale — done at chr1–3 (§3.8).** Germline and cancer pipelines
   re-run on a ~690 Mb chr1–3 subset of GRCh38 (~13× chr22): all metrics held or
   slightly improved (germline SNP 0.990 / indel 0.988; cancer somatic SNV 0.944 /
   indel 0.908), confirming nothing was overfit to chr22. Full whole-genome not
   pursued — no divergence to investigate, and Mutect2 T/N is the scale ceiling
   (≫48 h genome-wide without interval scatter).

2. **Multicore scaling / HPC tuning — COMPLETE (§3.6).** The thread regression
   was characterized (memory-bandwidth-bound, allocator- and NUMA-independent) and
   the strategy resolved and *measured*: multi-process region-sharding on exclusive
   nodes. The full GRCh38 sweep produced the shared-vs-exclusive contrast (3 h 41 m
   contended → 2 h 30 m clean at K=2, with a 6× tighter per-window ceiling), the
   superlinear packing curve (8→29→101 min/window at K=1/2/4), and the three-regime
   HPC-usage recommendation. The only follow-on is opportunistic: a low-K run on a
   *reservation* to demonstrate the ~30 min compute floor when whole nodes are
   actually available.

3. **Exercise the new cancer features deeply — SV validation DONE (§3.7).**
   Structural-variant realism was validated downstream with Manta + truvari: rneat's
   DEL/DUP/INV (to 1 Mb) recover at 0.88–1.00 recall, CNV confirmed by depth, BNDs
   emitted and partly reassembled by Manta, somatic specificity confirmed — the
   first end-to-end check of the SV machinery, with no simulator defect found.
   Remaining sub-items: exercise the **per-tissue cancer models** (BRCA / skin /
   lung) to confirm tissue-specific SV spectra downstream, and sweep the
   **purity-mixing** machinery across purity/coverage.

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
   sweep in item 2 is the first such measurement at whole-genome scale: a complete
   human genome at 30× in 2 h 30 m on ~30 exclusive nodes, ~30 min compute-bound.)

---

## 5. Phase 3 — stretch goals (time permitting)

- **Long-read simulation (ONT / PacBio-style)** (tracked: #319). Validate rneat against the
  long-read paradigm: simulate single-molecule reads (kb-scale lengths and the
  higher, indel-dominated, homopolymer-dependent error profile characteristic of
  nanopore / HiFi data, unpaired), align with **minimap2**, and score recall with
  long-read callers — **Clair3 / DeepVariant** for SNV/indel and **Sniffles2 /
  cuteSV** for SVs (truvari for SV scoring, as in §3.7). Because a single long read
  spans an SV breakpoint, this is the natural complement to the short-read SV
  validation: it stresses the SV machinery from the other side and would exercise
  the breakpoint-realism work in epic #311 directly. *Prerequisite:* a long-read
  mode in rneat (read-length distribution + long-read error model) — today rneat
  targets short paired-end Illumina data, so this is feature work gated ahead of
  the validation, mirroring the short-read harnesses already built
  (`germline_e2e` / `cancer_pipeline` / `sv_pipeline`).
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

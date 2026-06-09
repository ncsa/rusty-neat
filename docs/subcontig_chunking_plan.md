# Sub-contig chunking for multicore read generation (prototype)

> **Status: implemented, then DISABLED BY DEFAULT.** The machinery works and is
> correct, but benchmarking showed it does not help (read generation is
> memory-bandwidth bound — see **Conclusion** at the end). It ships behind an
> opt-in `chunk_size: <n>` config key, off by default. Notable findings along
> the way: (1) a real determinism bug was caught by testing — the R2 temp file
> was named with whole-contig coords, so chunks of a multi-chunk contig shared
> one R2 file and duplicated/raced R2 reads (fixed to use chunk coords); output
> is now verified byte-identical across thread counts. (2) Chunk-0 uses the
> plain per-contig RNG derivation, so with chunking off the read stream is
> identical to the pre-chunking (develop) behaviour.


## Motivation
`gen-reads` parallelizes per **contig** (one rayon task = one contig). On
few-contig genomes the longest chromosome pins a single core. Benchmark
(NEAT 4.5.3 vs rneat, 8 cores, c_elegans 97 MB, paired-end 10×):

| | 1 thread | 8 threads | speedup |
|---|---|---|---|
| NEAT 4 | 1162 s | 331 s | 3.5× |
| rneat (per-contig) | 379 s | 380 s | **~1.0× (none)** |

rneat is 3× faster single-threaded but gets **no multicore scaling** on
c_elegans (7 chromosomes, dominated by chr V) and ecoli (1 contig). NEAT wins
the default multicore case by splitting each contig into chunks.

## Approach (prototype)
Make a fixed-size **sub-contig chunk** the unit of parallel work instead of a
whole contig.

- **Work-list:** flatten every contig into `(contig_idx, chunk_idx, start, end)`
  items with a target size of `CHUNK_TARGET = 10 Mbp`, evenly split per contig.
  `par_bridge` over the flat list so both big-contig splitting *and*
  many-small-contig parallelism work.
- **Chunk size independent of `num_threads`** → output stays byte-identical
  across thread counts (NEAT ties chunking to thread count, so its output
  changes with `--threads`; rneat keeps its determinism guarantee — a
  differentiator).
- **Boundary fragments:** each chunk owns fragments whose **anchor (left
  coordinate)** falls in `[start, end)`, and may read past `end` into the full
  contig sequence (no read-content stitching). Disjoint anchor ranges ⇒ read
  names stay globally unique even though `frag_idx` resets per chunk.
- **Shared sequence:** an `Arc<SequenceBlock>` per contig, built lazily and
  shared by all its chunks (one sequence copy per contig — memory ≈ today, not
  multiplied by chunk count). `SequenceBlock` keeps `ref_start = 0` and the full
  sequence so `get_subseq` (absolute indexing) is unchanged.
- **Deterministic seeding:** per-chunk RNG = `base_rng.derive_child(contig_idx)
  .derive_child(chunk_idx)`.
- **Merge:** FASTQ temp files concatenated (order-independent, names already
  carry `ref_start`/`ref_end`); BAM bodies ordered by `(contig_idx, chunk_start)`
  = coordinate order; `AdCounter` (`HashMap<usize,(u32,u32)>`) summed per contig;
  `MutatedMap` read from the shared `Arc` for VCF (no per-chunk clone).

## Known prototype limitations (follow-ups)
- **Memory at hg38 scale:** the shared-sequence cache holds each built contig's
  sequence until end of run. Fine for ≤ few-hundred-MB genomes (benchmark);
  for hg38 use `Weak`-based eviction or Tier-2 sliced buffers (fix `get_subseq`
  to honor `ref_start != 0` so a chunk holds only `[start, end)` + overhang).
- **Junction double-count at chunk boundaries:** `suppress_junction_double_count`
  runs per chunk; an SV junction sitting exactly on a chunk boundary may have its
  crossing read-pairs split across two chunks, slightly perturbing suppression.
  Interior junctions (the vast majority) are unaffected.

## Target (hypothesis going in)
c_elegans 8-core: **380 s → ~150–200 s**, memory no worse than the current
~520 MB. Validate determinism (same seed ⇒ identical FASTQ regardless of
`num_threads`) and full test suite green.

## Conclusion (outcome — chunking does NOT help; disabled by default)

The hypothesis was wrong. Chunking delivered no multicore speedup, and the
investigation revealed why: **read generation is memory-bandwidth bound, not
CPU/parallelism-granularity bound.**

Benchmark (NEAT 4.5.3 vs rneat, 8-core desktop, c_elegans 97 Mbp, paired-end
10×, 3-rep medians; rneat timings cross-checked with a same-session develop-vs-
branch A/B):

| config | c_elegans 1-thread | c_elegans 8-thread |
| --- | --- | --- |
| rneat, chunking off (= develop) | ~365 s | ~400 s |
| rneat, chunking on (auto ~1 Mbp) | — | ~440 s |
| NEAT 4 | 1162 s | 331 s |

Evidence the bottleneck is bandwidth, not parallelism granularity:

1. **rneat 8-thread is *slower* than 1-thread regardless of chunking** (~400 s vs
   ~365 s), and develop shows the same — adding cores cannot help when there is
   no spare memory throughput. This is inherent, not a chunking artifact.
2. **Finer chunking monotonically worsened wall time** (25 Mbp → 396 s, 5 Mbp →
   431 s, 1 Mbp → 441 s): more, smaller work-items add overhead with no headroom
   to exploit.
3. **Swapping the global allocator to mimalloc changed nothing** (8-thread ~392 s
   vs ~400 s), ruling out allocator contention. The remaining cost is raw memory
   traffic — per-fragment `get_subseq().to_owned()` + `reverse_complement()`
   copies plus gzip, all bandwidth-heavy.

Why NEAT 4 "wins" multicore (1162 → 331 s, 3.5×) while rneat is flat: NEAT is
CPU-bound (Python interpreter overhead) and has abundant headroom to parallelize
into — but it only parallelizes far enough to *catch up* to rneat's already-fast
single thread (~365 s). rneat starts near the hardware bandwidth ceiling, so it
has little to gain.

**Decision:** keep the chunking machinery (correct, tested, deterministic) but
default it **off** behind `chunk_size`. Do not advertise multicore speed.

**When chunking might still help (the niche it's kept for).** Both conditions
should hold, otherwise leave it off:

1. **Few/one large contig** — the reference is dominated by one or a few big
   chromosomes, so per-contig parallelism leaves cores idle (single-chromosome
   assemblies; one chromosome dwarfing the rest). With many contigs, per-contig
   parallelism already fills the cores and chunking only adds overhead.
2. **Spare memory bandwidth relative to cores** — a multi-socket / many-channel
   HPC node (far higher aggregate bandwidth than the commodity desktop these
   benchmarks ran on), and/or compute-heavier settings (GC-bias-weighted
   coverage, long-read mode, very high coverage) where per-read CPU work
   dominates memory traffic, shifting the workload from bandwidth-bound toward
   CPU-bound.

Recommended starting point: `chunk_size ≈ longest_contig_bp / (4 × cores)`
(~4 chunks/core on the dominant contig). Always A/B against the default on the
target hardware and keep it only if wall time actually improves.

**The only real lever for more throughput** would be cutting per-read memory
traffic in the generation hot loop (reuse fragment buffers instead of allocating
+ copying per read; avoid the `reverse_complement` allocation). That reduces
bandwidth demand and would help *single-threaded* performance too — but it is a
substantial hot-loop rewrite with uncertain payoff, and single-thread is already
~3× faster than NEAT 4. Not pursued here.

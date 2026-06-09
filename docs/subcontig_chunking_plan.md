# Sub-contig chunking for multicore read generation (prototype)

> **Status (implemented):** done on `feature/subcontig-chunking`. Two deviations
> from the original plan below: (1) chunk size is **adaptive to genome size**
> (≈ genome/256, clamped to [1 Mbp, 25 Mbp]) rather than a fixed 10 Mbp, with a
> `chunk_size` config override (`0` disables); still independent of thread count.
> (2) A real bug was caught by determinism testing: the R2 temp file was named
> with whole-contig coords, so chunks of a multi-chunk contig shared one R2 file
> and duplicated/raced R2 reads — fixed to use chunk coords. Output is now
> verified byte-identical across thread counts (unit + e2e regression tests).


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

## Target
c_elegans 8-core: **380 s → ~150–200 s**, memory no worse than the current
~520 MB. Validate determinism (same seed ⇒ identical FASTQ regardless of
`num_threads`) and full test suite green.

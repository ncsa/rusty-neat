# Long-read simulation — epic scope

Status: **scoping draft** (2026-07-12). Tracks: #319 (long-read validation, ACCESS report §5)
and the breakpoint-realism epic #311. This is a *scope/architecture* doc — how long-read
support should be built into rneat — not an implementation plan for a single change.

## 1. Goal & motivation

rneat targets short paired-end Illumina today. A long-read mode (ONT / PacBio-style) would let
it simulate single-molecule reads — kb-scale lengths, unpaired, with the higher,
indel-dominated, homopolymer-dependent error profile of nanopore / HiFi data — and validate
against long-read callers (minimap2 → Clair3 / DeepVariant for SNV/indel; Sniffles2 / cuteSV
for SVs; truvari for SV scoring), per ACCESS report §5.

The payoff is structural: **a single long read spans an SV breakpoint**, so long reads are the
natural complement to the short-read SV work — they stress the SV / breakpoint machinery from
the other side and directly exercise epic #311. Concretely for the SCN track (see
`docs/scn_status.md`), long reads are what would resolve the **transposable-element insertions
and complex rearrangements that short-read Delly under-calls** (TEs live in repeats that short
reads can't span). This doc exists because that realization made long-read support a live
question, and it is genuinely **its own epic**, not a small toggle.

## 2. Current state (what exists, what's missing)

- **A stub toggle:** `config.long_reads: bool` (`rneat/src/gen_reads/utils/config.rs`) — but it
  only means "keep fragments shorter than `read_len` and emit truncated reads." It is a
  fragment-handling switch, **not** a long-read simulator.
- **Missing — read length:** `read_len` is a fixed scalar (default 151). Long reads need a
  read-**length distribution** (broad, right-skewed; ~log-normal), not a constant.
- **Missing — error model:** the seq-error model is Illumina per-position quality-score based
  (`common/src/models/sequencing_error_model.rs`). Long reads need an indel-dominated,
  homopolymer-length-dependent, position-along-read error/quality model.
- ACCESS report §5 already states the prerequisite: *"a long-read mode in rneat (read-length
  distribution + long-read error model) … this is feature work gated ahead of the validation."*

## 3. Architecture decision: a sequencing *profile* inside gen-reads — NOT a new subprocess

`gen-reads` splits into two halves with very different sharing:

- **The "what" — genome + variation + truth (SHARED; rneat's crown jewel):** reference reading,
  mutation / variant placement, the SV/CNV/BND machinery, ploidy, coverage targeting, and the
  **golden VCF/BAM truth set**. For long-read SV validation this must be **identical** to the
  short-read path — the entire point is benchmarking short- vs long-read callers against the
  **same injected truth**. Duplicating it in a separate subprocess would be wasteful and would
  drift out of sync.
- **The "how" — turning the mutated genome into reads (DIVERGES by platform):** read
  length + pairing (short: insert-size Normal + fixed `read_len` + paired; long: length
  distribution, single-molecule, unpaired) and the error/quality model.

**Decision:** keep **one** gen-reads engine and refactor the read-generation + error layers
behind a `SequencingProfile` strategy (`IlluminaShort | LongRead{ONT|HiFi}`), selected by a
platform toggle (the existing `long_reads` bool is the seed of this). Everything upstream
(variants / SV / truth) is untouched and shared.

A fully separate `gen-long-reads` **subcommand** is the wrong shape: it would either duplicate
the variant/truth engine (bad) or force that engine into a shared library — which *is* the
profile refactor, just with a second CLI entry point. A thin CLI alias that preselects the
long-read profile is acceptable, but the **engine stays one**.

## 4. Epic components (the real work, beyond the toggle)

1. **Read-length distribution model** — new. Fit a length distribution from a real long-read
   BAM; new/extended builder (e.g. generalize `gen-frag-length-model` into a read-structure
   model, or a `gen-read-length-model`) replacing the fixed `read_len`. Serialized like the
   other models.
2. **Long-read error model** — new. Indel-dominated, homopolymer-length-aware,
   position-along-read error + quality. New/extended builder trained on a real long-read BAM
   (CIGAR-based error-profile extraction). This is the largest piece.
3. **Read-generation layer refactor** — the `SequencingProfile` abstraction; single-molecule
   (unpaired) emission; **long reads spanning SV breakpoints** — a read crossing a DEL/INV/DUP
   junction must render the rearranged sequence correctly. This is where epic #311's
   breakpoint-realism work pays off directly, and it is exactly the TE / complex-SV realism
   short reads cannot reach.
4. **Validation harness** — `longread_pipeline` mirroring the short-read harnesses: minimap2
   alignment → Clair3 / DeepVariant (SNV/indel) + Sniffles2 / cuteSV (SV) → truvari scoring.
   Specified in ACCESS report §5; gated behind components 1–3.

## 5. Data prerequisite — and it exists

Components 1–2 require a **real long-read BAM** to train against. That data exists for our SCN
strains: the MM26 assembly (`GCA_040805935.1`, `UIUC_Hgly_MM26_v1`, BioProject **PRJNA852516**)
reports `sequencing_tech: "PacBio Sequel; Oxford Nanopore GridION; Illumina"`; PA3 is expected
to be analogous. So the SCN track double-serves: it both **needs** long reads (for TE/SV
realism) and **provides** the training data (PacBio + ONT reads on NCBI). Resolve the run
accessions via the ENA portal API (entrez-direct SSL-fails from the Delta login node):
`filereport?accession=PRJNA852516&result=read_run&fields=run_accession,instrument_platform,…`.

This makes finding the data a prerequisite that unblocks **both** near-term long-read SV
*calling* on real data and the future model *training* — the same reads serve both.

## 6. Relationship to other epics

- **#311 (breakpoint realism):** component 3 (reads across SV junctions) is shared work — do it
  once, both short- and long-read paths benefit.
- **#319 (long-read validation):** this doc is the feature-work scope that #319's validation is
  gated behind; component 4 *is* #319's harness.
- **Polyploidy / allele-dosage:** orthogonal, but the `SequencingProfile` refactor should not
  entangle ploidy — keep the "what" (which includes ploidy) cleanly separated from the "how".

## 7. Open questions

- Fit the length distribution and error model **per platform** (ONT vs HiFi differ a lot) — one
  profile enum with platform-specific parameters, or two profiles? Lean: one `LongRead` profile
  parameterized by a fitted model, so the platform difference lives in the *data*, not the code.
- How much of the Illumina error model (`QualityScoreModel`) is reusable vs. needs a parallel
  long-read type? Determines whether component 2 extends or replaces.
- Throughput: long-read runs are far fewer, much longer reads — check the per-read CPU cost and
  the chunking assumptions (`docs/subcontig_chunking_plan.md` already flags long-read mode as a
  stress case).

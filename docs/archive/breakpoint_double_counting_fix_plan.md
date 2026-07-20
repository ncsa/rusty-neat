# Fix plan: breakpoint double-counting at homozygous BND/INV junctions

> **✅ IMPLEMENTED & SHIPPED (v1.14.1, #236) — ARCHIVED.** The fix in this plan —
> per-fragment junction-crossing suppression (`collect_suppressible_junctions` /
> `suppress_junction_double_count`) — is merged and released. Retained as a design/decision
> record, not active work.

Status: **implemented** (v1.14.1, #236; plan drafted 2026-06-06).

## Problem (grounded in the code)

gen-reads produces SV junction reads in a **separate pass** from the regular
per-contig read generation, and the two passes don't coordinate:

- **Regular pass** (`process_contig` → `generate_fragments` → `write_block_fastq`):
  generates reads across the contig from the **unbroken forward reference**.
  Coverage is modulated per-base by `build_coverage_multipliers`
  (`runner.rs:2097`). For BND and INV, `coverage_multiplier_for`
  (`runner.rs:2205`) returns **`Some(1.0)`** — i.e. coverage-neutral, full
  reference coverage straight through the breakpoint.
- **Chimeric pass** (`process_chimeric_variants`, `runner.rs:838`): for each
  BND/INV junction it generates `num_frags = scale_coverage(coverage, mult)`
  junction-spanning fragments, where `mult = 1.0` (homozygous) or `1/ploidy`
  (heterozygous). INV emits this per junction × 2 junctions (`location` and
  `end`); BND emits it at `location`.

So at a junction `j`, the number of read-pairs crossing it is:
`regular (~coverage) + chimeric (coverage × broken_fraction)`.
- **Homozygous** (broken_fraction = 1.0): ~**2× coverage** at the junction.
- **Heterozygous** (1/ploidy = 0.5 at ploidy 2): ~1.5×.

The regular reads at a homozygous junction are also *wrong* (they transcribe the
unbroken reference, which doesn't exist on either allele), not just redundant.

## Why the existing coverage-multiplier path can't fix it

Coverage multipliers operate on **per-base coverage targets** over a span. But
the double-count is a property of **fragment geometry** — only read-pairs whose
sequenced read *crosses* the junction are doubled; pairs lying entirely on one
flank (within ±read_len of the junction) are legitimate on the broken allele and
must be kept. A per-base multiplier over a ±read_len window would suppress those
flank pairs too, punching coverage holes on both sides of every junction. The
fix must act at **fragment granularity**, not per-base. (This rules out the
"just add a reduction window to `build_coverage_multipliers`" approach.)

## Recommended design — per-fragment junction-crossing suppression

Model: *the chimeric pass owns `broken_fraction` of the read-pairs that cross
each junction.* So in the regular pass, **drop each junction-crossing pair with
probability `broken_fraction`**. The chimeric pass already adds
`coverage × broken_fraction` stitched pairs back. Net at junction `j`:
- homozygous: 0 regular + `coverage` chimeric = `coverage` ✓
- heterozygous: `0.5·coverage` regular (intact allele) + `0.5·coverage` chimeric = `coverage` ✓

### 1. Collect junctions (once per contig)
From `mutated_map.sv_records`, build a sorted `Vec<(usize /*pos_0based*/, f64 /*broken_fraction*/)>`:
- **BND**: one entry at `sv_rec.location`.
- **INV**: two entries, at `sv_rec.location` and the 0-based `end`
  (**coordinate care**: the chimeric INV branch derives `end` from `sv.end`
  (1-based END) or `location + span − 1`; the suppression position must be in
  the same 0-based block-relative space as fragment `(start,end)` — verify
  against `generate_inv_pair`'s junction coords so the two passes reference the
  identical base).
- `broken_fraction = 1.0` (Homozygous) or `1.0 / ploidy` (Heterozygous) — the
  exact `mult` the chimeric pass uses, so the books balance.

Size: O(#SVs) — kilobytes. Built from data already in memory.

### 2. Suppression predicate (per fragment)
A fragment `(start, end)` has a sequenced read crossing junction `j` when
`start ≤ j < start + read_len` (R1) **or**, paired-end,
`end − read_len ≤ j < end` (R2). (Junctions in the unsequenced insert gap don't
cross a read and aren't doubled — must NOT be suppressed; this matches where
`balanced_chimeric_offset` places the chimeric junction, i.e. inside a read.)
Use `partition_point` on the sorted junction list to find candidates in each
read window — O(log J + hits) per fragment.

If a crossing junction is found, draw one coin from the contig RNG and drop the
pair with probability `broken_fraction` (homozygous ⇒ always drop, can skip the
draw to save a little RNG churn).

### 3. Hook point
In `process_contig`, after `block_fragments` is assembled and **before**
`write_block_fastq` — filter `block_fragments` in place. Keeps the change in one
place, downstream (FASTQ/BAM/AD-counter) automatically sees the corrected set.

### 4. Determinism
The het coin-flip consumes the contig's `child_rng`, shifting its stream → output
changes (expected for a behavior fix). Keep the draw ordered deterministically
(iterate `block_fragments` in existing order; one draw per crossing fragment).
**Golden/determinism/parity baselines must be regenerated.** Homozygous-only
runs change too (fragments dropped), so baselines change regardless.

## Resource analysis (packages / memory / CPU)

| | Recommended (per-fragment suppression) | Coverage-window (rejected) | Materialize alt alleles (ideal, expensive) |
|---|---|---|---|
| **External pkgs** | **None** — std `sort`/`partition_point` + existing `NeatRng`/`rayon`. (#junctions is small; a `rust-lapper`/interval-tree dep is unnecessary — a sorted `Vec` + binary search is faster and dependency-free.) | none | none, but needs a segmented/"patched-reference" view to stay sane |
| **Memory** | O(#junctions) per contig (~KB). No large allocations. | O(#junctions) | **O(contig_len × ploidy)** if fully materialized — ~500 MB for a 250 Mb chr at ploidy 2 (1 byte/Nucleotide); whole-arm SVs make alt arms huge. Needs lazy/segmented representation to avoid. |
| **CPU** | O(fragments × log J) added to read-gen — **<1% overhead** (read-gen is already O(fragments × read_len)). | same, but incorrect | comparable to current read-gen, but replaces the chimeric special-case entirely |
| **Correctness** | Exact balance at fragment granularity | over-suppresses flanks | exact by construction (junctions emerge naturally) |
| **Risk/scope** | Localized, ~1 function + junction collection | n/a | large rewrite of the core variant-application + read-gen path |

→ **Recommended: per-fragment suppression.** No new dependency, negligible
memory/CPU, correct. The materialize-alt-alleles approach is the "right"
long-term architecture (it also deletes the entire chimeric special-case path
and the double-count can't exist), but its memory cost and blast radius make it a
separate redesign, not this fix.

## Scope: BND/INV + DEL/CNV-loss suppressed; DUP/CNV-gain intentionally not

**Implemented** (`collect_suppressible_junctions` / `suppress_junction_double_count`):

- **BND** — one point junction at POS, `broken_fraction` = genotype mult (1.0 / 1/ploidy).
- **INV** — two junctions (POS and END), genotype mult.
- **DEL** — one junction at POS (the anchor), genotype mult. A read spanning the
  whole deletion also covers POS, so the single POS junction suffices; the deleted
  interior is already coverage-zeroed by the multiplier, and suppression handles
  the flank reads crossing POS that leak on top. The two mechanisms target
  different reads, so they compose (verified by the full suite passing).
- **CNV-loss** (cn < ploidy) — DEL-like: one POS junction at `(ploidy − cn)/ploidy`
  (matches the chimeric CNV branch's mult exactly; CN=0 → 1.0).

**Intentionally NOT suppressed:**

- **DUP** and **CNV-gain** (cn ≥ ploidy). Their chimeric junction is the *tandem
  boundary* (END→POS) — a **novel adjacency** that the regular linear pass never
  reproduces, so it is not double-counted. (The extra-copy *depth* is handled by
  the coverage multiplier; a separate, subtler boundary-coverage question for
  tandem dups exists but is out of scope for the chimeric double-count fix.)
- **INS** — literal insertion, no junction double-count.

## Test plan

- **Coverage-depth regression** (the headline): a fixture with a homozygous BND
  and a homozygous INV; align the FASTQ; assert mean depth in a ±read_len window
  around each junction is ≈ `coverage`, not ≈ `2× coverage` (was the symptom).
  Add a heterozygous case asserting ≈ `coverage` (currently ~1.5×).
- **No flank holes**: assert depth on the flanks immediately outside the
  junction window stays ≈ `coverage` (guards against over-suppression).
- **Caller still recovers junctions**: re-run `cancer_sv_benchmark.sh`; BND/INV
  positional recall must not regress (the chimeric junction reads are unchanged;
  we only removed redundant reference reads).
- **Determinism**: same seed → identical output; regenerate the model/parity
  baselines in `tests/` and `scripts/regenerate_model_baselines.sh`.
- **Unit**: the junction-crossing predicate (R1/R2 windows, gap exclusion,
  coordinate space) and the suppression-fraction-per-genotype mapping.

## Risks / follow-ups

- **Coordinate mismatch** between the suppression junction positions and the
  chimeric pass's junction positions would unbalance the books — pin them to the
  exact same expressions (especially INV `end`) and assert in a unit test.
- **DEL/DUP/CNV interaction** with existing coverage multipliers — compose, test.
- **BND mate pairing**: each BND record suppresses at its own contig position;
  the mate is handled when its contig/record is processed. Confirm both ends get
  suppressed (no reliance on the chimeric canonical-ID dedup, which only
  dedupes chimeric emission, not the regular pass).
- The deeper "materialize alt alleles" redesign remains the eventual clean
  architecture; this fix is the low-risk interim that removes the symptom.

Relates to [[project_gen_reads_write_perf]] (the read-gen path this touches) and
the v1.12.0/v1.13.x chimeric work (`docs/cancer_simulator.md` known-limitations).
</content>

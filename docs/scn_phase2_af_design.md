# SCN Phase 2 — Pool allele-frequency reproduction (design)

> **✅ IMPLEMENTED in v1.21.0 (#398 / #399).** The per-variant `allele_fraction` enabling change
> shipped and Phase 2 is validated: on real MM26 Pool-seq, simulated vs real per-site AF reached
> **Pearson r = 0.983 at 150× coverage** (harness `scripts/delta/run_scn_af_validation.sh` +
> `scn_af_compare.py`). This document is retained as the design/rationale record; read it as such,
> not as pending work.

Status: **design draft** (2026-07-11; **code anchors re-verified 2026-07-17**; **implemented
2026-07-18, v1.21.0**). Track: realism epic #311, SCN data from João Gomes Viana. Phase 1 (build
real per-strain models from SCN Pool-seq data) is **done and validated**; see `docs/` history and
the memory note. This is the design for Phase 2.

> **2026-07-17 re-trace (supersedes the stale line numbers below).** The alt-fraction mechanism
> was traced end-to-end against current `develop`. Key corrections, which make this change
> **small and well-contained**:
> - **The primary lever is the per-read decision in `generate_read`
>   (`common/src/file_tools/fastq_tools.rs:462`):**
>   `if (variant.genotype == Homozygous) || (rng.random()? < 0.5) { …alt; entry.1 += 1 } else
>   { …ref; entry.0 += 1 }`. This one branch decides whether the read carries alt or ref **and**
>   increments the `AdCounter` — so it drives both the emitted reads and the measured AF. Input-VCF
>   variants reach it via `variant_map` (`flagged_positions.contains(pos)` → `variant_map[pos]`).
> - There is a **second** het flip in `MutatedMap::mutate_position`
>   (`common/src/structs/mutated_map.rs:106`), but it is reached only by `get_stitched_sequence`
>   → `get_mutated_subseq` (`runner.rs:2183`) for reads that **span an SV junction**, and it does
>   **not** touch the `AdCounter`. Honoring AF there is a secondary correctness item (minority of
>   reads, not AD-counted), not the main path.
> - The `genotype_fraction` closure at `runner.rs:2368` is **not** a lever — it only feeds
>   SV-junction double-count suppression. The coverage-multiplier path is SV-only.
> - Both flips are hardcoded **0.5**, not `1/ploidy` — the SNV read path ignores `ploidy`
>   entirely (they coincide only at diploid). So the reachable fractions today are `{0.5, 1.0}`.
> - **Open question RESOLVED:** input-VCF variants flow through the main path.
>   `read_vcf` (`runner.rs:145`, uses `Variant::from_file` — the full parser with INFO+FORMAT) →
>   `input_variants` → `generate_mutated_map` → literal variants land in `variant_map` →
>   `fastq_tools.rs:462`. One primary fraction-source to change (plus the secondary SV-junction one).
> - **Measuring simulated AF is already free:** the golden VCF emits `GT:AD:DP:AF` with
>   AF = alt_count/depth from the same `AdCounter` the `:462` branch increments (`vcf_tools.rs:184`).
>   Step 3 below (re-align to measure sim AF) is unnecessary — read the golden AF.

## 1. The goal (João's ask, verbatim intent)

> Use the PA3 and/or MM26 assembly to simulate reads based on the parameters observed in
> MM26A and/or MM-BD3, and check if the resulting **allele frequencies** are captured in the
> target variants.

Concretely: given a real Pool-seq sample (150 pooled *H. glycines* females → a **continuous
allele-frequency spectrum** over variant sites), produce simulated reads such that, when you
measure the alt-allele fraction at each site from the simulated BAM, it matches the **real
pool's observed AF** at that site.

### Scope decision (2026-07-11): reproductive replay, NOT a generative population model

This is deliberately **reproductive** — rneat is *fed* the pool's observed per-site AFs and
plays them back faithfully; success = simulated AF ≈ real AF. It does **not** model *why* the
frequencies are what they are (drift, selection, the 150-individual pool sampling, linkage) —
it takes the observed spectrum as input. That is exactly what João's question asks ("are the
resulting frequencies captured in the target variants"), and it makes rneat able to simulate
realistic pool / somatic-AF data, which it cannot today.

**Out of scope (possible future extension):** a *generative* AF model — fit an
allele-frequency distribution / site-frequency-spectrum from the pool and *sample* per-variant
frequencies from it during generation (like the trinucleotide / fragment models). That would
let rneat produce realistic spectra from first principles rather than replaying a fixed list.
It is strictly a superset: it reuses the same per-variant `allele_fraction` primitive below,
so nothing here is wasted if we pursue it later.

## 2. Why this does NOT map onto rneat today (the core constraint)

rneat is a single-sample simulator. Tracing the read path:

- A variant's realized alt fraction in the reads comes from `genotype_fraction`
  (`rneat/src/gen_reads/utils/runner.rs:2368`):
  ```
  Homozygous   => 1.0
  Heterozygous => 1.0 / ploidy
  ```
- `Genotype` (`common/src/structs/variants.rs:27`) is a **binary enum** — `Heterozygous |
  Homozygous`. There is no per-variant continuous fraction.
- `ploidy` is a **single global** `config.ploidy`, not per-variant.
- `generate_genotype` (`common/src/models/mutation_model.rs:206`) emits a per-ploid bitmask,
  but read generation collapses it to Het/Hom.

Consequence: the only alt fractions rneat can emit are **`{1/ploidy, 1.0}`** — two values,
the same for every het variant. A continuous spectrum (e.g. sites at 0.07, 0.13, 0.41, 0.88)
is unreachable. Neither lever rescues this:

- **Global ploidy** is shared by all variants → one het fraction for the whole run.
- **Cancer purity** (`gen-cancer-reads`) is a **global scalar** applied to all somatic
  variants (effective VAF = purity × genotype_fraction) → still one value per class, not a
  per-site spectrum.

So Phase 2 needs a **small, targeted enabling change**, not just a clever invocation.

## 3. Enabling change: optional per-variant allele fraction

Add an optional `allele_fraction: Option<f64>` carried on a variant, honored by the per-read
alt/ref decision when present. Confirmed to be a **small, contained change**:

1. `common/src/structs/variants.rs` — add `allele_fraction: Option<f64>` to `Variant`. This is
   the only churny part: ~32 struct-literal construction sites across 9 files each need
   `allele_fraction: None` added (mostly test fixtures + `compare_vcfs` writers where `None` is
   correct); only `from_file` sets it meaningfully.
2. `Variant::from_file` (`vcf_tools.rs`) — populate `allele_fraction` from the input VCF:
   `INFO/AF` if present, else compute from `FORMAT/AD` = alt_depth / (ref+alt). For pooled
   data use **AD, not GT** (§7). `bcftools +fill-tags -- -t AF` upstream also works.
   (`from_file_lean` only has INFO, no FORMAT/SAMPLE — leave it `None`; it feeds gen-mut-model
   training, not the gen-reads input path.)
3. **Primary:** `generate_read` in `common/src/file_tools/fastq_tools.rs:462` — generalize the
   branch. Today:
   ```
   if (variant.genotype == Homozygous) || (rng.random()? < 0.5) { …alt; entry.1 += 1 }
   else { …ref; entry.0 += 1 }
   ```
   becomes: when `variant.allele_fraction == Some(f)`, use `rng.random()? < f` as the alt
   predicate; otherwise the existing `Homozygous || <0.5`. The fixed 0.5 becomes any
   `f ∈ (0,1]`, per variant — and because this same branch increments the `AdCounter`, the
   golden AF tracks it automatically.
4. **Secondary (correctness):** `MutatedMap::mutate_position` (`mutated_map.rs:106`) — the same
   generalization for the SV-junction-spanning read path. Lower priority (minority of reads, not
   AD-counted); can land with (3) or as a follow-up.
5. Golden VCF output — **no change needed.** AF is already measured from the reads and emitted
   (`vcf_tools.rs:184`, `GT:AD:DP:AF`). The truth set already records the realized fraction.

This is **bounded and general**: "simulate reads matching an input VCF's allele frequencies"
is useful for any population / somatic-AF simulation, not just SCN — so it earns its place as
a real rneat feature (candidate issue under #311), not a one-off hack.

**Determinism note:** the `None` path must remain byte-identical — keep the exact same
`rng.random()?` draw and comparison so no existing fixture shifts. Only an input VCF that
actually carries AF changes output, so `model_parity`/output-fidelity baselines are unaffected
unless a fixture is deliberately given an AF-bearing input VCF.

**Overlap with parked polyploidy epic (#265–#270):** issue #266 ("dosage-aware genotype model +
read sampler") generalizes the *same* `:462` branch from Het/Hom to per-copy dosage. `allele_fraction`
and dosage are complementary levers on one choke point — design the AF field so a future dosage
model composes with it (e.g. dosage sets a default fraction, explicit `allele_fraction` overrides)
rather than the two fighting over the branch. The polyploidy epic is on hold (see memory), so AF
can land first; just leave the branch factored so #266 slots in cleanly.

## 4. Pipeline (per strain: MM26/MM26A, PA3/MM-BD3A)

1. **Extract the target AF spectrum** from the real pool BAM: per-site alt fraction from
   allele depths (`bcftools mpileup -a AD | bcftools call` then AF from AD, or `+fill-tags`).
   Output: `pool.af.vcf.gz` (positions + observed AF). This is both the **input** and the
   **truth** for the comparison.
2. **Simulate** reads from the strain assembly with gen-reads in the new honor-input-AF mode,
   consuming `pool.af.vcf.gz`. Use the Phase-1 built models (seq-error / frag / GC) for
   realistic reads. Coverage ≈ the real pool's, so AF estimation noise matches.
3. **Measure** the simulated AF: **read it straight from the golden VCF's `AF` field** — it is
   already alt_count/depth over the actual simulated reads (no re-alignment). Optionally
   re-align for an independent check, but the golden AF is the primary measurement.
4. **Compare** the golden VCF `AF` (sim) vs `pool.af` (target/input) at the shared sites.

## 5. Validation tool: `scripts/delta/scn_af_compare.py`

Dependency-free, mirroring `sbs96_compare.py`. Inputs: `--truth pool.af.vcf.gz` (input AF) and
`--sim golden.vcf.gz` (read the emitted `AF` field — no separate `sim.af.vcf.gz` step needed).
Reports:
- Pearson/Spearman correlation of per-site AF (truth vs sim).
- RMSE / mean abs error overall and per AF-decile bin (does it hold across rare→common?).
- A text scatter / binned table, and count of sites recovered vs dropped.
Success criterion (proposal): correlation ≥ 0.95 and per-bin MAE within AF-estimation noise
at that coverage. (Analogous to the SBS-96 cosine ≥ 0.9 bar from #372.)

## 6. Phase 2 Step 0 — quantify the gap first (NO code)

Before building anything: run gen-reads **today, zero code** with the pool VCF as `input_vcf`,
then read the golden VCF's `AF` field and run `scn_af_compare.py` against the input AF. Because
the golden VCF already measures realized AF from the reads, the gap is directly observable:
het sites pile at ~0.5, hom at 1.0, versus the real continuous spectrum. This gives the
baseline the enabling change must beat, motivates the feature, and shakes out the
AF-extraction + comparison tooling before any Rust change. (The het coin flip is a fixed 0.5
regardless of `ploidy`, so `ploidy` does not change this baseline for SNVs.)

## 7. Risks / open questions

- **AF from AD vs GT:** pooled genotypes (`GT`) are not diploid-meaningful (this is why
  Phase-1 `homozygous_frequency` differed by pool); use **allele depths (AD)**, not GT, for
  the pool AF.
- **Low-coverage / noisy sites:** filter to a min depth so AF targets aren't dominated by
  sampling noise; the comparison should gate on the same depth.
- **Indels:** AF from AD is cleaner for SNVs; decide whether Phase 2 scopes SNVs first.
- **Determinism / placement:** input variants are placed at fixed positions (not model-
  sampled), so #372 context-weighting does not apply to them — good, keeps AF the only
  variable under test.
- **Reference bias:** each pool → its fitting reference (MM26→MM26, BD3→PA3), consistent
  with Phase 1; a common-reference variant comparison is a separate, orthogonal question.

## 8. Sequenced deliverables

1. AF-extraction + `scn_af_compare.py` (input-AF vs golden-VCF-`AF`) + **Step 0 baseline**
   (no Rust change) — quantify the gap. Runnable today; the golden VCF already emits `AF`.
2. Enabling `allele_fraction` change (issue under #311) — the four edits in §3. The path is
   already traced (§ header note); no separate trace task needed.
3. Re-run pipeline with honor-input-AF; validate correlation ≥ target; A/B MM26 vs PA3/BD3.

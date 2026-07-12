# SCN Phase 2 — Pool allele-frequency reproduction (design)

Status: **design draft** (2026-07-11). Track: realism epic #311, SCN data from João Gomes
Viana. Phase 1 (build real per-strain models from SCN Pool-seq data) is **done and
validated**; see `docs/` history and the memory note. This is the design for Phase 2.

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

Add an optional `allele_fraction: Option<f64>` carried on a variant, honored by read
generation when present:

- `common/src/structs/variants.rs` — add `allele_fraction: Option<f64>` to `Variant`.
- Input-VCF ingestion (the path `gen-cancer-reads` already uses to carry germline into the
  tumor pass) — parse a per-site AF into it. Source of the AF:
  - `INFO/AF` if present, else compute from `FORMAT/AD` (alt_depth / total_depth) — a pooled
    VCF's per-site alt fraction. `bcftools +fill-tags -- -t AF` or a direct AD ratio.
- Read-gen coverage path — where `genotype_fraction` (runner.rs:2368) and the coverage
  multipliers (`coverage_multiplier_for`, runner.rs:2580; `scale_coverage`, 2694) decide how
  many alt-bearing fragments to emit: when `allele_fraction` is `Some(f)`, use `f` directly
  instead of the Het/Hom value. The machinery to scale fragment counts by a fraction already
  exists (it's how Het vs Hom already differ) — this generalizes the fraction from
  `{1/ploidy, 1.0}` to any `f ∈ (0,1]`.
- Golden VCF output — emit the intended AF so the truth set records it.

This is **bounded and general**: "simulate reads matching an input VCF's allele frequencies"
is useful for any population / somatic-AF simulation, not just SCN — so it earns its place as
a real rneat feature (candidate issue under #311), not a one-off hack.

**Scope note / open question:** confirm exactly how the input-VCF path currently instantiates
input variants into the read stream (do they flow through the same coverage/`genotype_fraction`
path as model-generated variants, or a separate one?). That determines whether the change is a
single fraction-source swap or two. First implementation task = trace that path.

## 4. Pipeline (per strain: MM26/MM26A, PA3/MM-BD3A)

1. **Extract the target AF spectrum** from the real pool BAM: per-site alt fraction from
   allele depths (`bcftools mpileup -a AD | bcftools call` then AF from AD, or `+fill-tags`).
   Output: `pool.af.vcf.gz` (positions + observed AF). This is both the **input** and the
   **truth** for the comparison.
2. **Simulate** reads from the strain assembly with gen-reads in the new honor-input-AF mode,
   consuming `pool.af.vcf.gz`. Use the Phase-1 built models (seq-error / frag / GC) for
   realistic reads. Coverage ≈ the real pool's, so AF estimation noise matches.
3. **Measure** the simulated AF: align the simulated reads back to the assembly (or use the
   golden BAM), compute per-site alt fraction the same way as step 1 → `sim.af.vcf.gz`.
4. **Compare** `sim.af` vs `pool.af` at the shared sites.

## 5. Validation tool: `scripts/delta/scn_af_compare.py`

Dependency-free, mirroring `sbs96_compare.py`. Inputs: `--truth pool.af.vcf.gz --sim
sim.af.vcf.gz`. Reports:
- Pearson/Spearman correlation of per-site AF (truth vs sim).
- RMSE / mean abs error overall and per AF-decile bin (does it hold across rare→common?).
- A text scatter / binned table, and count of sites recovered vs dropped.
Success criterion (proposal): correlation ≥ 0.95 and per-bin MAE within AF-estimation noise
at that coverage. (Analogous to the SBS-96 cosine ≥ 0.9 bar from #372.)

## 6. Phase 2 Step 0 — quantify the gap first (NO code)

Before building anything: run gen-reads with the pool VCF as input at `ploidy=2` (so
het→0.5, hom→1.0) and run `scn_af_compare.py`. This measures how badly the current binary
model reproduces a continuous spectrum — expected to be poor (everything piles at 0.5/1.0) —
and gives the baseline the enabling change must beat. Cheap, motivates the feature, and
shakes out the AF-extraction + comparison tooling before any Rust change.

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

1. AF-extraction + `scn_af_compare.py` + **Step 0 baseline** (no Rust change) — quantify gap.
2. Trace the input-VCF → read-stream path; write the enabling `allele_fraction` change
   (issue under #311).
3. Re-run pipeline with honor-input-AF; validate correlation ≥ target; A/B MM26 vs PA3/BD3.

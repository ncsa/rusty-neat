5/24/2026
=========
## rneat v1.10.1

### Patch release: SV-model refinements from real-data validation
- **`<CNV>` records sampled from a model with empty `cnv_copy_number_distribution` now carry a fallback CN.** Trained models that observed `<CNV>` records but had no `INFO/CN` (gnomAD-SV is the canonical case — per-sample CN lives in FORMAT) used to produce `<CNV>` outputs with `copy_number = None`, which `gen-reads` warn-and-skipped for coverage modulation, so every CNV in the output VCF was uninterpretable for downstream callers. `SvModel::sample_variants` now falls back to a hard-coded `FALLBACK_CN_DISTRIBUTION` (mirroring the bundled-default CN shape) when the trained distribution is empty but `Cnv` is in the type list. An `info!` log at sampler entry flags when the fallback is firing.
- **`gen-mut-model` now ingests whole-genome VCFs without OOM.** `read_vcf` loaded every `Variant` with its full `info`, `id`, `filter`, `genotype_str`, `format`, `sample`, and `quality_score` fields populated — none of which the trainer reads. On a full-genome gnomAD-SV ingest (~1.5M records after the standard INV/CPX/CTX/BND filter), the `info` column alone (population AFs, VEP annotations, etc.) carried ~2–3 GB of unused string data and OOM'd the trainer with the desktop session as collateral damage. Added `read_vcf_lean` / `Variant::from_file_lean` which leave the unused fields empty for callers that never touch them. Peak RSS on the full-genome run dropped from 3.6 GB (OOM'd) to 1.3 GB; the full-fat `read_vcf` stays unchanged for `gen-reads` and `compare-vcfs`, both of which still need the round-trip fields. Also demoted the per-record `Found genotype in vcf file` log from `debug` to `trace` — under `--log-level debug` it was emitting one line per VCF record (44 MB on the OOM'd attempt).
- **Bundled `sv_model_defaults` refit from the full-genome gnomAD-SV v4.1 validation run.** The lean reader unblocked the `full`-mode validation the v1.10.1 PR had originally deferred. The fitted values supersede both the v1.10.0 literature defaults and the interim chr22-only fits committed earlier in this release cycle. Three of the five parameters are now genome-wide MLEs: `per_base_rate = 4.841e-4`, type proportions DEL/DUP/CNV = 0.8171/0.1824/0.0005, and per-type log-normal lengths (DEL median 703 bp, DUP median 864 bp, CNV median 15.3 kb). The `per_base_rate` change from 8.9e-6 (chr22 interim) to 4.8e-4 (full genome) is **not** sampling drift — the trainer normalizes by `total_reflen_seen` (the full FASTA), so the chr22 interim default effectively reported "fraction of all-genome SVs that happen to fall on chr22," not a true per-base rate. The full-genome refit corrects that. `homozygous_frequency` and `cnv_copy_number_distribution` stay at literature-derived values: gnomAD-SV is sites-only with no GT, and per-sample CN lives in FORMAT — neither parameter is fittable from this corpus. `tools/validate_with_gnomad.sh full` reproduces the run that generated these numbers.

### Behavioral change — read before upgrading
- **Default SV generation rate increased ~54× relative to v1.10.0.** The refit above raises `MutationModel::default()`'s `per_base_rate` from `3.0e-6` (v1.10.0 literature default) to `4.841e-4` (v1.10.1 full-genome MLE). Anyone with a saved gen-reads YAML that sets `sv_rate_scale > 0` *and* doesn't supply an explicit `mutation_model:` will see roughly **54× more `<DEL>` / `<DUP>` / `<CNV>` records per base** after upgrading. The previous default was meaningfully below the true gnomAD-SV rate and produced unrealistically sparse SV output — the new default is correct for genome-wide simulation, but the change is large enough that downstream pipelines tuned against v1.10.0 numbers may need recalibration.
- **Migration paths:** (a) leave `sv_rate_scale` unset / `0.0` to disable de novo SV generation entirely (this remains the unchanged default and is unaffected); (b) scale your existing `sv_rate_scale` down by ~50× to approximate v1.10.0's output volume; or (c) pin the v1.10.0 defaults by passing an explicit `mutation_model:` file trained against your own corpus.
- `<CNV>` records also drop from ~10% of all generated SVs (v1.10.0 hand-rolled) to ~0.05% (full-genome MLE) — anyone exercising CNV-specific downstream paths will see ~200× fewer of them. If you need denser CNVs for testing, fit a model from a CNV-enriched VCF and override the default.

5/24/2026
=========
## rneat v1.10.0

### De novo CNV / SV generation from a learned model (#105)
- **New `SvModel` on `MutationModel`.** Optional field carrying five learned distributions: `per_base_rate` (Poisson λ per base), `type_probabilities` over `<DEL>` / `<DUP>` / `<CNV>`, per-type log-normal `length_log_normal` parameters, `cnv_copy_number_distribution`, and `homozygous_frequency`. Marked `#[serde(default)]` so v1.9 mutation-model JSON files load unchanged with `sv_model = None`.
- **`gen-mut-model` fits the SV component.** Symbolic ALTs in the training VCF feed a parallel accumulator; after the SNP/indel pass completes, `SvModel::fit_from_observations` returns `Some(SvModel)` if any type cleared the 2-observation per-type fit bar, else `None`. Filter chain (each stage logs a count): non-symbolic ALTs ignored → `<INV>` / breakends / unknown tags dropped (v1 only generates `<DEL>` / `<DUP>` / `<CNV>`) → records with no derivable span dropped → spans below 50bp dropped (they're indels, not SVs) → types with < 2 surviving samples dropped from the model.
- **`gen-reads` samples and emits SVs.** Per-contig Poisson draw → weighted type pick → log-normal length (rejected outside `[MIN_SV_LENGTH_BP, contig_len / 4]`) → uniform anchor with overlap and N-gap rejection → `INFO/CN` draw for `<CNV>` → Bernoulli for genotype using the trained homozygous fraction. De novo records merge with `input_vcf` SVs and flow through the existing depth-modulation path, so simulated coverage and the golden VCF round-trip behave identically for both sources.
- **`sv_rate_scale` config knob, default `0.0`.** De novo SV generation is **opt-in** — v1.9 users see no behavior change. `1.0` reproduces the rate of the trained corpus; larger values scale proportionally for stress testing.
- **Bundled default `SvModel`.** `MutationModel::default()` (the embedded "no model file supplied" path) now attaches a hand-rolled SV model with approximate gnomAD-SV v2.1 parameters (DEL/DUP/CNV proportions, per-type log-normal length, CN distribution, homozygous fraction). Parameter provenance is documented inline in `common/src/models/sv_model_defaults.rs`. Users doing serious benchmarking should retrain via `gen-mut-model` against the SV-rich VCF of their choice; the default exists for tire-kicking.
- **`gen-mut-model` now handles SV-only training corpora.** v1.9.1's `NoLiteralVariants` guard rejected VCFs with zero literal SNPs/indels (e.g. the gnomAD-SV sites-only release) — surfaced during the v1.10 release validation against real gnomAD-SV data. v1.10 instead produces an SV-only model with `mutation_rate = 0` and a populated `sv_model`. gen-reads loaded with this model generates zero SNPs/indels and relies entirely on de novo SVs (via `sv_rate_scale > 0`) or `input_vcf` records. The `NoLiteralVariants` error variant is retained in the enum but no longer fires from this path.
- **Bgzip-compressed VCFs and FASTAs now read correctly.** `read_gzip_lines` previously used `flate2::GzDecoder`, which only reads the first gzip stream and silently stops. Bgzip files — produced by every standard bioinformatics tool (`bgzip`, `tabix`, `bcftools view -Oz`, `samtools faidx`) — are a series of concatenated gzip streams (one per ~64 KB block), so all but the first block were dropped. Practical effect: 100% of records lost from every bgzipped VCF (including gnomAD-SV, dbSNP, anything tabix-indexed). Switched to `MultiGzDecoder`, which is a strict superset — it handles both single-stream plain gzip and multi-stream bgzip identically. Bundled test fixtures (`H1N1_read1.fq.gz`, `H1N1_read2.fq.gz`) were bgzipped, so the `seq_error.canonical.json.gz` baseline was re-blessed (the trainer now sees the full FASTQ instead of just the first block).
- **Sites-only VCFs (no FORMAT / SAMPLE columns) now parse instead of being silently dropped.** Real-world callsets like gnomAD-SV ship without per-sample genotype columns. The reader used to skip every record because GT was absent; it now defaults to heterozygous (`0/1`) for sites-only records and records with FORMAT lacking GT. The trainer's `homozygous_frequency` will report `0` for a sites-only corpus — an honest "we don't know" rather than a synthesized 50/50 split.
- **Default log level lowered from `trace` to `info`.** The `--log-level` default was `trace`, which meant the file logger captured per-base debug events from `gen-reads` (sequencing-error generation, read placement, etc.) — on the order of `coverage × read_len × reference_bp` log lines. A single `gen-reads` run against hg38 at coverage 10 produced a multi-GB `.neat.log` and burned most of its wall-clock on disk writes. Users who actually need verbose output now opt in with `--log-level debug` or `--log-level trace`.

### Known limitations
- The bundled default model is **literature-derived**, not refit on the actual gnomAD VCFs. Users who care about distribution fidelity should retrain.
- `<INV>` and breakends still aren't generated — neither the read-content nor the junction modeling for them has landed.
- The sampler's overlap rejection retries up to `10 × n` times; very high `sv_rate_scale` values may saturate retries and warn-and-stop with fewer SVs than the Poisson draw requested.

5/24/2026
=========
## rneat v1.9.1

### Patch release on top of v1.9.0
- **Cargo workspace version catchup.** v1.9.0 shipped with `Cargo.toml`'s `workspace.version` still at `1.8.1` — anything inspecting `cargo` metadata or the compiled binary's reported version read the old value. Bumped to `1.9.1` so the package metadata matches the release.
- **Hardened `.as_literal().unwrap()` sites against future symbolic-ALT leaks (#169).** Added `debug_assert!(v.alternate.is_literal(), ...)` guards at the four call sites in `gen_mut_model::runner`, `fastq_tools::generate_read`, `compare_vcfs::equivalence::apply_variants`, and `compare_vcfs::runner::VariantKey::from`. Also added an explicit `is_symbolic()` skip at the top of `gen_mut_model::runner`'s variant loop — symbolic records were already falling into the `_` arm of the `variant_type` match, but the loop still incremented `homozygous_count` for them, biasing the trained model.
- **Fixed `SvData::span` silently accepting `INFO/END < POS` (#170).** Used to compute `end.saturating_sub(location).saturating_add(1)`, so a malformed `END` strictly less than `POS` silently produced `Some(1)`. The `<DEL>` path masked this via a downstream empty-range check, but `<DUP>` / `<CNV>` / `<INV>` got a bogus 1-base coverage modulation instead of warn-and-skip. `span()` now returns `None` for `end < location`; `Variant::from_file` also clears the bad `END` at parse time with a `warn!` naming both `POS` and `END`.
- **Fixed `gen-mut-model` NaN-poisoning on training VCFs with no literal variants.** A training VCF containing only symbolic SVs (zero SNPs / insertions / deletions) NaN-divided the SNP-frequency math (`0/0`), which then serialized to JSON as `null` and failed to deserialize on re-read. The runner now returns `GenMutationModelError::NoLiteralVariants` with a clear message and writes no output file. Real training corpora always carry literal variants alongside SVs, so this is an edge-case hardening rather than a fix for an observed user incident.
- **Comment and test-file cleanup.**
    - Scrubbed all remaining iterative-development "Phase N" references in source and test comments (most of them landed alongside #169, with one straggler in `compare_vcfs::runner` cleared up in this release).
    - Renamed `compare_vcfs_phase{1,2,3,4}.rs` integration test files to descriptive names (`compare_vcfs_exact_match.rs`, `compare_vcfs_equivalence.rs`, `compare_vcfs_attribution.rs`, `compare_vcfs_reasons_vcf.rs`) matching the doc comments already inside them.

5/24/2026
=========
## rneat v1.9.0

### Symbolic / structural variant support: input → depth modulation → output round-trip (#105, #106, #107)
- **New variant data model.** `AlternateType` is now an enum — `Literal(Vec<Nucleotide>)` for SNPs / indels / multi-base substitutions and `Symbolic(SvData)` for structural ALTs. `SvData` carries the verbatim ALT string, an `SvType` tag (`Del`, `Dup`, `Cnv`, `Ins`, `Inv`, `Bnd`, `Unknown`), and optional `INFO/END`, `INFO/SVLEN`, `INFO/CN`. Helper methods `as_literal` / `as_symbolic` / `is_literal` / `is_symbolic` let callers dispatch without unwrapping. Both VCF writers (`common::file_tools::vcf_tools::write_vcf` and `compare_vcfs::utils::vcf_writer`) emit `SvData.raw_alt` verbatim so the original symbolic notation survives the round-trip.
- **Input VCF parser handles symbolic ALTs.** `parse_alternate` classifies `<DEL>` / `<DUP>` / `<CNV>` / `<INS>` / `<INS:ME:ALU>` / `<INV>`, breakend notation (`G]17:198982]`, `]13:123456]T`), and other `<TAG>` forms (recorded as `SvType::Unknown` with the raw ALT preserved). `Variant::from_file` populates `SvData` from `INFO/END` / `INFO/SVLEN` / `INFO/CN` via a new `parse_sv_info` helper; `INFO/SVTYPE` disagreement vs the ALT tag warns but trusts the ALT (it's the more canonical signal). Malformed INFO values, bare flags, and `.` / empty INFO all parse without error.
- **Per-region depth modulation in `gen-reads`.** Symbolic SVs from `input_vcf` are no longer dropped; they drive per-region coverage scaling instead. `<DEL>` zeroes out reads in its span (hom → ×0, het → ×(ploidy−1)/ploidy), `<DUP>` inflates (hom → ×2, het → ×(ploidy+1)/ploidy), and any SV with `INFO/CN` is scaled by CN / ploidy. `<CNV>` records without `INFO/CN` warn and pass through unmodulated. `<INS>`, `<INV>`, breakends, and unknown tags pass through at ×1 — round-tripped to the output VCF but not yet biologically simulated (orientation flip / junction reads are out of scope for this release). Overlapping SVs compose multiplicatively. The fragment loop splits each region of interest by multiplier boundaries and calls `generate_fragments` / `generate_weighted_fragments` once per sub-region with scaled coverage.
- **`MutatedMap` separates symbolic SVs from per-base mutation flags.** New `sv_records: Vec<Variant>` lives alongside `variant_map`; `from_interval` routes any `AlternateType::Symbolic` variant there so its position is never flagged for per-base mutation. `write_vcf` emits `sv_records` after `variant_map` and preserves the original `INFO` field verbatim so `SVLEN` / `END` / `CN` survive the round-trip. `mutate_position` falls back to the reference if a symbolic record somehow reaches `variant_map`, so per-base read generation can't panic on `as_literal().unwrap()`.
- **`compare-vcfs` handles symbolic ALTs cleanly.** Previously `filter_vcf` / `VariantKey::from` would panic on `as_literal().unwrap()` if the golden VCF contained a round-tripped `<DEL>` / `<DUP>` / `<CNV>`. They're now counted into a new `skipped.symbolic` bucket (golden / called) and excluded from classification — `compare-vcfs` is byte-wise REF/ALT comparison and has no notion of equivalent symbolic edits. Report **schema version bumped to 1.3.0**; consumers parsing the JSON should default missing `symbolic` to 0.
- **End-to-end coverage.** New `pipeline_e2e.rs` test runs `gen-reads` through the real binary with a homozygous `<DEL>` at `H1N1_HA:500-1500` and asserts (a) the symbolic record appears verbatim in the output golden VCF with `INFO/END` preserved, and (b) zero reads start inside the deleted span while the surrounding region is still covered. Unit tests cover the new helpers (`build_coverage_multipliers`, `coverage_multiplier_for`, `apply_coverage_factor`, `split_region_by_multipliers`, `scale_coverage`) plus the `MutatedMap` routing change and the `compare-vcfs` symbolic-skip path.

### Known limitations
- `<INV>` and breakend records round-trip but do not yet drive orientation flipping or junction-read modeling — reads from inversion / breakend spans still come from the forward-strand reference.
- The mutation-rate map excludes the anchor base of each input variant from de-novo SNP sampling but does not yet exclude full SV spans. A homozygous `<DEL>` produces zero coverage in its span, so generated SNPs inside the deletion are invisible in the FASTQ, but they will appear in the output golden VCF.
- Multi-allelic SV sites (`,`-separated ALTs) are still skipped at the reader, as they are for literal multi-allelic records.
- NEAT does not yet *generate* CNVs / SVs de novo from a model — the v1.9.0 work is round-trip + depth modulation for user-supplied symbolic variants only. The mutation-model integration called out in #105 remains open.

5/23/2026
=========
## rneat v1.8.1

### Internal refactors (no user-visible behavior change)
- **Split `gen-gc-bias-model` `RunConfiguration` into walker params + `GcBiasModelParams` (#148).** `run_from_coverage` now takes only the FASTA-sweep / model-write fields it actually uses (`reference`, `bed_table`, `output_file`, `overwrite_output`, `window_size`, `window_stride`, `min_windows_per_bin`). The standalone CLI keeps `bam_file` and `min_mapq` for its own walk and embeds the model params as `config.model`. The unified `gen-bam-models` runner now constructs `GcBiasModelParams` directly from its `GcBiasSection` — no more dummy `bam_file` clone or hard-coded `min_mapq: 0`.
- **Factor `check_overwrite` into `common::file_tools::file_io` (#149).** Replaces four duplicated "output already exists and `overwrite_output` is false" checks in `gen_bam_models` (frag_length + gc_bias sub-sections), `gen_frag_length_model`, `gen_seq_error_model`, and `gen_gc_bias_model`. Standalone error messages are uniform; each module wraps the returned `io::Error` into its own error type.

5/22/2026
=========
## rneat v1.8.0

### `compare-vcfs`: new subcommand for validating a downstream caller against a NEAT golden VCF
- New `rneat compare-vcfs` subcommand. Classifies every variant in the called VCF as TP, FN, or FP relative to a NEAT-simulated golden VCF; produces `comparison_summary.json`, `comparison_summary.txt`, `FN_with_reasons.vcf`, and (optional) `FP.vcf`. Pure-Rust port of NEAT 2.1's `vcf_compare.py` algorithm plus NEAT 4's FN-attribution layer — no `hap.py` / pysam / Python subprocess required.
- **Exact-match classification.** Per contig, variants are keyed by `(position, ref, alt)`. Intersection → TP; golden-only → FN; called-only → FP. Per-VCF skip counters surface in the report so multi-allelic / hom-ref / non-PASS / outside-target-bed / outside-simulated-contigs records aren't silently lost.
- **Equivalence sweep** (ported from NEAT 2.1, default on). For each FP, take a ±`equivalence_window` bp window of the reference and apply both the FP set and FN set within it. If the resulting byte sequences match (case-insensitive), the two sets are denotation-different spellings of the same edit — every consumed FN is promoted to TP, every consumed FP is dropped. Catches left/right-aligned indels, compound SNPs vs single multi-base substitution, and other denotation differences exact matching would mis-classify. Disable with `fast: true`.
- **NEAT-aware FN attribution** (ported from NEAT 4). Every surviving FN gets one or more reason tags: `outside_simulated_contigs` (reported alone; BED checks are skipped since they presuppose simulation), `outside_mutation_bed`, `outside_target_bed`, `unknown`. Counts roll up into `fn_attribution` in the report; `FN_with_reasons.vcf` annotates each FN with a `NEAT_REASON` INFO tag.
- **Chrom aliases.** Optional two-column TSV (`bed_chrom\treference_chrom`) remaps BED chrom names at load time so users can compare against a reference using one convention (e.g., `chr1`) when their BED uses another (e.g., `1`). Applied to both `mutation_bed` and `target_bed`.
- **Chrom-naming-mismatch warnings.** Surfaced in `comparison_summary.warnings` when any BED chrom (post-alias) is absent from the reference's FASTA contig list. Catches single-typo BEDs (`chr_MP` in an otherwise correct H1N1 BED) that would otherwise silently misattribute every variant on that contig. Distinct from NEAT 4's "zero overlap" semantics, which would swallow the single-typo case.
- Report shape is schema-versioned at `1.2.0`. Outputs section enumerates every artifact path; precision / recall / F1 derived from TP/FN/FP counts (with `null` returned for undefined ratios). New `equivalents_promoted` counter in totals so users can see how many TPs were rescued by the sweep.
- See README → "Comparing a downstream caller's VCF against the golden VCF" for the end-to-end runbook; `template_config/compare_vcfs_template.yml` is the canonical config.

### Model-parity regression suite + release shakeout
- New `rneat/tests/model_parity.rs` integration suite. Builds each of the four data-driven models (`gen-mut-model`, `gen-seq-error-model`, `gen-frag-length-model`, `gen-gc-bias-model`) from fixed inputs and compares the gzipped JSON output against checked-in canonical baselines under `rneat/test_data/baseline_models/`. Also covers `gen-bam-models` GC-bias byte-equality vs the standalone runner (complements the in-tree fragment-length parity test which only covered frag-length).
- Cross-process stability is achieved by canonicalizing the parsed JSON before comparison: object keys are sorted via `serde_json::Map`'s BTreeMap backing, and arrays of tuple-key pairs (the shape serde produces for `HashMap<(K1,K2), V>` fields like `MutationModel::statistical_models`) are sorted by serialized key. Neutralizes Rust's per-process HashMap hasher randomization.
- Baselines are gzipped (~120 KB total). Re-bless with `BLESS_BASELINES=1 cargo test --test model_parity` or via `scripts/regenerate_model_baselines.sh` — only when an intentional algorithmic change is being made.
- Added `scripts/release_shakeout.sh` — full-size hg38 chr22 manual shakeout runbook for pre-publish QA. Not run in CI. Trains each model on real chr22 data, asserts byte-parity between the unified `gen-bam-models` and the standalone runners (one config per output, since the unified `min_mapq` gates both observers), runs `gen-reads` end-to-end with the trained models, and reports peak RSS + elapsed time per step for cross-release comparison.
- Added `scripts/RELEASE_CHECKLIST.md` — tickable gate-list for the publish workflow (lint, full test suite, shakeout, cross-version determinism, tag + GitHub release).
- Cross-version determinism verified for v1.7.0 → v1.8.0 on H1N1 paired-end with a fixed seed: byte-identical FASTQ multiset hashes across binaries built from each branch. The new `compare-vcfs` work doesn't touch any `gen-reads` code paths, so this remains true.

### Minor
- `cargo fmt --check` drift in `gen-bam-models` + `main.rs` (carried in from a merge resolution) is fixed; workspace formatting is once again clean.
- `gen-bam-models` template config / README explicitly call out that the unified `min_mapq` gates *both* observers — surfaced during the v1.7.0 chr22 shakeout when the standalone `gen-frag-length-model` (which hard-codes MAPQ > 10) couldn't be matched with a unified `min_mapq: 20` config.
- Added HPC friendly binaries to the auto workflow for releases.
- 
5/22/2026
=========
## rneat v1.7.0

### `gen-bam-models`: unified single-pass model builder
- New `rneat gen-bam-models` subcommand that walks a single BAM once and builds whichever subset of BAM-derived models the config requests. Currently supports `frag_length` and `gc_bias` sections; both are optional, at least one must be present. Saves a BAM iteration per additional model relative to running the per-tool commands separately.
- Config is two optional sub-sections under one top-level `bam_file` / `min_mapq`:
  ```yaml
  bam_file: /path/to/aligned.bam
  min_mapq: 0
  frag_length:
    output_file: /path/to/frag.json.gz
    min_reads: 2
  gc_bias:
    reference: /path/to/reference.fa
    output_file: /path/to/gc.json.gz
    window_size: 100
    window_stride: 100
    min_windows_per_bin: 10
  ```
- `FragLengthObserver` is now self-filtering — it applies its own paired/first-in-pair/mate-mapped/same-ref/MAPQ checks internally regardless of the walker filter, so it remains correct when paired with `BamWalkFilter::for_coverage()` in the unified runner. The standalone `gen-frag-length-model` command continues to pass `for_frag_length()` to the walker; the observer's internal checks are then redundant no-ops, so behavior on that path is unchanged.
- Internal refactor: `gen_frag_length_model::utils::runner::run_from_tlens` and `gen_gc_bias_model::utils::runner::run_from_coverage` are now public so the unified runner can call them with observer outputs from the shared BAM walk. The standalone runners are thin wrappers over the same helpers.
- A parity test asserts the unified runner and the standalone `gen-frag-length-model` produce numerically identical fragment-length models against the same BAM.

### Maintenance / cleanup
- Workspace-wide `cargo fmt` pass. `cargo fmt --check` now exits 0.
- `cargo clippy --workspace --all-targets` now passes with no deny-level errors. Two pre-existing hard failures fixed: `absurd_extreme_comparisons` in `common::rng::range_i64` (`if max > i64::MAX` is always false — replaced with `if min > max` so the function actually validates its range) and `unused_io_amount` in `filter_reads::filter_lib` tests (switched `write(...)` to `write_all(...)`).
- Idiomatic cleanups via `cargo clippy --fix` and a final manual pass: dropped redundant `.clone()` on Copy types, redundant `.into()` on identity conversions, redundant `&` on already-borrowed args; collapsed nested `if X { if Y { ... } }` into combined `if X && Y { ... }`; rewrote `(a + b - 1) / b` → `a.div_ceil(b)`; modernized `repeat(x).take(n)` → `repeat_n(x, n)`, `x >= a && x < b` → `(a..b).contains(&x)`, `x % 2 == 0` → `x.is_multiple_of(2)`; replaced three `assert!(true)` placeholders with `.expect(...)` calls; added a `Default` impl for `NucleotideSelector`; converted `impl Into<usize> for Nucleotide` to the canonical `impl From<Nucleotide> for usize` (same for `TrinucFrame`).
- CLI dispatch arms in `main.rs` (filter-reads, gen-mut-model, gen-seq-error-model, gen-frag-length-model, gen-gc-bias-model, gen-bam-models + the previously-touched gen-reads) now use a single `if let Some((name, cmd)) = subcommand && cmd.contains_id("configuration_yaml")` form instead of two nested ifs. Drive-by: fixed a copy-paste log message in `gen-mut-model` that said "Running rneat filter-reads".
- `gen-reads::run_neat` now takes `&RunConfiguration` instead of `&Box<RunConfiguration>`. Call sites no longer need to `Box::new(config.clone())`.
- `gen-reads`'s test VCF reader uses `.map_while(|l| l.ok())` instead of `.filter_map(|l| l.ok())` so a misbehaving reader can't loop forever on repeated `Err`.
- `gen-reads` config parser accepts `rng_seed` as either a YAML string (`"42 hello"`) or a bare integer (`42`); previously bare integers panicked with `ConfigReadError("rng_seed", "String")`. Two new tests pin both forms.
- `BedRecord` now has an `is_empty()` companion to `len()`.
- `SnpTrinucModel::from_file` no longer takes a dummy `&self` (it never used it).
- `Nucleotide::is_masked` rewritten with `matches!` instead of a boolean-returning `match`.
- Added explicit `#[allow(...)]` annotations for three cases where clippy's suggestion doesn't fit: `should_implement_trait` on the seven `default()` methods that return `Result` (can't implement infallible `std::Default`), `field_reassign_with_default` in the `gen_reads/utils/config.rs` test module (clearer than struct-init for sparse field overrides on a 30-field config), and `same_item_push` on `fastq_tools::generate_read` (pushing `'D'` N times is the CIGAR encoding, not a mistake).

5/21/2026
=========
## rneat v1.6.0

### `gen-gc-bias-model`: BAM input replaces pre-computed coverage files
- `gen-gc-bias-model` now reads an aligned BAM directly. Coverage is accumulated in process per reference base using `samtools depth`-equivalent defaults (skip unmapped/secondary/supplementary, configurable `min_mapq`, default 0). The external `samtools depth` / `bedtools genomecov` preprocessing step is no longer required.
- **Breaking change.** Config fields `coverage_file`, `coverage_format`, and the `CoverageFormat` enum (`samtools-depth`, `bedtools-genomecov-d`, `bedtools-genomecov-dz`) have been removed. Replace `coverage_file` with `bam_file` pointing at the aligned BAM the depth file was originally generated from. `min_mapq` is the new optional knob (default 0; set to 20 to mirror typical `samtools depth -q 20` filtering).
- Memory profile changes: peak memory is now ~4 bytes × total reference bases that received at least one record (~13 GB for hg38 WGS at 30×). The plain-text coverage-file path's per-contig-on-demand behavior is gone; if memory is a constraint on whole-genome runs, split the BAM by chromosome.
- The legacy `CoverageReader`, `CoverageData`, and `coverage_reader.rs` module (459 lines) have been deleted along with the samtools-depth / bedtools-genomecov line parsers.

### Shared BAM-walking infrastructure
- Introduced `common::file_tools::bam_reader::walk_bam` — a single-pass BAM iterator with a configurable `BamWalkFilter` and a `RecordObserver` visitor trait. Multiple observers receive each kept record in one pass, avoiding redundant BAM opens across consumers.
- Concrete observers: `FragLengthObserver` (TLEN histogram, used by `gen-frag-length-model`), `TransitionObserver` (4×4 read-vs-reference mismatch matrix from MD tags, used by `gen-seq-error-model`), and `CoverageObserver` (per-base reference depth, used by `gen-gc-bias-model`).
- `read_fragment_lengths` and `read_bam_transitions` are now thin wrappers over `walk_bam`; their signatures are unchanged so existing call sites continue to compile.
- `BamWalkFilter::for_frag_length()` / `for_transitions()` / `for_coverage()` centralize the per-tool filter policy that was previously duplicated across modules.
- An `on_start(&Header)` hook lets observers cache reference_sequence_id → contig name/length lookup tables once per walk.

### Migration
- In `gen_gc_bias_model.yml`, replace:
  ```yaml
  coverage_file: /path/to/coverage.txt
  coverage_format: samtools-depth
  ```
  with:
  ```yaml
  bam_file: /path/to/aligned.bam
  min_mapq: 20  # optional; default 0 matches `samtools depth` defaults
  ```
  No other config fields changed.

5/20/2026
=========
## rneat v1.5.1
- `gen-reads` now correctly handles reference genomes that contain IUPAC ambiguity codes (R/Y/M/K/S/W/H/B/V/D). Previously these bases were silently mapped to N, producing excess-N reads indistinguishable from assembly gaps. Each IUPAC code is now stochastically resolved to one of its constituent bases at reference-load time using the per-contig simulation RNG, so results are fully reproducible given the same seed. N in the reference retains its existing gap semantics and is unaffected. A single `WARN` line is emitted per contig listing the count of resolved bases. `gen-mut-model` and `gen-gc-bias-model` continue to treat IUPAC codes as N (appropriate since they analyze existing VCF/coverage data where variant callers typically skip ambiguous positions).
- Added `resolve_iupac_bases(raw, rng)` to `common::file_tools::fasta_stream`, available to any future consumer that needs RNG-seeded IUPAC resolution.
- `FastaStream` now yields raw sequence strings rather than `Vec<Nucleotide>`, allowing each caller to choose its own conversion strategy.
- Unit tests: all-code round-trip (output contains only ACGTN), ACGTN passthrough, masked-base passthrough, R two-way uniform distribution (1 000 draws), H three-way uniform distribution (3 000 draws).
- Integration test: `gen_reads_with_iupac_reference_produces_no_iupac_in_output` — runs the full pipeline against a synthetic FASTA with every IUPAC code and asserts no IUPAC letter appears in the output FASTQ.
- Added a binned scoring feature. Activated by user input when forming the quality score model, this will force rneat to bin the quality scores of the input data (if they are already binned, this is trivial) to user-defined bins. When this model is used in gen_reads, it creates a binned output.
- Added a number of unit tests and integration tests.

## rneat v1.5.0

### Memory improvements in `gen-reads`
- Replaced the per-position `Vec<f64>` bias map (one `f64` per chromosome base) with a compact `Vec<(start, end, rate)>` segment list. Measured peak RSS for rice (373 Mbp, 12 chromosomes) dropped from 3.06 GB → 1.05 GB (~3× reduction); major page faults dropped from 109 K → 0. Genomes previously too large to run safely (miscanthus ~11.9 GB estimated) now fit comfortably (~5.2 GB estimated).
- `cover_dataset` (single-ended path): fragment pool reduced from `num_frags` entries (`O(coverage × contig_length)`) to a single-entry pool with a cyclic index. Eliminates up to ~23 MB of VecDeque allocation per large chromosome and millions of redundant RNG shuffle calls.

### `gen-mut-model` fixes for real-world VCF data
- VCF parser no longer panics on missing allele calls (`./1`, `./.`) in the GT field. Missing alleles are skipped during genotype determination; variants with no allele data are treated as homozygous alt.
- Multi-allelic ALT fields (comma-separated, e.g. `T,TAATGGAATGG`) are now skipped with a debug log rather than being stored as garbled Complex variants.
- Variants whose FORMAT or SAMPLE fields are unparseable now warn-and-skip rather than aborting the model build, making the parser robust to mixed-format real-world VCFs.
- Soft-masked FASTA bases (`Maskeda/c/g/t`) are now canonicalized to uppercase (`A/C/G/T`) before comparing against the VCF REF allele in the trinucleotide context pass. This previously caused a hard `BaseMismatch` error on any variant in a repeat-annotated region of a soft-masked reference (e.g. hg38). Genuine mismatches (different base, not just masking) now warn-and-skip instead of aborting.
- Optional `bed_file` key in `gen-mut-model` config now handled safely; previously a missing key triggered a HashMap index panic.
- Per-variant "Found genotype in vcf file" log message demoted from `INFO` to `DEBUG`, eliminating millions of log lines when building models from whole-genome VCFs (NA12878 WGS produced ~4 M `INFO` lines before this fix).
- Human NA12878 whole-genome VCF (hg38) now builds successfully in ~2:15 (wall clock), peak RSS 5.6 GB, yielding `mutation_rate: 0.001457`, `homozygous_frequency: 0.390`, variant distribution 86% SNP / 6.6% INS / 7.4% DEL.

5/9/2026
=========

## rneat v1.4.1
- `gen-reads` new feature: `long_reads` configuration flag (boolean, default `false`). When `true`, fragments sampled from the fragment length model that are shorter than `read_len` are accepted rather than discarded, and the resulting read is truncated to the actual fragment length. This reproduces the realistic read-length distribution seen in real long-read datasets (PacBio HiFi, Oxford Nanopore), where size selection is imperfect and a tail of shorter reads is always present. When `false` (the default, appropriate for Illumina short reads), fragments shorter than `read_len` continue to be rejected — for Illumina, short fragments produce adapter read-through rather than truncated reads, which is outside the scope of this simulator. The flag affects both the uniform-coverage fragment path and the GC-bias-weighted path, as well as the FASTQ and BAM writers, which compute an effective read length of `min(fragment_len, read_len)` per fragment when `long_reads: true`. Added to `template_config/gen_reads_template.yml` with inline documentation.

5/9/2026
=========

## rneat v1.4.0
- `gen-reads` bug fix: `generate_read` in `fastq_tools.rs` was rejecting fragments whose length equals `read_length` (`<=` guard); all single-ended reads were silently discarded. Changed to `<` so fragments of exactly `read_length` are accepted.
- `gen-reads` bug fix: `map_buffer` in `fasta_stream.rs` was treating soft-masked bases (`Maskeda/c/g/t`) as N-regions, excluding repeat-annotated but valid sequence from read generation. Only `N` and `X` now delimit true gap regions.
- `gen-reads` bug fix: `generate_fragments` was drawing fragment lengths from the fragment model for single-ended runs, causing severe under-coverage when the model mean exceeds `read_length`. Single-ended fragments are now fixed at exactly `read_length`; the fragment model applies only to paired-ended runs.
- `gen-reads` bug fix: `cover_dataset` accumulated `start` across loop iterations instead of resetting to `temp_end`, causing reads to pile up at the beginning of each region rather than spanning it evenly.
- `gen-reads` / `mutation_model` bug fix: soft-masked nucleotides (`Maskeda/c/g/t`) were leaking into VCF REF/ALT fields and FASTQ read sequences as lowercase letters (`a/c/g/t`). `generate_mutation` now applies `check_base()` to the extracted `ref_base` and to each base of the deletion reference slice, unmasking all bases before they are written to output.
- New integration tests in `gen-reads`: `test_paired_ended_bam_flags_correct` (SAM flags on every R1/R2 record), `test_paired_ended_insert_size_matches_model` (TLEN mean within ±20% of configured fragment mean), `test_input_vcf_snp_appears_in_bam_reads` (seeded input-VCF SNP visible in BAM reads at the correct position).
- `gen-seq-error-model`: FASTQ input is now streamed one record at a time instead of being fully loaded into memory. Peak memory is bounded by a single record rather than the entire file, which matters for large (multi-GB) FASTQ inputs. The `max_reads` early-exit now also fires without reading the rest of the file.
- `gen-mut-model`: Trinucleotide probability computation is now O(n) instead of O(n²). Transition counts are pre-grouped by reference frame once before the probability loop, eliminating a redundant linear scan per trinucleotide context.
- `filter-reads` (`filter_fastq` and `filter_vcf`): Replaced `format!("{}\n", line)` heap allocations on every written line with two `write_all` calls (bytes + newline). Eliminated double HashMap lookups (`contains_key` + index) in favour of a single `get` call. Contig name reconstruction changed from a manual push loop to `join("_")`. The read-name prefix check simplified from `find` + index comparison to `starts_with`.
- `filter-reads` config: Removed a redundant identity `match` on a bool.
- `gen-frag-length-model`: `filter_lengths` no longer clones its input `Vec` before sorting — the owned value is sorted in place, eliminating one allocation proportional to the fragment count.
- `gen-gc-bias-model`: `overall_mean` computation replaced two separate filtered iterations with a single `fold`, scanning the 101-bin array once instead of twice.
- Streaming FASTA reader (`FastaStream`) moved from `gen-gc-bias-model` to the `common` crate, making it available to all subcommands. `gen-reads` now streams the reference one contig at a time instead of writing temp JSON block files, eliminating up to ~570 MB of disk I/O for human-scale genomes and capping peak memory to a single contig rather than the full genome.
- Contig-level parallelism added to `gen-reads` via rayon `par_bridge`. All contigs are now processed concurrently using rayon's work-stealing thread pool regardless of whether BAM output is requested; the output ordering is preserved by sorting results by contig index before writing.
- BAM creation is now fully parallel. Each contig worker writes its alignment records to a private temp BAM body file (no header, plain BGZF blocks). After all contigs finish, a single concatenation pass strips the BGZF EOF block from each intermediate file and assembles them in reference order into the final coordinate-sorted BAM. No `samtools sort` step is required.
- New `num_threads` config key in `gen-reads`. Set to an integer to cap the rayon thread pool used for parallel contig processing. Omit (or set to `.`) to use all available cores (default). Set to `1` to disable parallelism entirely.
- `NeatRng::derive_child(idx)` added to the common RNG: each parallel contig task derives a deterministic per-contig seed from the parent RNG state and the contig index, so results are fully reproducible regardless of thread scheduling.
- Removed `fasta_reader` and `fasta_map` modules. `SequenceBlock`, `SequenceMap`, and `RegionType` now live in `common::structs::sequence_block`; `map_buffer` and `apply_n_substitution` moved to `common::file_tools::fasta_stream`. Eliminates the temp JSON block files that were written to disk for every contig.
- `gen-mut-model` migrated to a single-pass `FastaStream` loop: trinucleotide counting and variant processing now share one iteration over each contig's sequence instead of two separate passes.

5/3/2026
=========

## rneat v1.3.0
- Added gen_mut_model mimicking the functionality in python Neat. The purpose is to allow users to create custom models based on real data.
- Added `gen-seq-error-model` subcommand. Reads FASTQ (.fastq or .fastq.gz) input and outputs a serialized `SequencingErrorModel` + `QualityScoreModel` (JSON, gzip-compressed). Supports optional BAM file input to infer a SNP transition matrix from real aligned reads using MD tags, and an optional TSV file to supply a fully custom transition matrix. Priority order: TSV override > BAM inference > default matrix.
- BAM-based transition matrix inference reads MD tags from aligned reads (via noodles), walking CIGAR + MD tokens to identify reference/alt base pairs at each SNP position and accumulating per-context counts into the 4×4 substitution matrix.
- Added FASTQ read shuffle in gen-reads: output reads are now randomly shuffled so alignment-tool ordering assumptions are not baked in. Simplified read naming to reduce verbosity.
- Added `template_config/gen_seq_error_model_template.yml` documenting all config fields with inline comments.
- Updated README with a full "Generating a Sequencing Error Model" section covering config schema, BAM MD-tag requirement, TSV format, and priority order.
- Bug fixes: degenerate zero-row handling in transition matrix normalization; edge-position SNP handling in MD-tag walker.
- Added `gen-frag-length-model` subcommand. Reads a paired-end BAM or SAM file and outputs a `FragmentLengthModel::Normal` (gzipped JSON) fitted to the observed template lengths. Filters using a median + 10×MAD outlier ceiling and a configurable per-length `min_reads` count floor (default 2; set 0 to disable). Added `template_config/gen_frag_length_model_template.yml` and README section.
- Implemented `produce_bam` / `output_bam` in `gen-reads` config. Setting `produce_bam: true` now writes a golden BAM file whose reads exactly match the generated FASTQs (same sequence with variants and sequencing errors applied, same quality scores). CIGAR strings reflect both variant positions and sequencing error indels: deletion errors contribute `D` ops and insertion errors contribute `I` ops, so the BAM alignment accurately represents the simulated read relative to the reference. The BAM is written in coordinate-sorted order using a windowed flush strategy (inspired by NEAT 2.1): records are buffered per block, flushed in sorted order at each block boundary, and any cross-block paired-end R2 reads are carried forward to the correct position. No `samtools sort` is needed; run `samtools index` directly on the output. Note: heterozygous variants use the same probabilistic application as the FASTQ, so reads drawn to the reference allele will not show the variant in either file.
- FASTQ shuffle is now global across all contigs rather than per-contig, so no chromosome ordering is visible in the output. All per-block temp files from every contig are collected into memory together and shuffled in a single pass before writing. For large genomes (> 500 Mbp) a runtime warning is emitted recommending `seqkit shuffle` as a memory-bounded post-processing alternative.
- Added `target_bed` config key to `gen-reads`. When supplied, reads and variants are generated only over regions in the BED; contigs absent from the BED are skipped entirely (no blocks read, no temp files written). This is the recommended approach for exome or targeted-panel runs on large genomes, where post-hoc `filter-reads` would otherwise process the full genome first. Both `.bed` and `.bed.gz` inputs are accepted.
- Added `input_vcf` config key to `gen-reads`. A single-sample VCF (`.vcf` or `.vcf.gz`) can now be supplied to force specific variants into the simulation. Each record must carry a `GT` field; SNPs and indels are applied at the given position before random variant generation runs on the remaining mutable bases. Complex variants (multi-base REF and ALT) are skipped with a warning. Random variants still fire at `mutation_rate` for positions not covered by the input VCF; set `mutation_rate: 0` to suppress them entirely. Caveats: only the first ALT allele of multi-allelic records is used; REF alleles are not verified against the reference sequence.
- Added `gen-gc-bias-model` subcommand. Reads a reference FASTA and a per-base coverage file (samtools depth, bedtools genomecov -d, or bedtools genomecov -dz), tiles the reference in sliding windows, and accumulates mean coverage per integer GC% bin to produce a `GcBiasModel` (gzipped JSON). Bins with fewer than `min_windows_per_bin` observations receive a neutral weight of 1.0 so sparse GC% values do not distort the model. Designed to scale to large genomes: coverage is indexed once and loaded one contig at a time; the reference is streamed without temp files; both GC counting and coverage summing use O(1)-amortized sliding windows so peak memory is bounded by the longest contig rather than total genome size. Added `template_config/gen_gc_bias_model.yml` documenting all config fields.
- Added GC bias simulation to `gen-reads`. When a `gc_bias_model` is supplied, fragment start positions are sampled directly from a per-position probability distribution derived from GC content, replacing the previous rejection-sampling approach. The model's `window_size` is embedded in the model file so the application-time window always matches the training-time window. Set `gc_bias_normalize_coverage: true` (default) to inflate fragment count to compensate for low-weight regions; set `false` to preserve the raw requested coverage.

9/19/2025
=========

## rneat v1.2.0
- Executable name change to `rneat`. In order to differentiate rusty-neat from the main neat, we are changing the binary name to `rneat` and we have expanded subcommands and made them more robust.
- Added bed filtering functions. Now if you supply a "filter_output" bed file to your gen-reads run, it will produce the original files plus a filtered version showing only reads and variants that overlap with the regions in the bed.
- It's also possible to run this utility separately with `rneat filter-reads`

9/13/2025
=========

## rust-neat v1.1.2
- I missed an intermediate version. Basically, it was converting to the command `neat gen-reads` instead of `generate-reads`, which will allow us to add subcommands, similar to how python NEAT worked. I also fixed some bugs in the fastq generation and paired-ended creation.
- In this release, we are adding a bed reader that can limit the fasta regions of interest. Basically, it will take an input bed and only generate reads over the sections in the bed.

9/11/2025
=========

## rusty-neat v1.1.0
- We have refactored the fastq code. Now a temporary fastq is written to file per block read in by the fasta reader. This should allow us to parallelize the temporary fastqs. Then it's simply a matter of concatenating them together. Works great! I got through single-ended E coli depth 10 in 36 seconds. I got through paired-ended E coli in 78 seconds.
- Caveats: the fastqs are sorted by read start. Included in each name is the exact location where the sequence was drawn from. Now, when you run an alignment, the reads will be named in such a way that will have <chromosome name>_start-pos_length in the name, so you can shuffle the reads and run an alignment, and easily check to see if your alignment is in the correct contig and location.

9/9/2025
=========

## rusty-neat v1.0.0
- It's been a couple years since I had an update to rusty-neat. This new version is version 1.0.0! It is not yet fully functional, as it is missing the BAM creation, though it does work on some test data: H1N1.fa and ecoli.fa. It's very fast for H1N1, but that's not too surprising. It's slow right now on fastq creation. I already know what I need to do to fix that part, I put an [issue](https://github.com/ncsa/rusty-neat/issues/65) up here that I am working on in develop. But for now, rusty-neat functions and can be executed on some arbitrary data. Will update the readme as well. Once I get the fastq generation idea to work (will require changes to fastq_tools and runner scripts), next step will be parallelizing the main contig loop in runner. That should give us some speedup. Also, we went from using way too much memory to using almost no memory. I'm not sure which is better. It would be nice to find a balance. I still think that a database could work. I'm wondering, maybe we index on position. Just use like a json-based database: {"contig": "chr1", "sequence": {"1":A, "2": C, etc}} We'd want something that had fast retrieval for writing. Something to investigate, maybe.
- The code is significantly different now, I read in all the old models from the original NEAT2.0, and stored them as default models. I used a lot of the same algorithms that they used in NEAT 2.0, figuring it was closer to rust than NEAT 3.X.
- rusty-neat writes out a ton of files. It should clean them up afterward, but I'm not sure it does yet, may need to investigate if there is an added command required. But it stores them in temp and linux cleans that up. Basically, we are splitting up the input fasta into smaller chunks (tbd on if there is an optimal size) then writing those to file for later retrival. Need to test various chunk sizes. I am basically doing blocked gzip-type things, but with some twists. If someone really understood bgzip format, it could probably be used directly. I'm not there yet.
- My runtimes running on a WSL2 Ubuntu setup with the latest rust installed: 311 milliseconds for H1N1.fa, 13 minutes for ecoli. It's all the IO during that step. I knew it was ugly as I was writing it but it's a rough draft. Obviously all the repeat code in fastq_tools is not ideal.
- The output files, so far have looked good. I realized I haven't tested paired-ended yet.
- I decided to map out the N-regions then overwrite them with random valid bases. The main reason being that it kept dipping into N regions and trying to use them for snp creation and that isn't going to work. NEAT 2.0 had a data structure that kept track of how long N-regions were, then would skip them if they were longer than a read length or overwrite them if they were short. That might work. Potential room for improvement, but it's fast so far.
- Eliminated bam creation for now. The option is in there, but it warns you that it doesn't work yet. Once fastqs are sorted and maybe after parallelism, it'll be next up.

9/11/2023
=========

## rusty-neat v0.1.0

Moved the code from a branch under NEAT (https://github.com/ncsa/neat) to its own repo https://github.com/ncsa/rusty-neat. This project will be considered a separate, related project. We will use the same underlying design ideas as NEAT, but it will essentially be different and won't be guaranteed in any way to be able to reproduce what the original NEAT did (though we aim to get it as close as possible, the differences in Python and Rust would make byte-by-byte compatibility very difficult). I have implemented an idea I prototyped in a private repo (https://github.com/joshfactorial/basic_neat), of reading in the file in what I hope is a straightforward way. I was looking for stuff like running out of memory and speed. It can read in the `E. coli` fasta file in about 30 seconds on my laptop, which is pretty fast, compared to Python. What's more is it can store the whole thing in memory, due to the efficient way Rust uses memory. Instead of copying the entire data set every time in memory, quickly eating it up, Rust forces you to decide when to pass a read-only version, copy the dataset, or use up the dataset once and lose it forever (all valid choices depending on context). So I'm feeling hopeful that this will work. I am publishing an Alpha called v0.1.0 for people to test, if they feel so inclined. See the README for more details.

## rusty-neat v0.1

Goal: To release a version of NEAT in rust capable of reading an arbitrary fasta file,
and outputting a version of that file with 1% of the bases mutated as SNPs.

9/5/2023
=========
- Initialized Cargo package with hello world default script.
- Added function to select a random integer from a list (for selecting positions)

9/9/2023
=========
- Created a first pass program that mutates a hardcoded DNA string at one random non-N position.
- Added a logging function to print out time-stamped log items both to stdout and to file.
- Need to improve error handling all around, taking advantage of Rusts error handling system.
- Design decision: After toying with some options to create a standardized RNG to use throughout the program, it seems the best option is to use an on-the-fly rng per function for now. This eliminates reproducibility between runs, but should enable parallelism of processing.
- Redesign from the ground up. Trying to use KISS as much as possible, and only introduce complexity as it becomes necessary, or as I learn Rust tricks that may help this endeavor.

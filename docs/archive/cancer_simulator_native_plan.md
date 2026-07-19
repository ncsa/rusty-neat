# Native cancer simulator: `rneat gen-cancer-reads` implementation plan

> **✅ IMPLEMENTED & SHIPPED (v1.15.0, #239) — ARCHIVED.** `rneat gen-cancer-reads` is a
> first-class subcommand, and the same-seed parity test (`rneat/tests/cancer_parity.rs`) landed
> in v1.15.1 — so the retention condition for the `cancer_simulate.sh` wrapper below is met.
> Retained as a design/decision record.

Status: **implemented** (v1.15.0, #239; parity test v1.15.1). `rneat gen-cancer-reads` ships the
design below. Drafted 2026-06-04.

This is the v2 step the cancer-simulator design doc deferred
(`docs/cancer_simulator.md:15-25`): porting the `tools/cancer_simulate.sh`
shell orchestration into a native `rneat` subcommand.

## Decisions locked in

1. **Full native merge** — reimplement the FASTQ tag+concat and the golden-VCF
   origin merge in Rust, dropping the runtime dependency on `bcftools`/`bgzip`/`awk`.
2. **Inherit, don't fix, documented caveats** — junction-read double-counting
   (~2x coverage at homozygous BND/INV junctions) and symbolic-SV `AD/DP/AF`
   placeholders are real-pass `run_neat` behaviors, out of scope here.
3. **Strict tumor/normal API**, but shape internals so N-way subclonal support
   is a later loop swap rather than a rewrite.
4. **Keep `cancer_simulate.sh`** until the native path has parity tests, then
   deprecate.

## What's already first-class in Rust (the building blocks)

| Piece | Location | Why it matters |
|---|---|---|
| `RunConfiguration` | `gen_reads/utils/config.rs:16` | Public-field struct + `Default` — build two in memory, no YAML round-trip |
| `run_neat(&config, &mut rng) -> Result<Vec<PathBuf>>` | `gen_reads/utils/runner.rs:84` | Reusable core entry; returns written paths; already callable twice |
| `NeatRng::new_from_seed(&Vec<String>)` | `common/src/rng/mod.rs:42` | Per-pass seed derivation (append `-normal`/`-tumor`) |
| `INFO/NEAT_PROVENANCE` (`denovo`/`input`) | `common/src/file_tools/vcf_tools.rs:105,169` + `Provenance` enum in `structs/variants.rs` | The origin-merge's input signal already exists natively |
| `MutationModel` + `sv_model` | `common/src/models/mutation_model.rs` | Tumor model loads identically to germline; bundled COSMIC at `tools/cosmic_v104_pancancer_model.json.gz` |
| VCF reader/writer | `common/src/file_tools/vcf_tools.rs` | Read per-pass goldens back as `Vec<Variant>` for the merge |

Genuinely new code is only **FASTQ name-tag+concat** and **VCF origin-merge**.

## Module layout (`rneat/src/gen_cancer_reads/`)

Mirror the existing per-subcommand convention.

```
rneat/src/gen_cancer_reads/
├── mod.rs              # pub fn main(&PathBuf) -> Result<(), GenCancerReadsError>
├── errors.rs           # GenCancerReadsError
└── utils/
    ├── config.rs       # CancerConfig: parse cancer YAML → two RunConfigurations
    ├── runner.rs       # pass1 → pass2 → merge_fastqs → merge_goldens
    ├── fastq_merge.rs   # stream-tag read names (N_/T_), concat (replaces awk)
    └── vcf_merge.rs     # PROVENANCE → ORIGIN set-merge (replaces bcftools isec)
```

### 1. `mod.rs`
`pub fn main(config: &PathBuf) -> Result<(), GenCancerReadsError>`. Load
`CancerConfig::from_yaml_file`; overwrite-guard the three merged-output paths
(copy `gen_reads/mod.rs:24-40`); call `runner::run_cancer`.

### 2. `errors.rs`
`GenCancerReadsError` via `thiserror`: `#[from] GenerateReadsError`,
`#[from] std::io::Error`, `ConfigError(String)`, `MergeError(String)`,
`PurityOutOfRange(f64)`, `PerPassCoverageZero { normal: usize, tumor: usize }`.

### 3. `utils/config.rs` — `CancerConfig`

```rust
pub struct CancerConfig {
    // shared
    pub reference: PathBuf,
    pub output_dir: PathBuf,
    pub output_prefix: String,
    pub total_coverage: usize,
    pub purity: f64,                  // strict (0,1); endpoints -> error (use gen-reads)
    pub read_len: usize,
    pub paired_ended: bool,
    pub fragment_mean: Option<f64>,
    pub fragment_st_dev: Option<f64>,
    pub rng_seed_root: Option<String>,
    // per-pass
    pub normal_model: Option<PathBuf>,
    pub tumor_model: Option<PathBuf>,
    pub normal_mutation_rate: Option<f64>,
    pub tumor_mutation_rate: Option<f64>,
    pub germline_vcf: Option<PathBuf>,   // skip de-novo germline if supplied
    pub sv_rate_scale: f64,              // tumor pass only
    pub keep_per_pass: bool,             // keep normal/tumor files or just merged
}
```

- Parse with the same `serde_yml::from_reader` -> `HashMap<String,Value>` ->
  match-arm loop as `gen_reads/config.rs:149-470` (reuse the `.`-default convention).
- Validation (port `cancer_simulate.sh:162-194`): reference exists;
  `0 < purity < 1` else `PurityOutOfRange`; paired ⇒ both fragment params present;
  per-pass coverage rounds ≥1 each else `PerPassCoverageZero`.
- `fn shared_run_config(&self) -> RunConfiguration` — reference/read_len/paired/
  fragment/output_dir/overwrite/produce_fastq=true/produce_vcf=true.
- `fn normal_pass(&self)` / `fn tumor_pass(&self, germline_vcf: PathBuf)`:

```rust
fn normal_pass(&self) -> RunConfiguration {
    let mut c = RunConfiguration {
        coverage: ((1.0 - self.purity) * self.total_coverage as f64).round() as usize,
        mutation_rate: self.normal_mutation_rate,
        mutation_model: self.normal_model.clone(),
        input_vcf: self.germline_vcf.clone(),
        sv_rate_scale: 0.0,                       // no somatic SVs in normal
        output_filename: format!("{}_normal", self.output_prefix),
        rng_seed: Some(format!("{}-normal", self.seed_root())),
        ..self.shared_run_config()
    };
    RunConfiguration::check_and_log_config(&mut c).unwrap(); // populates output paths
    c
}
// tumor_pass: coverage = purity*total, tumor_model, tumor_mutation_rate,
// sv_rate_scale = self.sv_rate_scale, input_vcf = germline_vcf, _tumor suffix.
```

- `fn seed_root(&self)` -> `rng_seed_root` or `"cancer-simulate"`.

This eliminates the shell `write_config` temp-YAML dance (`cancer_simulate.sh:206`).

### 4. `utils/runner.rs` — `run_cancer`

```rust
pub fn run_cancer(cfg: &CancerConfig) -> Result<(), GenCancerReadsError> {
    let normal = cfg.normal_pass();
    let mut rng_n = NeatRng::new_from_seed(&normal.seed_vec)?;
    run_neat(&normal, &mut rng_n)?;

    let germline_vcf = cfg.germline_vcf.clone()
        .unwrap_or_else(|| normal.output_vcf.clone().unwrap());

    let tumor = cfg.tumor_pass(germline_vcf);
    let mut rng_t = NeatRng::new_from_seed(&tumor.seed_vec)?;
    run_neat(&tumor, &mut rng_t)?;

    merge_fastqs(&normal, &tumor, cfg)?;
    merge_goldens(&normal.output_vcf.clone().unwrap(),
                  &tumor.output_vcf.clone().unwrap(), cfg)?;
    // emit summary block (port cancer_simulate.sh:428-454)
    // if !cfg.keep_per_pass: unlink per-pass FASTQs (keep per-pass VCFs)
    Ok(())
}
```

### 5. `utils/fastq_merge.rs`
`fn merge_fastqs(normal, tumor, cfg)`. For r1 (and r2 if paired): open each
pass's `output_fastq_*` with `MultiGzDecoder` (pattern in `gen_reads/mod.rs`
tests:124), rewrite record line 1 (`NR%4==1`) `@`->`@N_`/`@T_`, re-encode through
the project bgzf/flate2 writer (`common/src/file_tools/fastq_tools.rs`), concat
into `<prefix>_merged_r{1,2}.fastq.gz`. Carry the read-name-collision rationale
comment from `cancer_simulate.sh:280` (collisions MarkDuplicates would drop).

### 6. `utils/vcf_merge.rs`
`fn merge_goldens(normal_vcf, tumor_vcf, cfg)`. Read both via `vcf_tools` reader
-> `Vec<Variant>`. Key by `(contig, pos, ref, alt)`. Partition:

| In normal | In tumor | `NEAT_ORIGIN` |
|---|---|---|
| yes | no | `germline` |
| no | yes (denovo) | `somatic` |
| yes | yes | `shared` |

Set `INFO/NEAT_ORIGIN`, sort by (contig-order, pos), write one
`<prefix>_merged_truth.vcf.gz`. The Rust writer declares SV INFO once, killing
the `cancer_simulate.sh:394-413` bcftools-header-injection workaround.

**Check first:** whether `compare_vcfs/utils/equivalence.rs:42 sweep()` should
drive matching (normalization / left-align) instead of exact-key.

### 7. Wiring (`rneat/src/main.rs`)
- `pub mod gen_cancer_reads;` (~line 16).
- `NeatErrors`: add `GenCancerReads(#[from] GenCancerReadsError)` (`:39`).
- `neat_commands()`: return type `[Command; 8]` -> `[Command; 9]`; add
  `Command::new("gen-cancer-reads")` (`:59`).
- Add `Some(("gen-cancer-reads", _))` match arm (`:199`), copy gen-reads shape.

### 8. Tests (`rneat/tests/`)
- `cancer_reads_e2e.rs` — H1N1, purity 0.5: merged FASTQ exists, record count ≈
  sum of passes; read names carry `N_`/`T_`; merged truth has all three
  `NEAT_ORIGIN` values; determinism (same seed root -> identical output, mirror
  `determinism.rs`).
- `cancer_reads_config.rs` — purity bounds, per-pass-zero, paired-without-fragment.
- Parity vs `cancer_simulate.sh` on a fixed seed (origin counts match) before
  deprecating the script.

### 9. Docs
- Flip `docs/cancer_simulator.md` status to "native subcommand shipped".
- Add a `template_config/` cancer YAML example.
- CHANGELOG entry.

## Effort
~6 new files, ~4 small `main.rs` edits, plus tests. Most of the subcommand is a
thin driver over `run_neat`; new logic is confined to `fastq_merge` + `vcf_merge`.

## Resolved decisions (2026-06-06; tracked in #239)
1. **`vcf_merge` uses exact-key `(contig, pos, ref, alt)` matching**, not
   `compare_vcfs::sweep()`. Both inputs are rneat-golden VCFs with identical
   representation — the tumor pass carries normal's germline verbatim via
   `input_vcf` (tagged `NEAT_PROVENANCE=input`), somatic are `denovo` — so there
   is no cross-tool representation drift requiring `sweep()`'s positional/
   normalized matching. Exact-key is correct and avoids coupling to compare_vcfs.
   Origin rules: tumor `denovo` → `somatic`; tumor `input` → `shared`; normal key
   absent from tumor → `germline`. Revisit only if inputs ever become
   external/called VCFs.
2. **Strict tumor/normal for v1; N-way subclonal deferred.** Matches the
   validated 2-population orchestration; N-way needs shared-mutation tracking +
   an N-way merge + its own validation. `CancerConfig` is shaped around `purity`
   so a future `subclones:` list is an additive, non-breaking extension.

## Orthogonal / out of scope
- Junction-read double-counting and symbolic-SV `AD/DP/AF` placeholders
  (real-pass `run_neat` behaviors).
- #218 PCAWG SV-model refit (what SVs the tumor model emits, not how passes
  are orchestrated).
</content>
</invoke>

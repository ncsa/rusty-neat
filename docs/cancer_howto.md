# Cancer simulation how-to

A practical, copy-paste guide to simulating tumor/normal sequencing data with
`rneat gen-cancer-reads`. For the design rationale, mutation-rate calibration
math, and the SV roadmap, see [`cancer_simulator.md`](cancer_simulator.md).

## What it does

A real tumor biopsy is a **mixture**: some fraction of cells are tumor (carrying
germline + somatic variants), the rest are normal-tissue contamination (germline
only). `rneat gen-cancer-reads` reproduces this in one command by running two
`gen-reads` passes over the same reference and merging them:

```
                      coverage          variants
  normal pass   →   (1 − purity)·C   →  germline only
  tumor  pass   →     purity·C       →  germline (shared) + somatic (de novo)
                                         │
            merge: tag reads N_/T_, concatenate, and
            combine the two truth VCFs with origin tags
                                         ▼
   <prefix>_merged_r1.fastq.gz   ← feed this to your aligner + somatic caller
   <prefix>_merged_truth.vcf.gz  ← score the caller against this
```

The tumor pass consumes the normal pass's golden VCF as its germline, so both
populations share the same patient germline — biologically correct, and what lets
a somatic caller separate germline from somatic. No `bcftools`/`awk` needed; the
merge is native Rust.

## Quick start

Smoke test on the bundled ~14 kb H1N1 reference — finishes in seconds and proves
your build works end-to-end:

```bash
cat > cancer.yml <<'YAML'
reference: rneat/test_data/references/H1N1.fa
output_dir: ./cancer_out
output_prefix: smoketest
total_coverage: 30
purity: 0.6
read_len: 70
paired_ended: true
fragment_mean: 250
fragment_st_dev: 30
rng_seed: demo
overwrite_output: true
YAML

rneat gen-cancer-reads -c cancer.yml
ls cancer_out/
# smoketest_normal.vcf.gz   smoketest_tumor.vcf.gz
# smoketest_normal_r1.fastq.gz  ...  smoketest_merged_r1.fastq.gz  smoketest_merged_r2.fastq.gz
# smoketest_merged_truth.vcf.gz
```

Copy `template_config/gen_cancer_reads_template.yml` as a fully-commented starting
point for your own config.

## Output files

| File | What it is |
|---|---|
| `<prefix>_merged_r1.fastq.gz` (+ `_r2`) | **The deliverable.** `N_`/`T_`-tagged, concatenated reads — the simulated tumor biopsy. Feed this to your aligner. |
| `<prefix>_merged_truth.vcf.gz` | Origin-tagged truth set. Every record carries `INFO/NEAT_ORIGIN ∈ {germline, somatic, shared}`. Score your caller against this. |
| `<prefix>_normal.vcf.gz` | Germline-only truth (the patient's genome). |
| `<prefix>_tumor.vcf.gz` | Germline + somatic truth. |
| `<prefix>_normal_r1.fastq.gz`, `<prefix>_tumor_r1.fastq.gz` | Per-pass reads (set `keep_per_pass: false` to delete after merge). |

The `N_`/`T_` read-name tags matter: without them a normal and a tumor read
sampled at the same coordinate get identical QNAMEs, which Picard MarkDuplicates
silently drops — halving effective coverage at those loci.

## Worked examples

All examples assume `rneat` is on your `PATH` and a reference at
`~/code/data/chr22.fa` (any `.fa`/`.fa.gz` works).

### 1. Realistic WGS tumor/normal with the bundled COSMIC model

70% purity, 60× combined coverage, somatic SNVs/indels from the pan-cancer COSMIC
model at a realistic per-tumor rate (the `1e-5` default):

```bash
cat > brca_wgs.yml <<'YAML'
reference: ~/code/data/chr22.fa
output_dir: ./brca_out
output_prefix: tumor70
total_coverage: 60
purity: 0.7
read_len: 151
paired_ended: true
fragment_mean: 350
fragment_st_dev: 50
tumor_model: tools/cosmic_v104_pancancer_model.json.gz
tumor_mutation_rate: 1e-5
rng_seed: brca-demo
overwrite_output: true
YAML

rneat gen-cancer-reads -c brca_wgs.yml
```

### 2. Tissue-specific model (breast) with structural variants enabled

Bundled **per-tissue models** stratify *both* halves by primary site
(`tools/cosmic_per_tissue_{BRCA,skin,lung}.json.gz`): a per-tissue COSMIC
SNP/indel spectrum (e.g. skin/melanoma is UV-C>T-heavy and SNV-dominated; breast
carries more indels) paired with a per-tissue `sv_model` (BRCA DUP-dominant, skin
BND-enriched, lung DEL-dominant). Turn on de novo SV generation with
`sv_rate_scale`:

```bash
cat > brca_sv.yml <<'YAML'
reference: ~/code/data/chr22.fa
output_dir: ./brca_sv_out
output_prefix: brca_sv
total_coverage: 60
purity: 0.7
read_len: 151
paired_ended: true
fragment_mean: 350
fragment_st_dev: 50
tumor_model: tools/cosmic_per_tissue_BRCA.json.gz
sv_rate_scale: 1.0          # 1.0 = the model's nominal rate; higher stress-tests SV callers
rng_seed: brca-sv-demo
overwrite_output: true
YAML

rneat gen-cancer-reads -c brca_sv.yml
```

### 3. One germline, many tumor scenarios

To compare callers across purities while holding the patient germline fixed,
generate the germline once, then point each run at it via `germline_vcf:`. Any
normal-pass golden VCF (or a real germline VCF) works:

```bash
# Reuse an existing germline truth as the shared germline:
for p in 0.3 0.5 0.8; do
  cat > scenario_$p.yml <<YAML
reference: ~/code/data/chr22.fa
output_dir: ./purity_sweep
output_prefix: purity_$p
total_coverage: 60
purity: $p
read_len: 151
paired_ended: true
fragment_mean: 350
fragment_st_dev: 50
germline_vcf: ./brca_out/tumor70_normal.vcf.gz
tumor_model: tools/cosmic_v104_pancancer_model.json.gz
rng_seed: sweep-$p
overwrite_output: true
YAML
  rneat gen-cancer-reads -c scenario_$p.yml
done
```

### 4. Low-purity sample

Poor-yield biopsies can be <30% tumor. Just lower `purity` — but keep
`total_coverage` high enough that the tumor pass (`purity·C`) doesn't round to a
useless depth:

```bash
# purity 0.2 at 100x → 80x normal / 20x tumor
rneat gen-cancer-reads -c <(sed 's/^purity:.*/purity: 0.2/; s/^total_coverage:.*/total_coverage: 100/' brca_wgs.yml)
```

> **Small / viral references:** the de novo SV sampler caps SV length at
> `contig_len/4` by default. On bacterial or viral contigs that starves the larger
> SV types. Raise it per the structural-variants section of
> `template_config/gen_reads_template.yml` (`sv_max_length_fraction`) — though note
> that knob lives on `gen-reads`; for cancer SV work on small contigs, supply the
> events through a `germline_vcf:`/input VCF instead.

## Benchmarking a somatic caller

The repo ships Docker-based scoring pipelines so you can go from simulated reads
to a scored caller in one command.

**SNV/indel** — BWA-MEM → Mutect2 → `som.py`:

```bash
tools/cancer_benchmark.sh \
    --reference   ~/code/data/chr22.fa \
    --normal-fastq ./brca_out/tumor70_normal_r1.fastq.gz \
    --tumor-fastq  ./brca_out/tumor70_merged_r1.fastq.gz \
    --truth-vcf    ./brca_out/tumor70_merged_truth.vcf.gz \
    --output-dir   ./benchmark_out
```

**Structural variants** — BWA-MEM → Manta → `truvari`:

```bash
tools/cancer_sv_benchmark.sh \
    --reference   ~/code/data/chr22.fa \
    --normal-fastq ./brca_sv_out/brca_sv_normal_r1.fastq.gz \
    --tumor-fastq  ./brca_sv_out/brca_sv_merged_r1.fastq.gz \
    --truth-vcf    ./brca_sv_out/brca_sv_merged_truth.vcf.gz \
    --output-dir   ./sv_benchmark_out
```

Because the truth VCF carries `INFO/NEAT_ORIGIN`, you can filter to somatic-only
records before scoring (the benchmark scripts expose `--truth-filter`), giving
clean somatic TP/FP/FN without germline contaminating the numbers.

## Training your own tumor model

The bundled COSMIC pan-cancer model is a convenience. To fit a model from a corpus
you control, the two shipped adapters chain into `rneat gen-mut-model`:

```bash
# Open-tier TCGA MC3 (no registration; SNPs only):
tools/fetch_tumor_corpus.sh   # → writes a normalized VCF
# Or COSMIC GenomeScreensMutant (academic registration; SNV + indel):
tools/fetch_cosmic_corpus.sh

# then train:
rneat gen-mut-model -c your_gen_mut_model_config.yml
```

Point `tumor_model:` at the resulting `.json.gz`.

## Calibration notes

- **Somatic rate.** `tumor_mutation_rate` defaults to `1e-5` (typical solid
  tumor; `~1e-4` for high-burden/MSI). Corpus-aggregated model rates (e.g. the
  COSMIC model's ~5.5e-3) measure "fraction of bp mutated *anywhere across the
  catalog*" and overstate per-tumor burden by 50–500×, so the de novo somatic rate
  is set here, not taken from the model. Pass `tumor_mutation_rate: model` only if
  you deliberately want the model's fitted rate.
- **Purity drives VAF.** Somatic variant allele fractions scale with `purity`;
  sweep it (example 3) to test a caller's low-VAF sensitivity.
- **Reproducibility.** `rng_seed` seeds both passes (suffixed `-normal`/`-tumor`).
  The same seed + config reproduces a run; the seed used is printed to the log.

## Relationship to `tools/cancer_simulate.sh`

The native subcommand is a drop-in replacement for the original shell
orchestrator. They are verified to produce equivalent output by the parity test
`rneat/tests/cancer_parity.rs` (identical merged-FASTQ record multisets, per-pass
golden VCFs, and origin classifications). The shell script is retained for now as
the reference implementation and for the Docker benchmark wiring; prefer the
native subcommand for new work.

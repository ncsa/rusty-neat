# Cancer simulator design

Status: **planning** (no implementation yet — see #129 for the umbrella issue and #183–#192 for the breakdown).

## Goal

Produce simulated short-read sequencing data with the statistical and structural properties of a real cancer sample, including:

1. A **mixture of tumor + normal DNA** at a configurable purity (real tumor biopsies are typically 30–80% tumor cells; the rest is normal-tissue contamination).
2. **Germline variation shared between the two populations** (the tumor and the normal both came from the same patient).
3. **Tumor-only somatic variation** in the tumor population — elevated SNV / indel rates, cancer-specific mutation signatures, and a richer SV catalog than the gnomAD-SV-derived germline default.

The output must be **benchmarkable** — a downstream somatic variant caller (Mutect2, Strelka, GRIDSS, Manta, …) should be able to run on rneat's simulated reads and have its calls scored against a ground-truth VCF that distinguishes germline from somatic.

## Approach: orchestration first, native subcommand later

The recommended v1 path is a **thin orchestration layer** in `tools/` that drives `rneat gen-reads` twice — once per population — and concatenates the resulting FASTQs. No core rneat changes for the MVP; we lean on the existing primitives (`input_vcf`, `mutation_model`, `sv_rate_scale`, paired-end fragment generation, depth modulation).

Why this approach:

- **Validated by precedent.** NEAT 2.1 shipped this exact workflow as its `genReadsTumorTutorial` (~30-line shell script, archived at `~/code/neat2/models/genReadsTumorTutorial.zip`). It worked then; it'll work now.
- **Minimum new code.** Until the simulation produces credible output end-to-end, every line of new core code is risk. The orchestration path gets us to a working artifact with zero rneat-core diffs.
- **Reversible.** If we later decide a native `rneat gen-cancer-reads` subcommand is worth the engineering, the orchestration script is the design spec for that subcommand's behavior — it codifies the data flow.

A native subcommand is **not** an explicit non-goal — it's a v2+ option pending validation that the orchestration approach produces samples that real somatic callers handle correctly.

**Update (#239): the native subcommand `rneat gen-cancer-reads` is now implemented** — the orchestration validated cleanly (Mutect2 SNV + Manta SV recall), so the two-pass flow is now a first-class subcommand: `rneat gen-cancer-reads -c <yaml>` runs both passes and merges them (tagged FASTQs + an origin-tagged truth VCF), with no `bcftools`/`awk` runtime dependency. Config template: `template_config/gen_cancer_reads_template.yml`; design + decisions: `docs/cancer_simulator_native_plan.md`. `tools/cancer_simulate.sh` is retained until a same-seed parity test against the native path lands.

## Data flow

```
                  ┌─────────────────┐
reference FASTA ──┤                 │
                  │  Pass 1: normal │── normal_reads.fq.gz    ─┐
                  │                 │── normal_golden.vcf.gz   │
                  │  cov = (1−p)·C  │      ↓                   │
                  └─────────────────┘      │ (germline truth)  │
                                           │                   │
                  ┌─────────────────┐      │                   ├──► concat ──► merged.fq.gz
reference FASTA ──┤                 │      │                   │
        normal_golden.vcf.gz   ─────┤  Pass 2: tumor  │── tumor_reads.fq.gz   ─┘
                                    │                 │── tumor_golden.vcf.gz
                                    │  cov = p·C      │      ↓
                                    │  + tumor model  │      │ (germline + somatic)
                                    └─────────────────┘      │
                                                             ↓
                                        normal_golden + tumor_golden ──► merge ──► truth.vcf.gz
                                                                          (origin tags:
                                                                           germline / somatic / shared)
```

Definitions:
- `p` = **purity**, fraction of merged reads coming from the tumor pass. Real tumor purity 0.3–0.8 typical; biopsies of poor-yield tissues can go lower.
- `C` = total per-base coverage of the merged output (e.g. 30× for WGS).
- The germline VCF from pass 1 feeds pass 2 as `input_vcf:`, guaranteeing the tumor cells carry the same germline variants as the normal cells (this is biologically correct — both populations descend from the same zygote).
- Somatic variants are added only in pass 2 via the tumor's `mutation_model:` and `sv_rate_scale:`. They appear in `tumor_golden.vcf.gz` but not in `normal_golden.vcf.gz`.

## Component breakdown

Each component maps to a sub-ticket of #129:

| Component | Ticket | Effort | Notes |
|---|---|---|---|
| Source tumor training corpus | **#183** | Medium | Two adapters shipped: TCGA MC3 MAF (open) and COSMIC GenomeScreensMutant VCF (academic) |
| Orchestration script (`tools/cancer_simulate.sh`) | **#184** | Low | Ports NEAT 2.1's `simulate.sh` idiom |
| Variant-origin annotation in golden VCF | **#185** | Medium | gen-reads emits `INFO/NEAT_PROVENANCE={denovo,input}`; `cancer_simulate.sh` resolves to `INFO/NEAT_ORIGIN={germline,somatic,shared}` via post-merge |
| Train v1 tumor models | **#186** | Low (per tissue) | Run `gen-mut-model` on the chosen corpus. Likely BRCA / SKCM / LUAD |

#184 can develop in parallel with #183 using a hand-rolled or NEAT-2.1-recovered model as a stand-in. #186 blocks on #183 (data) but not on anything else.

## Source-data decisions

### Chosen for v1

Two open-tier sources are supported, each with its own adapter in `tools/`. The two are complementary, not competing — pick by access constraints and indel needs.

| Source | Adapter | License | Records (v1 snapshot) | Notes |
|---|---|---|---|---|
| **TCGA MC3 PUBLIC MAF** | `tools/fetch_tumor_corpus.sh` | Open, no login | ~2.8M PASS SNPs (v0.2.8) | Lower-friction default — no registration. SNPs only; indels are dropped because MAF's `-` allele convention requires reference-anchor lookup that the adapter does not implement (v2 follow-up). |
| **COSMIC GenomeScreensMutant VCF** | `tools/fetch_cosmic_corpus.sh` | Free for academic use; one-time registration at cancer.sanger.ac.uk | ~16.7M unique SNV/INS/DEL records after dedup (v104) | Larger and indel-capable. Source is already a properly-anchored VCF, so indels round-trip into `gen-mut-model` directly. The adapter deduplicates the ~11.8M multi-transcript-annotation duplicates that COSMIC ships and remaps `MT` → `chrM` for UCSC compatibility. |

Both adapters can chain into `rneat gen-mut-model` via `--train`. The output mutation model is structurally identical regardless of which adapter produced the input VCF — "cancer-ness" is in the training data, not the schema.

A pre-trained pan-cancer model fit from the COSMIC adapter against v104 is checked into the repo at `tools/cosmic_v104_pancancer_model.json.gz` (9.2 KB). It's a temporary placement so the cancer simulator can be exercised end-to-end without requiring users to download the ~1 GB COSMIC corpus first. Final placement (bundled / hosted with a fetch helper / user-trains-their-own) is tracked at **#186**.

**A note on the fitted `mutation_rate`.** `gen-mut-model` normalizes variant count by total reference length, so a corpus aggregated across many donors (COSMIC, MC3) will produce a `mutation_rate` that reflects "fraction of bp with any observed mutation across the catalog," not per-tumor burden. The v104 COSMIC fit lands at ~5.5e-3 — high if mistaken for a per-tumor rate. Users simulating a specific tumor purity should treat this as a corpus-level constant and shape the per-sample burden via coverage and purity in the orchestration layer, not by rescaling the model's rate.

### Considered and ruled out (for v1)

| Source | Why not for v1 |
|---|---|
| **Legacy ICGC25K** (NEAT 2.1's source) | The original `dcc.icgc.org` portal was retired June 2024. Open-tier bulk data is reportedly mirrored on AWS S3, but URL stability is not guaranteed and the format (simple-somatic-mutation TSV) requires a custom parser that `gen-mut-model` doesn't have. |
| **PCAWG consensus calls** | Excellent quality, includes SV consensus, but smaller donor count than TCGA aggregate. Better suited to **#191** (chromothripsis / chromoplexy training) than to general somatic-SNV models. |
| **EGA controlled-access** (e.g. CLLE-ES, the Spanish CLL project NEAT 2.1 included) | Controlled tier — not viable as a default for an open-source tool. Users with EGA access can train their own. |

## Mutation model: no new schema needed

NEAT 2.1's three tumor models (BRCA, CLLE, SKCM) used the **same pickled-dict schema** as its germline NA12878 model. The "tumor-ness" came purely from the input training data — the model schema captures whatever statistics the input VCF supplies.

rneat's existing `common::models::mutation_model::MutationModel` (mutation_rate, snp_trinuc_model, insertion_lengths, deletion_lengths, sv_model, homozygous_frequency, variant_dist) is structurally sufficient to serve as a tumor model. **No new `CancerModel` type is needed for v1.** The cancer-ness is a property of the training corpus, the orchestration, and the SV gap closures — not of the model schema.

The one feature NEAT 2.1 had that rneat does not is `HIGH_MUT_REGIONS` — a per-region rate scaler (~1.7k entries in their BRCA model, used to model mutation-rate hotspot beds). This is a v2 extension; the v1 MVP does not need it.

## SV gap analysis

The cancer MVP's credibility depends on rneat being able to generate the SVs that define cancer karyotypes. v1.10.1 generated `<DEL>` / `<DUP>` / `<CNV>` (covering the coverage-altering axis); **v1.12.0 closes the structural axis** by adding `<INS>` / `<INV>` / `<BND>` — which is what most clinically-recognizable cancer SVs depend on.

| Tier | Ticket | Gap | Cancer relevance | Status |
|---|---|---|---|---|
| **1** | **#187** | `<BND>` translocations with junction reads | **Critical** — BCR-ABL, PML-RARA, EWSR1-FLI1, MYC-IGH, TMPRSS2-ERG, etc. | ✅ Shipped in v1.12.0 |
| **1** | **#188** | `<INV>` inversions with read-strand flipping | Moderate — inv(16) AML, inv(3) MDS, focal inversions in solid tumors | ✅ Shipped in v1.12.0 |
| **2** | **#189** | Verify whole-arm / whole-chromosome aneuploidy via existing DEL/DUP | **High** (universal in cancer) but mostly verification, not new code | ✅ Resolved in v1.10.3 |
| **2** | **#190** | De novo `<INS>` with novel-sequence generation | Moderate — L1 retrotransposition (lung, colorectal, HCC); viral integrations (HPV, HBV, EBV) | ✅ Shipped in v1.12.0 (random-novel-sequence MVP) |
| **3** | **#191** | Complex rearrangement patterns (chromothripsis, chromoplexy, BFB) | **High in specific cancers** (sarcoma, GBM, prostate) — but requires correlated multi-event sampler | Design pending |
| **3** | **#192** | ecDNA / double minutes | **High in some aggressive cancers** (~25–50% of GBM) but needs first-class extrachromosomal data-model concept | Design pending |

## Roadmap

### Stage 0 — foundation (already shipped, v1.10.1)
- ✅ `<DEL>` / `<DUP>` / `<CNV>` de novo generation
- ✅ Lean VCF reader (enables training on whole-genome corpora)
- ✅ Full-genome gnomAD-SV-derived germline SV defaults

### Stage 1 — cancer MVP (orchestration validation)
1. ✅ **#189** — aneuploidy verification (resolved). Depth-modulation math (`build_coverage_multipliers` + `coverage_multiplier_for`) verified to handle whole-contig `<DEL>` / `<DUP>` / `<CNV>` events correctly at chromosome scale via `input_vcf:` — see `pipeline_e2e::whole_contig_*` tests. **De-novo aneuploidy via `SvModel::sample_variants` is structurally limited:** the sampler rejects lengths above `contig_len / 4` (see `sv_model.rs:431`) to keep its overlap-rejection search space tractable. Whole-chromosome and most whole-arm events therefore can't be generated by random draws from the log-normal length distribution. **The supported workflow for aneuploidy is `input_vcf:`** — supply segmentation-call output (GATK-CNV, ASCAT, PURPLE, etc.) and the depth-modulation primitive handles it correctly. A dedicated "aneuploidy event" sampler (with per-arm frequencies rather than a global rate) is a v2 follow-up, intentionally deferred.
2. **#183** — survey + commit to a tumor training corpus. Shipped: TCGA MC3 MAF adapter (`tools/fetch_tumor_corpus.sh`) and COSMIC GenomeScreensMutant VCF adapter (`tools/fetch_cosmic_corpus.sh`).
3. **#184** — orchestration script. Can develop in parallel with #183 using a stand-in tumor model.
4. **#186** — train initial tumor models from the chosen corpus.
5. **#185** — variant-origin annotation in the merged truth VCF. **Shipped:** every golden VCF now carries `INFO/NEAT_PROVENANCE` (`denovo` for sampled variants, `input` for those carried through from `input_vcf:`). `tools/cancer_simulate.sh` resolves these to `INFO/NEAT_ORIGIN` (`germline`/`somatic`/`shared`) via a `bcftools isec` post-merge — see the script's "Merge golden VCFs" section.

**Exit criterion:** running the orchestration script with reasonable defaults produces a merged FASTQ that a current somatic SNV caller (Mutect2 or Strelka) calls into a VCF that scores reasonably well against the origin-tagged truth.

### Stage 2 — cancer SV credibility (shipped in v1.12.0)
6. ✅ **#187** — `<BND>` translocations with junction reads. Chimeric-pair generation via `process_chimeric_variants` → `generate_chimeric_pair` → `get_bnd_pieces` + `get_stitched_sequence`, correctly handling the four VCF 4.2 breakend orientation forms (`t[p[`, `t]p]`, `[p[t`, `]p]t`). Paired-BND deduplication via canonical-ID HashSet.
7. ✅ **#188** — `<INV>` inversions with strand flipping. Two-junction model (start of inversion + end of inversion), each independently sampling reads. `generate_inv_pair` + `get_inv_pieces` handle orientation flipping per junction. Smarter fragment-length sampling + graceful `TruncatedRead` handling shipped alongside, applied to BND too.
8. ✅ **#190** — de novo `<INS>` with novel sequence. Emitted as LITERAL Insertion records (`REF="A"`, `ALT="A<novel-bases>"`) — same shape Manta / DELLY produce when they resolve the inserted sequence. Routes through the existing literal-insertion machinery. Novel sequence drawn from single-base composition of a ±250 bp window around the anchor.

**Exit criterion:** a current somatic SV caller (GRIDSS or Manta) calls translocation breakpoints from the simulated data with high enough recall that the simulator could plausibly be used as a benchmark fixture. **Status:** functionally validated against a synthetic 5 kb fixture (see v1.12.0 CHANGELOG); a real-caller benchmark (Manta / DELLY swapped into `tools/cancer_benchmark.sh` for Mutect2) is the natural v1.12.1 follow-up.

#### Known limitations carried into v1.12.0
- **Breakpoint double-counting.** The regular per-contig pass still generates reads covering BND / INV breakpoint positions (it reads from the unbroken reference), so the breakpoint locus ends up with regular reads PLUS junction reads — roughly 2× coverage for homozygous variants, closer to correct for heterozygous. A proper fix would teach the regular pass to skip the broken-allele fraction at these positions; deferred to v2.
- **Default model has no BND / INV.** ~~The bundled gnomAD-SV-derived `default_sv_model()` filters BND / INV out at training time.~~ **Resolved in v1.12.1 (#217)** via `tools/inject_cancer_sv_model.py`: the bundled COSMIC tumor model (`tools/cosmic_v104_pancancer_model.json.gz`) now carries a literature-derived `sv_model` covering all six SV types. The germline-side `default_sv_model()` was already populated at heuristic rates. A data refit from real cancer-SV calls **shipped in v1.14.0 (#218)** — the bundled tumor model's `sv_model` is now counted from the PCAWG consensus SV/CNV callsets (+ gnomAD-SV for INS length); see the "Cancer-specific SvModel defaults" note below.
- **`FORMAT/AD/DP/AF` asymmetric across SV types.** Symbolic SVs (BND, INV, DEL, DUP, CNV) emit `.,.:.:.` placeholders since span-based depth doesn't fit a point-based counter. INS is now literal and gets real counts.
- **No mobile-element library for de novo INS.** Insertions are random novel sequence, not drawn from L1 / Alu / SVA / HERV consensus. Library-driven path is a #190 v2 extension; viral integration is a separate design problem.

### Stage 3 — advanced cancer patterns (defer until stage 2 lands)
9. **#191** — chromothripsis / chromoplexy / BFB design pass.
10. **#192** — ecDNA / double-minute design pass.

Stage 3 starts with **design notes**, not implementation. Either of these could prove out-of-scope for rneat and be better served by a sibling tool that consumes rneat's output.

## Open questions

These are documented here so future work has a record of what's been decided vs what's still up for grabs:

- **Subclonal heterogeneity.** The orchestration approach models two populations. Real tumors often have 3–5+ subclones at different fractions. v1 ignores this; v2 should think about whether to support N populations via the orchestration layer (cheap, N-way concat) or via a native sampler (requires shared-mutation tracking).
- **Mutation signatures (COSMIC SBS).** Cancer mutation patterns are dominated by tissue- and exposure-specific signatures (UV → SBS7 in melanoma, tobacco → SBS4 in lung, MMR deficiency → SBS6/15 in colorectal, etc.). rneat's trinucleotide-frequency model can in principle capture these per-corpus, but explicit signature support (load a COSMIC SBS profile, mix signatures at known proportions) is a separate enhancement.
- **Variant allele fraction (VAF) annotation in the truth VCF.** **Shipped via #176.** Golden VCFs now carry `FORMAT/GT:AD:DP:AF` on every literal SNP/indel record, populated from a per-variant allelic-depth counter incremented by the gen-reads fragment loop at the existing per-read coin-flip site (`fastq_tools.rs::generate_read`). Symbolic SVs emit `.` placeholders for AD/DP/AF since their depth semantics are span-based, not point-based — that's tracked as a separate v2 enhancement. Realistic `FILTER` values (`LowQual` etc.) were left out of v1; everything is still `PASS`.
- **Cancer-specific SvModel defaults.** First shipped in v1.12.1 (#217) as a *literature-derived* `sv_model` injection (`tools/inject_cancer_sv_model.py`, now deprecated). **Replaced in v1.14.0 (#218) with a data-derived refit**: the bundled tumor model's `sv_model` is counted from the PCAWG consensus SV/CNV callsets (Li et al. 2020, n=2,748 donors) with gnomAD-SV supplying INS length, via `tools/fetch_pcawg_sv_corpus.sh` → `tools/build_pcawg_sv_vcf.py` → `gen-mut-model` → `tools/normalize_pcawg_sv_model.py`. The refit corrected the heuristic type mix (notably INV 0.08→0.134, BND 0.35→0.232) to match the published PCAWG ranking (DEL>DUP>BND>INV); the fitted per-base rate 3.46e-8 confirmed the prior 3.8e-8 estimate. INS and focal-CNV *rates* remain literature estimates (PCAWG MEI calls are controlled-access; focal-CNA segment counts are segmentation, not discrete events). The germline-side `default_sv_model()` stays gnomAD-derived. See CHANGELOG v1.14.0.
- **Per-tissue SV models (#202).** **Shipped.** Three tissue-specific tumor models
  are bundled alongside the pan-cancer one:
  `tools/cosmic_pancancer_sv_{BRCA,skin,lung}.json.gz`. Each is the pan-cancer
  COSMIC SNP/indel base with a `sv_model` fitted from that tissue's PCAWG donors.
  `tools/build_pcawg_sv_vcf.py` gained `--sample-sheet` + `--projects` (filtering
  donors by `aliquot_id → dcc_project_code` via `pcawg_sample_sheet.tsv`), and
  `tools/build_per_tissue_sv_models.sh` drives build→fit→normalize per tissue.

  Tissue → PCAWG project codes (SV-donor counts): **BRCA** → `BRCA` (211);
  **skin** → `SKCM-US,MELA-AU` (106, melanoma); **lung** → `LUAD-US,LUSC-US` (85).
  PCAWG's codes diverge from TCGA's and SKCM/LUAD alone are underpowered (~37
  each), so skin/lung group the related projects. The fits show expected
  tissue-specific signal — BRCA is **DUP-dominant** (Dup 0.317, the tandem-
  duplicator phenotype) at ~2× SV burden; skin is **BND-enriched** (0.306);
  lung is **DEL-dominant** (0.349).

  Pick one at run time: `cancer_simulate.sh --tumor-model tools/cosmic_pancancer_sv_BRCA.json.gz`.

  Caveats: only the **SV component** is tissue-specific — the SNP/indel side stays
  pan-cancer COSMIC (per-tissue SNP/indel is the MC3/COSMIC-classification half of
  #202, still open; COSMIC's classification file is a separate academic download).
  INS length is shared from the gnomAD-derived base (germline, tissue-agnostic);
  INS and focal-CNV *rates* remain pan-cancer literature values, so a tissue's CNV
  fraction shifts mainly because its SV burden (the denominator) differs.

## References

- **Umbrella issue:** [#129 — Cancer modeling](https://github.com/ncsa/rusty-neat/issues/129)
- **Cancer MVP sub-tickets:** [#183](https://github.com/ncsa/rusty-neat/issues/183), [#184](https://github.com/ncsa/rusty-neat/issues/184), [#185](https://github.com/ncsa/rusty-neat/issues/185), [#186](https://github.com/ncsa/rusty-neat/issues/186)
- **SV gap tickets:** [#187](https://github.com/ncsa/rusty-neat/issues/187), [#188](https://github.com/ncsa/rusty-neat/issues/188), [#189](https://github.com/ncsa/rusty-neat/issues/189), [#190](https://github.com/ncsa/rusty-neat/issues/190), [#191](https://github.com/ncsa/rusty-neat/issues/191), [#192](https://github.com/ncsa/rusty-neat/issues/192)
- **NEAT 2.1 prior art:** `~/code/neat2/models/genReadsTumorTutorial.zip` — 30-line `simulate.sh` reference implementation; `~/code/neat2/utilities/genMutModel.py` — original tumor-model trainer that parsed ICGC simple-somatic-mutation TSVs.
- **gnomAD-SV-derived germline defaults** that the cancer side will complement: `common/src/models/sv_model_defaults.rs` (refit in v1.10.1).

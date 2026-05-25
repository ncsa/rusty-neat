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
| Source tumor training corpus | **#183** | Medium | Survey GDC / COSMIC / PCAWG / legacy ICGC. Working assumption: GDC TCGA MAFs |
| Orchestration script (`tools/cancer_simulate.sh`) | **#184** | Low | Ports NEAT 2.1's `simulate.sh` idiom |
| Variant-origin annotation in golden VCF | **#185** | Medium | `INFO/NEAT_ORIGIN={germline,somatic,shared}` |
| Train v1 tumor models | **#186** | Low (per tissue) | Run `gen-mut-model` on the chosen corpus. Likely BRCA / SKCM / LUAD |

#184 can develop in parallel with #183 using a hand-rolled or NEAT-2.1-recovered model as a stand-in. #186 blocks on #183 (data) but not on anything else.

## Source-data decisions

### Chosen for v1: NCI Genomic Data Commons (GDC) TCGA MAFs

- **Open access**, no DACO equivalent required.
- **Programmatic API** via `gdc-client`.
- **Project-stratified** so models can be trained per-tissue.
- **Modern format** (MAF v1.0+) with consistent column semantics.
- Conversion to VCF is a well-supported one-time preprocessing step.

### Considered and ruled out (for v1)

| Source | Why not for v1 |
|---|---|
| **Legacy ICGC25K** (NEAT 2.1's source) | The original `dcc.icgc.org` portal was retired June 2024. Open-tier bulk data is reportedly mirrored on AWS S3, but URL stability is not guaranteed and the format (simple-somatic-mutation TSV) requires a custom parser that `gen-mut-model` doesn't have. |
| **COSMIC** | Free for academic use but requires login + license click-through, which complicates an open-source default recipe. Reconsider for v2 if GDC coverage is insufficient. |
| **PCAWG consensus calls** | Excellent quality, includes SV consensus, but smaller donor count than TCGA aggregate. Better suited to **#191** (chromothripsis / chromoplexy training) than to general somatic-SNV models. |
| **EGA controlled-access** (e.g. CLLE-ES, the Spanish CLL project NEAT 2.1 included) | Controlled tier — not viable as a default for an open-source tool. Users with EGA access can train their own. |

## Mutation model: no new schema needed

NEAT 2.1's three tumor models (BRCA, CLLE, SKCM) used the **same pickled-dict schema** as its germline NA12878 model. The "tumor-ness" came purely from the input training data — the model schema captures whatever statistics the input VCF supplies.

rneat's existing `common::models::mutation_model::MutationModel` (mutation_rate, snp_trinuc_model, insertion_lengths, deletion_lengths, sv_model, homozygous_frequency, variant_dist) is structurally sufficient to serve as a tumor model. **No new `CancerModel` type is needed for v1.** The cancer-ness is a property of the training corpus, the orchestration, and the SV gap closures — not of the model schema.

The one feature NEAT 2.1 had that rneat does not is `HIGH_MUT_REGIONS` — a per-region rate scaler (~1.7k entries in their BRCA model, used to model mutation-rate hotspot beds). This is a v2 extension; the v1 MVP does not need it.

## SV gap analysis

The cancer MVP's credibility depends on rneat being able to generate the SVs that define cancer karyotypes. v1.10.1 generates `<DEL>` / `<DUP>` / `<CNV>` (covering the coverage-altering axis) but not `<INS>` / `<INV>` / `<BND>` (the structural axis — which is what most clinically-recognizable cancer SVs depend on).

| Tier | Ticket | Gap | Cancer relevance | Effort |
|---|---|---|---|---|
| **1** | **#187** | `<BND>` translocations with junction reads | **Critical** — BCR-ABL, PML-RARA, EWSR1-FLI1, MYC-IGH, TMPRSS2-ERG, etc. | Medium-high |
| **1** | **#188** | `<INV>` inversions with read-strand flipping | Moderate — inv(16) AML, inv(3) MDS, focal inversions in solid tumors | Medium |
| **2** | **#189** | Verify whole-arm / whole-chromosome aneuploidy via existing DEL/DUP | **High** (universal in cancer) but mostly verification, not new code | Low |
| **2** | **#190** | De novo `<INS>` with novel-sequence generation | Moderate — L1 retrotransposition (lung, colorectal, HCC); viral integrations (HPV, HBV, EBV) | Low |
| **3** | **#191** | Complex rearrangement patterns (chromothripsis, chromoplexy, BFB) | **High in specific cancers** (sarcoma, GBM, prostate) — but requires correlated multi-event sampler | Medium-high (design first) |
| **3** | **#192** | ecDNA / double minutes | **High in some aggressive cancers** (~25–50% of GBM) but needs first-class extrachromosomal data-model concept | High (design first) |

## Roadmap

### Stage 0 — foundation (already shipped, v1.10.1)
- ✅ `<DEL>` / `<DUP>` / `<CNV>` de novo generation
- ✅ Lean VCF reader (enables training on whole-genome corpora)
- ✅ Full-genome gnomAD-SV-derived germline SV defaults

### Stage 1 — cancer MVP (orchestration validation)
1. **#189** — verify aneuploidy via existing DEL/DUP works at chromosome scale. Pure verification; cheap warmup task.
2. **#183** — survey + commit to a tumor training corpus. Working assumption: GDC TCGA MAFs.
3. **#184** — orchestration script. Can develop in parallel with #183 using a stand-in tumor model.
4. **#186** — train initial tumor models from the chosen corpus.
5. **#185** — variant-origin annotation in the merged truth VCF.

**Exit criterion:** running the orchestration script with reasonable defaults produces a merged FASTQ that a current somatic SNV caller (Mutect2 or Strelka) calls into a VCF that scores reasonably well against the origin-tagged truth.

### Stage 2 — cancer SV credibility
6. **#187** — `<BND>` translocations with junction reads.
7. **#188** — `<INV>` inversions with strand flipping.
8. **#190** — de novo `<INS>` with novel sequence.

**Exit criterion:** a current somatic SV caller (GRIDSS or Manta) calls translocation breakpoints from the simulated data with high enough recall that the simulator could plausibly be used as a benchmark fixture.

### Stage 3 — advanced cancer patterns (defer until stage 2 lands)
9. **#191** — chromothripsis / chromoplexy / BFB design pass.
10. **#192** — ecDNA / double-minute design pass.

Stage 3 starts with **design notes**, not implementation. Either of these could prove out-of-scope for rneat and be better served by a sibling tool that consumes rneat's output.

## Open questions

These are documented here so future work has a record of what's been decided vs what's still up for grabs:

- **Subclonal heterogeneity.** The orchestration approach models two populations. Real tumors often have 3–5+ subclones at different fractions. v1 ignores this; v2 should think about whether to support N populations via the orchestration layer (cheap, N-way concat) or via a native sampler (requires shared-mutation tracking).
- **Mutation signatures (COSMIC SBS).** Cancer mutation patterns are dominated by tissue- and exposure-specific signatures (UV → SBS7 in melanoma, tobacco → SBS4 in lung, MMR deficiency → SBS6/15 in colorectal, etc.). rneat's trinucleotide-frequency model can in principle capture these per-corpus, but explicit signature support (load a COSMIC SBS profile, mix signatures at known proportions) is a separate enhancement.
- **Variant allele fraction (VAF) annotation in the truth VCF.** Critical for benchmarking VAF-aware somatic callers (Mutect2, Strelka), and required by `hap.py` / `som.py` for VAF-stratified TP/FP/FN scoring. Tracked at **#176** (with `FORMAT/AD`, `FORMAT/DP`, `FORMAT/AF`, and FILTER companions). Either lands before the cancer MVP or co-lands with **#184** orchestration — the cancer MVP isn't truly benchmarkable without it.
- **Cancer-specific SvModel defaults.** The bundled `default_sv_model()` is gnomAD-SV-derived (germline). A cancer-specific default trained against PCAWG SV consensus would be the natural complement once **#187** / **#188** land. Possibly bundle as `default_cancer_sv_model()` alongside the germline default.

## References

- **Umbrella issue:** [#129 — Cancer modeling](https://github.com/ncsa/rusty-neat/issues/129)
- **Cancer MVP sub-tickets:** [#183](https://github.com/ncsa/rusty-neat/issues/183), [#184](https://github.com/ncsa/rusty-neat/issues/184), [#185](https://github.com/ncsa/rusty-neat/issues/185), [#186](https://github.com/ncsa/rusty-neat/issues/186)
- **SV gap tickets:** [#187](https://github.com/ncsa/rusty-neat/issues/187), [#188](https://github.com/ncsa/rusty-neat/issues/188), [#189](https://github.com/ncsa/rusty-neat/issues/189), [#190](https://github.com/ncsa/rusty-neat/issues/190), [#191](https://github.com/ncsa/rusty-neat/issues/191), [#192](https://github.com/ncsa/rusty-neat/issues/192)
- **NEAT 2.1 prior art:** `~/code/neat2/models/genReadsTumorTutorial.zip` — 30-line `simulate.sh` reference implementation; `~/code/neat2/utilities/genMutModel.py` — original tumor-model trainer that parsed ICGC simple-somatic-mutation TSVs.
- **gnomAD-SV-derived germline defaults** that the cancer side will complement: `common/src/models/sv_model_defaults.rs` (refit in v1.10.1).

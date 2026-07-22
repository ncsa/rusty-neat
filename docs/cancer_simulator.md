# Cancer simulator design

Status: **shipped** — the cancer epic (#129) is closed. The two-pass orchestration, the native `gen-cancer-reads` subcommand, SV/CNV support, per-tissue models, and the PCAWG SV refit all landed across v1.12.0–v1.16.0 (see CHANGELOG); #183–#192 are closed. This document is retained as the design/rationale record — read the "Update" notes inline for current state.

## Goal

Produce simulated short-read sequencing data with the statistical and structural properties of a real cancer sample, including:

1. A **mixture of tumor + normal DNA** at a configurable purity (real tumor biopsies are typically 30–80% tumor cells; the rest is normal-tissue contamination).
2. **Germline variation shared between the two populations** (the tumor and the normal both came from the same patient).
3. **Tumor-only somatic variation** in the tumor population — elevated SNV / indel rates, cancer-specific mutation signatures, and a richer SV catalog than the gnomAD-SV-derived germline default.

The output must be **benchmarkable** — a downstream somatic variant caller (Mutect2, Strelka, GRIDSS, Manta, …) should be able to run on eidolon's simulated reads and have its calls scored against a ground-truth VCF that distinguishes germline from somatic.

## Approach: orchestration first, native subcommand later

The recommended v1 path is a **thin orchestration layer** in `tools/` that drives `eidolon gen-reads` twice — once per population — and concatenates the resulting FASTQs. No core eidolon changes for the MVP; we lean on the existing primitives (`input_vcf`, `mutation_model`, `sv_rate_scale`, paired-end fragment generation, depth modulation).

Why this approach:

- **Validated by precedent.** NEAT 2.1 shipped this exact workflow as its `genReadsTumorTutorial` (~30-line shell script, archived at `~/code/neat2/models/genReadsTumorTutorial.zip`). It worked then; it'll work now.
- **Minimum new code.** Until the simulation produces credible output end-to-end, every line of new core code is risk. The orchestration path gets us to a working artifact with zero eidolon-core diffs.
- **Reversible.** If we later decide a native `eidolon gen-cancer-reads` subcommand is worth the engineering, the orchestration script is the design spec for that subcommand's behavior — it codifies the data flow.

A native subcommand was **not** an explicit non-goal — it was a v2+ option pending validation that the orchestration produces samples real somatic callers handle correctly. That validation passed, and the native subcommand shipped (v1.15.0, #239) — see the Update below.

**Update (#239): the native subcommand `eidolon gen-cancer-reads` is now implemented** — the orchestration validated cleanly (Mutect2 SNV + Manta SV recall), so the two-pass flow is now a first-class subcommand: `eidolon gen-cancer-reads -c <yaml>` runs both passes and merges them (tagged FASTQs + an origin-tagged truth VCF), with no `bcftools`/`awk` runtime dependency. Config template: `template_config/gen_cancer_reads_template.yml`; design + decisions: `docs/archive/cancer_simulator_native_plan.md` (archived — shipped). `tools/cancer_simulate.sh` is retained; the same-seed parity test (`eidolon/tests/cancer_parity.rs`) landed in v1.15.1.

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

Both adapters can chain into `eidolon gen-mut-model` via `--train`. The output mutation model is structurally identical regardless of which adapter produced the input VCF — "cancer-ness" is in the training data, not the schema.

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

eidolon's existing `eidolon_core::models::mutation_model::MutationModel` (mutation_rate, snp_trinuc_model, insertion_lengths, deletion_lengths, sv_model, homozygous_frequency, variant_dist) is structurally sufficient to serve as a tumor model. **No new `CancerModel` type is needed for v1.** The cancer-ness is a property of the training corpus, the orchestration, and the SV gap closures — not of the model schema.

The one feature NEAT 2.1 had that eidolon does not is `HIGH_MUT_REGIONS` — a per-region rate scaler (~1.7k entries in their BRCA model, used to model mutation-rate hotspot beds). This is a v2 extension (tracked in #413); the v1 MVP does not need it.

## SV gap analysis

The cancer MVP's credibility depends on eidolon being able to generate the SVs that define cancer karyotypes. v1.10.1 generated `<DEL>` / `<DUP>` / `<CNV>` (covering the coverage-altering axis); **v1.12.0 closes the structural axis** by adding `<INS>` / `<INV>` / `<BND>` — which is what most clinically-recognizable cancer SVs depend on.

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
1. ✅ **#189** — aneuploidy verification (resolved). Depth-modulation math (`build_coverage_multipliers` + `coverage_multiplier_for`) verified to handle whole-contig `<DEL>` / `<DUP>` / `<CNV>` events correctly at chromosome scale via `input_vcf:` — see `pipeline_e2e::whole_contig_*` tests. **De-novo aneuploidy via `SvModel::sample_variants` is bounded by a configurable length cap:** the sampler rejects lengths above `contig_len * sv_max_length_fraction` to keep its overlap-rejection search space tractable. The default fraction is `0.25` (`contig_len / 4`), which reproduces the historical behavior bit-for-bit and is correct for human chromosomes — whole-chromosome and most whole-arm events still can't be generated by random draws under the default. **#229** made this cap a gen-reads config knob: set `sv_max_length_fraction:` toward `1.0` to allow large de-novo events (the sampler also scales its retry budget proportionally and warns when it starves). This is most useful on small contigs (bacteria, viruses), where the default `0.25` sits below the trained DEL median and starves DEL/DUP/INV into BND; for chromosome-scale aneuploidy with realistic per-arm frequencies, **the supported workflow remains `input_vcf:`** — supply segmentation-call output (GATK-CNV, ASCAT, PURPLE, etc.) and the depth-modulation primitive handles it correctly. A dedicated "aneuploidy event" sampler (with per-arm frequencies rather than a global rate) is a v2 follow-up, intentionally deferred.
2. **#183** — survey + commit to a tumor training corpus. Shipped: TCGA MC3 MAF adapter (`tools/fetch_tumor_corpus.sh`) and COSMIC GenomeScreensMutant VCF adapter (`tools/fetch_cosmic_corpus.sh`).
3. **#184** — orchestration script. Can develop in parallel with #183 using a stand-in tumor model.
4. **#186** — train initial tumor models from the chosen corpus.
5. **#185** — variant-origin annotation in the merged truth VCF. **Shipped:** every golden VCF now carries `INFO/NEAT_PROVENANCE` (`denovo` for sampled variants, `input` for those carried through from `input_vcf:`). `tools/cancer_simulate.sh` resolves these to `INFO/NEAT_ORIGIN` (`germline`/`somatic`/`shared`) via a `bcftools isec` post-merge — see the script's "Merge golden VCFs" section.

**Exit criterion:** running the orchestration script with reasonable defaults produces a merged FASTQ that a current somatic SNV caller (Mutect2 or Strelka) calls into a VCF that scores reasonably well against the origin-tagged truth.

### Stage 2 — cancer SV credibility (shipped in v1.12.0)
6. ✅ **#187** — `<BND>` translocations with junction reads. Chimeric-pair generation via `process_chimeric_variants` → `generate_chimeric_pair` → `get_bnd_pieces` + `get_stitched_sequence`, correctly handling the four VCF 4.2 breakend orientation forms (`t[p[`, `t]p]`, `[p[t`, `]p]t`). Paired-BND deduplication via canonical-ID HashSet.
7. ✅ **#188** — `<INV>` inversions with strand flipping. Two-junction model (start of inversion + end of inversion), each independently sampling reads. `generate_inv_pair` + `get_inv_pieces` handle orientation flipping per junction. Smarter fragment-length sampling + graceful `TruncatedRead` handling shipped alongside, applied to BND too.
8. ✅ **#190** — de novo `<INS>` with novel sequence. Emitted as LITERAL Insertion records (`REF="A"`, `ALT="A<novel-bases>"`) — same shape Manta / DELLY produce when they resolve the inserted sequence. Routes through the existing literal-insertion machinery. Novel sequence drawn from single-base composition of a ±250 bp window around the anchor.

**Exit criterion:** a current somatic SV caller (GRIDSS or Manta) calls translocation breakpoints from the simulated data with high enough recall that the simulator could plausibly be used as a benchmark fixture. **Status:** validated — functionally against a synthetic 5 kb fixture (v1.12.0 CHANGELOG), and against real callers on Delta (Manta + Delly per-type SV recall, ACCESS report §3.7; cross-caller coverage #317, closed).

#### Known limitations carried into v1.12.0
- **Breakpoint double-counting — resolved in v1.14.1 (#236).** The regular per-contig pass used to also generate reads across BND / INV breakpoint positions (from the unbroken reference), leaving those loci with regular + junction reads — ~2× coverage for homozygous variants. Fixed: the regular pass now suppresses the broken-allele fraction of junction-crossing read-pairs (`collect_suppressible_junctions` / `suppress_junction_double_count`; DEL/BND/INV/CNV-loss, with DUP/CNV-gain intentionally excluded).
- **Default model has no BND / INV.** ~~The bundled gnomAD-SV-derived `default_sv_model()` filters BND / INV out at training time.~~ **Resolved in v1.12.1 (#217)** via `tools/inject_cancer_sv_model.py`: the bundled COSMIC tumor model (`tools/cosmic_v104_pancancer_model.json.gz`) now carries a literature-derived `sv_model` covering all six SV types. The germline-side `default_sv_model()` was already populated at heuristic rates. A data refit from real cancer-SV calls **shipped in v1.14.0 (#218)** — the bundled tumor model's `sv_model` is now counted from the PCAWG consensus SV/CNV callsets (+ gnomAD-SV for INS length); see the "Cancer-specific SvModel defaults" note below.
- **`FORMAT/AD/DP/AF` asymmetric across SV types.** Symbolic SVs (BND, INV, DEL, DUP, CNV) emit `.,.:.:.` placeholders since span-based depth doesn't fit a point-based counter. INS is now literal and gets real counts. (Span-based depth for symbolic SVs is tracked in #411.)
- **No mobile-element library for de novo INS.** Insertions are random novel sequence, not drawn from L1 / Alu / SVA / HERV consensus. Library-driven path is a #190 v2 extension; viral integration is a separate design problem.

### Stage 3 — advanced cancer patterns (defer until stage 2 lands)
9. **#191** — chromothripsis / chromoplexy / BFB design pass.
10. **#192** — ecDNA / double-minute design pass.

Stage 3 starts with **design notes**, not implementation. Either of these could prove out-of-scope for eidolon and be better served by a sibling tool that consumes eidolon's output.

## Open questions

These are documented here so future work has a record of what's been decided vs what's still up for grabs:

- **Subclonal heterogeneity.** The orchestration approach models two populations. Real tumors often have 3–5+ subclones at different fractions. v1 ignores this; v2 should think about whether to support N populations via the orchestration layer (cheap, N-way concat) or via a native sampler (requires shared-mutation tracking).
- **Mutation signatures (COSMIC SBS).** Cancer mutation patterns are dominated by tissue- and exposure-specific signatures (UV → SBS7 in melanoma, tobacco → SBS4 in lung, MMR deficiency → SBS6/15 in colorectal, etc.). eidolon's trinucleotide-frequency model captures these per-corpus, and **context-weighted SNP placement (#372, v1.20.0) now reproduces context-specific signatures** — validated on HCC1395 (SBS-96 cosine 0.72 → 0.99) and via signature extraction on chr1–3 (SBS1 + APOBEC recovered; #320 closed). What remains optional is explicit signature *authoring* (load a named COSMIC SBS profile and mix signatures at chosen proportions), a separate enhancement.
- **Variant allele fraction (VAF) annotation in the truth VCF.** **Shipped via #176.** Golden VCFs now carry `FORMAT/GT:AD:DP:AF` on every literal SNP/indel record, populated from a per-variant allelic-depth counter incremented by the gen-reads fragment loop at the existing per-read coin-flip site (`fastq_tools.rs::generate_read`). Symbolic SVs emit `.` placeholders for AD/DP/AF since their depth semantics are span-based, not point-based — that's tracked in #411. Realistic `FILTER` values (`LowQual` etc.) were left out of v1; everything is still `PASS`.
- **Cancer-specific SvModel defaults.** First shipped in v1.12.1 (#217) as a *literature-derived* `sv_model` injection (`tools/inject_cancer_sv_model.py`, now deprecated). **Replaced in v1.14.0 (#218) with a data-derived refit**: the bundled tumor model's `sv_model` is counted from the PCAWG consensus SV/CNV callsets (Li et al. 2020, n=2,748 donors) with gnomAD-SV supplying INS length, via `tools/fetch_pcawg_sv_corpus.sh` → `tools/build_pcawg_sv_vcf.py` → `gen-mut-model` → `tools/normalize_pcawg_sv_model.py`. The refit corrected the heuristic type mix (notably INV 0.08→0.134, BND 0.35→0.232) to match the published PCAWG ranking (DEL>DUP>BND>INV); the fitted per-base rate 3.46e-8 confirmed the prior 3.8e-8 estimate. INS and focal-CNV *rates* remain literature estimates (PCAWG MEI calls are controlled-access; focal-CNA segment counts are segmentation, not discrete events). The germline-side `default_sv_model()` stays gnomAD-derived. See CHANGELOG v1.14.0.
- **Per-tissue tumor models (#202).** **Shipped — both halves now tissue-specific.**
  Three complete tissue models are bundled alongside the pan-cancer one:
  `tools/cosmic_per_tissue_{BRCA,skin,lung}.json.gz`. Each pairs a **per-tissue
  COSMIC SNP/indel spectrum** with a **per-tissue `sv_model`**.

  *SV half (first shipped in #237).* The `sv_model` is fitted from that tissue's
  PCAWG donors. `tools/build_pcawg_sv_vcf.py` gained `--sample-sheet` + `--projects`
  (filtering donors by `aliquot_id → dcc_project_code` via `pcawg_sample_sheet.tsv`),
  and `tools/build_per_tissue_sv_models.sh` drives build→fit→normalize per tissue,
  emitting `tools/cosmic_pancancer_sv_{BRCA,skin,lung}.json.gz` (pan-cancer SNP/indel
  base + per-tissue SV). Tissue → PCAWG project codes (SV-donor counts): **BRCA** →
  `BRCA` (211); **skin** → `SKCM-US,MELA-AU` (106, melanoma); **lung** →
  `LUAD-US,LUSC-US` (85). PCAWG's codes diverge from TCGA's and SKCM/LUAD alone are
  underpowered (~37 each), so skin/lung group related projects. The fits show
  expected signal — BRCA **DUP-dominant** (Dup 0.317, the tandem-duplicator
  phenotype) at ~2× SV burden; skin **BND-enriched** (0.306); lung **DEL-dominant**
  (0.349).

  *SNP/indel half (#202 completion).* Each tissue's SNP/indel spectrum is trained
  from its COSMIC GenomeScreensMutant subset. Because the COSMIC VCF is
  tissue-aggregated, the split uses the sample-level **TSV** export:
  `tools/fetch_cosmic_per_tissue_corpus.sh` reads the Classification file
  (`COSMIC_PHENOTYPE_ID → PRIMARY_SITE`), collects every COSV mutation ID seen in
  that tissue's phenotypes from the mutant TSV, and filters the GenomeScreensMutant
  VCF (which carries VCF-anchored alleles, so indels need no re-anchoring) down to
  those IDs. `tools/build_per_tissue_models.sh` orchestrates the full chain
  (fetch → `gen-mut-model` train → `tools/graft_sv_model.py` grafts the per-tissue
  `sv_model` on). The fits reproduce known tissue biology: **skin** is the most
  SNV-dominated (melanoma's UV-driven C>T; indels only ~1.6% of records), while
  **breast** carries the most indels (~5.2%); **lung** sits between (~3.6%).

  Pick one at run time: `gen-cancer-reads` with
  `tumor_model: tools/cosmic_per_tissue_BRCA.json.gz`.

  Caveats: the bundled `mutation_rate` is COSMIC-corpus-aggregated (fraction of bp
  mutated across that tissue's catalog), so it overstates single-tumor burden — set
  per-tumor burden via `tumor_mutation_rate` (defaults to 1e-5); what's genuinely
  per-tissue is the *spectrum* (trinucleotide context, indel ratio, SV mix). INS
  length is shared from the gnomAD-derived base (germline, tissue-agnostic); INS and
  focal-CNV *rates* remain pan-cancer literature values, so a tissue's CNV fraction
  shifts mainly because its SV burden (the denominator) differs. The
  `cosmic_pancancer_sv_*` files (pan-cancer SNP/indel + per-tissue SV) are retained
  as the `sv_model` donor for the graft step.

## References

- **Umbrella issue:** [#129 — Cancer modeling](https://github.com/ncsa/eidolon/issues/129)
- **Cancer MVP sub-tickets:** [#183](https://github.com/ncsa/eidolon/issues/183), [#184](https://github.com/ncsa/eidolon/issues/184), [#185](https://github.com/ncsa/eidolon/issues/185), [#186](https://github.com/ncsa/eidolon/issues/186)
- **SV gap tickets:** [#187](https://github.com/ncsa/eidolon/issues/187), [#188](https://github.com/ncsa/eidolon/issues/188), [#189](https://github.com/ncsa/eidolon/issues/189), [#190](https://github.com/ncsa/eidolon/issues/190), [#191](https://github.com/ncsa/eidolon/issues/191), [#192](https://github.com/ncsa/eidolon/issues/192)
- **NEAT 2.1 prior art:** `~/code/neat2/models/genReadsTumorTutorial.zip` — 30-line `simulate.sh` reference implementation; `~/code/neat2/utilities/genMutModel.py` — original tumor-model trainer that parsed ICGC simple-somatic-mutation TSVs.
- **gnomAD-SV-derived germline defaults** that the cancer side will complement: `eidolon-core/src/models/sv_model_defaults.rs` (refit in v1.10.1).

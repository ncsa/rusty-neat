# SCN track — status

Consolidated status for the soybean-cyst-nematode (SCN, *Heterodera glycines*) real-data
work under the realism epic (#311). Data from **João Gomes Viana** (Matt Hudson lab).
Keep this current as the SCN work progresses. Last updated: **2026-07-19**.

Goal: use real SCN Pool-seq data to (Phase 1) build genuine per-strain eidolon models covering
**all variant classes — SNP/indel, structural variants, and copy-number** — and (Phase 2)
reproduce the population allele-frequency spectrum in simulated reads.

**Status: the planned track is COMPLETE.** Both strains have full four-class models built from
real data and verified end-to-end (build → simulate → render); the virulent-vs-avirulent
contrast is done. Remaining items are enhancements (see *Open / next steps*).

## Data & accessions (BioProject PRJNA1055977)

**Strain reference assemblies** (Walden et al. 2025, doi:10.1186/s12864-025-12493-x; via NCBI `datasets`):

| strain | assembly | phenotype | genome |
|---|---|---|---|
| MM26 | `GCA_040805935.1` | virulent | 144.9 Mb |
| PA3  | `GCA_040805705.1` | avirulent ("picky eater") | 123.7 Mb |

**Population Pool-seq reads** (Illumina NextSeq 2000, paired, pools of 150 virgin females):

| pool | SRA run | aligned to |
|---|---|---|
| MM26A   | `SRR27329602` | MM26 (`GCA_040805935.1`) |
| MM-BD3A | `SRR27329600` | PA3 (`GCA_040805705.1`) |

**Reference strategy:** each pool is aligned to its best-fit strain reference (MM26→MM26,
BD3→PA3) to mitigate reference bias, per João. **Consequence for the contrast:** MM26A is a
*same-strain* pairing while BD3A→PA3 is *cross-population*, so cross-strain raw comparisons are
confounded (see Contrast). A clean virulence comparison needs a **common reference** (not yet done).

**Pool-seq caveat:** pooled samples (150 individuals → allele **frequencies**, not diploid
genotypes). Use allele **depths (AD)**, not `GT`, for anything AF-related.

## Phase 1 — per-strain SNP/indel models — **DONE ✅**

`scripts/delta/stage_scn.sh` stages a self-consistent reference + BAM + VCF (align pool reads →
call VCF *from* that BAM), then `model_builders.sbatch` builds all five models. Verified by the
actual statistics (not just PASS):

| statistic | MM26 / MM26A | PA3 / MM-BD3A |
|---|---|---|
| mutation_rate | 2.44×10⁻³ | 3.71×10⁻³ |
| homozygous_frequency | 0.126 | 0.298 |
| error_rate | 4.35×10⁻³ | 3.99×10⁻³ (consistent — sequencing property ✓) |
| frag length | Normal(516, 154) | Normal(448, 178) |

## Phase 1b — Structural variants + copy-number — **DONE ✅ (the "complex organism" payoff)**

bcftools only calls SNP/indel, so SV/CNV are added with Delly and folded into the model:

- **SVs** — `scripts/delta/call_scn_sv.sh` (Delly 2.x `sr`) → merge with the short VCF → rebuild.
  MM26A: **9,410 PASS SVs** (DEL 55% / INV 16% / BND 15% / DUP 11% / INS 2%). The rebuilt
  `mut_model` has a **populated `sv_model`** (per-base rate 4.3×10⁻⁵; length log-normals DEL ~6.6 kb,
  DUP ~16 kb, INV ~60 kb). INS/TEs are under-called (short reads can't span the repeats they live
  in) — the empirical case for long reads. (`type_probabilities` differ from raw counts because
  `SvModel::fit_from_observations` drops SVs < `MIN_SV_LENGTH_BP`=50.)
- **CNV** — `scripts/delta/call_scn_cnv.sh` (Delly `cnv`; auto-builds a `dicey` mappability map).
  MM26A: **433 PASS CNV segments**. ⚠️ **Format gap:** `gen-mut-model` reads copy number from
  `INFO/CN` (`variants.rs` `parse_sv_info`) but Delly writes `FORMAT/RDCN` — so the merge must
  **lift `RDCN`→`INFO/CN`** (and strip `RDCN`, which also resolves a `bcftools concat` header
  type-conflict). After the lift, `cnv_copy_number_distribution` populates (real, ≠ the bundled
  default): MM26 ~92% losses (CN0 0.42 / CN1 0.50), ~8% gains.
- **Rendering (loop closed)** — `scripts/delta/simulate_scn_sv.sh` runs gen-reads with
  `sv_rate_scale=1.0` and confirms SVs **and CNVs render** into the golden VCF. MM26 sim:
  6,701 SVs incl. **465 `SVTYPE=CNV`** (~model weight). So: real SCN SV/CNV → model → simulated output.

## Contrast — MM26 (virulent) vs PA3/BD3 (avirulent) — **DONE ✅**

Same pipeline on the avirulent line (`ACC=GCA_040805705.1 SRR=SRR27329600`): 10,495 PASS SVs,
320 CNV segments, full four-class model built. Model comparison:

| statistic | MM26 (virulent) | PA3/BD3 (avirulent) |
|---|---|---|
| mutation_rate | 2.44×10⁻³ | 3.71×10⁻³ |
| homozygous_frequency | 0.126 | 0.298 |
| SV type mix | Del-heavy (Del .32, Cnv .065) | rearrangement-heavy (Inv .27, Bnd .26) |
| CNV spectrum | ~92% losses, ~8% gains | ~78% losses, ~22% gains |

The two models are genuinely distinct — PA3/BD3 shows higher diversity, ~2.4× more near-fixation,
a more rearrangement-heavy SV profile, and ~3× more CNV gains. **But this is confounded by the
strain-fitting design** (same-strain vs cross-population reference pairing), so it is **not** a
clean virulence signal. What the contrast establishes: the pipeline reproduces per-strain
four-class models on a second real line, and the results differ. A causal virulence claim needs
a common-reference alignment.

## Phase 2 — reproduce the pool allele-frequency spectrum — **DONE ✅ (v1.21.0)**

Full design: `docs/scn_phase2_af_design.md` (implemented). Reproductive replay: feed the pool's
observed per-site AFs, emit reads that match them, validate sim-AF ≈ real-AF. The enabling change —
an optional per-variant `allele_fraction` on `Variant`, parsed from the input VCF (`INFO/AF`, else
`FORMAT/AD`) — shipped as **v1.21.0 (#398 / #399)**. Validated on real MM26 Pool-seq: simulated vs
real per-site allele frequency reached **Pearson r = 0.983 at 150× coverage** (0.936 at 50×),
unbiased, with the residual tracking the binomial coverage-noise floor. Harness:
`scripts/delta/run_scn_af_validation.sh` + `scripts/delta/scn_af_compare.py`.

## Merged tooling

| PR | what |
|---|---|
| #383 / #384 / #385 | `stage_scn.sh` + robust ENA download + `set -u` hotfix |
| #388 / #390 / #391 | `call_scn_sv.sh` (Delly SV) + static-binary/Boost override + Delly 2.x `call`→`sr` |
| #392 | `simulate_scn_sv.sh` (SV/CNV render) + `call_scn_cnv.sh` (CNV) + this status doc |
| #382 | cancer-pipeline re-run staleness fix (adjacent, found during this track) |
| #386 / #389 | Phase-2 AF design + long-read epic scope docs |
| #393 / #394 | **RNG out-of-range fix** (`mash()` u32 mask) → released as **v1.20.1** |
| #396 | conda recipe → 1.20.1 |

## Operational notes (Delta)

- **`$SCRATCH` must be the per-user path** `/scratch/bhrd/jallen17` (not `/scratch/bhrd`) — the
  eidolon binary lives at `$SCRATCH/cargo-target/...`. Share data with João via perms/ACL, not by
  repointing `$SCRATCH`.
- Tools: NCBI `datasets` (assemblies), Delly **static binary** `DELLY=$SCRATCH/bin/delly` (bioconda
  Delly fails a Boost load; Delta's `module load boost` is the wrong version), `dicey` (CNV map).
- Pass strain overrides explicitly to jobs — e.g. `sbatch --export=ALL,ACC=…,SRR=… …` — and
  **check the job banner** (`ref=`/`bam=`) before trusting output; a strain mix builds silently-wrong
  models (a PA3-ref + MM26-BAM build failed exactly this way).
- **Verify the statistic, not the PASS:** confirm non-zero variant/SV/CNV counts and that
  `cnv_copy_number_distribution` populates — a green build can still be hollow.
- Rebuild the binary (`setup.sh`) before whole-genome gen-reads so the RNG fix (v1.20.1) is in play.

## Open / next steps

The planned **validation** track is complete (per-strain sim-vs-truth fidelity, build → simulate →
render, on both lines). Remaining items split into tool/validation work that can proceed and
**biology extensions parked pending colleague feedback**.

**Validation / tooling (can proceed):**

1. **Phase 2 — DONE (v1.21.0).** Allele-frequency replay shipped and validated (r = 0.983 at
   150×); see the Phase 2 section above. No longer pending.
2. **Harness (issue #395)** — the correctness fixes (silent-contig-drop + pipefail guards) merged
   in **#409**; #395 remains open only for the de-dup refactor (hoist shared `stage_scn`/`stage_soy`
   helpers into `lib_report.sh`) and baking the `RDCN`→`INFO/CN` lift into `call_scn_cnv.sh`.

**Biology extensions — PARKED pending feedback from João / colleagues:**

- **Common-reference variant comparison** — align both pools to one reference to remove the
  strain-fitting confound. This is *not* a validation gap: per-strain sim-vs-truth fidelity stands
  on its own regardless of reference. It only matters for turning the strain contrast into a
  defensible virulent-vs-avirulent *biological* claim, which is a project riding on top of the
  validated tool. Do not start without colleague sign-off that the biological claim is wanted.
- **CNV biology sanity** — the ~78–92% loss-dominated CNV spectra on same/near-strain alignments
  are surprising; check with João (likely partly low-mappability/repeat regions).
- **Long reads** — `delly lr` on the PacBio/ONT data (PRJNA852516) to recover TEs short reads miss.

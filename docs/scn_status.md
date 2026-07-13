# SCN track — status

Consolidated status for the soybean-cyst-nematode (SCN, *Heterodera glycines*) real-data
work under the realism epic (#311). Data from **João Gomes Viana** (Matt Hudson lab).
Keep this current as the SCN work progresses. Last updated: **2026-07-12**.

Goal: use real SCN Pool-seq data to (Phase 1) build genuine per-strain rneat models incl.
**structural variants**, and (Phase 2) reproduce the population allele-frequency spectrum in
simulated reads.

## Data & accessions (BioProject PRJNA1055977)

**Strain reference assemblies** (Walden et al. 2025, doi:10.1186/s12864-025-12493-x;
download via NCBI `datasets`):

| strain | assembly | phenotype | genome |
|---|---|---|---|
| MM26 | `GCA_040805935.1` | virulent | 144.9 Mb |
| PA3  | `GCA_040805705.1` | avirulent ("picky eater") | 123.7 Mb |

**Population Pool-seq reads** (Illumina NextSeq 2000, paired, pools of 150 virgin females):

| pool | SRA run | aligned to | R1/R2 md5 (ENA) |
|---|---|---|---|
| MM26A   | `SRR27329602` | MM26 (`GCA_040805935.1`) | `8e233e2d…` / `6843e5e2…` |
| MM-BD3A | `SRR27329600` | PA3 (`GCA_040805705.1`) | `41e8850a…` / `7b9a8b9f…` |

**Reference strategy:** SCN strains differ in virulence, so each pool is aligned to its
best-fit strain reference to mitigate reference bias (MM26→MM26, BD3→PA3), per João. Other
assemblies exist (TN10 by Masonbrink; X12) but PA3/MM26 fit their categories best.

**Pool-seq caveat:** these are *pooled* samples (150 individuals → allele **frequencies**,
not diploid genotypes). Use allele **depths (AD)**, not `GT`, for anything AF-related.

## Phase 1 — build real per-strain models — **DONE ✅**

`scripts/delta/stage_scn.sh` (SCN analogue of `stage_soy.sh`) stages a self-consistent
reference + BAM + VCF (align pool reads → call VCF *from* that BAM), then
`model_builders.sbatch` builds all five models and round-trips them through gen-reads.

Both strains built and **verified by the actual statistics** (not just the PASS):

| statistic | MM26 / MM26A (job 20039284) | PA3 / MM-BD3A (job 20046157) | note |
|---|---|---|---|
| mutation_rate | 2.44×10⁻³ | 3.68×10⁻³ | pooled diversity vs own ref |
| homozygous_frequency | 0.126 | 0.298 | finite/real (pooled near-fixation fraction) |
| error_rate | 4.35×10⁻³ | 3.99×10⁻³ | **consistent** — sequencing property ✓ |
| frag length | Normal(516, 154) | Normal(448, 178) | real inserts |
| called variants | 475,587 | 611,155 | vs *own* reference |
| model rate ÷ naive (vars/genome) | 74.3% | 74.5% | **identical ratio → builder is consistent** ✓ |
| round-trip reads | 219,404 pairs | 243,464 pairs | models usable |

Artifacts archived to `/projects/bhrd/jallen17/rneat-access-results/modelbuild/job_{20039284,20046157}`.

**Interpretation:** both builds are sound (error_rate consistent; identical rate/naive ratio).
BD3 shows higher diversity + near-fixation, but this is **confounded** — each pool is measured
against a *different* reference, so raw variant counts aren't a clean cross-strain comparison.
A true virulence-linked variant comparison needs both pools on a **common** reference (a
separate, orthogonal task, not yet done).

## Phase 1b — Structural variants — **DONE ✅ (the "complex organism" payoff)**

The short-variant staging yields only SNPs + small indels, so the model's `sv_model` was empty
and rneat's flagship SV machinery was untested on SCN. `scripts/delta/call_scn_sv.sh` fixes that:
Delly (`sr` on 2.x) on the staged BAM → merge with the short VCF → rebuild the model.

MM26A (Delly v2.3.1, PASS SVs): **9,410 SVs** — DEL 5,169 (55%), INV 1,541 (16%), BND 1,449 (15%),
DUP 1,081 (11%), INS 170 (2%). The rebuilt `mut_model` (job 20101361) has a **populated
`sv_model`**: per-base SV rate 4.28×10⁻⁵; type mix Del .34 / Inv .25 / Bnd .23 / Dup .17; length
log-normals DEL ~6.6 kb, DUP ~16 kb, INV ~60 kb (BND lengthless).

- **Biology matches predictions:** DEL-dominant, **DUP present** (SCN's copy-number/duplication
  character), and **INS/TEs suppressed** (2%) — the empirical demonstration that short reads
  can't span the repeats where TEs live (→ the long-read motivation).
- **Model type mix ≠ raw counts** (DEL .34 vs .56) is *expected*: `SvModel::fit_from_observations`
  drops SVs < `MIN_SV_LENGTH_BP` = 50 (`common/src/structs/sv_model.rs:32`), filtering short DELs
  and all 170 small INS. Not a bug.
- **`cnv_copy_number_distribution` is empty** because `delly sr` gives DEL/DUP/INV/BND, not
  explicit copy-number — that's what `delly cnv` (below) adds.

Follow-ups (this PR):
- **`scripts/delta/simulate_scn_sv.sh`** — simulate reads *from* the SV-inclusive model with
  `sv_rate_scale=1.0` and verify SVs actually **render** in the golden VCF (closes the loop).
- **`scripts/delta/call_scn_cnv.sh`** — `delly cnv` read-depth copy-number (the biologically
  central SCN signal). Needs a mappability map (auto-built with `dicey`; `dicey` CLI is
  version-sensitive — pass `MAP=<prebuilt>` to skip). Merge its calls to populate
  `cnv_copy_number_distribution`.

Delly 2.x note: `call` → `sr` (short-read SV); `cnv` (read-depth CNV); **`lr` = native long-read
SV discovery** — directly usable on the PRJNA852516 PacBio/ONT data (no rneat long-read
generation needed to *call* SVs on real long reads).

## Phase 2 — reproduce the pool allele-frequency spectrum — **designed**

Full design: **`docs/scn_phase2_af_design.md`** (PR #386). Summary:

- **Scope: reproductive replay** (decided 2026-07-11) — feed rneat the pool's observed
  per-site AFs, emit reads that match them, validate sim-AF ≈ real-AF (João's question). NOT a
  generative population/SFS model (a future superset that reuses the same primitive).
- **Core constraint:** rneat can't do this today — `genotype_fraction` (`runner.rs:2368`) emits
  only `{1/ploidy, 1.0}`, `Genotype` is a binary Het/Hom enum, and `ploidy` is global. Cancer
  `purity` is a global scalar, so it doesn't help either.
- **Enabling change:** an optional per-variant `allele_fraction: Option<f64>` on `Variant`,
  read from the input VCF (from AD) and honored by the read-gen coverage path. Bounded and
  general (useful for any population/somatic-AF simulation).
- **Validation:** a dependency-free `scripts/delta/scn_af_compare.py` (mirrors
  `sbs96_compare.py`) — per-site correlation + per-decile MAE, target ≥ 0.95.
- **Step 0 (no code, do first):** run at `ploidy=2` and measure how far the current binary
  model is from the real spectrum — quantifies the gap, motivates the change, and builds the
  AF-extraction + comparison tooling before any Rust.

## Merged tooling

| PR | what |
|---|---|
| #383 | `stage_scn.sh` — SCN staging (Phase 1) |
| #384 | robust, ENA-verified read download (retry + size/md5) |
| #385 | `set -u` hotfix (unbound `mate` in `fetch_read`) |
| #382 | cancer-pipeline re-run staleness fix (adjacent, found during this track) |
| #386 | Phase 2 (allele-frequency) design doc |
| #388 | `call_scn_sv.sh` — Delly structural-variant calling → SV model |
| #389 | long-read epic scope doc |
| #390 / #391 | Delly static-binary override + Boost fix; Delly 2.x `call`→`sr` |
| #387 | this status doc |
| *(this PR)* | `simulate_scn_sv.sh` (SV render check) + `call_scn_cnv.sh` (CNV) |

## Operational notes (Delta)

- **`$SCRATCH` must be the per-user path** `/scratch/bhrd/jallen17` (not the project root
  `/scratch/bhrd`) — the rneat binary lives at `$SCRATCH/cargo-target/...`. Sharing data with
  João = perms/ACL on the data dir, **not** repointing `$SCRATCH`.
- `stage_scn.sh` needs NCBI `datasets` on PATH (for the assembly) and writes to
  `$SCRATCH/neat_data/scn`. Downloads ~45 GB/pool then subsamples to ~30× (`MAX_PAIRS`).
- Always **verify the staging summary's variant count is non-zero** before building models —
  a contig-naming mismatch silently produces hollow-but-valid models.

## Next steps

1. **Run the SV follow-ups (in progress):** `simulate_scn_sv.sh` (confirm SVs render from the
   trained model) + `call_scn_cnv.sh` (populate the CNV distribution).
2. **PA3 / MM-BD3 SV contrast** — run the SV + CNV pipeline on the avirulent line and compare SV
   spectra virulent-vs-avirulent (reference-bias caveat applies).
3. **Long reads:** `delly lr` on the PRJNA852516 PacBio/ONT data to recover the TEs short reads
   under-call (no rneat long-read *generation* needed). See `docs/longread_epic_scope.md`.
4. **Phase 2 Step 0** — AF extraction + `scn_af_compare.py` + ploidy=2 gap baseline (no code),
   then the `allele_fraction` enabling change (issue under #311).

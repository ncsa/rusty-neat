# SCN track — status

Consolidated status for the soybean-cyst-nematode (SCN, *Heterodera glycines*) real-data
work under the realism epic (#311). Data from **João Gomes Viana** (Matt Hudson lab).
Keep this current as the SCN work progresses. Last updated: **2026-07-11**.

Goal: use real SCN Pool-seq data to (Phase 1) build genuine per-strain rneat models, and
(Phase 2) reproduce the population allele-frequency spectrum in simulated reads.

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
| #386 | Phase 2 design doc *(open)* |

## Operational notes (Delta)

- **`$SCRATCH` must be the per-user path** `/scratch/bhrd/jallen17` (not the project root
  `/scratch/bhrd`) — the rneat binary lives at `$SCRATCH/cargo-target/...`. Sharing data with
  João = perms/ACL on the data dir, **not** repointing `$SCRATCH`.
- `stage_scn.sh` needs NCBI `datasets` on PATH (for the assembly) and writes to
  `$SCRATCH/neat_data/scn`. Downloads ~45 GB/pool then subsamples to ~30× (`MAX_PAIRS`).
- Always **verify the staging summary's variant count is non-zero** before building models —
  a contig-naming mismatch silently produces hollow-but-valid models.

## Next steps

1. **Phase 2 Step 0** — AF extraction + `scn_af_compare.py` + ploidy=2 gap baseline (no code).
2. Trace the input-VCF → read-stream path; implement the `allele_fraction` enabling change
   (issue under #311); re-run + validate; A/B MM26 vs PA3/BD3.
3. *(Orthogonal, optional)* common-reference variant comparison for a clean virulence signal.

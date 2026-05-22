# Release Checklist

Pre-publish gate for `rneat` releases. Tick every box before tagging.

## 0. Branch hygiene

- [ ] All PRs targeted at the release are merged into `develop`.
- [ ] `develop` → `main` PR is open with the release notes pasted into the body.
- [ ] CHANGELOG.md has a top entry for the new version with bullets describing
      behavioral and CLI-visible changes (not just refactors).
- [ ] `Cargo.toml` `workspace.package.version` matches the CHANGELOG header.
- [ ] `Cargo.lock` has been refreshed (`cargo update -p common -p rneat`).

## 1. Workspace-local checks

- [ ] `cargo fmt --check` exits 0.
- [ ] `cargo clippy --workspace --all-targets` reports no deny-level errors.
- [ ] `cargo test --workspace --release` is green, including:
  - `model_parity::*` (model-output regression — load-bearing for this gate)
  - `determinism::*` (same-seed → same multiset)
  - `pipeline_e2e::*` (gen-reads end-to-end with default models)
  - `cli_smoke::*` (binary surface)
- [ ] If `model_parity::*` fails, **do not** re-bless the baselines unless the
      diff is intentional. The point of the test is to catch silent drift.

## 2. Full-size shakeout

For the items below you need a working directory with:

| File | Source |
|------|--------|
| `hg38_chr22.fa`     | GRCh38 chr22, single-contig FASTA. ~50 MB. |
| `hg38_chr22.bam`    | Aligned reads against `hg38_chr22.fa`, MD tags present (`samtools calmd` if not). |
| `hg38_chr22.vcf.gz` | Multi-sample or single-sample VCF on chr22, `GT` populated. |

Then run:

```bash
DATA_DIR=/scratch/rneat-shakeout ./scripts/release_shakeout.sh
```

The script will:

1. Re-run `cargo fmt`, `clippy`, full test suite in `--release`.
2. Train `mut`, `frag-length`, `gc-bias`, and unified `gen-bam-models` on chr22.
3. Assert byte-equality of unified vs standalone frag-length and gc-bias outputs.
4. Run `gen-reads` end-to-end with the trained models.
5. Report peak RSS + elapsed time for each step under `$WORK/*.time`.

- [ ] Shakeout exits 0.
- [ ] Peak RSS and elapsed time for each step are within ±20% of the previous
      release's numbers (eyeballed comparison against the prior shakeout log,
      which lives in `$DATA_DIR/shakeout-*/`).
- [ ] Output FASTQs from step 5 are well-formed (the script already asserts
      non-empty, but spot-check with `seqkit stats` if you suspect anything).

## 3. Cross-version determinism (optional but recommended)

If the release touches anything in the RNG, fragment generation, mutation
sampling, or output assembly paths:

- [ ] Check out the previous release tag (e.g. `git checkout v1.6.0`).
- [ ] Run `gen-reads` with a fixed seed on H1N1 (in-tree, ~14 KB reference).
- [ ] Save the output FASTQ multiset hash:
      `zcat *.fastq.gz | sort | sha256sum`.
- [ ] Check out the new release branch and re-run with the same seed.
- [ ] Compare hashes. If they diverge and the release was not intended to
      change output, **stop** and investigate before proceeding.

## 4. Publication

- [ ] Squash-merge `develop` → `main` (or merge commit, per project convention).
- [ ] Tag: `git tag -a vX.Y.Z -m 'rneat vX.Y.Z' && git push origin vX.Y.Z`.
- [ ] Create GitHub Release on the tag with the CHANGELOG section pasted in.
- [ ] Open a follow-up PR on `develop` bumping to `vX.Y.(Z+1)-dev` (or next
      planned minor), so subsequent work doesn't clobber the released version.

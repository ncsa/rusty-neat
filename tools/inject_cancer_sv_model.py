#!/usr/bin/env python3
# ============================================================================
# DEPRECATED as of v1.14.0 (#218). This script injects HAND-ROLLED literature
# values. The bundled tumor model's sv_model is now DATA-DERIVED from the PCAWG
# consensus SV/CNV callsets (+ gnomAD-SV for INS length) via this pipeline:
#
#   tools/fetch_pcawg_sv_corpus.sh        # download open-tier PCAWG SV + CNV
#   tools/build_pcawg_sv_vcf.py           # -> symbolic-SV VCF (+ counts sidecar)
#   eidolon gen-mut-model -c <fit.yml>      # fit length / CN distributions
#   tools/normalize_pcawg_sv_model.py     # per-tumor rate/type-mix + splice in
#
# The data refit corrected this script's heuristics substantially — e.g. INV
# 0.08 -> 0.134 and BND 0.35 -> 0.232 (this file's "~35% BND" was a misreading
# of the PCAWG breakdown; inversions are represented as paired junctions, which
# the data pipeline de-duplicates). Kept only for historical reference and for
# bootstrapping a null-sv_model model from scratch; do NOT use it to regenerate
# the bundled model.
# ============================================================================
#
# Inject a literature-derived cancer SvModel component into an eidolon
# MutationModel JSON (gzipped or plain). The bundled COSMIC tumor model
# (tools/cosmic_v104_pancancer_model.json.gz) is trained from COSMIC
# GenomeScreensMutant, which is SNV/indel-only — so its top-level
# sv_model field is null and the tumor pass of cancer_simulate.sh emits
# zero symbolic SVs, even with --sv-rate-scale > 0.
#
# This script populates that field with pan-cancer somatic SV defaults
# derived from PCAWG (Li et al., "Patterns of somatic structural
# variation in human cancer genomes", Nature 578, 112-121, 2020).
# The parameters are heuristic literature values, not a refit from raw
# data — see issue #218 for the v1.13 plan to retrain from a real
# cancer-SV corpus (PCAWG consensus calls or COSMIC StructVar).
#
# Usage:
#   tools/inject_cancer_sv_model.py <model.json.gz>          # in-place
#   tools/inject_cancer_sv_model.py <in.json.gz> <out.json.gz>
#
# The script refuses to overwrite a non-null sv_model unless --force.

import argparse
import gzip
import json
import sys
from pathlib import Path

# ─── Pan-cancer somatic SV parameters (PCAWG, Li et al. 2020) ──────────────
#
# Per-base rate: median ~110 somatic SVs per tumor across the PCAWG
# pan-cancer cohort (n=2,658), normalized over ~2.9 Gb of GRCh38 primary
# autosomes+sex chromosomes. gen-reads multiplies this by contig_len ×
# sv_rate_scale, so chr22 (~50 Mb) at scale=1.0 expects ~2 somatic SVs
# per pass — a sensible scale for a synthetic tumor.
PER_BASE_RATE = 3.8e-8

# Type breakdown: pan-cancer averages from PCAWG Fig 1c / Table 1.
# Cancer over-represents translocations vs germline (gnomAD: ~19% BND;
# PCAWG: ~35% BND) and under-represents small deletions in favor of
# large oncogene/tumor-suppressor deletions.
TYPE_PROBABILITIES = {
    "Del": 0.35,
    "Dup": 0.17,
    "Bnd": 0.35,
    "Inv": 0.08,
    "Ins": 0.04,
    "Cnv": 0.01,
}

# Length log-normals: medians from PCAWG, sigmas reflect the heavy
# right tail (multi-Mb deletions, chromosome-scale duplications). Quick
# gut-check on the medians (= exp(mu)):
#   DEL: ~5.0 kb   (germline ~700 bp; cancer DELs are larger)
#   DUP: ~30 kb    (tandem duplications dominate)
#   INV: ~22 kb    (typical somatic INV)
#   CNV: ~1 Mb     (large copy-number segments)
#   INS: ~300 bp   (MEI-typical, matching germline default)
#   BND: 1 bp      (point notation, length unused)
LENGTH_LOG_NORMAL = {
    "Del": [8.52, 2.0],
    "Dup": [10.31, 2.5],
    "Inv": [9.90, 2.3],
    "Cnv": [13.82, 1.5],
    "Ins": [5.7, 1.0],
    "Bnd": [0.0, 0.0],
}

# CNV copy-number distribution. Cancer copy-number profiles are bimodal:
# losses (CN=0, CN=1) from tumor-suppressor inactivation and gains
# (CN=3,4,5+) from oncogene amplification. Compressing the high-level
# amplification tail into CN=5 keeps the table small without changing
# the simulated coverage signal much. Excludes CN=2 (normal diploid).
CNV_COPY_NUMBER_DISTRIBUTION = {
    "0": 0.10,
    "1": 0.20,
    "3": 0.30,
    "4": 0.20,
    "5": 0.20,
}

# Homozygous fraction: somatic SVs are predominantly heterozygous in
# the tumor genome unless followed by loss-of-heterozygosity. Lower
# than germline (gnomAD-derived default 0.20).
HOMOZYGOUS_FREQUENCY = 0.10


def cancer_sv_model() -> dict:
    return {
        "per_base_rate": PER_BASE_RATE,
        "type_probabilities": TYPE_PROBABILITIES,
        "length_log_normal": LENGTH_LOG_NORMAL,
        "cnv_copy_number_distribution": CNV_COPY_NUMBER_DISTRIBUTION,
        "homozygous_frequency": HOMOZYGOUS_FREQUENCY,
    }


def load(path: Path) -> dict:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as f:
        return json.load(f)


def save(path: Path, data: dict) -> None:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "wt") as f:
        json.dump(data, f)


def main() -> int:
    p = argparse.ArgumentParser(
        description="Inject pan-cancer SvModel defaults (PCAWG) into a "
        "MutationModel JSON."
    )
    p.add_argument("input", type=Path, help="Input MutationModel JSON (.json or .json.gz)")
    p.add_argument(
        "output",
        type=Path,
        nargs="?",
        help="Output path (default: overwrite input)",
    )
    p.add_argument(
        "--force",
        action="store_true",
        help="Overwrite an existing non-null sv_model field",
    )
    args = p.parse_args()

    if not args.input.exists():
        print(f"error: {args.input} not found", file=sys.stderr)
        return 2
    out_path = args.output or args.input

    model = load(args.input)
    if "sv_model" not in model:
        print(
            f"warning: {args.input} has no sv_model field — adding one. "
            "If this model predates the sv_model addition, gen-reads will "
            "still consume it (the field is #[serde(default)] = None).",
            file=sys.stderr,
        )
    if model.get("sv_model") is not None and not args.force:
        print(
            f"error: {args.input} already has a non-null sv_model. "
            "Pass --force to overwrite.",
            file=sys.stderr,
        )
        return 3

    sv = cancer_sv_model()
    total = sum(sv["type_probabilities"].values())
    if abs(total - 1.0) > 1e-9:
        print(
            f"internal error: type_probabilities sum to {total}, not 1.0",
            file=sys.stderr,
        )
        return 4
    cn_total = sum(sv["cnv_copy_number_distribution"].values())
    if abs(cn_total - 1.0) > 1e-9:
        print(
            f"internal error: cnv_copy_number_distribution sums to "
            f"{cn_total}, not 1.0",
            file=sys.stderr,
        )
        return 4

    model["sv_model"] = sv
    save(out_path, model)
    print(
        f"wrote {out_path} with PCAWG-derived sv_model "
        f"(per_base_rate={sv['per_base_rate']:.2e}, "
        f"BND={sv['type_probabilities']['Bnd']:.0%}, "
        f"DEL={sv['type_probabilities']['Del']:.0%})"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())

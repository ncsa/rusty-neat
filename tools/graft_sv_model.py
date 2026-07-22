#!/usr/bin/env python3
"""Graft an sv_model from one eidolon MutationModel JSON onto another.

Completes the per-tissue tumor models (#202): a model trained by gen-mut-model
from a per-tissue COSMIC SNV/indel corpus carries the per-tissue SNP/indel
spectrum (mutation_rate, variant_dist, statistical_models) but no sv_model. The
per-tissue SV models from #237 carry the right per-tissue sv_model but a
pan-cancer SNP/indel base. This takes the SNP/indel half from --snv-model and the
sv_model from --sv-source, producing a fully per-tissue model.

    tools/graft_sv_model.py \
        --snv-model cosmic_breast_model.json.gz \
        --sv-source tools/cosmic_pancancer_sv_BRCA.json.gz \
        --out       tools/cosmic_per_tissue_BRCA.json.gz
"""
import argparse
import gzip
import json
import sys


def load_gz_json(path):
    with gzip.open(path, "rt") as fh:
        return json.load(fh)


def write_gz_json(path, obj):
    with gzip.open(path, "wt") as fh:
        json.dump(obj, fh)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--snv-model", required=True,
                    help="model providing the SNP/indel half (the new per-tissue base)")
    ap.add_argument("--sv-source", required=True,
                    help="model providing the sv_model to graft on")
    ap.add_argument("--out", required=True, help="output model JSON.gz")
    args = ap.parse_args()

    base = load_gz_json(args.snv_model)
    src = load_gz_json(args.sv_source)

    sv = src.get("sv_model")
    if sv is None:
        sys.exit(f"--sv-source {args.sv_source} has no sv_model to graft")

    base["sv_model"] = sv
    write_gz_json(args.out, base)

    # Brief provenance so the build log shows what came from where.
    rate = base.get("mutation_rate")
    types = sorted((sv.get("type_probabilities") or {}).keys())
    print(f"  grafted sv_model ({', '.join(types) or 'no types'}) onto "
          f"SNP/indel base (mutation_rate={rate}) -> {args.out}")


if __name__ == "__main__":
    main()

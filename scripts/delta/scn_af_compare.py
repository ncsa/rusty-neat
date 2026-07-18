#!/usr/bin/env python3
"""Compare per-site alt-allele fractions between two VCFs.

Purpose: validate SCN Phase 2 (issue #398 / realism epic #311). Feed the REAL
pool's observed per-site AFs (A, the input to gen-reads) and the SIMULATED golden
VCF gen-reads emitted (B). A high correlation + low per-site error means the
per-variant `allele_fraction` feature replayed the observed spectrum faithfully —
i.e. simulated AF ~= real AF, which is exactly Joao's question.

Sites are matched by (CHROM, POS, REF, ALT). For each side the AF is taken from,
in order: FORMAT/AF (the golden VCF's measured field), INFO/AF, or FORMAT/AD
(alt_depth / total_depth). Use AD-derived AF for pooled truth sets, where GT is
not diploid-meaningful.

Dependency-free (Python 3 stdlib only): no numpy/scipy/pysam. Pearson and
Spearman are computed directly.

Usage:
  python3 scripts/delta/scn_af_compare.py \
      --truth pool.af.vcf.gz --sim af_run.vcf.gz [--min-depth 20]
"""
import argparse
import gzip
import math
import sys
from collections import OrderedDict


def _open(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def _field_index(fmt, key):
    parts = fmt.split(":")
    return parts.index(key) if key in parts else None


def _af_from_record(info, fmt, sample):
    """Extract an alt fraction from one VCF record's INFO/FORMAT/SAMPLE.

    Order: FORMAT/AF, then INFO/AF, then FORMAT/AD (sum alt / total). Returns
    (af, depth) — depth is total AD when available (for min-depth gating), else
    None. Returns (None, None) when no fraction can be derived.
    """
    # FORMAT/AF
    if fmt and sample:
        af_i = _field_index(fmt, "AF")
        if af_i is not None:
            vals = sample.split(":")
            if af_i < len(vals):
                try:
                    return float(vals[af_i].split(",")[0]), _depth_from_ad(fmt, sample)
                except ValueError:
                    pass
    # INFO/AF
    for kv in info.split(";"):
        if kv.startswith("AF="):
            try:
                return float(kv[3:].split(",")[0]), _depth_from_ad(fmt, sample)
            except ValueError:
                pass
    # FORMAT/AD
    if fmt and sample:
        ad_i = _field_index(fmt, "AD")
        if ad_i is not None:
            vals = sample.split(":")
            if ad_i < len(vals):
                try:
                    counts = [float(x) for x in vals[ad_i].split(",")]
                    total = sum(counts)
                    if len(counts) >= 2 and total > 0:
                        return sum(counts[1:]) / total, total
                except ValueError:
                    pass
    return None, None


def _depth_from_ad(fmt, sample):
    if not fmt or not sample:
        return None
    ad_i = _field_index(fmt, "AD")
    if ad_i is None:
        return None
    vals = sample.split(":")
    if ad_i >= len(vals):
        return None
    try:
        return sum(float(x) for x in vals[ad_i].split(","))
    except ValueError:
        return None


def load_af(path):
    """Map (chrom, pos, ref, alt) -> (af, depth) for every usable record."""
    out = OrderedDict()
    with _open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            chrom, pos, ref, alt = f[0], f[1], f[3], f[4]
            if alt.startswith("<") or "[" in alt or "]" in alt:
                continue  # symbolic / breakend — not a literal-fraction site
            info = f[7] if len(f) > 7 else "."
            fmt = f[8] if len(f) > 8 else ""
            sample = f[9] if len(f) > 9 else ""
            af, depth = _af_from_record(info, fmt, sample)
            if af is not None:
                out[(chrom, pos, ref, alt)] = (af, depth)
    return out


def pearson(xs, ys):
    n = len(xs)
    if n < 2:
        return float("nan")
    mx, my = sum(xs) / n, sum(ys) / n
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    if sxx == 0 or syy == 0:
        return float("nan")
    return sxy / math.sqrt(sxx * syy)


def _ranks(vals):
    order = sorted(range(len(vals)), key=lambda i: vals[i])
    ranks = [0.0] * len(vals)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and vals[order[j + 1]] == vals[order[i]]:
            j += 1
        avg = (i + j) / 2.0
        for k in range(i, j + 1):
            ranks[order[k]] = avg
        i = j + 1
    return ranks


def spearman(xs, ys):
    if len(xs) < 2:
        return float("nan")
    return pearson(_ranks(xs), _ranks(ys))


def main():
    ap = argparse.ArgumentParser(description="Per-site AF correlation of two VCFs (truth vs sim)")
    ap.add_argument("--truth", required=True, help="VCF A — real/input per-site AFs")
    ap.add_argument("--sim", required=True, help="VCF B — simulated golden VCF")
    ap.add_argument("--min-depth", type=float, default=0.0,
                    help="skip sites whose total AD is below this on either side (default 0)")
    args = ap.parse_args()

    a, b = load_af(args.truth), load_af(args.sim)
    shared = [k for k in a if k in b]
    only_a = len(a) - len(shared)
    only_b = len(b) - len(shared)

    xs, ys, gated = [], [], 0
    for k in shared:
        (af_a, dp_a), (af_b, dp_b) = a[k], b[k]
        if args.min_depth > 0:
            if (dp_a is not None and dp_a < args.min_depth) or (
                dp_b is not None and dp_b < args.min_depth
            ):
                gated += 1
                continue
        xs.append(af_a)
        ys.append(af_b)

    print(f"truth sites={len(a)}  sim sites={len(b)}  shared={len(shared)}"
          f"  (only-truth={only_a}, only-sim={only_b})")
    if args.min_depth > 0:
        print(f"min-depth {args.min_depth:g}: {gated} shared sites gated, {len(xs)} compared")
    if len(xs) < 2:
        sys.exit("fewer than 2 comparable sites — check inputs / --min-depth")

    diffs = [y - x for x, y in zip(xs, ys)]
    mae = sum(abs(d) for d in diffs) / len(diffs)
    rmse = math.sqrt(sum(d * d for d in diffs) / len(diffs))
    print(f"n={len(xs)}  Pearson r={pearson(xs, ys):.4f}  Spearman rho={spearman(xs, ys):.4f}")
    print(f"MAE={mae:.4f}  RMSE={rmse:.4f}  mean(sim-truth)={sum(diffs)/len(diffs):+.4f}")
    print("  target: r>=0.95 and per-bin MAE within AF-estimation noise at that coverage")

    print("per-AF-decile MAE (truth bin):")
    for lo in [i / 10 for i in range(10)]:
        hi = lo + 0.1
        bd = [abs(y - x) for x, y in zip(xs, ys) if (lo <= x < hi or (hi >= 1.0 and x == 1.0))]
        if bd:
            print(f"  [{lo:.1f},{hi:.1f})  n={len(bd):6d}  MAE={sum(bd)/len(bd):.4f}")


if __name__ == "__main__":
    main()

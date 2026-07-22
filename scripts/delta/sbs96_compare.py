#!/usr/bin/env python3
"""Compare the SBS-96 trinucleotide spectra of two somatic-SNV VCFs.

Purpose: validate that an eidolon cancer simulation reproduces the mutational
signature its tumor_model was trained on. Feed the REAL somatic truth set the
model learned from (A) and the SIMULATED somatic SNVs the sim produced (B); a
high cosine similarity means the trinucleotide spectrum — i.e. the signature —
carried through the build->simulate path, not just the overall mutation rate.

SBS-96 = the 96 pyrimidine-strand trinucleotide contexts (6 substitution types x
16 flanking-base combinations), the standard COSMIC/SigProfiler mutation-context
representation. This computes each VCF's 96-vector and their cosine similarity.

Dependency-free (Python 3 stdlib only): reads a plain or bgzipped/gzipped FASTA
for one contig into memory and looks up context by position — no pysam/faidx.

Usage:
  python3 scripts/delta/sbs96_compare.py --ref chr1.fa --chrom chr1 \
      --a real_somatic.vcf.gz --b simulated_somatic.vcf.gz
"""
import argparse
import gzip
import math
import sys
from collections import Counter

BASES = set("ACGT")
COMP = str.maketrans("ACGT", "TGCA")
PYR = {"C", "T"}
# 96 contexts in a fixed order (pyrimidine ref): 6 sub types x 16 flanks.
CONTEXTS = [
    f"{a}[{r}>{s}]{b}"
    for r in "CT"
    for s in "ACGT"
    if s != r
    for a in "ACGT"
    for b in "ACGT"
]


def _open(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def rc(s):
    return s.translate(COMP)[::-1]


def load_contig(fa, want):
    """Return the uppercase sequence of contig `want` from FASTA `fa`."""
    seq, on = [], False
    with _open(fa) as fh:
        for line in fh:
            if line.startswith(">"):
                on = line[1:].split()[0] == want
            elif on:
                seq.append(line.strip())
    if not seq:
        sys.exit(f"contig {want!r} not found in {fa}")
    return "".join(seq).upper()


def spectrum(vcf, seq, chrom):
    """Tally SBS-96 contexts for single-base substitutions on `chrom` in `vcf`."""
    ctr = Counter()
    skipped_refmismatch = 0
    with _open(vcf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 5 or f[0] != chrom:
                continue
            ref, alt = f[3].upper(), f[4].upper()
            if len(ref) != 1 or len(alt) != 1 or ref not in BASES or alt not in BASES or ref == alt:
                continue  # not a clean SNV (indel / multi-allelic / N)
            pos = int(f[1])  # 1-based
            if pos < 2 or pos > len(seq) - 1:
                continue
            tri = seq[pos - 2 : pos + 1]  # [5', ref, 3']
            if len(tri) != 3 or any(b not in BASES for b in tri):
                continue
            if tri[1] != ref:
                skipped_refmismatch += 1
                continue  # VCF ref disagrees with the reference base
            if ref not in PYR:  # collapse onto the pyrimidine strand
                tri, ref, alt = rc(tri), rc(ref), rc(alt)
            ctr[f"{tri[0]}[{ref}>{alt}]{tri[2]}"] += 1
    if skipped_refmismatch:
        print(f"  ({vcf}: {skipped_refmismatch} SNVs skipped — ref base != reference)", file=sys.stderr)
    return ctr


def cosine(u, v):
    du = math.sqrt(sum(x * x for x in u))
    dv = math.sqrt(sum(x * x for x in v))
    if du == 0 or dv == 0:
        return float("nan")
    return sum(a * b for a, b in zip(u, v)) / (du * dv)


def main():
    ap = argparse.ArgumentParser(description="SBS-96 spectrum cosine similarity of two somatic VCFs")
    ap.add_argument("--ref", required=True, help="reference FASTA (plain or .gz) for context lookup")
    ap.add_argument("--chrom", default="chr1", help="contig to compare (default: chr1)")
    ap.add_argument("--a", required=True, help="VCF A — the real/input somatic truth set")
    ap.add_argument("--b", required=True, help="VCF B — the simulated somatic SNVs")
    args = ap.parse_args()

    seq = load_contig(args.ref, args.chrom)
    ca, cb = spectrum(args.a, seq, args.chrom), spectrum(args.b, seq, args.chrom)
    na, nb = sum(ca.values()), sum(cb.values())
    if na == 0 or nb == 0:
        sys.exit(f"no SNVs tallied on {args.chrom} (A={na}, B={nb}) — check --chrom / inputs")

    va = [ca.get(k, 0) for k in CONTEXTS]
    vb = [cb.get(k, 0) for k in CONTEXTS]
    cos = cosine(va, vb)

    print(f"chrom={args.chrom}   A(real) SNVs={na}   B(sim) SNVs={nb}")
    print(f"SBS-96 cosine similarity(A,B) = {cos:.4f}")
    print("  >=0.95 spectra match (signature reproduced) | 0.85-0.95 close | <0.85 divergent")
    print("top contexts   A% / B%:")
    for k in sorted(CONTEXTS, key=lambda k: -(ca.get(k, 0) + cb.get(k, 0)))[:12]:
        print(f"  {k}   {100*ca.get(k,0)/na:5.2f}% / {100*cb.get(k,0)/nb:5.2f}%")


if __name__ == "__main__":
    main()

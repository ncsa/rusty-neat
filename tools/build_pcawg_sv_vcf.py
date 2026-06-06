#!/usr/bin/env python3
"""Convert the PCAWG open-tier SV/CNV corpus (+ gnomAD-SV INS) into a single
symbolic-SV VCF that `rneat gen-mut-model` can fit into an SvModel.

This is the conversion half of the #218 cancer-SvModel refit. Run
tools/fetch_pcawg_sv_corpus.sh first to download/extract the inputs.

Type sourcing (see #218 / docs/cancer_simulator.md):
  DEL, DUP, BND, INV   <- PCAWG consensus_sv BEDPE (somatic, per-donor)
  CNV                  <- PCAWG consensus_cnv CNA segments (somatic, focal,
                          non-neutral relative to per-sample ploidy)
  INS (length only)    <- gnomAD-SV germline callset (its one strongly
                          populated type; the same mobile-element sizes apply
                          germline vs somatic). INS *rate* is NOT taken from
                          gnomAD -- it stays a somatic literature estimate set
                          in the downstream normalization step.

Why a sidecar counts file:
  gen-mut-model computes per_base_rate = n_observations / reference_length and
  type_probabilities from raw record counts. Our records span DIFFERENT donor
  counts per source (SV ~2748 donors, CNA ~2778) and gnomAD INS is a germline
  site catalogue, not per-donor somatic counts. So the fit's per_base_rate and
  type_probabilities are NOT directly usable -- only its length / CN
  distributions are. The sidecar <out>.counts.json records per-type totals and
  donor counts so the normalization step can compute per-tumor type_probs and
  per_base_rate with the correct per-type denominators.

Output VCF record shapes (REF=N anchor; symbolic ALT):
  DEL: <DEL>  INFO=SVTYPE=DEL;END=<end>
  DUP: <DUP>  INFO=SVTYPE=DUP;END=<end>
  INV: <INV>  INFO=SVTYPE=INV;END=<end>
  BND: <BND>  INFO=SVTYPE=BND
  CNV: <CNV>  INFO=SVTYPE=CNV;END=<end>;CN=<total_cn>
  INS: <INS>  INFO=SVTYPE=INS;SVLEN=<len>

Usage:
  tools/build_pcawg_sv_vcf.py \
      --bedpe-dir      ~/code/data/pcawg/consensus_sv \
      --cna-dir        ~/code/data/pcawg/consensus_cnv \
      --purity-ploidy  ~/code/data/pcawg/consensus_cnv/consensus.20170217.purity.ploidy.txt.gz \
      --gnomad-ins-vcf ~/code/data/gnomad/gnomad_v4_sv.sites.vcf.gz \
      --out pcawg_sv_corpus.vcf.gz
"""
import argparse
import glob
import gzip
import json
import os
import re
import sys

MIN_SV_LENGTH_BP = 50  # mirrors common::structs::sv_model::MIN_SV_LENGTH_BP


def log(msg):
    print(msg, file=sys.stderr, flush=True)


def norm_chrom(c, chr_prefix):
    """Normalize a chromosome name to the chosen convention."""
    c = c.strip()
    if c.startswith("chr"):
        bare = c[3:]
    else:
        bare = c
    if bare == "MT":
        bare = "M"
    if not chr_prefix:
        return bare
    return "chr" + bare


def open_maybe_gz(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def aliquot_of(fp):
    """BEDPE/CNA per-donor files are named <aliquot_id>.<...>; return the UUID."""
    return os.path.basename(fp).split(".")[0]


def load_tissue_aliquots(sample_sheet, projects):
    """Set of tumor aliquot_ids whose dcc_project_code starts with any of
    `projects`. Matching is by prefix, so 'BRCA' picks up BRCA-US/UK/EU while an
    explicit group like ['SKCM-US', 'MELA-AU'] (skin melanoma) also works.
    Sourced from PCAWG's pcawg_sample_sheet.tsv (aliquot_id, dcc_project_code)."""
    prefixes = tuple(projects)
    allowed = set()
    with open_maybe_gz(sample_sheet) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            i_aliquot = header.index("aliquot_id")
            i_project = header.index("dcc_project_code")
        except ValueError:
            log("!! sample sheet missing aliquot_id / dcc_project_code columns")
            return allowed
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= max(i_aliquot, i_project):
                continue
            if f[i_project].startswith(prefixes):
                allowed.add(f[i_aliquot])
    log(f"   tissue filter: {len(allowed)} aliquots match {list(projects)}")
    return allowed


# ── PCAWG BEDPE → DEL/DUP/INV/BND ────────────────────────────────────────
def parse_bedpe(bedpe_dir, chr_prefix, records, counts, allowed=None):
    files = sorted(glob.glob(os.path.join(bedpe_dir, "**", "*.bedpe.gz"),
                             recursive=True))
    if not files:
        log(f"!! no *.bedpe.gz under {bedpe_dir}")
    donors = 0
    for fp in files:
        if allowed is not None and aliquot_of(fp) not in allowed:
            continue
        donors += 1
        with gzip.open(fp, "rt") as fh:
            for line in fh:
                if line.startswith("chrom1") or not line.strip():
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 11:
                    continue
                chrom1, start1, chrom2, start2, svclass = (
                    f[0], f[1], f[3], f[4], f[10])
                # Map svclass -> SvType
                if svclass == "DEL":
                    svt = "DEL"
                elif svclass == "DUP":
                    svt = "DUP"
                elif svclass in ("h2hINV", "INV"):
                    # A balanced inversion appears as TWO BEDPE rows: one
                    # head-to-head (h2hINV) + one tail-to-tail (t2tINV) junction
                    # (per-donor h2h/t2t counts are near-equal — they pair up).
                    # rneat models an inversion as ONE event, so count it once,
                    # off the h2hINV row, and skip t2tINV below to avoid the ~2x
                    # double-count. Both junctions bracket the same footprint, so
                    # the h2h span is a valid length observation.
                    svt = "INV"
                elif svclass == "t2tINV":
                    continue  # partner junction of h2hINV; already counted
                elif svclass == "TRA":
                    svt = "BND"
                else:
                    continue
                try:
                    s1 = int(start1)
                    s2 = int(start2)
                except ValueError:
                    continue
                c1 = norm_chrom(chrom1, chr_prefix)
                if svt == "BND":
                    # No length; one symbolic <BND> record at breakpoint 1.
                    records.append((c1, s1 + 1, "<BND>", "SVTYPE=BND", "0/1"))
                    counts["BND"] += 1
                    continue
                # Intrachromosomal span event. Require same chrom + valid span.
                if norm_chrom(chrom2, chr_prefix) != c1:
                    continue
                pos = min(s1, s2) + 1
                end = max(s1, s2)
                if end - pos + 1 < MIN_SV_LENGTH_BP:
                    continue
                records.append(
                    (c1, pos, f"<{svt}>", f"SVTYPE={svt};END={end}", "0/1"))
                counts[svt] += 1
    counts["_donors_sv"] = donors
    log(f"   BEDPE: {donors} donors -> "
        f"DEL={counts['DEL']} DUP={counts['DUP']} "
        f"INV={counts['INV']} BND={counts['BND']}")


# ── purity/ploidy table → per-sample neutral copy number ─────────────────
def load_ploidy(path):
    ploidy = {}
    with open_maybe_gz(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            i_name = header.index("samplename")
            i_ploidy = header.index("ploidy")
        except ValueError:
            log("!! purity/ploidy header missing samplename/ploidy columns")
            return ploidy
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= max(i_name, i_ploidy):
                continue
            try:
                ploidy[f[i_name]] = float(f[i_ploidy])
            except ValueError:
                continue
    log(f"   ploidy table: {len(ploidy)} samples")
    return ploidy


def sample_id_from_cna_path(fp):
    """CNA files are named <aliquot_id>.consensus.<date>.somatic.cna.txt."""
    base = os.path.basename(fp)
    return base.split(".")[0]


# ── PCAWG CNA segments → focal non-neutral CNV ───────────────────────────
def parse_cna(cna_dir, ploidy, chr_prefix, max_focal_bp, records, counts, allowed=None):
    files = sorted(glob.glob(os.path.join(cna_dir, "**", "*.somatic.cna.txt"),
                             recursive=True))
    if not files:
        files = sorted(glob.glob(os.path.join(cna_dir, "**", "*.cna.txt"),
                                 recursive=True))
    if not files:
        log(f"!! no *.cna.txt under {cna_dir}")
    donors = 0
    missing_ploidy = 0
    for fp in files:
        if allowed is not None and aliquot_of(fp) not in allowed:
            continue
        donors += 1
        sid = sample_id_from_cna_path(fp)
        sample_ploidy = ploidy.get(sid)
        if sample_ploidy is None:
            # Fall back to diploid neutral when the sample isn't in the table.
            neutral = 2
            missing_ploidy += 1
        else:
            neutral = int(round(sample_ploidy))
        with open(fp, "r") as fh:
            header = fh.readline()  # chromosome start end total_cn major_cn ...
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 4:
                    continue
                try:
                    start = int(f[1])
                    end = int(f[2])
                    total_cn = int(f[3])
                except ValueError:
                    continue
                if total_cn == neutral:
                    continue  # copy-number-neutral segment
                span = end - start + 1
                if span < MIN_SV_LENGTH_BP:
                    continue
                if max_focal_bp > 0 and span > max_focal_bp:
                    continue  # arm/chromosome-scale -> depth-modulation path (#189)
                c = norm_chrom(f[0], chr_prefix)
                pos = start + 1
                # Homozygous deletion (CN=0) is genuinely homozygous; other
                # CN changes have no zygosity info in this table -> het.
                gt = "1/1" if total_cn == 0 else "0/1"
                records.append(
                    (c, pos, "<CNV>",
                     f"SVTYPE=CNV;END={end};CN={total_cn}", gt))
                counts["CNV"] += 1
                counts["_cn_hist"][total_cn] = \
                    counts["_cn_hist"].get(total_cn, 0) + 1
    counts["_donors_cnv"] = donors
    if missing_ploidy:
        log(f"   CNA: {missing_ploidy} donors missing from ploidy table "
            f"(used diploid neutral)")
    log(f"   CNA: {donors} donors -> CNV(focal,non-neutral)={counts['CNV']}")


# ── gnomAD-SV → INS (length only) ────────────────────────────────────────
def parse_gnomad_ins(vcf, chr_prefix, records, counts):
    svtype_re = re.compile(r"SVTYPE=([^;]+)")
    svlen_re = re.compile(r"SVLEN=(-?\d+)")
    n = 0
    with gzip.open(vcf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            info = f[7]
            m = svtype_re.search(info)
            if not m or m.group(1) != "INS":
                continue
            ml = svlen_re.search(info)
            if not ml:
                continue
            svlen = abs(int(ml.group(1)))
            if svlen < MIN_SV_LENGTH_BP:
                continue
            c = norm_chrom(f[0], chr_prefix)
            pos = int(f[1])
            records.append(
                (c, pos, "<INS>", f"SVTYPE=INS;SVLEN={svlen}", "0/1"))
            n += 1
    counts["INS"] = n
    log(f"   gnomAD INS: {n} records (length source only)")


# ── VCF writer ───────────────────────────────────────────────────────────
HEADER = """##fileformat=VCFv4.2
##source=tools/build_pcawg_sv_vcf.py
##reference=PCAWG-consensus-GRCh38 (SV/CNV) + gnomAD-SV (INS length)
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=INV,Description="Inversion">
##ALT=<ID=CNV,Description="Copy number variant">
##ALT=<ID=BND,Description="Breakend">
##ALT=<ID=INS,Description="Insertion">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT">
##INFO=<ID=CN,Number=1,Type=Integer,Description="Copy number of segment">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype (synthetic; see #218)">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNEAT_pcawg_corpus
"""


def chrom_key(c):
    bare = c[3:] if c.startswith("chr") else c
    special = {"X": 100, "Y": 101, "M": 102}
    if bare in special:
        return (special[bare],)
    if bare.isdigit():
        return (int(bare),)
    return (200,)


def write_vcf(out_path, records):
    records.sort(key=lambda r: (chrom_key(r[0]), r[1]))
    op = gzip.open(out_path, "wt") if out_path.endswith(".gz") else \
        open(out_path, "w")
    with op as out:
        out.write(HEADER)
        for c, pos, alt, info, gt in records:
            out.write(f"{c}\t{pos}\t.\tN\t{alt}\t.\tPASS\t{info}\tGT\t{gt}\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bedpe-dir", required=True,
                    help="consensus_sv dir (has icgc/ + tcga/ BEDPE trees)")
    ap.add_argument("--cna-dir",
                    help="consensus_cnv dir (CNA segment trees). Omit to skip CNV.")
    ap.add_argument("--purity-ploidy",
                    help="purity/ploidy table (.txt.gz). Needed for CNV neutral CN.")
    ap.add_argument("--gnomad-ins-vcf",
                    help="gnomAD-SV sites VCF for INS length. Omit to skip INS.")
    ap.add_argument("--out", required=True, help="output VCF[.gz]")
    ap.add_argument("--cnv-max-focal-bp", type=int, default=3_000_000,
                    help="drop CNA segments longer than this as arm/chromosome "
                         "aneuploidy (handled via depth-modulation, #189). "
                         "0 = no cap. Default 3,000,000.")
    ap.add_argument("--no-chr-prefix", action="store_true",
                    help="emit bare chrom names (1/X/M) instead of chr-prefixed")
    ap.add_argument("--sample-sheet",
                    help="pcawg_sample_sheet.tsv (aliquot_id -> dcc_project_code). "
                         "Required with --projects for per-tissue filtering.")
    ap.add_argument("--projects",
                    help="comma-separated dcc_project_code prefixes to include "
                         "(e.g. 'BRCA' for BRCA-US/UK/EU, or 'SKCM-US,MELA-AU' for "
                         "skin melanoma). Omit for pan-cancer (all donors).")
    ap.add_argument("--tissue-label",
                    help="label recorded in the sidecar (e.g. BRCA, skin, lung)")
    args = ap.parse_args()

    chr_prefix = not args.no_chr_prefix
    records = []
    counts = {"DEL": 0, "DUP": 0, "INV": 0, "BND": 0, "CNV": 0, "INS": 0,
              "_cn_hist": {}}

    # Per-tissue filter: restrict the somatic BEDPE/CNA donors to a tissue's
    # aliquots. gnomAD INS is germline (tissue-agnostic) so it is never filtered.
    allowed = None
    if args.projects:
        if not args.sample_sheet:
            ap.error("--projects requires --sample-sheet")
        projects = [p.strip() for p in args.projects.split(",") if p.strip()]
        log(f">> Per-tissue filter: projects {projects}")
        allowed = load_tissue_aliquots(args.sample_sheet, projects)
        if not allowed:
            ap.error(f"no aliquots matched --projects {projects}")

    log(">> Parsing PCAWG BEDPE (DEL/DUP/INV/BND)...")
    parse_bedpe(args.bedpe_dir, chr_prefix, records, counts, allowed)

    if args.cna_dir:
        if not args.purity_ploidy:
            log("!! --cna-dir given without --purity-ploidy; "
                "using diploid neutral for all samples")
            ploidy = {}
        else:
            log(">> Loading purity/ploidy table...")
            ploidy = load_ploidy(args.purity_ploidy)
        log(">> Parsing PCAWG CNA (focal non-neutral CNV)...")
        parse_cna(args.cna_dir, ploidy, chr_prefix, args.cnv_max_focal_bp,
                  records, counts, allowed)

    if args.gnomad_ins_vcf:
        log(">> Parsing gnomAD-SV INS (length source)...")
        parse_gnomad_ins(args.gnomad_ins_vcf, chr_prefix, records, counts)

    log(f">> Writing {len(records)} records -> {args.out}")
    write_vcf(args.out, records)

    # Sidecar: per-type totals + per-source donor counts for the downstream
    # per-tumor rate/type-probability normalization (NOT derivable from the
    # fit because sources have different denominators).
    sidecar = {
        "type_totals": {k: counts[k] for k in
                        ("DEL", "DUP", "INV", "BND", "CNV", "INS")},
        "donors_sv": counts.get("_donors_sv", 0),
        "donors_cnv": counts.get("_donors_cnv", 0),
        "cnv_copy_number_hist": counts["_cn_hist"],
        "notes": {
            "tissue_label": args.tissue_label,
            "projects": (args.projects.split(",") if args.projects else None),
            "ins_source": "gnomad-sv germline (length only; rate is literature)",
            "cnv_focal_max_bp": args.cnv_max_focal_bp,
            "per_donor_means": {
                k: (counts[k] / counts["_donors_sv"]
                    if k in ("DEL", "DUP", "INV", "BND")
                    and counts.get("_donors_sv") else
                    counts[k] / counts["_donors_cnv"]
                    if k == "CNV" and counts.get("_donors_cnv") else None)
                for k in ("DEL", "DUP", "INV", "BND", "CNV")
            },
        },
    }
    sidecar_path = (args.out[:-7] if args.out.endswith(".vcf.gz")
                    else os.path.splitext(args.out)[0]) + ".counts.json"
    with open(sidecar_path, "w") as fh:
        json.dump(sidecar, fh, indent=2)
    log(f">> Sidecar counts -> {sidecar_path}")
    log("   (DEL/DUP/INV/BND per-donor means + CNV; feed these to the "
        "normalization step to set type_probabilities + per_base_rate.)")


if __name__ == "__main__":
    main()

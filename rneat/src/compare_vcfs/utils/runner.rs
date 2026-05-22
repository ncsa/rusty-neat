//! Phase 1 driver for `compare-vcfs`: exact-match classification.
//!
//! Algorithm:
//!   1. Read golden + called VCFs via `common::file_tools::vcf_tools::read_vcf`.
//!   2. For each VCF, drop variants the user has filtered out: outside the
//!      simulated contigs, outside the target BED, hom-ref, non-PASS — each
//!      counted into a [`SkipCounts`] for the report.
//!   3. Per contig, key surviving variants by `(position, ref, alt)` into a
//!      `HashSet<VariantKey>`. TP = golden ∩ called; FN = golden \ called;
//!      FP = called \ golden.
//!   4. Aggregate per-contig and total counts, compute precision/recall/F1,
//!      write `comparison_summary.{json,txt}`.
//!
//! Equivalence detection (denotation-different indels) and NEAT-aware FN
//! attribution are out of scope; see issue #127.
use crate::compare_vcfs::{
    errors::CompareVcfsError,
    utils::{
        config::RunConfiguration,
        report::{
            ComparisonSummary, ContigCounts, Metrics, SCHEMA_VERSION, SkipCounts, SummaryInputs,
            write_json, write_txt,
        },
    },
};
use common::{
    file_tools::{bed_reader::read_bed, vcf_tools::read_vcf},
    structs::{bed_record::BedRecord, nucleotides::Nucleotide, variants::Variant},
};
use log::*;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs;

/// Comparison key: the three columns that decide a true positive. Hashing on
/// `Variant` directly is wrong because its derived `Hash` includes
/// genotype/filter/info, none of which the truth-vs-called comparison cares about.
#[derive(Debug, Clone, Hash, PartialEq, Eq)]
struct VariantKey {
    location: usize,
    reference: Vec<Nucleotide>,
    alternate: Vec<Nucleotide>,
}

impl From<&Variant> for VariantKey {
    fn from(v: &Variant) -> Self {
        VariantKey {
            location: v.location,
            reference: v.reference.clone(),
            alternate: v.alternate.clone(),
        }
    }
}

pub fn runner(config: &RunConfiguration) -> Result<(), CompareVcfsError> {
    fs::create_dir_all(&config.output_dir)?;

    let target_bed = match &config.target_bed {
        Some(p) => Some(read_bed(p, false)?),
        None => None,
    };
    let simulated_set = config
        .contigs_simulated
        .as_ref()
        .map(|v| v.iter().cloned().collect::<HashSet<_>>());

    info!("Reading golden VCF: {}", config.golden_vcf.display());
    let golden_raw = read_vcf(config.golden_vcf.clone())?;
    info!("Reading called VCF: {}", config.called_vcf.display());
    let called_raw = read_vcf(config.called_vcf.clone())?;

    let (golden_by_contig, skipped_golden) = filter_vcf(
        &golden_raw,
        target_bed.as_ref(),
        simulated_set.as_ref(),
        config.include_homs,
        config.include_filtered,
    );
    let (called_by_contig, skipped_called) = filter_vcf(
        &called_raw,
        target_bed.as_ref(),
        simulated_set.as_ref(),
        config.include_homs,
        config.include_filtered,
    );

    let (per_contig, totals) = classify(&golden_by_contig, &called_by_contig);
    let metrics = Metrics::from_counts(totals.tp, totals.fn_, totals.fp);
    info!(
        "Classification: TP={} FN={} FP={}",
        totals.tp, totals.fn_, totals.fp
    );

    let summary = ComparisonSummary {
        schema_version: SCHEMA_VERSION,
        inputs: SummaryInputs {
            golden_vcf: config.golden_vcf.clone(),
            called_vcf: config.called_vcf.clone(),
            reference: config.reference.clone(),
            target_bed: config.target_bed.clone(),
            contigs_simulated: config.contigs_simulated.clone(),
            include_homs: config.include_homs,
            include_filtered: config.include_filtered,
        },
        totals,
        metrics,
        per_contig,
        skipped_golden,
        skipped_called,
    };

    let json_path = write_json(&summary, &config.output_dir, config.overwrite_output)?;
    let txt_path = write_txt(&summary, &config.output_dir, config.overwrite_output)?;
    info!("Wrote {}", json_path.display());
    info!("Wrote {}", txt_path.display());
    Ok(())
}

/// Drops disqualified variants and tallies the reasons. The returned map keeps
/// only contigs whose surviving variant list is non-empty so downstream
/// classification doesn't iterate empty buckets.
fn filter_vcf(
    raw: &HashMap<String, Vec<Variant>>,
    target_bed: Option<&HashMap<String, Vec<BedRecord>>>,
    simulated_set: Option<&HashSet<String>>,
    include_homs: bool,
    include_filtered: bool,
) -> (HashMap<String, Vec<Variant>>, SkipCounts) {
    let mut kept: HashMap<String, Vec<Variant>> = HashMap::new();
    let mut skipped = SkipCounts::default();
    for (chrom, variants) in raw {
        if let Some(sim) = simulated_set
            && !sim.contains(chrom)
        {
            skipped.outside_simulated_contigs += variants.len();
            continue;
        }
        for v in variants {
            if !include_homs && v.reference == v.alternate {
                skipped.homozygous_ref += 1;
                continue;
            }
            if !include_filtered && !filter_is_pass(v) {
                skipped.filtered += 1;
                continue;
            }
            if let Some(bed) = target_bed {
                let regions = bed.get(chrom);
                let in_target = regions
                    .map(|rs| rs.iter().any(|r| r.contains(chrom, v.location)))
                    .unwrap_or(false);
                if !in_target {
                    skipped.outside_target_bed += 1;
                    continue;
                }
            }
            kept.entry(chrom.clone()).or_default().push(v.clone());
        }
    }
    (kept, skipped)
}

fn filter_is_pass(v: &Variant) -> bool {
    match v.filter.as_deref() {
        None | Some("PASS") | Some(".") | Some("") => true,
        Some(_) => false,
    }
}

fn classify(
    golden: &HashMap<String, Vec<Variant>>,
    called: &HashMap<String, Vec<Variant>>,
) -> (BTreeMap<String, ContigCounts>, ContigCounts) {
    let mut per_contig: BTreeMap<String, ContigCounts> = BTreeMap::new();
    let mut totals = ContigCounts {
        tp: 0,
        fn_: 0,
        fp: 0,
    };

    let all_contigs: HashSet<&String> = golden.keys().chain(called.keys()).collect();
    for chrom in all_contigs {
        let g_keys: HashSet<VariantKey> = golden
            .get(chrom)
            .map(|vs| vs.iter().map(VariantKey::from).collect())
            .unwrap_or_default();
        let c_keys: HashSet<VariantKey> = called
            .get(chrom)
            .map(|vs| vs.iter().map(VariantKey::from).collect())
            .unwrap_or_default();
        let tp = g_keys.intersection(&c_keys).count();
        let fn_ = g_keys.len() - tp;
        let fp = c_keys.len() - tp;
        totals.tp += tp;
        totals.fn_ += fn_;
        totals.fp += fp;
        per_contig.insert(chrom.clone(), ContigCounts { tp, fn_, fp });
    }
    (per_contig, totals)
}

#[cfg(test)]
mod tests {
    use super::*;
    use common::structs::variants::{Genotype, VariantType};

    fn snp(loc: usize, ref_b: Nucleotide, alt_b: Nucleotide, filter: Option<&str>) -> Variant {
        Variant {
            variant_type: VariantType::SNP,
            location: loc,
            reference: vec![ref_b],
            alternate: vec![alt_b],
            genotype_str: "0/1".to_string(),
            genotype: Genotype::Heterozygous,
            id: None,
            quality_score: None,
            filter: filter.map(str::to_string),
            info: None,
            format: Vec::new(),
            sample: Vec::new(),
        }
    }

    #[test]
    fn classify_all_tp() {
        let v = snp(10, Nucleotide::A, Nucleotide::C, Some("PASS"));
        let mut g = HashMap::new();
        g.insert("chr1".into(), vec![v.clone()]);
        let mut c = HashMap::new();
        c.insert("chr1".into(), vec![v]);
        let (per, tot) = classify(&g, &c);
        assert_eq!(tot.tp, 1);
        assert_eq!(tot.fn_, 0);
        assert_eq!(tot.fp, 0);
        assert_eq!(per["chr1"].tp, 1);
    }

    #[test]
    fn classify_mixed() {
        // golden has (10,A,C) and (20,G,T); called has (10,A,C) and (30,T,A).
        // → 1 TP, 1 FN, 1 FP.
        let mut g = HashMap::new();
        g.insert(
            "chr1".into(),
            vec![
                snp(10, Nucleotide::A, Nucleotide::C, Some("PASS")),
                snp(20, Nucleotide::G, Nucleotide::T, Some("PASS")),
            ],
        );
        let mut c = HashMap::new();
        c.insert(
            "chr1".into(),
            vec![
                snp(10, Nucleotide::A, Nucleotide::C, Some("PASS")),
                snp(30, Nucleotide::T, Nucleotide::A, Some("PASS")),
            ],
        );
        let (_per, tot) = classify(&g, &c);
        assert_eq!(tot.tp, 1);
        assert_eq!(tot.fn_, 1);
        assert_eq!(tot.fp, 1);
    }

    #[test]
    fn classify_contig_only_in_one_side() {
        // Golden has variants on chr1 only; called has variants on chr2 only.
        // Per-contig: chr1 → 1 FN, chr2 → 1 FP. Totals: 0 TP, 1 FN, 1 FP.
        let mut g = HashMap::new();
        g.insert(
            "chr1".into(),
            vec![snp(10, Nucleotide::A, Nucleotide::C, Some("PASS"))],
        );
        let mut c = HashMap::new();
        c.insert(
            "chr2".into(),
            vec![snp(50, Nucleotide::G, Nucleotide::A, Some("PASS"))],
        );
        let (per, tot) = classify(&g, &c);
        assert_eq!(tot.tp, 0);
        assert_eq!(tot.fn_, 1);
        assert_eq!(tot.fp, 1);
        assert_eq!(per["chr1"].fn_, 1);
        assert_eq!(per["chr2"].fp, 1);
    }

    #[test]
    fn filter_skips_homref_when_not_included() {
        let mut raw = HashMap::new();
        raw.insert(
            "chr1".into(),
            vec![snp(10, Nucleotide::A, Nucleotide::A, Some("PASS"))],
        );
        let (kept, skipped) = filter_vcf(&raw, None, None, false, false);
        assert!(kept.is_empty());
        assert_eq!(skipped.homozygous_ref, 1);
    }

    #[test]
    fn filter_skips_non_pass_when_not_included() {
        let mut raw = HashMap::new();
        raw.insert(
            "chr1".into(),
            vec![snp(10, Nucleotide::A, Nucleotide::C, Some("LowQual"))],
        );
        let (kept, skipped) = filter_vcf(&raw, None, None, false, false);
        assert!(kept.is_empty());
        assert_eq!(skipped.filtered, 1);
    }

    #[test]
    fn filter_includes_non_pass_when_flag_set() {
        let mut raw = HashMap::new();
        raw.insert(
            "chr1".into(),
            vec![snp(10, Nucleotide::A, Nucleotide::C, Some("LowQual"))],
        );
        let (kept, skipped) = filter_vcf(&raw, None, None, false, true);
        assert_eq!(kept.get("chr1").unwrap().len(), 1);
        assert_eq!(skipped.filtered, 0);
    }

    #[test]
    fn filter_drops_unsimulated_contigs() {
        let mut raw = HashMap::new();
        raw.insert(
            "chr1".into(),
            vec![snp(10, Nucleotide::A, Nucleotide::C, Some("PASS"))],
        );
        raw.insert(
            "chr2".into(),
            vec![snp(20, Nucleotide::G, Nucleotide::T, Some("PASS"))],
        );
        let sim: HashSet<String> = ["chr1".to_string()].into_iter().collect();
        let (kept, skipped) = filter_vcf(&raw, None, Some(&sim), false, false);
        assert!(kept.contains_key("chr1"));
        assert!(!kept.contains_key("chr2"));
        assert_eq!(skipped.outside_simulated_contigs, 1);
    }
}

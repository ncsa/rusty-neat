//! Output artifacts for `compare-vcfs`: a schema-versioned JSON rollup and a
//! plain-text rendering of the same data.
//!
//! Schema is versioned so downstream tooling can detect format drift; bump
//! [`SCHEMA_VERSION`] only when a backward-incompatible field change lands.
use crate::compare_vcfs::errors::CompareVcfsError;
use crate::compare_vcfs::utils::attribution::{ChromNamingWarning, Reason};
use serde::Serialize;
use std::collections::BTreeMap;
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};

pub const SCHEMA_VERSION: &str = "1.2.0";

#[derive(Debug, Clone, Serialize)]
pub struct ContigCounts {
    pub tp: usize,
    pub fn_: usize,
    pub fp: usize,
    /// FNs promoted to TP by the equivalence sweep. Already included in `tp`
    /// and excluded from `fn_`; surfaced separately so users can see how many
    /// matches were rescued by denotation-aware comparison.
    #[serde(default)]
    pub equivalents_promoted: usize,
}

#[derive(Debug, Clone, Serialize)]
pub struct Metrics {
    /// `None` when TP + FP == 0.
    pub precision: Option<f64>,
    /// `None` when TP + FN == 0.
    pub recall: Option<f64>,
    /// `None` when precision or recall is `None`, or when both are 0.
    pub f1: Option<f64>,
}

impl Metrics {
    pub fn from_counts(tp: usize, fn_: usize, fp: usize) -> Self {
        let precision = if tp + fp > 0 {
            Some(tp as f64 / (tp + fp) as f64)
        } else {
            None
        };
        let recall = if tp + fn_ > 0 {
            Some(tp as f64 / (tp + fn_) as f64)
        } else {
            None
        };
        let f1 = match (precision, recall) {
            (Some(p), Some(r)) if p + r > 0.0 => Some(2.0 * p * r / (p + r)),
            _ => None,
        };
        Metrics {
            precision,
            recall,
            f1,
        }
    }
}

/// Per-VCF parser-skip counters surfaced to the report so users aren't
/// surprised that some records didn't make it into the classification.
#[derive(Debug, Clone, Default, Serialize)]
pub struct SkipCounts {
    pub multiallelic: usize,
    pub homozygous_ref: usize,
    pub filtered: usize,
    pub outside_target_bed: usize,
    pub outside_simulated_contigs: usize,
}

#[derive(Debug, Clone, Serialize)]
pub struct ComparisonSummary {
    pub schema_version: &'static str,
    pub inputs: SummaryInputs,
    pub totals: ContigCounts,
    pub metrics: Metrics,
    /// Keyed by contig name. BTreeMap so JSON / TXT output is sorted.
    pub per_contig: BTreeMap<String, ContigCounts>,
    pub skipped_golden: SkipCounts,
    pub skipped_called: SkipCounts,
    /// Rolled-up false-negative reason counts. Each FN contributes once per
    /// tagged reason — see `attribution::Reason` for the taxonomy.
    #[serde(default)]
    pub fn_attribution: BTreeMap<Reason, usize>,
    /// Warnings surfaced during attribution. Currently: BED-vs-reference
    /// chrom-naming mismatches.
    #[serde(default)]
    pub warnings: Vec<ChromNamingWarning>,
    /// Absolute paths to every artifact this run wrote. Lets downstream
    /// tooling locate sibling files without re-deriving them from the
    /// output_dir layout.
    pub outputs: SummaryOutputs,
}

#[derive(Debug, Clone, Serialize)]
pub struct SummaryOutputs {
    pub comparison_summary_json: PathBuf,
    pub comparison_summary_txt: PathBuf,
    pub fn_with_reasons_vcf: PathBuf,
    /// `None` when `write_fp_vcf` was false.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub fp_vcf: Option<PathBuf>,
}

#[derive(Debug, Clone, Serialize)]
pub struct SummaryInputs {
    pub golden_vcf: PathBuf,
    pub called_vcf: PathBuf,
    pub reference: PathBuf,
    pub target_bed: Option<PathBuf>,
    pub contigs_simulated: Option<Vec<String>>,
    pub include_homs: bool,
    pub include_filtered: bool,
    pub equivalence_window: usize,
    pub fast: bool,
    pub mutation_bed: Option<PathBuf>,
    pub chrom_aliases: Option<PathBuf>,
}

pub fn write_json(
    summary: &ComparisonSummary,
    output_dir: &Path,
    overwrite_output: bool,
) -> Result<PathBuf, CompareVcfsError> {
    let path = output_dir.join("comparison_summary.json");
    check_overwrite(&path, overwrite_output)?;
    let file = fs::File::create(&path)?;
    serde_json::to_writer_pretty(file, summary)?;
    Ok(path)
}

pub fn write_txt(
    summary: &ComparisonSummary,
    output_dir: &Path,
    overwrite_output: bool,
) -> Result<PathBuf, CompareVcfsError> {
    let path = output_dir.join("comparison_summary.txt");
    check_overwrite(&path, overwrite_output)?;
    let mut file = fs::File::create(&path)?;
    file.write_all(render_txt(summary).as_bytes())?;
    Ok(path)
}

fn check_overwrite(path: &Path, allow: bool) -> Result<(), CompareVcfsError> {
    if path.is_file() && !allow {
        return Err(CompareVcfsError::OverwriteFileError(
            path.display().to_string(),
        ));
    }
    Ok(())
}

fn render_txt(s: &ComparisonSummary) -> String {
    let mut out = String::new();
    out.push_str("rneat compare-vcfs report\n");
    out.push_str("=========================\n\n");
    out.push_str(&format!("Schema version: {}\n\n", s.schema_version));
    out.push_str("Inputs\n------\n");
    out.push_str(&format!(
        "  Golden VCF: {}\n",
        s.inputs.golden_vcf.display()
    ));
    out.push_str(&format!(
        "  Called VCF: {}\n",
        s.inputs.called_vcf.display()
    ));
    out.push_str(&format!("  Reference:  {}\n", s.inputs.reference.display()));
    if let Some(p) = &s.inputs.target_bed {
        out.push_str(&format!("  Target BED: {}\n", p.display()));
    }
    if let Some(c) = &s.inputs.contigs_simulated {
        out.push_str(&format!("  Contigs simulated: {}\n", c.join(", ")));
    }
    out.push_str(&format!("  Include hom-refs:  {}\n", s.inputs.include_homs));
    out.push_str(&format!(
        "  Include filtered:  {}\n",
        s.inputs.include_filtered
    ));
    if s.inputs.fast {
        out.push_str("  Equivalence sweep: SKIPPED (fast mode)\n");
    } else {
        out.push_str(&format!(
            "  Equivalence window: ±{} bp\n",
            s.inputs.equivalence_window
        ));
    }
    out.push('\n');

    out.push_str("Totals\n------\n");
    out.push_str(&format!("  True positives  (TP): {}\n", s.totals.tp));
    out.push_str(&format!(
        "    ↳ promoted by equivalence sweep: {}\n",
        s.totals.equivalents_promoted
    ));
    out.push_str(&format!("  False negatives (FN): {}\n", s.totals.fn_));
    out.push_str(&format!("  False positives (FP): {}\n\n", s.totals.fp));

    out.push_str("Metrics\n-------\n");
    out.push_str(&format!(
        "  Precision: {}  (TP / (TP + FP))\n",
        fmt_metric(s.metrics.precision)
    ));
    out.push_str(&format!(
        "  Recall:    {}  (TP / (TP + FN))\n",
        fmt_metric(s.metrics.recall)
    ));
    out.push_str(&format!("  F1:        {}\n\n", fmt_metric(s.metrics.f1)));

    out.push_str("Per-contig counts\n-----------------\n");
    if s.per_contig.is_empty() {
        out.push_str("  (none)\n\n");
    } else {
        let name_width = s
            .per_contig
            .keys()
            .map(|k| k.len())
            .max()
            .unwrap_or(0)
            .max(6);
        out.push_str(&format!(
            "  {:<width$}    TP      FN      FP\n",
            "Contig",
            width = name_width
        ));
        for (chrom, c) in &s.per_contig {
            out.push_str(&format!(
                "  {:<width$}  {:>6}  {:>6}  {:>6}\n",
                chrom,
                c.tp,
                c.fn_,
                c.fp,
                width = name_width
            ));
        }
        out.push('\n');
    }

    out.push_str("Skipped records (golden / called)\n---------------------------------\n");
    out.push_str(&format!(
        "  Multi-allelic:               {} / {}\n",
        s.skipped_golden.multiallelic, s.skipped_called.multiallelic
    ));
    out.push_str(&format!(
        "  Homozygous reference:        {} / {}\n",
        s.skipped_golden.homozygous_ref, s.skipped_called.homozygous_ref
    ));
    out.push_str(&format!(
        "  Filtered (non-PASS):         {} / {}\n",
        s.skipped_golden.filtered, s.skipped_called.filtered
    ));
    out.push_str(&format!(
        "  Outside target BED:          {} / {}\n",
        s.skipped_golden.outside_target_bed, s.skipped_called.outside_target_bed
    ));
    out.push_str(&format!(
        "  Outside simulated contigs:   {} / {}\n",
        s.skipped_golden.outside_simulated_contigs, s.skipped_called.outside_simulated_contigs
    ));

    out.push_str("\nFN attribution\n--------------\n");
    if s.fn_attribution.is_empty() {
        out.push_str("  (no false negatives)\n");
    } else {
        let label_width = s
            .fn_attribution
            .keys()
            .map(|r| r.as_str().len())
            .max()
            .unwrap_or(0)
            + 2;
        for (reason, count) in &s.fn_attribution {
            out.push_str(&format!(
                "  {:<width$} {}\n",
                reason.as_str(),
                count,
                width = label_width
            ));
        }
    }

    if !s.warnings.is_empty() {
        out.push_str("\nWarnings\n--------\n");
        for w in &s.warnings {
            out.push_str(&format!("  - {}\n", w.message));
        }
    }

    out.push_str("\nOutputs\n-------\n");
    out.push_str(&format!(
        "  Summary (JSON): {}\n",
        s.outputs.comparison_summary_json.display()
    ));
    out.push_str(&format!(
        "  Summary (TXT):  {}\n",
        s.outputs.comparison_summary_txt.display()
    ));
    out.push_str(&format!(
        "  Annotated FNs:  {}\n",
        s.outputs.fn_with_reasons_vcf.display()
    ));
    if let Some(p) = &s.outputs.fp_vcf {
        out.push_str(&format!("  False positives:{}\n", p.display()));
    }

    out
}

fn fmt_metric(v: Option<f64>) -> String {
    match v {
        Some(x) => format!("{x:.4}"),
        None => "N/A".to_string(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_metrics_basic() {
        let m = Metrics::from_counts(80, 20, 5);
        assert!((m.precision.unwrap() - (80.0 / 85.0)).abs() < 1e-9);
        assert!((m.recall.unwrap() - 0.8).abs() < 1e-9);
        assert!(m.f1.unwrap() > 0.0);
    }

    #[test]
    fn test_metrics_zero_truth() {
        let m = Metrics::from_counts(0, 0, 3);
        assert!(m.recall.is_none());
        assert_eq!(m.precision, Some(0.0));
        assert!(m.f1.is_none());
    }

    #[test]
    fn test_render_txt_includes_totals() {
        let s = ComparisonSummary {
            schema_version: SCHEMA_VERSION,
            inputs: SummaryInputs {
                golden_vcf: PathBuf::from("/g.vcf"),
                called_vcf: PathBuf::from("/c.vcf"),
                reference: PathBuf::from("/r.fa"),
                target_bed: None,
                contigs_simulated: None,
                include_homs: false,
                include_filtered: false,
                equivalence_window: 50,
                fast: false,
                mutation_bed: None,
                chrom_aliases: None,
            },
            totals: ContigCounts {
                tp: 7,
                fn_: 1,
                fp: 2,
                equivalents_promoted: 0,
            },
            metrics: Metrics::from_counts(7, 1, 2),
            per_contig: BTreeMap::new(),
            skipped_golden: SkipCounts::default(),
            skipped_called: SkipCounts::default(),
            fn_attribution: BTreeMap::new(),
            warnings: Vec::new(),
            outputs: SummaryOutputs {
                comparison_summary_json: PathBuf::from("/out.json"),
                comparison_summary_txt: PathBuf::from("/out.txt"),
                fn_with_reasons_vcf: PathBuf::from("/out_fn.vcf"),
                fp_vcf: None,
            },
        };
        let txt = render_txt(&s);
        assert!(txt.contains("True positives  (TP): 7"));
        assert!(txt.contains("False negatives (FN): 1"));
        assert!(txt.contains("False positives (FP): 2"));
    }

    #[test]
    fn test_write_json_round_trips() {
        let dir = tempfile::tempdir().unwrap();
        let s = ComparisonSummary {
            schema_version: SCHEMA_VERSION,
            inputs: SummaryInputs {
                golden_vcf: PathBuf::from("/g.vcf"),
                called_vcf: PathBuf::from("/c.vcf"),
                reference: PathBuf::from("/r.fa"),
                target_bed: None,
                contigs_simulated: None,
                include_homs: false,
                include_filtered: false,
                equivalence_window: 50,
                fast: false,
                mutation_bed: None,
                chrom_aliases: None,
            },
            totals: ContigCounts {
                tp: 3,
                fn_: 0,
                fp: 0,
                equivalents_promoted: 0,
            },
            metrics: Metrics::from_counts(3, 0, 0),
            per_contig: BTreeMap::new(),
            skipped_golden: SkipCounts::default(),
            skipped_called: SkipCounts::default(),
            fn_attribution: BTreeMap::new(),
            warnings: Vec::new(),
            outputs: SummaryOutputs {
                comparison_summary_json: PathBuf::from("/out.json"),
                comparison_summary_txt: PathBuf::from("/out.txt"),
                fn_with_reasons_vcf: PathBuf::from("/out_fn.vcf"),
                fp_vcf: None,
            },
        };
        let path = write_json(&s, dir.path(), false).unwrap();
        let raw = fs::read_to_string(&path).unwrap();
        let val: serde_json::Value = serde_json::from_str(&raw).unwrap();
        assert_eq!(val["schema_version"], "1.2.0");
        assert_eq!(val["totals"]["tp"], 3);
    }

    #[test]
    fn test_overwrite_refused() {
        let dir = tempfile::tempdir().unwrap();
        let existing = dir.path().join("comparison_summary.json");
        fs::write(&existing, "stale").unwrap();
        let s = ComparisonSummary {
            schema_version: SCHEMA_VERSION,
            inputs: SummaryInputs {
                golden_vcf: PathBuf::from("/g"),
                called_vcf: PathBuf::from("/c"),
                reference: PathBuf::from("/r"),
                target_bed: None,
                contigs_simulated: None,
                include_homs: false,
                include_filtered: false,
                equivalence_window: 50,
                fast: false,
                mutation_bed: None,
                chrom_aliases: None,
            },
            totals: ContigCounts {
                tp: 0,
                fn_: 0,
                fp: 0,
                equivalents_promoted: 0,
            },
            metrics: Metrics::from_counts(0, 0, 0),
            per_contig: BTreeMap::new(),
            skipped_golden: SkipCounts::default(),
            skipped_called: SkipCounts::default(),
            fn_attribution: BTreeMap::new(),
            warnings: Vec::new(),
            outputs: SummaryOutputs {
                comparison_summary_json: PathBuf::from("/out.json"),
                comparison_summary_txt: PathBuf::from("/out.txt"),
                fn_with_reasons_vcf: PathBuf::from("/out_fn.vcf"),
                fp_vcf: None,
            },
        };
        let err = write_json(&s, dir.path(), false).unwrap_err();
        assert!(matches!(err, CompareVcfsError::OverwriteFileError(_)));
    }
}

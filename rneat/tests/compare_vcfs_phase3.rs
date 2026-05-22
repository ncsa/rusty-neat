//! Phase 3 integration tests: NEAT-aware FN attribution + chrom aliases.
//!
//! Each test arranges golden + called VCFs so the surviving FN(s) hit a
//! specific reason path, then asserts the JSON report's `fn_attribution`
//! map has the expected counts.

mod common;

use common::rneat;
use std::io::Write;
use std::path::Path;

fn write_fasta(path: &Path, contig: &str, sequence: &str) {
    let mut f = std::fs::File::create(path).unwrap();
    writeln!(f, ">{contig}").unwrap();
    writeln!(f, "{sequence}").unwrap();
}

fn write_vcf(path: &Path, body_lines: &[&str]) {
    let mut f = std::fs::File::create(path).unwrap();
    writeln!(f, "##fileformat=VCFv4.2").unwrap();
    writeln!(
        f,
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    )
    .unwrap();
    writeln!(
        f,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"
    )
    .unwrap();
    for line in body_lines {
        writeln!(f, "{line}").unwrap();
    }
}

fn write_bed(path: &Path, regions: &[(&str, usize, usize)]) {
    let mut f = std::fs::File::create(path).unwrap();
    for (chrom, start, end) in regions {
        writeln!(f, "{chrom}\t{start}\t{end}").unwrap();
    }
}

fn write_text(path: &Path, content: &str) {
    let mut f = std::fs::File::create(path).unwrap();
    write!(f, "{content}").unwrap();
}

fn load_summary(dir: &Path) -> serde_json::Value {
    let raw = std::fs::read_to_string(dir.join("comparison_summary.json")).unwrap();
    serde_json::from_str(&raw).unwrap()
}

/// An FN at position 500 of chr1, with `mutation_bed: chr1 0 100` (covers
/// only the first 100bp), should be tagged `outside_mutation_bed`.
#[test]
fn fn_outside_mutation_bed_is_attributed() {
    let tmp = tempfile::tempdir().unwrap();
    let fa = tmp.path().join("ref.fa");
    let golden = tmp.path().join("golden.vcf");
    let called = tmp.path().join("called.vcf");
    let mbed = tmp.path().join("mutation.bed");
    let out = tmp.path().join("report");
    // 1000bp reference filled with `A`.
    write_fasta(&fa, "chr1", &"A".repeat(1000));
    // Golden has a SNP at position 500 — outside the mutation BED region.
    write_vcf(&golden, &["chr1\t500\t.\tA\tC\t60\tPASS\t.\tGT\t0/1"]);
    // Called VCF is empty (just the header).
    write_vcf(&called, &[]);
    // Mutation BED covers only positions 0-99.
    write_bed(&mbed, &[("chr1", 0, 100)]);

    let yaml = tmp.path().join("cfg.yml");
    write_text(
        &yaml,
        &format!(
            "golden_vcf: {}\ncalled_vcf: {}\nreference: {}\noutput_dir: {}\n\
             overwrite_output: true\nmutation_bed: {}\n",
            golden.display(),
            called.display(),
            fa.display(),
            out.display(),
            mbed.display(),
        ),
    );

    rneat()
        .args(["compare-vcfs", "-c"])
        .arg(&yaml)
        .assert()
        .success();

    let summary = load_summary(&out);
    assert_eq!(summary["totals"]["fn_"], 1);
    assert_eq!(summary["fn_attribution"]["outside_mutation_bed"], 1);
    assert!(
        summary["fn_attribution"].get("unknown").is_none()
            || summary["fn_attribution"]["unknown"] == 0
    );
}

/// An FN on a contig that isn't in `contigs_simulated` should be tagged
/// `outside_simulated_contigs` even if it falls inside the mutation BED.
/// (The unsimulated-contig short-circuit fires before BED checks.)
#[test]
fn fn_outside_simulated_contigs_short_circuits_bed_checks() {
    let tmp = tempfile::tempdir().unwrap();
    let fa = tmp.path().join("ref.fa");
    let golden = tmp.path().join("golden.vcf");
    let called = tmp.path().join("called.vcf");
    let mbed = tmp.path().join("mutation.bed");
    let out = tmp.path().join("report");
    write_fasta(&fa, "chr2", &"A".repeat(1000));
    // Golden has an FN on chr2 at pos 50 — inside the mutation BED — but
    // the simulator only ran on chr1, so this should be reported as
    // `outside_simulated_contigs` alone.
    write_vcf(&golden, &["chr2\t50\t.\tA\tC\t60\tPASS\t.\tGT\t0/1"]);
    write_vcf(&called, &[]);
    write_bed(&mbed, &[("chr2", 0, 100)]);

    let yaml = tmp.path().join("cfg.yml");
    write_text(
        &yaml,
        &format!(
            "golden_vcf: {}\ncalled_vcf: {}\nreference: {}\noutput_dir: {}\n\
             overwrite_output: true\ncontigs_simulated: chr1\nmutation_bed: {}\n",
            golden.display(),
            called.display(),
            fa.display(),
            out.display(),
            mbed.display(),
        ),
    );

    rneat()
        .args(["compare-vcfs", "-c"])
        .arg(&yaml)
        .assert()
        .success();

    let summary = load_summary(&out);
    // The variant on chr2 is filtered out of the comparison entirely because
    // chr2 isn't in contigs_simulated — and contig-filter happens BEFORE
    // classification. So there's no FN to attribute, and the report should
    // show 0 FN and 0 attribution counts. Confirm the skip counter caught
    // it instead.
    assert_eq!(summary["totals"]["fn_"], 0);
    assert_eq!(summary["skipped_golden"]["outside_simulated_contigs"], 1);
}

/// FN inside both BEDs and on a simulated contig has no attribution reason
/// and should be tagged `unknown`.
#[test]
fn fn_with_no_simulator_reason_is_unknown() {
    let tmp = tempfile::tempdir().unwrap();
    let fa = tmp.path().join("ref.fa");
    let golden = tmp.path().join("golden.vcf");
    let called = tmp.path().join("called.vcf");
    let out = tmp.path().join("report");
    write_fasta(&fa, "chr1", &"A".repeat(1000));
    write_vcf(&golden, &["chr1\t500\t.\tA\tC\t60\tPASS\t.\tGT\t0/1"]);
    write_vcf(&called, &[]);

    let yaml = tmp.path().join("cfg.yml");
    write_text(
        &yaml,
        &format!(
            "golden_vcf: {}\ncalled_vcf: {}\nreference: {}\noutput_dir: {}\noverwrite_output: true\n",
            golden.display(),
            called.display(),
            fa.display(),
            out.display(),
        ),
    );
    rneat()
        .args(["compare-vcfs", "-c"])
        .arg(&yaml)
        .assert()
        .success();
    let summary = load_summary(&out);
    assert_eq!(summary["totals"]["fn_"], 1);
    assert_eq!(summary["fn_attribution"]["unknown"], 1);
}

/// A mutation BED whose chrom names use the `1` convention should still
/// attribute correctly when the VCFs use `chr1`, via a chrom_aliases TSV.
#[test]
fn chrom_aliases_remap_bed_names_for_attribution() {
    let tmp = tempfile::tempdir().unwrap();
    let fa = tmp.path().join("ref.fa");
    let golden = tmp.path().join("golden.vcf");
    let called = tmp.path().join("called.vcf");
    let mbed = tmp.path().join("mutation.bed");
    let aliases = tmp.path().join("aliases.tsv");
    let out = tmp.path().join("report");
    write_fasta(&fa, "chr1", &"A".repeat(1000));
    // VCFs use `chr1` (matching the reference). Mutation BED uses `1`.
    write_vcf(&golden, &["chr1\t500\t.\tA\tC\t60\tPASS\t.\tGT\t0/1"]);
    write_vcf(&called, &[]);
    write_bed(&mbed, &[("1", 0, 100)]);
    write_text(&aliases, "1\tchr1\n");

    let yaml = tmp.path().join("cfg.yml");
    write_text(
        &yaml,
        &format!(
            "golden_vcf: {}\ncalled_vcf: {}\nreference: {}\noutput_dir: {}\n\
             overwrite_output: true\nmutation_bed: {}\nchrom_aliases: {}\n",
            golden.display(),
            called.display(),
            fa.display(),
            out.display(),
            mbed.display(),
            aliases.display(),
        ),
    );

    rneat()
        .args(["compare-vcfs", "-c"])
        .arg(&yaml)
        .assert()
        .success();

    let summary = load_summary(&out);
    assert_eq!(summary["totals"]["fn_"], 1);
    // With aliases applied, the `1` BED becomes `chr1`, the FN at chr1:500
    // falls outside the [0, 100) region, and we get the expected
    // `outside_mutation_bed` attribution.
    assert_eq!(summary["fn_attribution"]["outside_mutation_bed"], 1);
    // No naming-mismatch warning expected because aliases successfully
    // bridged the two conventions.
    assert!(
        summary["warnings"]
            .as_array()
            .map(|a| a.is_empty())
            .unwrap_or(true)
    );
}

/// Without aliases, a BED whose chroms don't overlap the VCFs should
/// surface a chrom_naming_mismatch warning in the report.
#[test]
fn chrom_naming_mismatch_warning_appears_in_report() {
    let tmp = tempfile::tempdir().unwrap();
    let fa = tmp.path().join("ref.fa");
    let golden = tmp.path().join("golden.vcf");
    let called = tmp.path().join("called.vcf");
    let mbed = tmp.path().join("mutation.bed");
    let out = tmp.path().join("report");
    write_fasta(&fa, "chr1", &"A".repeat(1000));
    write_vcf(&golden, &["chr1\t500\t.\tA\tC\t60\tPASS\t.\tGT\t0/1"]);
    write_vcf(&called, &[]);
    // BED uses the `1` convention; VCFs use `chr1`; no aliases configured.
    write_bed(&mbed, &[("1", 0, 100)]);

    let yaml = tmp.path().join("cfg.yml");
    write_text(
        &yaml,
        &format!(
            "golden_vcf: {}\ncalled_vcf: {}\nreference: {}\noutput_dir: {}\n\
             overwrite_output: true\nmutation_bed: {}\n",
            golden.display(),
            called.display(),
            fa.display(),
            out.display(),
            mbed.display(),
        ),
    );
    rneat()
        .args(["compare-vcfs", "-c"])
        .arg(&yaml)
        .assert()
        .success();
    let summary = load_summary(&out);
    let warnings = summary["warnings"].as_array().expect("warnings array");
    assert!(
        warnings.iter().any(|w| w["bed_label"] == "mutation_bed"),
        "expected mutation_bed warning, got: {warnings:?}"
    );
}

/// The TXT report should include an "FN attribution" section and a
/// "Warnings" section when warnings are present.
#[test]
fn txt_report_renders_attribution_and_warnings() {
    let tmp = tempfile::tempdir().unwrap();
    let fa = tmp.path().join("ref.fa");
    let golden = tmp.path().join("golden.vcf");
    let called = tmp.path().join("called.vcf");
    let mbed = tmp.path().join("mutation.bed");
    let out = tmp.path().join("report");
    write_fasta(&fa, "chr1", &"A".repeat(1000));
    write_vcf(&golden, &["chr1\t500\t.\tA\tC\t60\tPASS\t.\tGT\t0/1"]);
    write_vcf(&called, &[]);
    // Use non-matching BED chrom to trigger both an attribution count and a
    // naming-mismatch warning.
    write_bed(&mbed, &[("not_a_real_contig", 0, 100)]);

    let yaml = tmp.path().join("cfg.yml");
    write_text(
        &yaml,
        &format!(
            "golden_vcf: {}\ncalled_vcf: {}\nreference: {}\noutput_dir: {}\n\
             overwrite_output: true\nmutation_bed: {}\n",
            golden.display(),
            called.display(),
            fa.display(),
            out.display(),
            mbed.display(),
        ),
    );
    rneat()
        .args(["compare-vcfs", "-c"])
        .arg(&yaml)
        .assert()
        .success();
    let txt = std::fs::read_to_string(out.join("comparison_summary.txt")).unwrap();
    assert!(
        txt.contains("FN attribution"),
        "TXT missing FN attribution section:\n{txt}"
    );
    assert!(
        txt.contains("Warnings"),
        "TXT missing Warnings section:\n{txt}"
    );
}

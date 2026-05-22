//! Phase 2 integration tests: equivalence sweep on synthetic indels.
//!
//! Each test writes a tiny custom FASTA and a golden + called VCF pair where
//! the called VCF intentionally uses a different denotation of the same edit.
//! The compare-vcfs binary should classify these as TP after the sweep.

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

fn write_yaml(path: &Path, content: &str) {
    let mut f = std::fs::File::create(path).unwrap();
    write!(f, "{content}").unwrap();
}

fn load_summary(dir: &Path) -> serde_json::Value {
    let raw = std::fs::read_to_string(dir.join("comparison_summary.json")).unwrap();
    serde_json::from_str(&raw).unwrap()
}

fn cfg_yaml(out: &Path, golden: &Path, called: &Path, reference: &Path, extra: &str) -> String {
    format!(
        "golden_vcf: {}\ncalled_vcf: {}\nreference: {}\noutput_dir: {}\noverwrite_output: true\n{extra}",
        golden.display(),
        called.display(),
        reference.display(),
        out.display(),
    )
}

/// Reference contains `ACGTACGTACGT`. A `CGTA` deletion can be denoted at
/// either pos 1 (REF=ACGTA, ALT=A) or pos 5 (REF=ACGTA, ALT=A). With the
/// equivalence sweep enabled, the two should match as TP.
#[test]
fn equivalence_catches_left_vs_right_aligned_deletion() {
    let tmp = tempfile::tempdir().unwrap();
    let fa = tmp.path().join("ref.fa");
    let golden = tmp.path().join("golden.vcf");
    let called = tmp.path().join("called.vcf");
    let out = tmp.path().join("report");
    write_fasta(&fa, "chrTest", "ACGTACGTACGT");
    // Golden VCF: deletion described starting at pos 5
    write_vcf(&golden, &["chrTest\t5\t.\tACGTA\tA\t60\tPASS\t.\tGT\t0/1"]);
    // Called VCF: same deletion left-aligned at pos 1
    write_vcf(&called, &["chrTest\t1\t.\tACGTA\tA\t60\tPASS\t.\tGT\t0/1"]);

    let yaml = tmp.path().join("cfg.yml");
    write_yaml(&yaml, &cfg_yaml(&out, &golden, &called, &fa, ""));

    rneat()
        .args(["compare-vcfs", "-c"])
        .arg(&yaml)
        .assert()
        .success();

    let summary = load_summary(&out);
    assert_eq!(summary["totals"]["tp"], 1, "{:?}", summary);
    assert_eq!(summary["totals"]["fn_"], 0);
    assert_eq!(summary["totals"]["fp"], 0);
    assert_eq!(summary["totals"]["equivalents_promoted"], 1);
}

/// With `fast: true`, the same input should classify as 1 FN + 1 FP (no
/// sweep). This confirms the fast knob disables the equivalence pass and
/// confirms the sweep was the thing closing the gap above.
#[test]
fn fast_mode_skips_sweep_so_denotation_differs_as_fp_fn() {
    let tmp = tempfile::tempdir().unwrap();
    let fa = tmp.path().join("ref.fa");
    let golden = tmp.path().join("golden.vcf");
    let called = tmp.path().join("called.vcf");
    let out = tmp.path().join("report");
    write_fasta(&fa, "chrTest", "ACGTACGTACGT");
    write_vcf(&golden, &["chrTest\t5\t.\tACGTA\tA\t60\tPASS\t.\tGT\t0/1"]);
    write_vcf(&called, &["chrTest\t1\t.\tACGTA\tA\t60\tPASS\t.\tGT\t0/1"]);

    let yaml = tmp.path().join("cfg.yml");
    write_yaml(
        &yaml,
        &cfg_yaml(&out, &golden, &called, &fa, "fast: true\n"),
    );

    rneat()
        .args(["compare-vcfs", "-c"])
        .arg(&yaml)
        .assert()
        .success();

    let summary = load_summary(&out);
    assert_eq!(summary["totals"]["tp"], 0);
    assert_eq!(summary["totals"]["fn_"], 1);
    assert_eq!(summary["totals"]["fp"], 1);
    assert_eq!(summary["totals"]["equivalents_promoted"], 0);
    assert_eq!(summary["inputs"]["fast"], true);
}

/// A compound SNP pair (two FPs, one at each of two adjacent positions)
/// producing the same alt sequence as a single multi-base substitution FN
/// should be caught by the sweep.
#[test]
fn equivalence_catches_compound_snps_vs_multibase_substitution() {
    let tmp = tempfile::tempdir().unwrap();
    let fa = tmp.path().join("ref.fa");
    let golden = tmp.path().join("golden.vcf");
    let called = tmp.path().join("called.vcf");
    let out = tmp.path().join("report");
    write_fasta(&fa, "chrTest", "ACGTACGTACGT");
    // Golden: pos 3 GT → AA (one MNP-style record)
    write_vcf(&golden, &["chrTest\t3\t.\tGT\tAA\t60\tPASS\t.\tGT\t0/1"]);
    // Called: pos 3 G→A AND pos 4 T→A (two SNPs)
    write_vcf(
        &called,
        &[
            "chrTest\t3\t.\tG\tA\t60\tPASS\t.\tGT\t0/1",
            "chrTest\t4\t.\tT\tA\t60\tPASS\t.\tGT\t0/1",
        ],
    );

    let yaml = tmp.path().join("cfg.yml");
    write_yaml(&yaml, &cfg_yaml(&out, &golden, &called, &fa, ""));

    rneat()
        .args(["compare-vcfs", "-c"])
        .arg(&yaml)
        .assert()
        .success();

    let summary = load_summary(&out);
    // 1 FN got promoted; the 2 FPs got consumed.
    assert_eq!(summary["totals"]["tp"], 1, "{summary:#?}");
    assert_eq!(summary["totals"]["fn_"], 0);
    assert_eq!(summary["totals"]["fp"], 0);
    assert_eq!(summary["totals"]["equivalents_promoted"], 1);
}

/// Genuinely different variants must not be equivalence-matched.
#[test]
fn unrelated_variants_remain_fp_and_fn() {
    let tmp = tempfile::tempdir().unwrap();
    let fa = tmp.path().join("ref.fa");
    let golden = tmp.path().join("golden.vcf");
    let called = tmp.path().join("called.vcf");
    let out = tmp.path().join("report");
    write_fasta(&fa, "chrTest", "ACGTACGTACGT");
    write_vcf(&golden, &["chrTest\t3\t.\tG\tA\t60\tPASS\t.\tGT\t0/1"]);
    write_vcf(&called, &["chrTest\t8\t.\tT\tC\t60\tPASS\t.\tGT\t0/1"]);

    let yaml = tmp.path().join("cfg.yml");
    write_yaml(&yaml, &cfg_yaml(&out, &golden, &called, &fa, ""));

    rneat()
        .args(["compare-vcfs", "-c"])
        .arg(&yaml)
        .assert()
        .success();

    let summary = load_summary(&out);
    assert_eq!(summary["totals"]["tp"], 0);
    assert_eq!(summary["totals"]["fn_"], 1);
    assert_eq!(summary["totals"]["fp"], 1);
    assert_eq!(summary["totals"]["equivalents_promoted"], 0);
}

/// The TXT report should announce the equivalence window in its header and
/// the promoted-count line. (Smoke test on string contents.)
#[test]
fn txt_report_mentions_equivalence_window_and_promoted_count() {
    let tmp = tempfile::tempdir().unwrap();
    let fa = tmp.path().join("ref.fa");
    let golden = tmp.path().join("golden.vcf");
    let called = tmp.path().join("called.vcf");
    let out = tmp.path().join("report");
    write_fasta(&fa, "chrTest", "ACGTACGTACGT");
    write_vcf(&golden, &["chrTest\t5\t.\tACGTA\tA\t60\tPASS\t.\tGT\t0/1"]);
    write_vcf(&called, &["chrTest\t1\t.\tACGTA\tA\t60\tPASS\t.\tGT\t0/1"]);
    let yaml = tmp.path().join("cfg.yml");
    write_yaml(&yaml, &cfg_yaml(&out, &golden, &called, &fa, ""));
    rneat()
        .args(["compare-vcfs", "-c"])
        .arg(&yaml)
        .assert()
        .success();
    let txt = std::fs::read_to_string(out.join("comparison_summary.txt")).unwrap();
    assert!(
        txt.contains("Equivalence window: ±50 bp"),
        "txt missing window header:\n{txt}"
    );
    assert!(
        txt.contains("promoted by equivalence sweep: 1"),
        "txt missing promoted line:\n{txt}"
    );
}

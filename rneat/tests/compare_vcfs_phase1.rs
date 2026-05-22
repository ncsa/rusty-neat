//! Integration tests for `rneat compare-vcfs` Phase 1.
//!
//! Exercises the binary against tiny hand-crafted golden + called VCF pairs
//! and asserts the parsed `comparison_summary.json` contents match expected
//! TP / FN / FP counts.

mod common;

use common::rneat;
use std::io::Write;
use std::path::{Path, PathBuf};

fn h1n1_reference() -> PathBuf {
    PathBuf::from(format!(
        "{}/test_data/references/H1N1.fa",
        env!("CARGO_MANIFEST_DIR"),
    ))
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

#[test]
fn compare_vcfs_perfect_match_yields_all_tp() {
    let tmp = tempfile::tempdir().unwrap();
    let golden = tmp.path().join("golden.vcf");
    let called = tmp.path().join("called.vcf");
    let out_dir = tmp.path().join("report");
    let body = &[
        "H1N1_HA\t100\t.\tA\tC\t60\tPASS\t.\tGT\t0/1",
        "H1N1_HA\t200\t.\tG\tT\t60\tPASS\t.\tGT\t1/1",
    ];
    write_vcf(&golden, body);
    write_vcf(&called, body);

    let yaml = tmp.path().join("cfg.yml");
    write_yaml(
        &yaml,
        &format!(
            "golden_vcf: {}\ncalled_vcf: {}\nreference: {}\noutput_dir: {}\noverwrite_output: true\n",
            golden.display(),
            called.display(),
            h1n1_reference().display(),
            out_dir.display(),
        ),
    );

    rneat()
        .args(["compare-vcfs", "-c"])
        .arg(&yaml)
        .assert()
        .success();

    let summary = load_summary(&out_dir);
    assert_eq!(summary["schema_version"], "1.0.0");
    assert_eq!(summary["totals"]["tp"], 2);
    assert_eq!(summary["totals"]["fn_"], 0);
    assert_eq!(summary["totals"]["fp"], 0);
    assert_eq!(summary["metrics"]["precision"], 1.0);
    assert_eq!(summary["metrics"]["recall"], 1.0);
    assert_eq!(summary["metrics"]["f1"], 1.0);
}

#[test]
fn compare_vcfs_mixed_classification() {
    // Golden: (100,A,C), (200,G,T)
    // Called: (100,A,C), (300,T,A)
    // → 1 TP, 1 FN (200), 1 FP (300).
    let tmp = tempfile::tempdir().unwrap();
    let golden = tmp.path().join("golden.vcf");
    let called = tmp.path().join("called.vcf");
    let out_dir = tmp.path().join("report");
    write_vcf(
        &golden,
        &[
            "H1N1_HA\t100\t.\tA\tC\t60\tPASS\t.\tGT\t0/1",
            "H1N1_HA\t200\t.\tG\tT\t60\tPASS\t.\tGT\t1/1",
        ],
    );
    write_vcf(
        &called,
        &[
            "H1N1_HA\t100\t.\tA\tC\t60\tPASS\t.\tGT\t0/1",
            "H1N1_HA\t300\t.\tT\tA\t60\tPASS\t.\tGT\t0/1",
        ],
    );

    let yaml = tmp.path().join("cfg.yml");
    write_yaml(
        &yaml,
        &format!(
            "golden_vcf: {}\ncalled_vcf: {}\nreference: {}\noutput_dir: {}\noverwrite_output: true\n",
            golden.display(),
            called.display(),
            h1n1_reference().display(),
            out_dir.display(),
        ),
    );

    rneat()
        .args(["compare-vcfs", "-c"])
        .arg(&yaml)
        .assert()
        .success();

    let summary = load_summary(&out_dir);
    assert_eq!(summary["totals"]["tp"], 1);
    assert_eq!(summary["totals"]["fn_"], 1);
    assert_eq!(summary["totals"]["fp"], 1);
    assert_eq!(summary["per_contig"]["H1N1_HA"]["tp"], 1);
    assert_eq!(summary["per_contig"]["H1N1_HA"]["fn_"], 1);
    assert_eq!(summary["per_contig"]["H1N1_HA"]["fp"], 1);
}

#[test]
fn compare_vcfs_writes_txt_report_alongside_json() {
    let tmp = tempfile::tempdir().unwrap();
    let golden = tmp.path().join("golden.vcf");
    let called = tmp.path().join("called.vcf");
    let out_dir = tmp.path().join("report");
    write_vcf(&golden, &["H1N1_HA\t100\t.\tA\tC\t60\tPASS\t.\tGT\t0/1"]);
    write_vcf(&called, &["H1N1_HA\t100\t.\tA\tC\t60\tPASS\t.\tGT\t0/1"]);

    let yaml = tmp.path().join("cfg.yml");
    write_yaml(
        &yaml,
        &format!(
            "golden_vcf: {}\ncalled_vcf: {}\nreference: {}\noutput_dir: {}\noverwrite_output: true\n",
            golden.display(),
            called.display(),
            h1n1_reference().display(),
            out_dir.display(),
        ),
    );

    rneat()
        .args(["compare-vcfs", "-c"])
        .arg(&yaml)
        .assert()
        .success();

    let txt = std::fs::read_to_string(out_dir.join("comparison_summary.txt")).unwrap();
    assert!(txt.contains("True positives  (TP): 1"));
    assert!(txt.contains("Schema version: 1.0.0"));
}

#[test]
fn compare_vcfs_contigs_simulated_filter_drops_unrelated_calls() {
    // Golden has H1N1_HA only; called has both H1N1_HA and H1N1_NA.
    // With contigs_simulated=H1N1_HA, the H1N1_NA "FP" should not be counted —
    // it's outside the simulated contig set.
    let tmp = tempfile::tempdir().unwrap();
    let golden = tmp.path().join("golden.vcf");
    let called = tmp.path().join("called.vcf");
    let out_dir = tmp.path().join("report");
    write_vcf(&golden, &["H1N1_HA\t100\t.\tA\tC\t60\tPASS\t.\tGT\t0/1"]);
    write_vcf(
        &called,
        &[
            "H1N1_HA\t100\t.\tA\tC\t60\tPASS\t.\tGT\t0/1",
            "H1N1_NA\t50\t.\tG\tA\t60\tPASS\t.\tGT\t0/1",
        ],
    );

    let yaml = tmp.path().join("cfg.yml");
    write_yaml(
        &yaml,
        &format!(
            "golden_vcf: {}\ncalled_vcf: {}\nreference: {}\noutput_dir: {}\n\
             overwrite_output: true\ncontigs_simulated: H1N1_HA\n",
            golden.display(),
            called.display(),
            h1n1_reference().display(),
            out_dir.display(),
        ),
    );

    rneat()
        .args(["compare-vcfs", "-c"])
        .arg(&yaml)
        .assert()
        .success();

    let summary = load_summary(&out_dir);
    assert_eq!(summary["totals"]["tp"], 1);
    assert_eq!(summary["totals"]["fn_"], 0);
    assert_eq!(summary["totals"]["fp"], 0);
    assert_eq!(summary["skipped_called"]["outside_simulated_contigs"], 1);
}

#[test]
fn compare_vcfs_non_pass_filtered_by_default() {
    // Both VCFs contain a LowQual record on top of a PASS record. With
    // include_filtered=false (default), the LowQual records are dropped from
    // both — leaving a 1 TP comparison plus 1 + 1 skipped/filtered.
    let tmp = tempfile::tempdir().unwrap();
    let golden = tmp.path().join("golden.vcf");
    let called = tmp.path().join("called.vcf");
    let out_dir = tmp.path().join("report");
    write_vcf(
        &golden,
        &[
            "H1N1_HA\t100\t.\tA\tC\t60\tPASS\t.\tGT\t0/1",
            "H1N1_HA\t200\t.\tG\tT\t10\tLowQual\t.\tGT\t0/1",
        ],
    );
    write_vcf(
        &called,
        &[
            "H1N1_HA\t100\t.\tA\tC\t60\tPASS\t.\tGT\t0/1",
            "H1N1_HA\t200\t.\tG\tT\t10\tLowQual\t.\tGT\t0/1",
        ],
    );

    let yaml = tmp.path().join("cfg.yml");
    write_yaml(
        &yaml,
        &format!(
            "golden_vcf: {}\ncalled_vcf: {}\nreference: {}\noutput_dir: {}\noverwrite_output: true\n",
            golden.display(),
            called.display(),
            h1n1_reference().display(),
            out_dir.display(),
        ),
    );

    rneat()
        .args(["compare-vcfs", "-c"])
        .arg(&yaml)
        .assert()
        .success();

    let summary = load_summary(&out_dir);
    assert_eq!(summary["totals"]["tp"], 1);
    assert_eq!(summary["totals"]["fn_"], 0);
    assert_eq!(summary["totals"]["fp"], 0);
    assert_eq!(summary["skipped_golden"]["filtered"], 1);
    assert_eq!(summary["skipped_called"]["filtered"], 1);
}

#[test]
fn compare_vcfs_missing_config_field_fails_nonzero() {
    let tmp = tempfile::tempdir().unwrap();
    let yaml = tmp.path().join("cfg.yml");
    write_yaml(&yaml, "golden_vcf: /nonexistent/path.vcf\n");
    rneat()
        .args(["compare-vcfs", "-c"])
        .arg(&yaml)
        .assert()
        .failure();
}

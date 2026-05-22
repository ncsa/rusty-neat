//! End-to-end pipeline tests: train a model with one rneat subcommand, then consume it
//! from another. These exercise the real binary boundary (CLI → config parser →
//! model file → consumer) — the kind of wiring breakage that unit tests will miss
//! because they call internal functions directly.

mod common;

use common::{
    GenReadsConfig, fresh_workdir, h1n1_reference, read_gzip_fastq_lines, rneat,
    write_gen_seq_error_model_config, write_iupac_fasta, write_synthetic_fastq_gz,
};
use std::collections::HashSet;

/// Convert Phred-encoded quality bytes back to raw Q-scores under Phred+33.
fn qual_bytes_to_scores(line: &str) -> Vec<u8> {
    line.bytes().map(|b| b.saturating_sub(33)).collect()
}

/// Extract every quality line from a multi-stream gzipped FASTQ. Quality lines are the
/// 4th of every 4-line record group.
fn quality_lines(fastq_path: &std::path::Path) -> Vec<String> {
    let lines = read_gzip_fastq_lines(fastq_path);
    assert!(
        lines.len() % 4 == 0,
        "FASTQ line count must be a multiple of 4, got {}",
        lines.len(),
    );
    lines
        .into_iter()
        .enumerate()
        .filter_map(|(i, l)| if i % 4 == 3 { Some(l) } else { None })
        .collect()
}

#[test]
fn gen_reads_with_default_model_produces_well_formed_fastq() {
    // Smoke test: run gen-reads end-to-end through the real binary with all defaults
    // (no custom error model, no input VCF) and verify the output is well-formed
    // (multiple of 4 lines, seq length == qual length, names start with '@').
    let (_dir, out) = fresh_workdir();
    let config = GenReadsConfig::new(h1n1_reference(), out.clone(), "default_model");
    let yaml = config.write_yaml();

    rneat()
        .args(["gen-reads", "-c"])
        .arg(yaml.path())
        .assert()
        .success();

    let fastq = out.join("default_model_r1.fastq.gz");
    assert!(fastq.exists(), "FASTQ not produced at {fastq:?}");

    let lines = read_gzip_fastq_lines(&fastq);
    assert!(!lines.is_empty(), "FASTQ has no records");
    assert_eq!(lines.len() % 4, 0, "incomplete record at FASTQ tail");

    for chunk in lines.chunks(4) {
        assert!(
            chunk[0].starts_with('@'),
            "name line must start with @: {:?}",
            chunk[0]
        );
        assert_eq!(chunk[2], "+", "third line must be '+': {:?}", chunk[2]);
        assert_eq!(
            chunk[1].len(),
            chunk[3].len(),
            "seq/qual length mismatch in {chunk:?}",
        );
    }
}

#[test]
fn train_binned_model_then_gen_reads_emits_only_bin_qualities() {
    // Cross-binary integration: train a binned quality model with gen-seq-error-model,
    // hand it to gen-reads, then prove gen-reads honors the binning by inspecting the
    // FASTQ output. This is the user-visible promise of binned-quality support.
    let (_dir, work) = fresh_workdir();

    // Step 1: build a tiny FASTQ for training. Quality scores cycle 33..=42 — all of
    // which snap to bin 37 under the bin set we'll request.
    let train_fastq = work.join("train.fastq.gz");
    write_synthetic_fastq_gz(&train_fastq, 100, 50);

    // Step 2: train a binned sequencing error model via the real binary.
    let model_path = work.join("binned_model.json.gz");
    let bins = [2usize, 12, 23, 37];
    let model_cfg = write_gen_seq_error_model_config(&train_fastq, &model_path, Some(&bins));
    rneat()
        .args(["gen-seq-error-model", "-c"])
        .arg(model_cfg.path())
        .assert()
        .success();
    assert!(model_path.exists(), "model not produced at {model_path:?}");

    // Step 3: generate reads using the binned model.
    let mut config = GenReadsConfig::new(h1n1_reference(), work.clone(), "binned_pipeline");
    config.sequence_error_model = Some(model_path);
    let reads_cfg = config.write_yaml();
    rneat()
        .args(["gen-reads", "-c"])
        .arg(reads_cfg.path())
        .assert()
        .success();

    // Step 4: verify every emitted quality score is in the bin set.
    let fastq = work.join("binned_pipeline_r1.fastq.gz");
    let bin_set: HashSet<u8> = bins.iter().map(|&b| b as u8).collect();
    let mut total_qualities = 0usize;
    for ql in quality_lines(&fastq) {
        for score in qual_bytes_to_scores(&ql) {
            assert!(
                bin_set.contains(&score),
                "gen-reads emitted non-bin quality score {score} (bins = {bins:?})",
            );
            total_qualities += 1;
        }
    }
    assert!(total_qualities > 0, "no quality scores observed");
}

#[test]
fn gen_reads_with_iupac_reference_produces_no_iupac_in_output() {
    // Regression test: a reference containing IUPAC ambiguity codes (R/Y/M/K/S/W/H/B/V/D)
    // must be processed without crashing and must not leak any IUPAC letter into output
    // FASTQ sequence lines — every base emitted must be in [ACGTNacgtn].
    let (_dir, work) = fresh_workdir();
    let ref_path = work.join("iupac_ref.fa");
    write_iupac_fasta(&ref_path);

    let config = GenReadsConfig::new(ref_path, work.clone(), "iupac_run");
    let yaml = config.write_yaml();

    rneat()
        .args(["gen-reads", "-c"])
        .arg(yaml.path())
        .assert()
        .success();

    let fastq = work.join("iupac_run_r1.fastq.gz");
    assert!(fastq.exists(), "FASTQ not produced at {fastq:?}");

    let lines = read_gzip_fastq_lines(&fastq);
    assert!(!lines.is_empty(), "FASTQ has no records");
    assert_eq!(lines.len() % 4, 0);

    let iupac: HashSet<char> = "RYMKSWHBVDrymkswhbvd".chars().collect();
    for chunk in lines.chunks(4) {
        for c in chunk[1].chars() {
            assert!(
                !iupac.contains(&c),
                "IUPAC character {:?} leaked into read sequence: {:?}",
                c,
                chunk[1]
            );
        }
    }
}

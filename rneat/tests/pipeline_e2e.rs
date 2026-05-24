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
        lines.len().is_multiple_of(4),
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
fn gen_reads_with_symbolic_del_modulates_depth_and_round_trips_to_vcf() {
    // End-to-end contract: a symbolic <DEL> in the input VCF must (a) be
    // preserved verbatim in the output golden VCF, and (b) drive per-region
    // depth modulation in gen_reads. We park a homozygous <DEL> at
    // H1N1_HA:500-1500 (1-based, inclusive), which corresponds to 0-based
    // half-open [499, 1500). At a hom DEL the multiplier is 0, so zero reads
    // should start in that span; reads outside it should still be generated
    // at full coverage.
    use std::io::Write as _;

    let (_dir, work) = fresh_workdir();

    // Coverage modulation needs the SV's span — supply it via INFO/END.
    let input_vcf = work.join("input_sv.vcf");
    {
        let mut f = std::fs::File::create(&input_vcf).unwrap();
        writeln!(f, "##fileformat=VCFv4.2").unwrap();
        writeln!(f, "##ALT=<ID=DEL,Description=\"Deletion\">").unwrap();
        writeln!(
            f,
            "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position\">"
        )
        .unwrap();
        writeln!(
            f,
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
        )
        .unwrap();
        writeln!(f, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS").unwrap();
        // Homozygous <DEL> at H1N1_HA, 1-based POS=500, END=1500.
        writeln!(
            f,
            "H1N1_HA\t500\t.\tA\t<DEL>\t60\tPASS\tEND=1500\tGT\t1/1"
        )
        .unwrap();
    }

    let mut config = GenReadsConfig::new(h1n1_reference(), work.clone(), "sv_run");
    config.coverage = 50;
    config.produce_vcf = true;
    config.input_vcf = Some(input_vcf);
    // Drive enough de-novo mutations across the contig that, pre-fix, several
    // would have been sampled inside the DEL span. Post-fix: zero — the rate
    // is zeroed over hom-DEL spans before generate_variants runs.
    config.mutation_rate = Some(0.05);
    let yaml = config.write_yaml();
    rneat()
        .args(["gen-reads", "-c"])
        .arg(yaml.path())
        .assert()
        .success();

    // (a) Output VCF round-trip — the <DEL> should appear verbatim, with the
    // original INFO field preserved (END=1500).
    let out_vcf_gz = work.join("sv_run.vcf.gz");
    assert!(out_vcf_gz.exists(), "output VCF not produced: {out_vcf_gz:?}");
    let vcf_lines: Vec<String> = {
        use flate2::read::MultiGzDecoder;
        use std::io::{BufRead, BufReader};
        let r = BufReader::new(MultiGzDecoder::new(std::fs::File::open(&out_vcf_gz).unwrap()));
        r.lines().map(|l| l.unwrap()).collect()
    };
    let del_line = vcf_lines
        .iter()
        .find(|l| l.starts_with("H1N1_HA\t500\t") && l.contains("<DEL>"))
        .unwrap_or_else(|| panic!("expected <DEL> record in output VCF; got: {vcf_lines:?}"));
    assert!(
        del_line.contains("END=1500"),
        "expected INFO/END to round-trip; got: {del_line}"
    );
    // No de-novo SNPs should be sampled inside the deleted span: the mutation
    // rate over [POS+1, END] in 1-based (= [501, 1500]) is zeroed for hom DEL.
    let snps_inside_del: Vec<&String> = vcf_lines
        .iter()
        .filter(|l| !l.starts_with('#'))
        .filter(|l| l.starts_with("H1N1_HA\t"))
        .filter(|l| !l.contains("<DEL>"))
        .filter(|l| {
            let cols: Vec<&str> = l.split('\t').collect();
            cols.get(1)
                .and_then(|p| p.parse::<usize>().ok())
                .map(|p| (501..=1500).contains(&p))
                .unwrap_or(false)
        })
        .collect();
    assert!(
        snps_inside_del.is_empty(),
        "expected no de-novo SNPs inside the hom <DEL> span [501, 1500]; got: {snps_inside_del:?}"
    );

    // (b) Depth modulation — count read start positions on H1N1_HA. Names look
    // like `@RNEAT_generated_{contig}_{abs_start:010}_{abs_end:010}`.
    let fastq = work.join("sv_run_r1.fastq.gz");
    let lines = read_gzip_fastq_lines(&fastq);
    let mut inside_del = 0usize;
    let mut outside_del = 0usize;
    for chunk in lines.chunks(4) {
        let name = &chunk[0];
        if !name.starts_with("@RNEAT_generated_H1N1_HA_") {
            continue;
        }
        // Strip leading '@' and prefix; the remainder is "{start}_{end}".
        let rest = name.trim_start_matches("@RNEAT_generated_H1N1_HA_");
        let abs_start: usize = rest
            .split('_')
            .next()
            .and_then(|s| s.parse().ok())
            .unwrap_or_else(|| panic!("unparseable read name: {name}"));
        // DEL span (0-based half-open): [499, 1500). Use inset margins so a
        // single-base off-by-one in coverage modulation wouldn't pass.
        if (500..1450).contains(&abs_start) {
            inside_del += 1;
        } else if abs_start < 449 || abs_start >= 1550 {
            outside_del += 1;
        }
    }
    assert_eq!(
        inside_del, 0,
        "expected ZERO reads starting inside the hom <DEL> span; got {inside_del}"
    );
    assert!(
        outside_del > 200,
        "expected plenty of reads outside the DEL span at coverage=50 over ~675 unaffected \
         bases of H1N1_HA; got only {outside_del}"
    );
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

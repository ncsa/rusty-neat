use log::{info, warn};
#[cfg(test)]
use std::path::PathBuf;
use common::{
    file_tools::{
        bam_reader::read_bam_transitions,
        file_io::{read_gzip_lines, read_lines},
    },
    models::{
        quality_scores::QualityScoreModel,
        sequencing_error_model::SequencingErrorModel,
    },
    structs::transition_matrix::TransitionMatrix,
};
use common::file_tools::file_io::is_gzipped_file;
use crate::gen_seq_error_model::{
    errors::GenSeqErrorModelError,
    utils::config::RunConfiguration,
};

const MAX_SCORE: usize = 94;


/// Normalizes a raw 4×4 mismatch count matrix into a `TransitionMatrix`.
///
/// Each row is normalized independently. Rows with no observed mismatches get
/// equal probability distributed across the three off-diagonal positions.
fn build_transition_matrix_from_counts(
    counts: [[usize; 4]; 4],
) -> Result<TransitionMatrix, GenSeqErrorModelError> {
    let mut weights = [[0.0f64; 4]; 4];
    for i in 0..4 {
        let total: f64 = counts[i].iter().sum::<usize>() as f64;
        if total == 0.0 {
            for j in 0..4 {
                if i != j {
                    weights[i][j] = 1.0 / 3.0;
                }
            }
        } else {
            for j in 0..4 {
                weights[i][j] = counts[i][j] as f64 / total;
            }
            weights[i][i] = 0.0;
        }
    }
    Ok(TransitionMatrix::from(weights[0], weights[1], weights[2], weights[3])?)
}

pub fn runner(config: &RunConfiguration) -> Result<(), GenSeqErrorModelError> {
    let lines: Vec<String> = if is_gzipped_file(&config.fastq_file)? {
        read_gzip_lines(&config.fastq_file)?
            .collect::<Result<Vec<_>, _>>()?
    } else {
        read_lines(&config.fastq_file)?
            .collect::<Result<Vec<_>, _>>()?
    };

    if lines.len() < 4 {
        return Err(GenSeqErrorModelError::MalformedFastq(
            "FASTQ file has fewer than 4 lines".to_string(),
        ));
    }

    // Determine read length from the first record's quality line (index 3)
    let read_length = lines[3].len();
    if read_length == 0 {
        return Err(GenSeqErrorModelError::MalformedFastq(
            "Quality line in first record is empty".to_string(),
        ));
    }

    let qual_offset = config.qual_offset;

    let mut seed_counts = vec![0usize; MAX_SCORE];
    let mut transition_counts = vec![vec![vec![0usize; MAX_SCORE]; MAX_SCORE]; read_length - 1];
    let mut global_counts = vec![0usize; MAX_SCORE];
    let mut total_bases: usize = 0;
    let mut reads_processed: usize = 0;

    // Process all records: groups of 4 lines
    let mut i = 0;
    while i + 3 < lines.len() {
        let qual_line = &lines[i + 3];
        let scores: Vec<usize> = qual_line
            .bytes()
            .map(|b| (b as usize).saturating_sub(qual_offset).min(MAX_SCORE - 1))
            .collect();

        if scores.is_empty() {
            i += 4;
            continue;
        }

        seed_counts[scores[0]] += 1;
        for j in 1..scores.len().min(read_length) {
            transition_counts[j - 1][scores[j - 1]][scores[j]] += 1;
        }
        for &score in &scores {
            global_counts[score] += 1;
        }
        total_bases += scores.len();
        reads_processed += 1;

        if config.max_reads > 0 && reads_processed >= config.max_reads {
            break;
        }
        i += 4;
    }

    info!("Processed {} reads ({} bases)", reads_processed, total_bases);

    if total_bases == 0 {
        return Err(GenSeqErrorModelError::MalformedFastq(
            "No bases found in FASTQ file".to_string(),
        ));
    }

    // Compute error_rate
    let error_rate: f64 = (0..MAX_SCORE)
        .map(|q| 10f64.powf(-(q as f64) / 10.0) * global_counts[q] as f64)
        .sum::<f64>()
        / total_bases as f64;

    info!("Computed error_rate: {:.6}", error_rate);

    // Collect observed quality scores
    let quality_score_options: Vec<usize> = (0..MAX_SCORE)
        .filter(|&q| global_counts[q] > 0)
        .collect();

    let n_scores = quality_score_options.len();

    // Re-index seed counts
    let seed_weights: Vec<f64> = quality_score_options
        .iter()
        .map(|&q| seed_counts[q] as f64)
        .collect();

    // Re-index transition counts
    let trans_weights: Vec<Vec<Vec<f64>>> = (0..read_length - 1)
        .map(|pos| {
            (0..n_scores)
                .map(|p| {
                    let prev_raw = quality_score_options[p];
                    (0..n_scores)
                        .map(|c| {
                            let curr_raw = quality_score_options[c];
                            transition_counts[pos][prev_raw][curr_raw] as f64
                        })
                        .collect()
                })
                .collect()
        })
        .collect();

    let quality_score_model = QualityScoreModel::from_counts(
        quality_score_options,
        read_length,
        seed_weights,
        trans_weights,
    )?;

    // Determine transition matrix: TSV > BAM inference > defaults from original Python NEAT
    // (see https://github.com/ncsa/NEAT/blob/main/neat/model_sequencing_error/runner.py)
    let snp_transition_matrix: Option<TransitionMatrix> =
        if let Some(path) = &config.transition_matrix_file {
            info!("Loading SNP transition matrix from TSV: {:?}", path);
            Some(TransitionMatrix::from_tsv(path)?)
        } else if let Some(path) = &config.bam_file {
            info!("Inferring SNP transition matrix from BAM: {:?}", path);
            let counts = read_bam_transitions(path)?;
            let total_mismatches: usize = counts.iter().flatten().sum();
            if total_mismatches == 0 {
                warn!(
                    "BAM file {:?} had no MD-tagged records; using default transition matrix",
                    path
                );
                None
            } else {
                info!("Observed {} SNP mismatches across all records", total_mismatches);
                Some(build_transition_matrix_from_counts(counts)?)
            }
        } else {
            None
        };

    let model = SequencingErrorModel::from_raw_data(error_rate, quality_score_model, snp_transition_matrix)?;
    model.write_model(&config.output_file)?;

    info!("Wrote sequencing error model to {:?}", config.output_file);
    Ok(())
}

/// Write a minimal BAM file with `n_records` mapped records.
/// Each record uses `ref_bases` for the MD mismatch source and `read_bases` for the SEQ field.
/// Slices must be the same length. Pass `with_md = false` to omit the MD tag entirely.
#[cfg(test)]
fn write_test_bam(path: &PathBuf, n_records: usize, ref_bases: &[u8], read_bases: &[u8], with_md: bool) {
    use noodles::bam;
    use noodles::sam::{
        self as sam,
        alignment::{
            RecordBuf,
            io::Write as _,
            record::{
                cigar::{op::Kind, Op},
                data::field::Tag,
                Flags,
            },
            record_buf::{Cigar, Sequence, data::field::Value as BufValue},
        },
    };

    let n = read_bases.len();
    let header = sam::Header::default();

    // Build MD string: match counts interspersed with ref bases at mismatches
    let md_str: Option<String> = if with_md {
        let mut md = String::new();
        let mut match_count = 0usize;
        for (&r, &q) in ref_bases.iter().zip(read_bases.iter()) {
            if r.to_ascii_uppercase() == q.to_ascii_uppercase() {
                match_count += 1;
            } else {
                md.push_str(&match_count.to_string());
                md.push(r as char);
                match_count = 0;
            }
        }
        md.push_str(&match_count.to_string());
        Some(md)
    } else {
        None
    };

    let file = std::fs::File::create(path).unwrap();
    let mut writer = bam::io::Writer::new(file);
    writer.write_header(&header).unwrap();

    for _ in 0..n_records {
        let cigar: Cigar = [Op::new(Kind::Match, n)].into_iter().collect();
        let mut record = RecordBuf::default();
        *record.flags_mut() = Flags::empty();
        *record.cigar_mut() = cigar;
        *record.sequence_mut() = Sequence::from(read_bases);

        if let Some(ref md) = md_str {
            record
                .data_mut()
                .insert(Tag::MISMATCHED_POSITIONS, BufValue::from(md.as_str()));
        }
        writer.write_alignment_record(&header, &record).unwrap();
    }
}

#[cfg(test)]
fn make_test_fastq(path: &PathBuf, n_reads: usize, read_length: usize) {
    use std::io::Write;
    let seq: String = "ACGT".chars().cycle().take(read_length).collect();
    // Quality scores cycling through a range: I(40), J(41), K(42), etc.
    let qual: String = (0..read_length)
        .map(|i| char::from_u32(('!' as u32) + 33 + (i % 10) as u32).unwrap())
        .collect();
    let mut f = std::fs::File::create(path).unwrap();
    for i in 0..n_reads {
        writeln!(f, "@read{}\n{}\n+\n{}", i, seq, qual).unwrap();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gen_seq_error_model::utils::config::RunConfiguration;

    fn make_config(fastq: PathBuf, output: PathBuf) -> RunConfiguration {
        RunConfiguration {
            fastq_file: fastq,
            output_file: output,
            overwrite_output: true,
            max_reads: 0,
            qual_offset: 33,
            bam_file: None,
            transition_matrix_file: None,
        }
    }

    #[test]
    fn test_runner_basic() {
        let temp = tempfile::tempdir().unwrap();
        let fastq_path = temp.path().join("test.fastq");
        make_test_fastq(&fastq_path, 20, 50);
        let output_path = temp.path().join("model.json.gz");
        let config = make_config(fastq_path, output_path.clone());
        runner(&config).unwrap();
        assert!(output_path.exists());
        let model = SequencingErrorModel::from_file(&output_path).unwrap();
        let scores = model.generate_quality_scores(50, &mut common::rng::NeatRng::new_from_seed(
            &vec!["test".to_string()]
        ).unwrap()).unwrap();
        assert_eq!(scores.len(), 50);
    }

    #[test]
    fn test_runner_max_reads() {
        let temp = tempfile::tempdir().unwrap();
        let fastq_path = temp.path().join("test.fastq");
        make_test_fastq(&fastq_path, 100, 50);
        let output_path = temp.path().join("model.json.gz");
        let mut config = make_config(fastq_path, output_path.clone());
        config.max_reads = 10;
        runner(&config).unwrap();
        assert!(output_path.exists());
    }

    #[test]
    fn test_runner_empty_fastq_errors() {
        let temp = tempfile::tempdir().unwrap();
        let fastq_path = temp.path().join("empty.fastq");
        std::fs::write(&fastq_path, "").unwrap();
        let output_path = temp.path().join("model.json.gz");
        let config = make_config(fastq_path, output_path);
        let result = runner(&config);
        assert!(result.is_err());
    }

    #[test]
    fn test_runner_error_rate_accuracy() {
        // All bases at Q30: '?' = ASCII 63, score = 63 - 33 = 30
        // Expected error_rate = 10^(-30/10) = 0.001 exactly
        let temp = tempfile::tempdir().unwrap();
        let fastq_path = temp.path().join("q30.fastq");
        let seq = "A".repeat(40);
        let qual = "?".repeat(40);
        let content: String = (0..20)
            .map(|i| format!("@r{i}\n{seq}\n+\n{qual}\n"))
            .collect();
        std::fs::write(&fastq_path, content).unwrap();
        let output_path = temp.path().join("model.json.gz");
        let config = make_config(fastq_path, output_path.clone());
        runner(&config).unwrap();
        let model = SequencingErrorModel::from_file(&output_path).unwrap();
        let expected = 10f64.powf(-30.0 / 10.0);
        assert!(
            (model.error_rate() - expected).abs() < 1e-10,
            "expected error_rate {expected:.6}, got {:.6}",
            model.error_rate()
        );
    }

    #[test]
    fn test_runner_error_rate_survives_serialization() {
        let temp = tempfile::tempdir().unwrap();
        let fastq_path = temp.path().join("q20.fastq");
        let qual = "5".repeat(30);
        let seq = "C".repeat(30);
        let content: String = (0..10)
            .map(|i| format!("@r{i}\n{seq}\n+\n{qual}\n"))
            .collect();
        std::fs::write(&fastq_path, content).unwrap();
        let output_path = temp.path().join("model.json.gz");
        runner(&make_config(fastq_path, output_path.clone())).unwrap();
        let model = SequencingErrorModel::from_file(&output_path).unwrap();
        let expected = 10f64.powf(-20.0 / 10.0);
        assert!(
            (model.error_rate() - expected).abs() < 1e-10,
            "expected {expected:.6}, got {:.6}",
            model.error_rate()
        );
    }

    #[test]
    fn test_runner_zero_transition_row_handled() {
        let temp = tempfile::tempdir().unwrap();
        let fastq_path = temp.path().join("sparse.fastq");
        let content: String = (0..20)
            .map(|i| format!("@r{i}\nAAA\n+\n?=J\n"))
            .collect();
        std::fs::write(&fastq_path, content).unwrap();
        let output_path = temp.path().join("model.json.gz");
        runner(&make_config(fastq_path, output_path.clone())).unwrap();
        let model = SequencingErrorModel::from_file(&output_path).unwrap();
        let mut rng = common::rng::NeatRng::new_from_seed(&vec!["s".to_string()]).unwrap();
        let scores: Vec<usize> = (0..100)
            .map(|_| model.generate_quality_scores(3, &mut rng).unwrap()[2])
            .collect();
        assert!(
            scores.iter().any(|&s| s == 41),
            "expected Q41 to appear in last-position scores with uniform fallback; got {scores:?}"
        );
    }

    #[test]
    fn test_runner_qual_offset() {
        let temp = tempfile::tempdir().unwrap();
        let fastq_path = temp.path().join("phred64.fastq");
        let qual = "a".repeat(20);
        let seq = "G".repeat(20);
        let content: String = (0..15)
            .map(|i| format!("@r{i}\n{seq}\n+\n{qual}\n"))
            .collect();
        std::fs::write(&fastq_path, content).unwrap();
        let output_path = temp.path().join("model.json.gz");
        let mut config = make_config(fastq_path, output_path.clone());
        config.qual_offset = 64;
        runner(&config).unwrap();
        let model = SequencingErrorModel::from_file(&output_path).unwrap();
        let expected = 10f64.powf(-33.0 / 10.0);
        assert!(
            (model.error_rate() - expected).abs() < 1e-10,
            "expected {expected:.6e}, got {:.6e}",
            model.error_rate()
        );
    }

    #[test]
    fn test_runner_gzip_h1n1_r1() {
        let manifest_dir = env!("CARGO_MANIFEST_DIR");
        let fastq_path = PathBuf::from(format!("{}/test_data/H1N1_read1.fq.gz", manifest_dir));
        let temp = tempfile::tempdir().unwrap();
        let output_path = temp.path().join("h1n1_r1_model.json.gz");
        let config = make_config(fastq_path, output_path.clone());
        runner(&config).unwrap();
        assert!(output_path.exists());

        let model = SequencingErrorModel::from_file(&output_path).unwrap();
        let error_rate = model.error_rate();
        assert!(
            error_rate > 0.0001 && error_rate < 0.05,
            "error_rate {error_rate} outside expected Illumina range [0.0001, 0.05]"
        );
        let mut rng = common::rng::NeatRng::new_from_seed(&vec!["h1n1".to_string()]).unwrap();
        let scores = model.generate_quality_scores(151, &mut rng).unwrap();
        assert_eq!(scores.len(), 151);
        assert!(scores.iter().all(|&s| s <= 94), "all scores should be ≤ MAX_SCORE");
    }

    #[test]
    fn test_runner_gzip_h1n1_r2() {
        let manifest_dir = env!("CARGO_MANIFEST_DIR");
        let fastq_path = PathBuf::from(format!("{}/test_data/H1N1_read2.fq.gz", manifest_dir));
        let temp = tempfile::tempdir().unwrap();
        let output_path = temp.path().join("h1n1_r2_model.json.gz");
        let config = make_config(fastq_path, output_path.clone());
        runner(&config).unwrap();
        assert!(output_path.exists());

        let model = SequencingErrorModel::from_file(&output_path).unwrap();
        let error_rate = model.error_rate();
        assert!(
            error_rate > 0.0001 && error_rate < 0.05,
            "error_rate {error_rate} outside expected Illumina range [0.0001, 0.05]"
        );
    }

    #[test]
    fn test_runner_max_reads_reduces_processed_count() {
        let manifest_dir = env!("CARGO_MANIFEST_DIR");
        let fastq_path = PathBuf::from(format!("{}/test_data/H1N1_read1.fq.gz", manifest_dir));
        let temp = tempfile::tempdir().unwrap();
        let output_path = temp.path().join("h1n1_maxreads_model.json.gz");
        let mut config = make_config(fastq_path, output_path.clone());
        config.max_reads = 50;
        runner(&config).unwrap();
        assert!(output_path.exists());

        let model = SequencingErrorModel::from_file(&output_path).unwrap();
        let mut rng = common::rng::NeatRng::new_from_seed(&vec!["seed".to_string()]).unwrap();
        let scores = model.generate_quality_scores(151, &mut rng).unwrap();
        assert_eq!(scores.len(), 151);
    }

    #[test]
    fn test_transition_matrix_from_tsv() {
        let temp = tempfile::tempdir().unwrap();
        let fastq_path = temp.path().join("test.fastq");
        make_test_fastq(&fastq_path, 20, 50);
        let output_path = temp.path().join("model.json.gz");

        // Write a valid 4×4 TSV with a header line
        let tsv_path = temp.path().join("matrix.tsv");
        std::fs::write(
            &tsv_path,
            "from\\to\tA\tC\tG\tT\n\
             0.0\t0.5\t0.3\t0.2\n\
             0.5\t0.0\t0.3\t0.2\n\
             0.4\t0.3\t0.0\t0.3\n\
             0.3\t0.3\t0.4\t0.0\n",
        ).unwrap();

        let config = RunConfiguration {
            fastq_file: fastq_path,
            output_file: output_path.clone(),
            overwrite_output: true,
            max_reads: 0,
            qual_offset: 33,
            bam_file: None,
            transition_matrix_file: Some(tsv_path),
        };
        runner(&config).unwrap();
        assert!(output_path.exists());
    }

    #[test]
    fn test_runner_with_bam_md_tags() {
        // ref=AAAA, read=CCCC → every base is an A→C mismatch.
        // Verifies that runner accepts a bam_file and completes successfully.
        let temp = tempfile::tempdir().unwrap();
        let fastq_path = temp.path().join("test.fastq");
        make_test_fastq(&fastq_path, 20, 4);
        let bam_path = temp.path().join("test.bam");
        write_test_bam(&bam_path, 10, b"AAAA", b"CCCC", true);
        let output_path = temp.path().join("model.json.gz");

        let config = RunConfiguration {
            fastq_file: fastq_path,
            output_file: output_path.clone(),
            overwrite_output: true,
            max_reads: 0,
            qual_offset: 33,
            bam_file: Some(bam_path),
            transition_matrix_file: None,
        };
        runner(&config).unwrap();
        assert!(output_path.exists());
        SequencingErrorModel::from_file(&output_path).unwrap();
    }

    #[test]
    fn test_runner_bam_no_md_tags_falls_back() {
        // BAM records without MD tags → runner succeeds and uses the default matrix.
        let temp = tempfile::tempdir().unwrap();
        let fastq_path = temp.path().join("test.fastq");
        make_test_fastq(&fastq_path, 20, 4);
        let bam_path = temp.path().join("nomd.bam");
        write_test_bam(&bam_path, 5, b"ACGT", b"TGCA", false);
        let output_path = temp.path().join("model.json.gz");

        let config = RunConfiguration {
            fastq_file: fastq_path,
            output_file: output_path.clone(),
            overwrite_output: true,
            max_reads: 0,
            qual_offset: 33,
            bam_file: Some(bam_path),
            transition_matrix_file: None,
        };
        runner(&config).unwrap();
        assert!(output_path.exists());
    }

    #[test]
    fn test_tsv_takes_precedence_over_bam() {
        // When both transition_matrix_file and bam_file are set, the TSV wins.
        // The BAM has A→C biased mismatches; the TSV has a uniform matrix.
        // We can't easily inspect the matrix after serialization, so just verify
        // that the runner completes and writes a valid model.
        let temp = tempfile::tempdir().unwrap();
        let fastq_path = temp.path().join("test.fastq");
        make_test_fastq(&fastq_path, 20, 4);
        let bam_path = temp.path().join("test.bam");
        write_test_bam(&bam_path, 10, b"AAAA", b"CCCC", true);
        let tsv_path = temp.path().join("matrix.tsv");
        std::fs::write(
            &tsv_path,
            "A\tC\tG\tT\n\
             0.0\t0.5\t0.3\t0.2\n\
             0.5\t0.0\t0.3\t0.2\n\
             0.4\t0.3\t0.0\t0.3\n\
             0.3\t0.3\t0.4\t0.0\n",
        )
        .unwrap();
        let output_path = temp.path().join("model.json.gz");

        let config = RunConfiguration {
            fastq_file: fastq_path,
            output_file: output_path.clone(),
            overwrite_output: true,
            max_reads: 0,
            qual_offset: 33,
            bam_file: Some(bam_path),
            transition_matrix_file: Some(tsv_path),
        };
        runner(&config).unwrap();
        assert!(output_path.exists());
        SequencingErrorModel::from_file(&output_path).unwrap();
    }

    #[test]
    fn test_build_transition_matrix_from_counts_uniform_fallback() {
        use common::structs::nucleotides::Nucleotide;
        // A row with all zeros should produce uniform off-diagonal distribution —
        // verify by sampling and checking that A is never the result and all three
        // other bases are reachable.
        let mut counts = [[0usize; 4]; 4];
        counts[1][0] = 10; // C→A
        counts[1][2] = 5;  // C→G
        counts[1][3] = 5;  // C→T
        let tm = build_transition_matrix_from_counts(counts).unwrap();

        // Sample the A row at several evenly-spaced random values
        let a_dist = &tm[&Nucleotide::A];
        let mut seen = std::collections::HashSet::new();
        for k in 1..=99 {
            let r = k as f64 / 100.0;
            let result = a_dist.sample(r).unwrap();
            assert_ne!(result, Nucleotide::A, "A→A self-transition should be impossible");
            seen.insert(result as usize);
        }
        assert_eq!(seen.len(), 3, "all three off-diagonal bases should be reachable");
    }
}
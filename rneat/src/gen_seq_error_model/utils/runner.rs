use log::info;
#[cfg(test)]
use std::path::PathBuf;
use common::{
    file_tools::file_io::{read_gzip_lines, read_lines},
    models::{
        quality_scores::QualityScoreModel,
        sequencing_error_model::SequencingErrorModel,
    },
};
use crate::gen_seq_error_model::{
    errors::GenSeqErrorModelError,
    utils::config::RunConfiguration,
};

const MAX_SCORE: usize = 94;

pub fn runner(config: &RunConfiguration) -> Result<(), GenSeqErrorModelError> {
    let is_gz = config
        .fastq_file
        .extension()
        .map(|e| e == "gz")
        .unwrap_or(false);

    let lines: Vec<String> = if is_gz {
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

    let model = SequencingErrorModel::from_raw_data(error_rate, quality_score_model)?;
    model.write_model(&config.output_file)?;

    info!("Wrote sequencing error model to {:?}", config.output_file);
    Ok(())
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
        // just check it round-trips — error_rate must be finite and positive
        let scores = model.generate_quality_scores(50, &mut simple_rng::NeatRng::new_from_seed(
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
        // Verify error_rate is preserved exactly after write + reload.
        let temp = tempfile::tempdir().unwrap();
        let fastq_path = temp.path().join("q20.fastq");
        // Q20: 'A' + 20 = chr(53) = '5', score = 53 - 33 = 20
        let qual = "5".repeat(30);
        let seq = "C".repeat(30);
        let content: String = (0..10)
            .map(|i| format!("@r{i}\n{seq}\n+\n{qual}\n"))
            .collect();
        std::fs::write(&fastq_path, content).unwrap();
        let output_path = temp.path().join("model.json.gz");
        runner(&make_config(fastq_path, output_path.clone())).unwrap();
        let model = SequencingErrorModel::from_file(&output_path).unwrap();
        let expected = 10f64.powf(-20.0 / 10.0); // 0.01
        assert!(
            (model.error_rate() - expected).abs() < 1e-10,
            "expected {expected:.6}, got {:.6}",
            model.error_rate()
        );
    }

    #[test]
    fn test_runner_zero_transition_row_handled() {
        // Q41 ('J' = ASCII 74, score 41) only appears at position 2 of each read.
        // It never acts as a prev_score, so its transition rows would be all-zero
        // without the uniform fallback fix. Model must build and generate scores cleanly.
        let temp = tempfile::tempdir().unwrap();
        let fastq_path = temp.path().join("sparse.fastq");
        // Scores per read: [30, 28, 41] → chars [?, =, J]
        let content: String = (0..20)
            .map(|i| format!("@r{i}\nAAA\n+\n?=J\n"))
            .collect();
        std::fs::write(&fastq_path, content).unwrap();
        let output_path = temp.path().join("model.json.gz");
        runner(&make_config(fastq_path, output_path.clone())).unwrap();
        let model = SequencingErrorModel::from_file(&output_path).unwrap();
        let mut rng = simple_rng::NeatRng::new_from_seed(&vec!["s".to_string()]).unwrap();
        // 100 samples; without the fix all would be ≤ 28 (lowest score); with fix Q41 can appear
        let scores: Vec<usize> = (0..100)
            .map(|_| model.generate_quality_scores(3, &mut rng).unwrap()[2])
            .collect();
        // At least some last-position scores should be Q41 (uniform fallback is working)
        assert!(
            scores.iter().any(|&s| s == 41),
            "expected Q41 to appear in last-position scores with uniform fallback; got {scores:?}"
        );
    }

    #[test]
    fn test_runner_qual_offset() {
        // Phred+64 encoding: 'a' = ASCII 97, score = 97 - 64 = 33
        // Expected error_rate = 10^(-33/10) ≈ 5.012e-4
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
        // H1N1 FASTQ is real Illumina data: error rate should be in [0.0001, 0.05]
        let error_rate = model.error_rate();
        assert!(
            error_rate > 0.0001 && error_rate < 0.05,
            "error_rate {error_rate} outside expected Illumina range [0.0001, 0.05]"
        );
        // Quality scores must be generable at the FASTQ read length (151)
        let mut rng = simple_rng::NeatRng::new_from_seed(&vec!["h1n1".to_string()]).unwrap();
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
        // With max_reads=50, the model should still write successfully and produce
        // a valid model from a subset of the H1N1 reads (653 total).
        let manifest_dir = env!("CARGO_MANIFEST_DIR");
        let fastq_path = PathBuf::from(format!("{}/test_data/H1N1_read1.fq.gz", manifest_dir));
        let temp = tempfile::tempdir().unwrap();
        let output_path = temp.path().join("h1n1_maxreads_model.json.gz");
        let mut config = make_config(fastq_path, output_path.clone());
        config.max_reads = 50;
        runner(&config).unwrap();
        assert!(output_path.exists());

        // Model should still be valid and usable
        let model = SequencingErrorModel::from_file(&output_path).unwrap();
        let mut rng = simple_rng::NeatRng::new_from_seed(&vec!["seed".to_string()]).unwrap();
        let scores = model.generate_quality_scores(151, &mut rng).unwrap();
        assert_eq!(scores.len(), 151);
    }
}

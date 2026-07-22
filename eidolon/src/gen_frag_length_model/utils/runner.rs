use crate::gen_frag_length_model::{
    errors::GenFragLengthModelError, utils::config::RunConfiguration,
};
use eidolon_core::{
    file_tools::bam_reader::read_fragment_lengths, models::fragment_length::FragmentLengthModel,
};
use log::info;
use std::collections::HashMap;

// Only consider fragment lengths up to this many median absolute deviations above the median.
// Mirrors the FILTER_MEDDEV_M constant from Python NEAT.
const FILTER_MEDDEV_M: f64 = 10.0;

pub fn runner(config: &RunConfiguration) -> Result<(), GenFragLengthModelError> {
    info!("Reading fragment lengths from {:?}", config.input_file);
    let tlens = read_fragment_lengths(&config.input_file)?;
    run_from_tlens(tlens, config.min_reads, &config.output_file)
}

/// Builds and writes a fragment length model from a pre-collected list of
/// template lengths (e.g. produced by `FragLengthObserver` during a shared
/// BAM walk in the unified `gen-bam-models` runner). Applies MAD-based
/// outlier filtering and rare-length pruning, fits a normal distribution,
/// and writes the gzipped JSON model.
pub fn run_from_tlens(
    tlens: Vec<usize>,
    min_reads: usize,
    output_file: &std::path::PathBuf,
) -> Result<(), GenFragLengthModelError> {
    if tlens.is_empty() {
        return Err(GenFragLengthModelError::EmptyData);
    }
    info!("Collected {} raw fragment lengths", tlens.len());

    let filtered = filter_lengths(tlens, min_reads);

    if filtered.is_empty() {
        return Err(GenFragLengthModelError::FilteredToEmpty);
    }
    info!(
        "Retained {} fragment lengths after filtering",
        filtered.len()
    );

    let mean = compute_mean(&filtered);
    let std_dev = compute_std_dev(&filtered, mean);
    info!(
        "Fragment length model: mean={:.1}, std_dev={:.1}",
        mean, std_dev
    );

    let model = FragmentLengthModel::new_normal(mean, std_dev)?;
    model.write_file(output_file)?;
    info!("Wrote fragment length model to {:?}", output_file);
    Ok(())
}

/// Filters a raw list of template lengths, removing outliers and rare lengths.
///
/// If `min_reads` is 0, the original list is returned unchanged.
/// Otherwise, any fragment length l is retained only when:
///   - l > 0
///   - l <= median + FILTER_MEDDEV_M * MAD  (outlier ceiling)
///   - count(l) >= min_reads
fn filter_lengths(mut tlens: Vec<usize>, min_reads: usize) -> Vec<usize> {
    if min_reads == 0 {
        return tlens;
    }

    tlens.sort_unstable();

    let median = median_f64(&tlens);

    // MAD = median(|x - median|)
    let mut abs_devs: Vec<usize> = tlens
        .iter()
        .map(|&x| (x as f64 - median).abs().round() as usize)
        .collect();
    abs_devs.sort_unstable();
    let mad = median_f64(&abs_devs);

    let ceiling = median + FILTER_MEDDEV_M * mad;

    // Count occurrences of each length
    let mut counts: HashMap<usize, usize> = HashMap::new();
    for &l in &tlens {
        *counts.entry(l).or_default() += 1;
    }

    let mut output = Vec::new();
    let mut unique: Vec<usize> = counts.keys().copied().collect();
    unique.sort_unstable();
    for l in unique {
        if l > 0 && (l as f64) <= ceiling && counts[&l] >= min_reads {
            output.extend(std::iter::repeat_n(l, counts[&l]));
        }
    }
    output
}

fn median_f64(sorted: &[usize]) -> f64 {
    let n = sorted.len();
    if n == 0 {
        return 0.0;
    }
    if n.is_multiple_of(2) {
        (sorted[n / 2 - 1] + sorted[n / 2]) as f64 / 2.0
    } else {
        sorted[n / 2] as f64
    }
}

fn compute_mean(data: &[usize]) -> f64 {
    data.iter().map(|&x| x as f64).sum::<f64>() / data.len() as f64
}

fn compute_std_dev(data: &[usize], mean: f64) -> f64 {
    let variance = data
        .iter()
        .map(|&x| {
            let d = x as f64 - mean;
            d * d
        })
        .sum::<f64>()
        / data.len() as f64;
    variance.sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use eidolon_core::models::fragment_length::FragmentLengthModel;

    // ── filter_lengths ───────────────────────────────────────────────────────

    #[test]
    fn test_filter_lengths_zero_min_reads_passthrough() {
        let data = vec![100, 200, 100, 300, 50];
        let result = filter_lengths(data.clone(), 0);
        assert_eq!(result, data);
    }

    #[test]
    fn test_filter_lengths_removes_rare() {
        // 100 appears 3× (>= min_reads=2), 999 appears 1× (< min_reads=2)
        let mut data: Vec<usize> = vec![999];
        data.extend(std::iter::repeat_n(100, 3));
        let result = filter_lengths(data, 2);
        assert!(
            result.iter().all(|&x| x == 100),
            "rare length 999 should be removed"
        );
        assert_eq!(result.len(), 3);
    }

    #[test]
    fn test_filter_lengths_removes_outliers() {
        // Tight cluster around 200, plus one extreme outlier at 100_000
        let mut data: Vec<usize> = (180..=220).flat_map(|v| vec![v; 5]).collect();
        data.extend(std::iter::repeat_n(100_000usize, 5));
        let result = filter_lengths(data, 2);
        assert!(
            result.iter().all(|&x| x < 100_000),
            "outlier 100_000 should be filtered out"
        );
        assert!(!result.is_empty());
    }

    #[test]
    fn test_filter_lengths_empty_input() {
        let result = filter_lengths(vec![], 2);
        assert!(result.is_empty());
    }

    // ── statistics ───────────────────────────────────────────────────────────

    #[test]
    fn test_compute_mean() {
        let data = vec![100usize, 200, 300];
        let mean = compute_mean(&data);
        assert!((mean - 200.0).abs() < 1e-10);
    }

    #[test]
    fn test_compute_std_dev() {
        // Population std of [100, 200, 300] = sqrt(6666.66...) ≈ 81.65
        let data = vec![100usize, 200, 300];
        let mean = compute_mean(&data);
        let std = compute_std_dev(&data, mean);
        assert!((std - (20000.0f64 / 3.0).sqrt()).abs() < 1e-6);
    }

    #[test]
    fn test_median_even() {
        let data = vec![1usize, 2, 3, 4];
        assert!((median_f64(&data) - 2.5).abs() < 1e-10);
    }

    #[test]
    fn test_median_odd() {
        let data = vec![1usize, 2, 3, 4, 5];
        assert!((median_f64(&data) - 3.0).abs() < 1e-10);
    }

    // ── runner integration ────────────────────────────────────────────────────

    #[test]
    fn test_runner_with_bam() {
        let temp = tempfile::tempdir().unwrap();
        let bam_path = temp.path().join("frags.bam");
        write_test_frag_bam(
            &bam_path,
            &[150usize, 151, 152, 150, 151, 152, 150, 151, 152, 150],
        );
        let output = temp.path().join("model.json.gz");
        let config = RunConfiguration {
            input_file: bam_path,
            output_file: output.clone(),
            overwrite_output: true,
            min_reads: 2,
        };
        runner(&config).unwrap();
        assert!(output.exists());
        let model = FragmentLengthModel::discrete_from_file(&output).unwrap();
        match model {
            FragmentLengthModel::Normal { mean, st_dev } => {
                assert!(mean > 100.0 && mean < 200.0, "mean={mean}");
                assert!(st_dev >= 0.0, "st_dev={st_dev}");
            }
            _ => panic!("Expected Normal model"),
        }
    }

    #[test]
    fn test_runner_min_reads_zero_skips_filter() {
        let temp = tempfile::tempdir().unwrap();
        let bam_path = temp.path().join("frags.bam");
        // Every length appears exactly once → min_reads=2 would remove all; min_reads=0 keeps all
        write_test_frag_bam(
            &bam_path,
            &[100, 200, 300, 400, 500, 150, 250, 350, 450, 160],
        );
        let output = temp.path().join("model.json.gz");
        let config = RunConfiguration {
            input_file: bam_path,
            output_file: output.clone(),
            overwrite_output: true,
            min_reads: 0,
        };
        runner(&config).unwrap();
        assert!(output.exists());
    }

    #[test]
    fn test_runner_zero_tlen_filtered_to_empty_data() {
        // Records with TLEN=0 must be filtered out by the BAM reader (it requires tlen > 0).
        // With every record at TLEN=0 the runner sees zero usable fragments and returns
        // EmptyData — same error variant as a literally empty BAM.
        let temp = tempfile::tempdir().unwrap();
        let bam_path = temp.path().join("zero_tlen.bam");
        write_test_frag_bam(&bam_path, &[0usize; 8]);
        let output = temp.path().join("model.json.gz");
        let config = RunConfiguration {
            input_file: bam_path,
            output_file: output,
            overwrite_output: true,
            min_reads: 2,
        };
        let err = runner(&config).unwrap_err();
        assert!(
            matches!(err, GenFragLengthModelError::EmptyData),
            "expected EmptyData for all-zero-TLEN BAM, got {err:?}",
        );
    }

    #[test]
    fn test_runner_single_ended_bam_filtered_to_empty_data() {
        // A BAM where reads are not segmented (single-ended sequencing) yields no fragment
        // lengths — the reader's flag filter drops every record. Lock in EmptyData.
        let temp = tempfile::tempdir().unwrap();
        let bam_path = temp.path().join("single_ended.bam");
        write_single_ended_frag_bam(&bam_path, &[150usize, 200, 250]);
        let output = temp.path().join("model.json.gz");
        let config = RunConfiguration {
            input_file: bam_path,
            output_file: output,
            overwrite_output: true,
            min_reads: 2,
        };
        let err = runner(&config).unwrap_err();
        assert!(
            matches!(err, GenFragLengthModelError::EmptyData),
            "expected EmptyData for single-ended BAM, got {err:?}",
        );
    }

    #[test]
    fn test_runner_low_mapq_filtered_to_empty_data() {
        // Reads with MAPQ <= FRAG_FILTER_MAPQUAL (10) are dropped. With every record at
        // MAPQ=5 we expect EmptyData.
        let temp = tempfile::tempdir().unwrap();
        let bam_path = temp.path().join("low_mapq.bam");
        write_frag_bam_with_mapq(&bam_path, &[150usize, 200, 250], 5);
        let output = temp.path().join("model.json.gz");
        let config = RunConfiguration {
            input_file: bam_path,
            output_file: output,
            overwrite_output: true,
            min_reads: 2,
        };
        let err = runner(&config).unwrap_err();
        assert!(
            matches!(err, GenFragLengthModelError::EmptyData),
            "expected EmptyData for low-MAPQ BAM, got {err:?}",
        );
    }

    /// Single-ended variant of `write_test_frag_bam`: no SEGMENTED/FIRST_SEGMENT flags set,
    /// otherwise identical. Used to verify the read filter rejects single-ended BAMs.
    #[cfg(test)]
    fn write_single_ended_frag_bam(path: &std::path::PathBuf, tlens: &[usize]) {
        use noodles::bam;
        use noodles::sam::{
            self as sam,
            alignment::{
                RecordBuf,
                io::Write as _,
                record::{
                    Flags, MappingQuality,
                    cigar::{Op, op::Kind},
                },
                record_buf::{Cigar, Sequence},
            },
            header::record::value::{Map, map::ReferenceSequence},
        };
        let header = sam::Header::builder()
            .add_reference_sequence(
                b"chr1".to_vec(),
                Map::<ReferenceSequence>::new(std::num::NonZero::<usize>::new(1_000_000).unwrap()),
            )
            .build();
        let file = std::fs::File::create(path).unwrap();
        let mut writer = bam::io::Writer::new(file);
        writer.write_header(&header).unwrap();
        let seq = b"ACGT";
        for &tlen in tlens {
            let cigar: Cigar = [Op::new(Kind::Match, 4)].into_iter().collect();
            let mut record = RecordBuf::default();
            *record.flags_mut() = Flags::empty(); // single-ended: no SEGMENTED bit
            *record.cigar_mut() = cigar;
            *record.sequence_mut() = Sequence::from(seq.as_ref());
            *record.mapping_quality_mut() = Some(MappingQuality::try_from(30u8).unwrap());
            *record.reference_sequence_id_mut() = Some(0);
            *record.template_length_mut() = tlen as i32;
            writer.write_alignment_record(&header, &record).unwrap();
        }
    }

    /// Variant of `write_test_frag_bam` that lets the test choose MAPQ. Used to verify the
    /// MAPQ filter in read_fragment_lengths drops low-confidence alignments.
    #[cfg(test)]
    fn write_frag_bam_with_mapq(path: &std::path::PathBuf, tlens: &[usize], mapq: u8) {
        use noodles::bam;
        use noodles::sam::{
            self as sam,
            alignment::{
                RecordBuf,
                io::Write as _,
                record::{
                    Flags, MappingQuality,
                    cigar::{Op, op::Kind},
                },
                record_buf::{Cigar, Sequence},
            },
            header::record::value::{Map, map::ReferenceSequence},
        };
        let header = sam::Header::builder()
            .add_reference_sequence(
                b"chr1".to_vec(),
                Map::<ReferenceSequence>::new(std::num::NonZero::<usize>::new(1_000_000).unwrap()),
            )
            .build();
        let file = std::fs::File::create(path).unwrap();
        let mut writer = bam::io::Writer::new(file);
        writer.write_header(&header).unwrap();
        let seq = b"ACGT";
        for &tlen in tlens {
            let cigar: Cigar = [Op::new(Kind::Match, 4)].into_iter().collect();
            let mut record = RecordBuf::default();
            *record.flags_mut() = Flags::SEGMENTED | Flags::FIRST_SEGMENT;
            *record.cigar_mut() = cigar;
            *record.sequence_mut() = Sequence::from(seq.as_ref());
            *record.mapping_quality_mut() = Some(MappingQuality::try_from(mapq).unwrap());
            *record.reference_sequence_id_mut() = Some(0);
            *record.mate_reference_sequence_id_mut() = Some(0);
            *record.template_length_mut() = tlen as i32;
            writer.write_alignment_record(&header, &record).unwrap();
        }
    }

    #[test]
    fn test_runner_empty_bam_errors() {
        let temp = tempfile::tempdir().unwrap();
        let bam_path = temp.path().join("empty.bam");
        write_test_frag_bam(&bam_path, &[]);
        let output = temp.path().join("model.json.gz");
        let config = RunConfiguration {
            input_file: bam_path,
            output_file: output,
            overwrite_output: true,
            min_reads: 2,
        };
        assert!(runner(&config).is_err());
    }

    /// Writes a minimal BGZF BAM with paired, first-in-pair, mapq=30 records,
    /// one record per entry in `tlens`. All records are placed on reference 0.
    #[cfg(test)]
    fn write_test_frag_bam(path: &std::path::PathBuf, tlens: &[usize]) {
        use noodles::bam;
        use noodles::sam::{
            self as sam,
            alignment::{
                RecordBuf,
                io::Write as _,
                record::{
                    Flags, MappingQuality,
                    cigar::{Op, op::Kind},
                },
                record_buf::{Cigar, Sequence},
            },
            header::record::value::{Map, map::ReferenceSequence},
        };

        // Build a header with one reference sequence so refID=0 is valid
        let header = sam::Header::builder()
            .add_reference_sequence(
                b"chr1".to_vec(),
                Map::<ReferenceSequence>::new(std::num::NonZero::<usize>::new(1_000_000).unwrap()),
            )
            .build();

        let file = std::fs::File::create(path).unwrap();
        let mut writer = bam::io::Writer::new(file);
        writer.write_header(&header).unwrap();

        let seq = b"ACGT";
        for &tlen in tlens {
            let cigar: Cigar = [Op::new(Kind::Match, 4)].into_iter().collect();
            let mut record = RecordBuf::default();
            *record.flags_mut() = Flags::SEGMENTED | Flags::FIRST_SEGMENT;
            *record.cigar_mut() = cigar;
            *record.sequence_mut() = Sequence::from(seq.as_ref());
            *record.mapping_quality_mut() = Some(MappingQuality::try_from(30u8).unwrap());
            *record.reference_sequence_id_mut() = Some(0);
            *record.mate_reference_sequence_id_mut() = Some(0);
            *record.template_length_mut() = tlen as i32;
            writer.write_alignment_record(&header, &record).unwrap();
        }
    }
}

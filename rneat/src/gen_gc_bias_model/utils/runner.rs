use log::*;
use std::{collections::HashMap, path::PathBuf};

use crate::gen_gc_bias_model::{
    errors::GenGcBiasModelError,
    utils::config::{GcBiasModelParams, RunConfiguration},
};
use common::{
    file_tools::{
        bam_reader::{BamWalkFilter, CoverageObserver, walk_bam},
        fasta_stream::{FastaStream, non_n_regions},
    },
    models::gc_bias_model::GcBiasModel,
    structs::nucleotides::Nucleotide,
};

pub fn runner(path: &PathBuf) -> Result<(), GenGcBiasModelError> {
    let config = RunConfiguration::from(path)?;

    info!("Accumulating coverage from BAM {:?}", config.bam_file);
    let mut obs = CoverageObserver::default();
    let mut filter = BamWalkFilter::for_coverage();
    filter.min_mapq = config.min_mapq;
    walk_bam(&config.bam_file, &filter, &mut [&mut obs])?;
    run_from_coverage(&config.model, obs.into_by_contig())
}

/// Builds and writes a GC bias model from a pre-computed per-contig coverage
/// map (e.g. produced by `CoverageObserver` during a shared BAM walk in the
/// unified `gen-bam-models` runner). Reads the reference FASTA, sweeps
/// windows, accumulates GC% vs mean-coverage, and writes the gzipped JSON
/// model to `config.output_file`.
pub fn run_from_coverage(
    config: &GcBiasModelParams,
    mut cov_by_contig: HashMap<String, Vec<u32>>,
) -> Result<(), GenGcBiasModelError> {
    let mut gc_weight_sum = [0.0f64; 101];
    let mut gc_window_count = [0usize; 101];

    info!("Processing reference {:?}", config.reference);
    let fasta = FastaStream::open(&config.reference)?;

    for result in fasta {
        let (contig_name, raw) = result?;
        // IUPAC codes map to N here intentionally — GC bias model training works on
        // observed coverage and doesn't require stochastic base resolution.
        let sequence: Vec<Nucleotide> = raw.chars().map(Nucleotide::from).collect();
        let contig_len = sequence.len();

        if contig_len < config.window_size {
            debug!(
                "Skipping {} (length {} < window_size {})",
                contig_name, contig_len, config.window_size
            );
            continue;
        }

        let cov = match cov_by_contig.remove(&contig_name) {
            Some(d) => d,
            None => {
                debug!("No coverage data for {}, skipping", contig_name);
                continue;
            }
        };

        let regions: Vec<(usize, usize)> =
            if let Some(bed_regions) = config.bed_table.get(&contig_name) {
                bed_regions
                    .iter()
                    .map(|r| (r.start, r.end.min(contig_len)))
                    .collect()
            } else {
                non_n_regions(&sequence)
            };

        for (region_start, region_end) in regions {
            if region_end.saturating_sub(region_start) < config.window_size {
                continue;
            }
            accumulate_region(
                &sequence,
                &cov,
                region_start,
                region_end,
                config.window_size,
                config.window_stride,
                &mut gc_weight_sum,
                &mut gc_window_count,
            );
        }

        debug!("Processed {}", contig_name);
    }

    let total_windows: usize = gc_window_count.iter().sum();
    if total_windows == 0 {
        return Err(GenGcBiasModelError::ConfigError(
            "No windows were processed — verify that the BAM and reference share contig names"
                .to_string(),
        ));
    }

    // Only well-populated bins contribute to the reference mean. Including sparse bins
    // risks a single anomalous window (e.g. a repetitive region) inflating the mean and
    // pushing all other weights below 1.0.
    let overall_mean = {
        let (weight_sum, count_sum) = (0..=100usize)
            .filter(|&i| gc_window_count[i] >= config.min_windows_per_bin)
            .fold((0.0f64, 0usize), |(ws, cs), i| {
                (ws + gc_weight_sum[i], cs + gc_window_count[i])
            });
        if count_sum == 0 {
            warn!(
                "No GC bins met min_windows_per_bin ({}); all weights will be neutral (1.0). \
                 Try lowering min_windows_per_bin or using a larger genome/region.",
                config.min_windows_per_bin
            );
            1.0
        } else {
            weight_sum / count_sum as f64
        }
    };

    let weights: Vec<f64> = (0..=100usize)
        .map(|i| {
            if gc_window_count[i] >= config.min_windows_per_bin && overall_mean > 0.0 {
                (gc_weight_sum[i] / gc_window_count[i] as f64) / overall_mean
            } else {
                1.0
            }
        })
        .collect();

    let model = GcBiasModel::from_weights(weights, config.window_size)?;
    model.write_to_file(&config.output_file)?;
    info!("GC bias model written to {:?}", config.output_file);

    Ok(())
}

/// Accumulate GC% vs mean-coverage data for all windows in `[region_start, region_end)`.
///
/// Uses a sliding window for both GC counting and coverage summing — O(region_len)
/// regardless of window_size or window_stride.
fn accumulate_region(
    sequence: &[Nucleotide],
    cov: &[u32],
    region_start: usize,
    region_end: usize,
    window_size: usize,
    window_stride: usize,
    gc_weight_sum: &mut [f64; 101],
    gc_window_count: &mut [usize; 101],
) {
    let first = &sequence[region_start..region_start + window_size];
    let mut gc_count: usize = first
        .iter()
        .filter(|&&n| n == Nucleotide::G || n == Nucleotide::C)
        .count();
    let mut n_count: usize = first
        .iter()
        .filter(|&&n| n == Nucleotide::N || n == Nucleotide::X)
        .count();

    let mut cov_sum: u64 = (0..window_size)
        .map(|i| cov.get(region_start + i).copied().unwrap_or(0) as u64)
        .sum();

    let mut w = region_start;
    loop {
        let called = window_size - n_count;
        if called > 0 {
            let gc_pct = ((gc_count as f64 / called as f64) * 100.0).round() as usize;
            let gc_pct = gc_pct.min(100);
            gc_weight_sum[gc_pct] += cov_sum as f64 / window_size as f64;
            gc_window_count[gc_pct] += 1;
        }

        let next_w = w + window_stride;
        if next_w + window_size > region_end {
            break;
        }

        // Advance both sliding windows by stride steps.
        for j in 0..window_stride {
            match sequence[w + j] {
                Nucleotide::G | Nucleotide::C => gc_count -= 1,
                Nucleotide::N | Nucleotide::X => n_count -= 1,
                _ => {}
            }
            match sequence[w + window_size + j] {
                Nucleotide::G | Nucleotide::C => gc_count += 1,
                Nucleotide::N | Nucleotide::X => n_count += 1,
                _ => {}
            }
            cov_sum -= cov.get(w + j).copied().unwrap_or(0) as u64;
            cov_sum += cov.get(w + window_size + j).copied().unwrap_or(0) as u64;
        }

        w = next_w;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn write_temp(content: &str) -> NamedTempFile {
        let mut f = NamedTempFile::new().unwrap();
        write!(f, "{}", content).unwrap();
        f
    }

    fn write_bam_config(
        reference: &PathBuf,
        bam_file: &PathBuf,
        output_file: &PathBuf,
        window_size: usize,
        window_stride: usize,
        min_windows: usize,
    ) -> NamedTempFile {
        let yaml = format!(
            "reference: {}\nbam_file: {}\noutput_file: {}\noverwrite_output: true\n\
             window_size: {}\nwindow_stride: {}\nmin_windows_per_bin: {}\n",
            reference.display(),
            bam_file.display(),
            output_file.display(),
            window_size,
            window_stride,
            min_windows,
        );
        write_temp(&yaml)
    }

    /// Writes a BGZF BAM at `path` that produces per-base reference coverage equal
    /// to the sum of `depth` over all overlapping `(start_0based, end_0based, depth)`
    /// segments. Each segment is realized by stacking `depth` reads of length
    /// `end - start` at position `start + 1`.
    fn write_coverage_bam(
        path: &std::path::PathBuf,
        contigs: &[(&[u8], usize)],
        // (ref_id, start_0based, end_0based, depth)
        segments: &[(usize, usize, usize, usize)],
    ) {
        use noodles::bam;
        use noodles::core::Position;
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

        let mut builder = sam::Header::builder();
        for &(name, len) in contigs {
            builder = builder.add_reference_sequence(
                name.to_vec(),
                Map::<ReferenceSequence>::new(std::num::NonZero::<usize>::new(len).unwrap()),
            );
        }
        let header = builder.build();
        let file = std::fs::File::create(path).unwrap();
        let mut writer = bam::io::Writer::new(file);
        writer.write_header(&header).unwrap();

        for &(ref_id, start_0, end_0, depth) in segments {
            let read_len = end_0 - start_0;
            if read_len == 0 || depth == 0 {
                continue;
            }
            let cigar: Cigar = [Op::new(Kind::Match, read_len)].into_iter().collect();
            let seq = vec![b'A'; read_len];
            for _ in 0..depth {
                let mut record = RecordBuf::default();
                *record.flags_mut() = Flags::empty();
                *record.cigar_mut() = cigar.clone();
                *record.reference_sequence_id_mut() = Some(ref_id);
                *record.alignment_start_mut() = Position::new(start_0 + 1);
                *record.sequence_mut() = Sequence::from(seq.as_slice());
                *record.mapping_quality_mut() = Some(MappingQuality::try_from(30u8).unwrap());
                writer.write_alignment_record(&header, &record).unwrap();
            }
        }
    }

    #[test]
    fn test_uniform_coverage_produces_neutral_weights() {
        // Reference: 300 bp of ACGT repeating (50% GC), uniform depth 10.
        // All populated bins should have relative weight ~1.0 (no bias).
        let seq: String = "ACGT".repeat(75);
        let fasta = format!(">chr1\n{}\n", seq);
        let ref_file = write_temp(&fasta);

        let temp = tempfile::tempdir().unwrap();
        let bam_path = temp.path().join("uniform.bam");
        write_coverage_bam(&bam_path, &[(b"chr1", 300)], &[(0, 0, 300, 10)]);

        let out = NamedTempFile::new().unwrap();
        let out_path = out.path().to_path_buf();
        drop(out);

        let cfg = write_bam_config(
            &ref_file.path().to_path_buf(),
            &bam_path,
            &out_path,
            100,
            100,
            1,
        );
        runner(&cfg.path().to_path_buf()).unwrap();

        let model = GcBiasModel::from_file(&out_path).unwrap();
        for gc_pct in 0..=100 {
            let w = model.weight_for_gc_fraction(gc_pct as f64 / 100.0);
            assert!(
                (w - 1.0).abs() < 1e-6,
                "Expected weight ~1.0 at GC {}%, got {}",
                gc_pct,
                w
            );
        }
    }

    #[test]
    fn test_gc_correlated_coverage_produces_correct_relative_weights() {
        // Window A: 200 bp AT (0% GC) at depth 10; Window B: 200 bp GC (100% GC)
        // at depth 30. Overall mean coverage = 20, so weight[0%] ~ 0.5 and
        // weight[100%] ~ 1.5.
        let seq = format!("{}{}", "AT".repeat(100), "GC".repeat(100));
        let fasta = format!(">chr1\n{}\n", seq);
        let ref_file = write_temp(&fasta);

        let temp = tempfile::tempdir().unwrap();
        let bam_path = temp.path().join("two_window.bam");
        write_coverage_bam(
            &bam_path,
            &[(b"chr1", 400)],
            &[(0, 0, 200, 10), (0, 200, 400, 30)],
        );

        let out = NamedTempFile::new().unwrap();
        let out_path = out.path().to_path_buf();
        drop(out);

        let cfg = write_bam_config(
            &ref_file.path().to_path_buf(),
            &bam_path,
            &out_path,
            200,
            200,
            1,
        );
        runner(&cfg.path().to_path_buf()).unwrap();

        let model = GcBiasModel::from_file(&out_path).unwrap();
        assert!(
            (model.weight_for_gc_fraction(0.0) - 0.5).abs() < 1e-6,
            "Expected low-GC weight ~0.5, got {}",
            model.weight_for_gc_fraction(0.0)
        );
        assert!(
            (model.weight_for_gc_fraction(1.0) - 1.5).abs() < 1e-6,
            "Expected high-GC weight ~1.5, got {}",
            model.weight_for_gc_fraction(1.0)
        );
    }

    #[test]
    fn test_sparse_bin_gets_neutral_weight() {
        // One 100bp window at GC=50%, min_windows_per_bin=10 → all weights stay 1.0.
        let seq: String = "ACGT".repeat(25);
        let fasta = format!(">chr1\n{}\n", seq);
        let ref_file = write_temp(&fasta);

        let temp = tempfile::tempdir().unwrap();
        let bam_path = temp.path().join("sparse.bam");
        write_coverage_bam(&bam_path, &[(b"chr1", 100)], &[(0, 0, 100, 20)]);

        let out = NamedTempFile::new().unwrap();
        let out_path = out.path().to_path_buf();
        drop(out);

        let cfg = write_bam_config(
            &ref_file.path().to_path_buf(),
            &bam_path,
            &out_path,
            100,
            100,
            10,
        );
        runner(&cfg.path().to_path_buf()).unwrap();

        let model = GcBiasModel::from_file(&out_path).unwrap();
        assert_eq!(model.weight_for_gc_fraction(0.5), 1.0);
    }

    #[test]
    fn test_all_n_contig_is_skipped() {
        // chr1 is all-N (skipped); chr2 has real coverage. Runner must not panic
        // or error on the all-N contig.
        let fasta = ">chr1\nNNNNNNNNNN\n>chr2\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n";
        let ref_file = write_temp(fasta);

        let temp = tempfile::tempdir().unwrap();
        let bam_path = temp.path().join("n_skip.bam");
        write_coverage_bam(
            &bam_path,
            &[(b"chr1", 10), (b"chr2", 104)],
            &[(0, 0, 10, 10), (1, 0, 104, 10)],
        );

        let out = NamedTempFile::new().unwrap();
        let out_path = out.path().to_path_buf();
        drop(out);

        let cfg = write_bam_config(
            &ref_file.path().to_path_buf(),
            &bam_path,
            &out_path,
            100,
            100,
            1,
        );
        runner(&cfg.path().to_path_buf()).unwrap();
    }

    #[test]
    fn test_sparse_bin_excluded_from_overall_mean() {
        // 2 low-GC windows at depth 10 (well-populated, min_windows=2).
        // 1 high-GC window at depth 1000 (sparse — exactly 1, below min_windows=2).
        // If the sparse bin were folded into overall_mean it would inflate the mean
        // and crush low-GC weights well below 1.0. Excluding it as the code does
        // leaves low-GC weight at exactly 1.0.
        let low_gc: String = "AT".repeat(50);
        let high_gc: String = "GC".repeat(50);
        let seq = format!("{}{}{}", low_gc, low_gc, high_gc);
        let fasta = format!(">chr1\n{}\n", seq);
        let ref_file = write_temp(&fasta);

        let temp = tempfile::tempdir().unwrap();
        let bam_path = temp.path().join("sparse_high.bam");
        write_coverage_bam(
            &bam_path,
            &[(b"chr1", 300)],
            &[(0, 0, 200, 10), (0, 200, 300, 1000)],
        );

        let out = NamedTempFile::new().unwrap();
        let out_path = out.path().to_path_buf();
        drop(out);

        let cfg = write_bam_config(
            &ref_file.path().to_path_buf(),
            &bam_path,
            &out_path,
            100,
            100,
            2,
        );
        runner(&cfg.path().to_path_buf()).unwrap();

        let model = GcBiasModel::from_file(&out_path).unwrap();
        let w_low = model.weight_for_gc_fraction(0.0);
        assert!(
            (w_low - 1.0).abs() < 1e-6,
            "Expected low-GC weight ~1.0, got {} (sparse high-GC bin may have inflated mean)",
            w_low
        );
        assert_eq!(
            model.weight_for_gc_fraction(1.0),
            1.0,
            "Expected sparse high-GC bin to have neutral weight 1.0"
        );
    }

    #[test]
    fn test_overlapping_windows_correct_gc_counts() {
        // window_size=4, window_stride=2 → overlapping windows on "AAAAGGGGAAAA".
        // Windows: [0,4)=0%GC, [2,6)=50%GC, [4,8)=100%GC, [6,10)=50%GC, [8,12)=0%GC.
        // Uniform depth → all populated bins should have weight ~1.0.
        let seq = "AAAAGGGGAAAA";
        let fasta = format!(">chr1\n{}\n", seq);
        let ref_file = write_temp(&fasta);

        let temp = tempfile::tempdir().unwrap();
        let bam_path = temp.path().join("overlap.bam");
        write_coverage_bam(&bam_path, &[(b"chr1", 12)], &[(0, 0, 12, 10)]);

        let out = NamedTempFile::new().unwrap();
        let out_path = out.path().to_path_buf();
        drop(out);

        let cfg = write_bam_config(
            &ref_file.path().to_path_buf(),
            &bam_path,
            &out_path,
            4,
            2,
            1,
        );
        runner(&cfg.path().to_path_buf()).unwrap();

        let model = GcBiasModel::from_file(&out_path).unwrap();
        for &gc_frac in &[0.0f64, 0.5, 1.0] {
            let w = model.weight_for_gc_fraction(gc_frac);
            assert!(
                (w - 1.0).abs() < 1e-6,
                "Expected weight ~1.0 at GC {:.0}%, got {}",
                gc_frac * 100.0,
                w
            );
        }
    }
}

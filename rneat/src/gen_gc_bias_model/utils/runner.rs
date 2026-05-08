use std::path::PathBuf;
use log::*;

use common::{
    models::gc_bias_model::GcBiasModel,
    structs::nucleotides::Nucleotide,
};
use crate::gen_gc_bias_model::{
    errors::GenGcBiasModelError,
    utils::{
        config::RunConfiguration,
        coverage_reader::{CoverageData, CoverageIndex},
        fasta_stream::FastaStream,
    },
};

pub fn runner(path: &PathBuf) -> Result<(), GenGcBiasModelError> {
    let config = RunConfiguration::from(path)?;

    info!("Building coverage index from {:?}", config.coverage_file);
    let cov_index = CoverageIndex::build(&config.coverage_file)?;

    let mut gc_weight_sum = [0.0f64; 101];
    let mut gc_window_count = [0usize; 101];

    info!("Processing reference {:?}", config.reference);
    let fasta = FastaStream::open(&config.reference)?;

    for result in fasta {
        let (contig_name, sequence) = result?;
        let contig_len = sequence.len();

        if contig_len < config.window_size {
            debug!(
                "Skipping {} (length {} < window_size {})",
                contig_name, contig_len, config.window_size
            );
            continue;
        }

        let cov = match CoverageData::load(
            &config.coverage_file,
            &cov_index,
            &contig_name,
            config.coverage_format,
        )? {
            Some(c) => c,
            None => {
                debug!("No coverage data for {}, skipping", contig_name);
                continue;
            }
        };

        let regions: Vec<(usize, usize)> = if let Some(bed_regions) =
            config.bed_table.get(&contig_name)
        {
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
            "No windows were processed — verify that coverage file and reference share contig names".to_string(),
        ));
    }

    let overall_mean = {
        let weight_sum: f64 = (0..=100usize)
            .filter(|&i| gc_window_count[i] > 0)
            .map(|i| gc_weight_sum[i])
            .sum();
        let count_sum: usize = (0..=100usize)
            .filter(|&i| gc_window_count[i] > 0)
            .map(|i| gc_window_count[i])
            .sum();
        weight_sum / count_sum as f64
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
    cov: &CoverageData,
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
    let mut n_count: usize = first.iter().filter(|&&n| n == Nucleotide::N).count();

    let mut cov_sum: u64 = (0..window_size)
        .map(|i| cov.depth_at(region_start + i) as u64)
        .sum();

    let mut w = region_start;
    loop {
        let called = window_size - n_count;
        if called > 0 {
            let gc_pct = ((gc_count as f64 / called as f64) * 100.0)
                .round() as usize;
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
                Nucleotide::N => n_count -= 1,
                _ => {}
            }
            match sequence[w + window_size + j] {
                Nucleotide::G | Nucleotide::C => gc_count += 1,
                Nucleotide::N => n_count += 1,
                _ => {}
            }
            cov_sum -= cov.depth_at(w + j) as u64;
            cov_sum += cov.depth_at(w + window_size + j) as u64;
        }

        w = next_w;
    }
}

/// Returns the contiguous non-N regions of a sequence as `(start, end)` pairs.
fn non_n_regions(sequence: &[Nucleotide]) -> Vec<(usize, usize)> {
    let mut regions = Vec::new();
    let mut region_start: Option<usize> = None;

    for (i, &nuc) in sequence.iter().enumerate() {
        match (nuc, region_start) {
            (Nucleotide::N, Some(s)) => {
                regions.push((s, i));
                region_start = None;
            }
            (Nucleotide::N, None) | (Nucleotide::X, _) => {}
            (_, None) => region_start = Some(i),
            (_, Some(_)) => {}
        }
    }
    if let Some(s) = region_start {
        regions.push((s, sequence.len()));
    }
    regions
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn make_sequence(bases: &str) -> Vec<Nucleotide> {
        bases.chars().map(Nucleotide::from).collect()
    }

    // Tests for non_n_regions
    #[test]
    fn test_non_n_regions_no_ns() {
        let seq = make_sequence("ACGTACGT");
        let regions = non_n_regions(&seq);
        assert_eq!(regions, vec![(0, 8)]);
    }

    #[test]
    fn test_non_n_regions_leading_n() {
        let seq = make_sequence("NNACGT");
        let regions = non_n_regions(&seq);
        assert_eq!(regions, vec![(2, 6)]);
    }

    #[test]
    fn test_non_n_regions_interior_n() {
        let seq = make_sequence("ACGTNACGT");
        let regions = non_n_regions(&seq);
        assert_eq!(regions, vec![(0, 4), (5, 9)]);
    }

    #[test]
    fn test_non_n_regions_all_n() {
        let seq = make_sequence("NNNN");
        let regions = non_n_regions(&seq);
        assert!(regions.is_empty());
    }

    // Tests for accumulate_region — tested via runner integration below.
    // The end-to-end runner tests cover accumulation correctness.

    fn write_temp(content: &str) -> NamedTempFile {
        let mut f = NamedTempFile::new().unwrap();
        write!(f, "{}", content).unwrap();
        f
    }

    fn write_config(
        reference: &PathBuf,
        coverage_file: &PathBuf,
        output_file: &PathBuf,
        window_size: usize,
        window_stride: usize,
        min_windows: usize,
        overwrite: bool,
    ) -> NamedTempFile {
        let yaml = format!(
            "\
reference: {}
coverage_file: {}
coverage_format: samtools-depth
bed_file: .
output_file: {}
overwrite_output: {}
window_size: {}
window_stride: {}
min_windows_per_bin: {}
",
            reference.display(),
            coverage_file.display(),
            output_file.display(),
            overwrite,
            window_size,
            window_stride,
            min_windows,
        );
        write_temp(&yaml)
    }

    #[test]
    fn test_uniform_coverage_produces_neutral_weights() {
        // Reference: 300 bp of ACGT repeating (25% GC), uniform coverage
        let seq: String = "ACGT".repeat(75); // 300 bp, GC = 50%
        let fasta_content = format!(">chr1\n{}\n", seq);
        let ref_file = write_temp(&fasta_content);

        // samtools depth: 1-based positions, depth 10 everywhere
        let cov_lines: String = (1..=300).map(|i| format!("chr1\t{}\t10\n", i)).collect();
        let cov_file = write_temp(&cov_lines);

        let out = NamedTempFile::new().unwrap();
        let out_path = out.path().to_path_buf();
        drop(out);

        let config_file = write_config(
            &ref_file.path().to_path_buf(),
            &cov_file.path().to_path_buf(),
            &out_path,
            100,  // window_size
            100,  // window_stride (non-overlapping)
            1,    // min_windows_per_bin
            true,
        );

        runner(&config_file.path().to_path_buf()).unwrap();

        let model = GcBiasModel::from_file(&out_path).unwrap();
        // All populated bins should be ~1.0 (uniform coverage → no bias)
        for gc_pct in 0..=100 {
            let w = model.weight_for_gc_fraction(gc_pct as f64 / 100.0);
            // Only the GC% value present in the sequence (50%) is populated
            // All others fall back to 1.0. Either way the result should be ~1.0.
            assert!(
                (w - 1.0).abs() < 1e-6,
                "Expected weight ~1.0 at GC {}%, got {}",
                gc_pct, w
            );
        }
    }

    #[test]
    fn test_gc_correlated_coverage_produces_correct_relative_weights() {
        // Two windows: low-GC gets depth 10, high-GC gets depth 30.
        // Expect weight[low_gc] ~ 0.5 and weight[high_gc] ~ 1.5
        // (relative to overall mean of 20).
        //
        // Window A: 200 bp of AT (GC = 0%)  → depth 10
        // Window B: 200 bp of GC (GC = 100%) → depth 30
        let seq = format!("{}{}", "AT".repeat(100), "GC".repeat(100));
        let fasta_content = format!(">chr1\n{}\n", seq);
        let ref_file = write_temp(&fasta_content);

        let cov_lines: String = (1..=400)
            .map(|i| {
                let depth = if i <= 200 { 10 } else { 30 };
                format!("chr1\t{}\t{}\n", i, depth)
            })
            .collect();
        let cov_file = write_temp(&cov_lines);

        let out = NamedTempFile::new().unwrap();
        let out_path = out.path().to_path_buf();
        drop(out);

        let config_file = write_config(
            &ref_file.path().to_path_buf(),
            &cov_file.path().to_path_buf(),
            &out_path,
            200,  // window_size = full window
            200,  // window_stride
            1,
            true,
        );

        runner(&config_file.path().to_path_buf()).unwrap();

        let model = GcBiasModel::from_file(&out_path).unwrap();
        let w_low_gc = model.weight_for_gc_fraction(0.0);
        let w_high_gc = model.weight_for_gc_fraction(1.0);

        assert!(
            (w_low_gc - 0.5).abs() < 1e-6,
            "Expected low-GC weight ~0.5, got {}",
            w_low_gc
        );
        assert!(
            (w_high_gc - 1.5).abs() < 1e-6,
            "Expected high-GC weight ~1.5, got {}",
            w_high_gc
        );
    }

    #[test]
    fn test_sparse_bin_gets_neutral_weight() {
        // One window at GC=50%; min_windows_per_bin = 10 → weight should be 1.0
        let seq: String = "ACGT".repeat(25); // 100 bp, GC = 50%
        let fasta_content = format!(">chr1\n{}\n", seq);
        let ref_file = write_temp(&fasta_content);

        let cov_lines: String = (1..=100).map(|i| format!("chr1\t{}\t20\n", i)).collect();
        let cov_file = write_temp(&cov_lines);

        let out = NamedTempFile::new().unwrap();
        let out_path = out.path().to_path_buf();
        drop(out);

        let config_file = write_config(
            &ref_file.path().to_path_buf(),
            &cov_file.path().to_path_buf(),
            &out_path,
            100,
            100,
            10, // min_windows_per_bin = 10 (only 1 window available)
            true,
        );

        runner(&config_file.path().to_path_buf()).unwrap();

        let model = GcBiasModel::from_file(&out_path).unwrap();
        // Only 1 window < min_windows_per_bin, so all weights fall back to 1.0
        assert_eq!(model.weight_for_gc_fraction(0.5), 1.0);
    }

    #[test]
    fn test_all_n_contig_is_skipped() {
        let fasta_content = ">chr1\nNNNNNNNNNN\n";
        let ref_file = write_temp(fasta_content);
        let cov_lines: String = (1..=10).map(|i| format!("chr1\t{}\t10\n", i)).collect();
        let cov_file = write_temp(&cov_lines);

        let out = NamedTempFile::new().unwrap();
        let out_path = out.path().to_path_buf();
        drop(out);

        // Also add a second valid contig so runner doesn't error on zero windows
        let fasta_content2 = ">chr1\nNNNNNNNNNN\n>chr2\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n";
        let ref_file2 = write_temp(fasta_content2);

        let cov2_lines: String = (1..=10).map(|i| format!("chr1\t{}\t10\n", i))
            .chain((1..=104).map(|i| format!("chr2\t{}\t10\n", i)))
            .collect();
        let cov_file2 = write_temp(&cov2_lines);

        let config_file = write_config(
            &ref_file2.path().to_path_buf(),
            &cov_file2.path().to_path_buf(),
            &out_path,
            100,
            100,
            1,
            true,
        );

        // Should not panic or error — all-N window is simply skipped
        runner(&config_file.path().to_path_buf()).unwrap();
    }
}

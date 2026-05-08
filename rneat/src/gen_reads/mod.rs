pub mod errors;
pub mod utils;
use log::*;
use errors::GenerateReadsError;

#[cfg(test)]
mod tests {
    use std::io::{BufRead, BufReader, Write};
    use std::path::{Path, PathBuf};
    use flate2::read::MultiGzDecoder;
    use noodles::bam;
    use tempfile::{tempdir, NamedTempFile};
    use common::rng::NeatRng;
    use crate::gen_reads::utils::{config::RunConfiguration, runner::run_neat};

    // ── helpers ──────────────────────────────────────────────────────────────

    fn make_config(reference: &Path, output_dir: &Path, paired_ended: bool) -> RunConfiguration {
        let mut yaml = NamedTempFile::new().unwrap();
        let pair_section = if paired_ended {
            "paired_ended: true\nfragment_mean: 200.0\nfragment_st_dev: 30.0\n"
        } else {
            "paired_ended: false\n"
        };
        write!(yaml,
            "reference: {ref}\nread_len: 50\ncoverage: 2\n\
             produce_fastq: true\nproduce_bam: true\nproduce_vcf: false\n\
             output_dir: {out}\noutput_filename: bam_int_test\noverwrite_output: true\n\
             rng_seed: integration bam test\n{pair}",
            ref = reference.display(),
            out = output_dir.display(),
            pair = pair_section,
        ).unwrap();
        RunConfiguration::from_yaml_file(&yaml.path().to_path_buf()).unwrap()
    }

    fn run(config: RunConfiguration) {
        let mut rng = NeatRng::new_from_seed(&config.seed_vec).unwrap();
        run_neat(&Box::new(config), &mut rng).unwrap();
    }

    /// Returns (ref_id, 0-based alignment start) for every record in the BAM.
    fn bam_positions(path: &Path) -> Vec<(usize, usize)> {
        let mut reader = bam::io::Reader::new(std::fs::File::open(path).unwrap());
        reader.read_header().unwrap();
        reader.records()
            .map(|r| {
                let rec = r.unwrap();
                let ref_id = rec.reference_sequence_id()
                    .unwrap()  // Option → Some(Result<usize>)
                    .unwrap(); // Result  → usize
                let pos = rec.alignment_start()
                    .unwrap()  // Option → Some(Result<Position>)
                    .unwrap()  // Result → Position
                    .get() - 1;
                (ref_id, pos)
            })
            .collect()
    }

    fn fastq_record_count(path: &Path) -> usize {
        // The output FASTQ is multiple concatenated gzip streams (one per contig);
        // MultiGzDecoder reads all streams in sequence.
        let reader = BufReader::new(MultiGzDecoder::new(std::fs::File::open(path).unwrap()));
        reader.lines().count() / 4
    }

    fn h1n1_reference() -> PathBuf {
        PathBuf::from(format!("{}/test_data/references/H1N1.fa", env!("CARGO_MANIFEST_DIR")))
    }

    // ── tests ─────────────────────────────────────────────────────────────────

    #[test]
    fn test_bam_single_ended_produced_and_sorted() {
        let out = tempdir().unwrap();
        let config = make_config(&h1n1_reference(), out.path(), false);
        let bam_path = config.output_bam.clone().unwrap();
        run(config);

        let positions = bam_positions(&bam_path);
        assert!(!positions.is_empty(), "BAM should contain records");

        for window in positions.windows(2) {
            assert!(
                window[0] <= window[1],
                "BAM not coordinate-sorted: {:?} followed by {:?}",
                window[0], window[1]
            );
        }
    }

    #[test]
    fn test_bam_paired_ended_produced_and_sorted() {
        let out = tempdir().unwrap();
        let config = make_config(&h1n1_reference(), out.path(), true);
        let bam_path = config.output_bam.clone().unwrap();
        run(config);

        let positions = bam_positions(&bam_path);
        assert!(!positions.is_empty(), "BAM should contain records");

        for window in positions.windows(2) {
            assert!(
                window[0] <= window[1],
                "BAM not coordinate-sorted: {:?} followed by {:?}",
                window[0], window[1]
            );
        }
    }

    #[test]
    fn test_bam_record_count_matches_fastq_single_ended() {
        let out = tempdir().unwrap();
        let config = make_config(&h1n1_reference(), out.path(), false);
        let bam_path = config.output_bam.clone().unwrap();
        let fastq_path = config.output_fastq_1.clone().unwrap();
        run(config);

        let bam_count = bam_positions(&bam_path).len();
        let fastq_count = fastq_record_count(&fastq_path);
        assert_eq!(
            bam_count, fastq_count,
            "BAM record count ({bam_count}) must equal FASTQ record count ({fastq_count})"
        );
    }

    /// Writes a BED covering only H1N1_HA (the full segment) and verifies that:
    /// 1. Reads are generated (the covered contig produces output).
    /// 2. Fewer reads are produced than the no-BED run (7 of 8 segments are excluded).
    #[test]
    fn test_target_bed_limits_reads() {
        let out = tempdir().unwrap();

        // BED covering only H1N1_HA; the other 7 segments are absent → skipped entirely.
        let bed_path = out.path().join("targets.bed");
        std::fs::write(&bed_path, "H1N1_HA\t0\t1726\n").unwrap();

        let mut yaml = NamedTempFile::new().unwrap();
        write!(yaml,
            "reference: {ref}\nread_len: 50\ncoverage: 2\n\
             produce_fastq: true\nproduce_bam: false\nproduce_vcf: false\n\
             output_dir: {out}\noutput_filename: bed_test\noverwrite_output: true\n\
             rng_seed: bed test\ntarget_bed: {bed}\n",
            ref = h1n1_reference().display(),
            out = out.path().display(),
            bed = bed_path.display(),
        ).unwrap();
        let config = RunConfiguration::from_yaml_file(&yaml.path().to_path_buf()).unwrap();
        let fastq_path = config.output_fastq_1.clone().unwrap();
        run(config);

        let bed_count = fastq_record_count(&fastq_path);

        // Without a BED, all 8 segments produce reads.
        let out_full = tempdir().unwrap();
        let config_full = make_config(&h1n1_reference(), out_full.path(), false);
        let fastq_full = config_full.output_fastq_1.clone().unwrap();
        run(config_full);
        let full_count = fastq_record_count(&fastq_full);

        assert!(bed_count > 0, "BED-filtered run produced no reads");
        assert!(
            bed_count < full_count,
            "BED-filtered run ({bed_count} reads) should be less than full run ({full_count} reads)"
        );
    }

    /// A BED that names a contig not in the reference should produce an empty FASTQ.
    #[test]
    fn test_target_bed_unknown_contig_produces_no_reads() {
        let out = tempdir().unwrap();

        let bed_path = out.path().join("empty_targets.bed");
        std::fs::write(&bed_path, "nonexistent_contig\t0\t1000\n").unwrap();

        let mut yaml = NamedTempFile::new().unwrap();
        write!(yaml,
            "reference: {ref}\nread_len: 50\ncoverage: 2\n\
             produce_fastq: true\nproduce_bam: false\nproduce_vcf: false\n\
             output_dir: {out}\noutput_filename: empty_bed_test\noverwrite_output: true\n\
             rng_seed: empty bed test\ntarget_bed: {bed}\n",
            ref = h1n1_reference().display(),
            out = out.path().display(),
            bed = bed_path.display(),
        ).unwrap();
        let config = RunConfiguration::from_yaml_file(&yaml.path().to_path_buf()).unwrap();
        let fastq_path = config.output_fastq_1.clone().unwrap();
        run(config);

        let count = fastq_record_count(&fastq_path);
        assert_eq!(count, 0, "Expected no reads when BED references only unknown contigs, got {count}");
    }

    /// Supply a VCF with a single known SNP on H1N1_HA and verify it appears in the
    /// output VCF. The SNP is at 1-based position 100 (A→T); H1N1_HA position 99
    /// (0-based) is confirmed non-N so the variant will be applied.
    #[test]
    fn test_input_vcf_snp_appears_in_output() {
        let out = tempdir().unwrap();

        // Minimal single-sample VCF: position 100 (1-based), A→T, homozygous
        let vcf_path = out.path().join("input.vcf");
        std::fs::write(&vcf_path,
            "##fileformat=VCFv4.2\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
             H1N1_HA\t100\t.\tA\tT\t.\tPASS\t.\tGT\t1/1\n"
        ).unwrap();

        let mut yaml = NamedTempFile::new().unwrap();
        write!(yaml,
            "reference: {ref}\nread_len: 50\ncoverage: 2\n\
             produce_fastq: true\nproduce_bam: false\nproduce_vcf: true\n\
             output_dir: {out}\noutput_filename: vcf_input_test\noverwrite_output: true\n\
             rng_seed: vcf input test\ninput_vcf: {vcf}\n",
            ref = h1n1_reference().display(),
            out = out.path().display(),
            vcf = vcf_path.display(),
        ).unwrap();
        let config = RunConfiguration::from_yaml_file(&yaml.path().to_path_buf()).unwrap();
        let output_vcf = config.output_vcf.clone().unwrap();
        run(config);

        // Read the output VCF and look for a record at position 100 on H1N1_HA
        use std::io::{BufRead, BufReader};
        use flate2::read::GzDecoder;
        let reader = BufReader::new(GzDecoder::new(std::fs::File::open(&output_vcf).unwrap()));
        let found = reader.lines()
            .filter_map(|l| l.ok())
            .filter(|l| !l.starts_with('#'))
            .any(|l| {
                let cols: Vec<&str> = l.split('\t').collect();
                cols.len() >= 5 && cols[0] == "H1N1_HA" && cols[1] == "100"
            });
        assert!(found, "Expected SNP at H1N1_HA:100 in output VCF but did not find it");
    }
}

use std::path::PathBuf;
use crate::{
    gen_reads::utils::{
        config::RunConfiguration, 
        runner::run_neat
    }
};
use common::rng::NeatRng;

/// gen-reads is the primary read generation function of rneat. It reads a fasta file and generates a set of fastqs and/or a set of variants. It can now also filter reads by bed file.
pub fn main(config: &PathBuf) -> Result<(), GenerateReadsError> {   
    info!("////////////// Welcome to rusty-neat read generator! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\");
    // set up the config struct based on whether there was an input config. Input config
    // overrides any other inputs.
    let config = if &config.display().to_string() != "" {
        info!("Using Configuration file input: {:?}", &config);
        RunConfiguration::from_yaml_file(config)
            .expect("Error generating run configuration from input yaml")
    } else {
        panic!("Failed to supply config file!");
    };

    // Check that we are clear to write outputs
    if !config.overwrite_output {
        if let Some(filename) = &config.output_fastq_1 {
            if filename.is_file() {
                panic!("Attempting to overwrite an existing file: {:?}", filename);
            }
        }
        if let Some(filename) = &config.output_fastq_2 {
            if filename.is_file() {
                panic!("Attempting to overwrite an existing file: {:?}", filename);
            }
        }
        if let Some(filename) = &config.output_vcf {
            if filename.is_file() {
                panic!("Attempting to overwrite an existing file: {:?}", filename);
            }
        }
    }
    
    info!("////////////// Configuration successuful! Ready to run! \\\\\\\\\\\\\\\\\\\\\\\\\\");
    // Generate the RNG used for this run. If one was given in the config file, use that, or else
    // use thread_rng to generate a random seed, then seed using a SeedableRng based on StdRng
    let mut rng: NeatRng = NeatRng::new_from_seed(&config.seed_vec)
        .expect("Neat failed during rng creation!");
    // run the generate reads main script
    let result = run_neat(&Box::new(config.clone()), &mut rng);
    match result {
        Ok(_) => {
            // Continue on for bed filtering
            Ok(())
        },
        Err(error) => {
            error!("runner returned an error {:?}", error);
            Err(GenerateReadsError::RunnerError)
        },
    }
}
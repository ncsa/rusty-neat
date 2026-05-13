// This is the run configuration for this particular run, which holds the parameters needed by the
// various side functions. It is build with a ConfigurationBuilder, which can take either a
// config yaml file or command line arguments and turn them into the configuration.
use serde_yml::Value;
use std::collections::HashMap;
use std::path::PathBuf;
use std::fs;
use common::{
    file_tools::bed_reader::read_bed,
    structs::bed_record::BedRecord,
};
use crate::gen_gc_bias_model::errors::GenGcBiasModelError;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CoverageFormat {
    SamtoolsDepth,
    BedtoolsGenomecovD,
    BedtoolsGenomecovDz,
}

impl CoverageFormat {
    pub fn from_str(input: &str) -> Result<Self, GenGcBiasModelError> {
        match input {
            "samtools-depth" => Ok(Self::SamtoolsDepth),
            "bedtools-genomecov-d" => Ok(Self::BedtoolsGenomecovD),
            "bedtools-genomecov-dz" => Ok(Self::BedtoolsGenomecovDz),
            other => Err(GenGcBiasModelError::ConfigError(
                format!("Unsupported coverage_format: {}", other)
            )),
        }
    }

    pub fn is_zero_based(&self) -> bool {
        matches!(self, Self::BedtoolsGenomecovDz)
    }
}

#[derive(Debug, Clone)]
pub struct RunConfiguration {
    // This struct holds all the parameters for running GC Bias modeling.
    // It is built from user-supplied input in the form of a configuration yaml file
    pub reference: PathBuf,
    pub coverage_file: PathBuf,
    pub coverage_format: CoverageFormat,
    pub bed_table: HashMap<String, Vec<BedRecord>>,
    pub output_file: PathBuf,
    pub overwrite_output: bool,
    pub window_size: usize,
    pub window_stride: usize,
    pub min_windows_per_bin: usize
}

impl RunConfiguration {
    pub fn from(yml_file: &PathBuf) -> Result<Self, GenGcBiasModelError> {
        let file = fs::File::open(yml_file)?;
        let scrape_config: HashMap<String, Value> = serde_yml::from_reader(file)?;

        let reference = PathBuf::from(
            scrape_config["reference"]
                .as_str()
                .ok_or_else(|| GenGcBiasModelError::ConfigError("Missing reference".to_string()))?
        );
        if !reference.is_file() {
            return Err(GenGcBiasModelError::ConfigError(
                format!("Invalid reference file {:?}", reference)
            ));
        }

        let coverage_file = PathBuf::from(
            scrape_config["coverage_file"]
                .as_str()
                .ok_or_else(|| GenGcBiasModelError::ConfigError("Missing coverage_file".to_string()))?
        );
        if !coverage_file.is_file() {
            return Err(GenGcBiasModelError::ConfigError(
                format!("Invalid coverage file {:?}", coverage_file)
            ));
        }

        let coverage_format = CoverageFormat::from_str(
            scrape_config
                .get("coverage_format")
                .and_then(|v| v.as_str())
                .unwrap_or("samtools-depth")
        )?;

        let bed_file_raw = scrape_config
            .get("bed_file")
            .and_then(|v| v.as_str())
            .unwrap_or(".");

        let bed_table = if bed_file_raw == "." {
            HashMap::new()
        } else {
            let bed_file = PathBuf::from(bed_file_raw);
            if !bed_file.is_file() {
                return Err(GenGcBiasModelError::ConfigError(
                    format!("Invalid BED file {:?}", bed_file)
                ));
            }
            read_bed(&bed_file, false)?
        };

        let overwrite_output = scrape_config
            .get("overwrite_output")
            .and_then(|v| v.as_bool())
            .unwrap_or(false);

        let output_file = PathBuf::from(
            scrape_config["output_file"]
                .as_str()
                .ok_or_else(|| GenGcBiasModelError::ConfigError("Missing output_file".to_string()))?
        );

        if !overwrite_output && output_file.is_file() {
            return Err(GenGcBiasModelError::ConfigError(
                format!("Attempting to overwrite existing file {:?}", output_file)
            ));
        }

        let window_size = scrape_config
            .get("window_size")
            .and_then(|v| v.as_u64())
            .unwrap_or(100) as usize;

        if window_size == 0 {
            return Err(GenGcBiasModelError::ConfigError(
                "window_size must be > 0".to_string()
            ));
        }

        let window_stride = scrape_config
            .get("window_stride")
            .and_then(|v| v.as_u64())
            .unwrap_or(window_size as u64) as usize;

        if window_stride == 0 {
            return Err(GenGcBiasModelError::ConfigError(
                "window_stride must be > 0".to_string()
            ));
        }

        if window_stride > window_size {
            return Err(GenGcBiasModelError::ConfigError(format!(
                "window_stride ({}) must not exceed window_size ({}); \
                 for non-overlapping windows use window_stride = window_size",
                window_stride, window_size
            )));
        }

        let min_windows_per_bin = scrape_config
            .get("min_windows_per_bin")
            .and_then(|v| v.as_u64())
            .unwrap_or(10) as usize;

        Ok(RunConfiguration {
            reference,
            coverage_file,
            coverage_format,
            bed_table,
            output_file,
            overwrite_output,
            window_size,
            window_stride,
            min_windows_per_bin,
        })
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn write_config(yaml: &str) -> tempfile::NamedTempFile {
        let mut tmp = tempfile::NamedTempFile::new().unwrap();
        write!(tmp, "{}", yaml).unwrap();
        tmp
    }

    fn write_temp_file(contents: &str) -> tempfile::NamedTempFile {
        let mut tmp = tempfile::NamedTempFile::new().unwrap();
        write!(tmp, "{}", contents).unwrap();
        tmp
    }

    fn base_yaml(
        reference: &PathBuf,
        coverage_file: &PathBuf,
        output_file: &PathBuf,
    ) -> String {
        format!(
            "\
reference: {}
coverage_file: {}
coverage_format: samtools-depth
bed_file: .
output_file: {}
overwrite_output: true
window_size: 100
window_stride: 100
min_windows_per_bin: 10
",
            reference.display(),
            coverage_file.display(),
            output_file.display(),
        )
    }

    #[test]
    fn test_coverage_format_from_str_samtools_depth() {
        let format = CoverageFormat::from_str("samtools-depth").unwrap();
        assert_eq!(format, CoverageFormat::SamtoolsDepth);
        assert!(!format.is_zero_based());
    }

    #[test]
    fn test_coverage_format_from_str_bedtools_genomecov_d() {
        let format = CoverageFormat::from_str("bedtools-genomecov-d").unwrap();
        assert_eq!(format, CoverageFormat::BedtoolsGenomecovD);
        assert!(!format.is_zero_based());
    }

    #[test]
    fn test_coverage_format_from_str_bedtools_genomecov_dz() {
        let format = CoverageFormat::from_str("bedtools-genomecov-dz").unwrap();
        assert_eq!(format, CoverageFormat::BedtoolsGenomecovDz);
        assert!(format.is_zero_based());
    }

    #[test]
    fn test_coverage_format_from_str_rejects_unknown_format() {
        let result = CoverageFormat::from_str("not-a-real-format");

        assert!(
            result.is_err(),
            "Expected unknown coverage format to error"
        );
    }

    #[test]
    fn test_run_configuration_parses_required_fields() {
        let reference = write_temp_file(">chr1\nACGTACGT\n");
        let coverage = write_temp_file("chr1\t1\t10\nchr1\t2\t11\n");
        let output = tempfile::NamedTempFile::new().unwrap();
        let output_path = output.path().to_path_buf();
        drop(output);

        let yaml = base_yaml(
            &reference.path().to_path_buf(),
            &coverage.path().to_path_buf(),
            &output_path,
        );
        let config_file = write_config(&yaml);

        let config = RunConfiguration::from(&config_file.path().to_path_buf()).unwrap();

        assert_eq!(config.reference, reference.path());
        assert_eq!(config.coverage_file, coverage.path());
        assert_eq!(config.coverage_format, CoverageFormat::SamtoolsDepth);
        assert!(config.bed_table.is_empty());
        assert_eq!(config.output_file, output_path);
        assert!(config.overwrite_output);
        assert_eq!(config.window_size, 100);
        assert_eq!(config.window_stride, 100);
        assert_eq!(config.min_windows_per_bin, 10);
    }

    #[test]
    fn test_run_configuration_defaults_optional_fields() {
        let reference = write_temp_file(">chr1\nACGTACGT\n");
        let coverage = write_temp_file("chr1\t1\t10\nchr1\t2\t11\n");
        let output = tempfile::NamedTempFile::new().unwrap();
        let output_path = output.path().to_path_buf();
        drop(output);

        let yaml = format!(
            "\
reference: {}
coverage_file: {}
output_file: {}
overwrite_output: true
",
            reference.path().display(),
            coverage.path().display(),
            output_path.display(),
        );
        let config_file = write_config(&yaml);

        let config = RunConfiguration::from(&config_file.path().to_path_buf()).unwrap();

        assert_eq!(config.coverage_format, CoverageFormat::SamtoolsDepth);
        assert!(config.bed_table.is_empty());
        assert_eq!(config.window_size, 100);
        assert_eq!(config.window_stride, 100);
        assert_eq!(config.min_windows_per_bin, 10);
    }

    #[test]
    fn test_run_configuration_parses_bedtools_genomecov_d_format() {
        let reference = write_temp_file(">chr1\nACGTACGT\n");
        let coverage = write_temp_file("chr1\t1\t10\nchr1\t2\t11\n");
        let output = tempfile::NamedTempFile::new().unwrap();
        let output_path = output.path().to_path_buf();
        drop(output);

        let yaml = format!(
            "\
reference: {}
coverage_file: {}
coverage_format: bedtools-genomecov-d
bed_file: .
output_file: {}
overwrite_output: true
window_size: 50
window_stride: 25
min_windows_per_bin: 5
",
            reference.path().display(),
            coverage.path().display(),
            output_path.display(),
        );
        let config_file = write_config(&yaml);

        let config = RunConfiguration::from(&config_file.path().to_path_buf()).unwrap();

        assert_eq!(config.coverage_format, CoverageFormat::BedtoolsGenomecovD);
        assert!(!config.coverage_format.is_zero_based());
        assert_eq!(config.window_size, 50);
        assert_eq!(config.window_stride, 25);
        assert_eq!(config.min_windows_per_bin, 5);
    }

    #[test]
    fn test_run_configuration_parses_bedtools_genomecov_dz_format() {
        let reference = write_temp_file(">chr1\nACGTACGT\n");
        let coverage = write_temp_file("chr1\t0\t10\nchr1\t1\t11\n");
        let output = tempfile::NamedTempFile::new().unwrap();
        let output_path = output.path().to_path_buf();
        drop(output);

        let yaml = format!(
            "\
reference: {}
coverage_file: {}
coverage_format: bedtools-genomecov-dz
bed_file: .
output_file: {}
overwrite_output: true
",
            reference.path().display(),
            coverage.path().display(),
            output_path.display(),
        );
        let config_file = write_config(&yaml);

        let config = RunConfiguration::from(&config_file.path().to_path_buf()).unwrap();

        assert_eq!(config.coverage_format, CoverageFormat::BedtoolsGenomecovDz);
        assert!(config.coverage_format.is_zero_based());
    }

    #[test]
    fn test_run_configuration_parses_bed_file() {
        let reference = write_temp_file(">chr1\nACGTACGT\n");
        let coverage = write_temp_file("chr1\t1\t10\nchr1\t2\t11\n");
        let bed = write_temp_file("chr1\t0\t4\nchr1\t4\t8\n");
        let output = tempfile::NamedTempFile::new().unwrap();
        let output_path = output.path().to_path_buf();
        drop(output);

        let yaml = format!(
            "\
reference: {}
coverage_file: {}
coverage_format: samtools-depth
bed_file: {}
output_file: {}
overwrite_output: true
",
            reference.path().display(),
            coverage.path().display(),
            bed.path().display(),
            output_path.display(),
        );
        let config_file = write_config(&yaml);

        let config = RunConfiguration::from(&config_file.path().to_path_buf()).unwrap();

        let chr1_records = config
            .bed_table
            .get("chr1")
            .expect("Expected chr1 BED records");

        assert_eq!(chr1_records.len(), 2);
    }

    #[test]
    fn test_run_configuration_rejects_missing_reference_file() {
        let coverage = write_temp_file("chr1\t1\t10\n");
        let output = tempfile::NamedTempFile::new().unwrap();
        let output_path = output.path().to_path_buf();
        drop(output);

        let missing_reference = PathBuf::from("/tmp/rneat_missing_reference_for_gc_bias.fa");

        let yaml = format!(
            "\
reference: {}
coverage_file: {}
coverage_format: samtools-depth
bed_file: .
output_file: {}
overwrite_output: true
",
            missing_reference.display(),
            coverage.path().display(),
            output_path.display(),
        );
        let config_file = write_config(&yaml);

        let result = RunConfiguration::from(&config_file.path().to_path_buf());

        assert!(
            result.is_err(),
            "Expected missing reference file to return an error"
        );
    }

    #[test]
    fn test_run_configuration_rejects_missing_coverage_file() {
        let reference = write_temp_file(">chr1\nACGTACGT\n");
        let output = tempfile::NamedTempFile::new().unwrap();
        let output_path = output.path().to_path_buf();
        drop(output);

        let missing_coverage = PathBuf::from("/tmp/rneat_missing_coverage_for_gc_bias.depth");

        let yaml = format!(
            "\
reference: {}
coverage_file: {}
coverage_format: samtools-depth
bed_file: .
output_file: {}
overwrite_output: true
",
            reference.path().display(),
            missing_coverage.display(),
            output_path.display(),
        );
        let config_file = write_config(&yaml);

        let result = RunConfiguration::from(&config_file.path().to_path_buf());

        assert!(
            result.is_err(),
            "Expected missing coverage file to return an error"
        );
    }

    #[test]
    fn test_run_configuration_rejects_missing_bed_file() {
        let reference = write_temp_file(">chr1\nACGTACGT\n");
        let coverage = write_temp_file("chr1\t1\t10\n");
        let output = tempfile::NamedTempFile::new().unwrap();
        let output_path = output.path().to_path_buf();
        drop(output);

        let missing_bed = PathBuf::from("/tmp/rneat_missing_bed_for_gc_bias.bed");

        let yaml = format!(
            "\
reference: {}
coverage_file: {}
coverage_format: samtools-depth
bed_file: {}
output_file: {}
overwrite_output: true
",
            reference.path().display(),
            coverage.path().display(),
            missing_bed.display(),
            output_path.display(),
        );
        let config_file = write_config(&yaml);

        let result = RunConfiguration::from(&config_file.path().to_path_buf());

        assert!(
            result.is_err(),
            "Expected missing BED file to return an error"
        );
    }

    #[test]
    fn test_run_configuration_rejects_unknown_coverage_format() {
        let reference = write_temp_file(">chr1\nACGTACGT\n");
        let coverage = write_temp_file("chr1\t1\t10\n");
        let output = tempfile::NamedTempFile::new().unwrap();
        let output_path = output.path().to_path_buf();
        drop(output);

        let yaml = format!(
            "\
reference: {}
coverage_file: {}
coverage_format: definitely-not-real
bed_file: .
output_file: {}
overwrite_output: true
",
            reference.path().display(),
            coverage.path().display(),
            output_path.display(),
        );
        let config_file = write_config(&yaml);

        let result = RunConfiguration::from(&config_file.path().to_path_buf());

        assert!(
            result.is_err(),
            "Expected unknown coverage format to return an error"
        );
    }

    #[test]
    fn test_run_configuration_rejects_zero_window_size() {
        let reference = write_temp_file(">chr1\nACGTACGT\n");
        let coverage = write_temp_file("chr1\t1\t10\n");
        let output = tempfile::NamedTempFile::new().unwrap();
        let output_path = output.path().to_path_buf();
        drop(output);

        let yaml = format!(
            "\
reference: {}
coverage_file: {}
coverage_format: samtools-depth
bed_file: .
output_file: {}
overwrite_output: true
window_size: 0
window_stride: 100
",
            reference.path().display(),
            coverage.path().display(),
            output_path.display(),
        );
        let config_file = write_config(&yaml);

        let result = RunConfiguration::from(&config_file.path().to_path_buf());

        assert!(
            result.is_err(),
            "Expected zero window_size to return an error"
        );
    }

    #[test]
    fn test_run_configuration_rejects_zero_window_stride() {
        let reference = write_temp_file(">chr1\nACGTACGT\n");
        let coverage = write_temp_file("chr1\t1\t10\n");
        let output = tempfile::NamedTempFile::new().unwrap();
        let output_path = output.path().to_path_buf();
        drop(output);

        let yaml = format!(
            "\
reference: {}
coverage_file: {}
coverage_format: samtools-depth
bed_file: .
output_file: {}
overwrite_output: true
window_size: 100
window_stride: 0
",
            reference.path().display(),
            coverage.path().display(),
            output_path.display(),
        );
        let config_file = write_config(&yaml);

        let result = RunConfiguration::from(&config_file.path().to_path_buf());

        assert!(
            result.is_err(),
            "Expected zero window_stride to return an error"
        );
    }

    #[test]
    fn test_run_configuration_rejects_existing_output_when_overwrite_false() {
        let reference = write_temp_file(">chr1\nACGTACGT\n");
        let coverage = write_temp_file("chr1\t1\t10\n");
        let output = write_temp_file("already exists\n");

        let yaml = format!(
            "\
reference: {}
coverage_file: {}
coverage_format: samtools-depth
bed_file: .
output_file: {}
overwrite_output: false
",
            reference.path().display(),
            coverage.path().display(),
            output.path().display(),
        );
        let config_file = write_config(&yaml);

        let result = RunConfiguration::from(&config_file.path().to_path_buf());

        assert!(
            result.is_err(),
            "Expected existing output with overwrite_output=false to return an error"
        );
    }

    #[test]
    fn test_run_configuration_allows_existing_output_when_overwrite_true() {
        let reference = write_temp_file(">chr1\nACGTACGT\n");
        let coverage = write_temp_file("chr1\t1\t10\n");
        let output = write_temp_file("already exists\n");

        let yaml = format!(
            "\
reference: {}
coverage_file: {}
coverage_format: samtools-depth
bed_file: .
output_file: {}
overwrite_output: true
",
            reference.path().display(),
            coverage.path().display(),
            output.path().display(),
        );
        let config_file = write_config(&yaml);

        let config = RunConfiguration::from(&config_file.path().to_path_buf()).unwrap();

        assert_eq!(config.output_file, output.path());
        assert!(config.overwrite_output);
    }

    #[test]
    fn test_run_configuration_rejects_window_stride_greater_than_window_size() {
        let reference = write_temp_file(">chr1\nACGTACGT\n");
        let coverage = write_temp_file("chr1\t1\t10\n");
        let output = tempfile::NamedTempFile::new().unwrap();
        let output_path = output.path().to_path_buf();
        drop(output);

        let yaml = format!(
            "\
reference: {}
coverage_file: {}
coverage_format: samtools-depth
bed_file: .
output_file: {}
overwrite_output: true
window_size: 100
window_stride: 200
",
            reference.path().display(),
            coverage.path().display(),
            output_path.display(),
        );
        let config_file = write_config(&yaml);

        let result = RunConfiguration::from(&config_file.path().to_path_buf());

        assert!(
            result.is_err(),
            "Expected window_stride > window_size to return an error"
        );
    }

    #[test]
    fn test_run_configuration_allows_window_stride_equal_to_window_size() {
        let reference = write_temp_file(">chr1\nACGTACGT\n");
        let coverage = write_temp_file("chr1\t1\t10\n");
        let output = tempfile::NamedTempFile::new().unwrap();
        let output_path = output.path().to_path_buf();
        drop(output);

        let yaml = format!(
            "\
reference: {}
coverage_file: {}
coverage_format: samtools-depth
bed_file: .
output_file: {}
overwrite_output: true
window_size: 100
window_stride: 100
",
            reference.path().display(),
            coverage.path().display(),
            output_path.display(),
        );
        let config_file = write_config(&yaml);

        let config = RunConfiguration::from(&config_file.path().to_path_buf()).unwrap();
        assert_eq!(config.window_stride, 100);
    }
}

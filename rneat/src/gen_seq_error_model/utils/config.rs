use serde_yml::Value;
use std::collections::HashMap;
use std::fs;
use std::path::PathBuf;
use crate::gen_seq_error_model::errors::GenSeqErrorModelError;

#[derive(Debug, Clone)]
pub struct RunConfiguration {
    pub fastq_file: PathBuf,
    pub output_file: PathBuf,
    pub overwrite_output: bool,
    pub max_reads: usize,
    pub qual_offset: usize,
}

impl RunConfiguration {
    pub fn from(yml_file: &PathBuf) -> Result<Self, GenSeqErrorModelError> {
        let f = fs::File::open(yml_file).map_err(|e| {
            GenSeqErrorModelError::ConfigurationError(format!(
                "Cannot open config file {:?}: {}",
                yml_file, e
            ))
        })?;
        let scrape_config: HashMap<String, Value> = serde_yml::from_reader(f)
            .map_err(|e| GenSeqErrorModelError::ConfigurationError(e.to_string()))?;

        let fastq_file = PathBuf::from(
            scrape_config
                .get("fastq_file")
                .and_then(|v| v.as_str())
                .ok_or_else(|| {
                    GenSeqErrorModelError::ConfigurationError(
                        "Missing required field: fastq_file".to_string(),
                    )
                })?,
        );
        if !fastq_file.is_file() {
            return Err(GenSeqErrorModelError::ConfigurationError(format!(
                "fastq_file not found: {:?}",
                fastq_file
            )));
        }

        let output_file = PathBuf::from(
            scrape_config
                .get("output_file")
                .and_then(|v| v.as_str())
                .ok_or_else(|| {
                    GenSeqErrorModelError::ConfigurationError(
                        "Missing required field: output_file".to_string(),
                    )
                })?,
        );

        let overwrite_output = scrape_config
            .get("overwrite_output")
            .and_then(|v| v.as_bool())
            .unwrap_or(false);

        if !overwrite_output && output_file.is_file() {
            return Err(GenSeqErrorModelError::ConfigurationError(format!(
                "output_file already exists and overwrite_output is false: {:?}",
                output_file
            )));
        }

        let max_reads = scrape_config
            .get("max_reads")
            .and_then(|v| v.as_u64())
            .unwrap_or(0) as usize;

        let qual_offset = scrape_config
            .get("qual_offset")
            .and_then(|v| v.as_u64())
            .unwrap_or(33) as usize;

        Ok(RunConfiguration {
            fastq_file,
            output_file,
            overwrite_output,
            max_reads,
            qual_offset,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn test_run_configuration() {
        let temp_dir = tempfile::tempdir().unwrap();
        let fastq_path = temp_dir.path().join("test.fastq");
        std::fs::write(&fastq_path, "@read1\nACGT\n+\nIIII\n").unwrap();
        let output_path = temp_dir.path().join("model.json.gz");

        let yaml = format!(
            "fastq_file: {}\noutput_file: {}\noverwrite_output: true\nmax_reads: 1000\nqual_offset: 33\n",
            fastq_path.display(),
            output_path.display()
        );
        let mut tmp = tempfile::NamedTempFile::new().unwrap();
        write!(tmp, "{}", yaml).unwrap();

        let config = RunConfiguration::from(&tmp.path().to_path_buf()).unwrap();
        assert_eq!(config.max_reads, 1000);
        assert_eq!(config.qual_offset, 33);
        assert!(config.overwrite_output);
    }

    #[test]
    fn test_missing_fastq_file_errors() {
        let temp_dir = tempfile::tempdir().unwrap();
        let output_path = temp_dir.path().join("model.json.gz");
        let yaml = format!(
            "fastq_file: /nonexistent/path.fastq\noutput_file: {}\n",
            output_path.display()
        );
        let mut tmp = tempfile::NamedTempFile::new().unwrap();
        write!(tmp, "{}", yaml).unwrap();
        let result = RunConfiguration::from(&tmp.path().to_path_buf());
        assert!(result.is_err());
    }
}

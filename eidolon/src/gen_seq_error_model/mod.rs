pub mod errors;
pub mod utils;

use crate::gen_seq_error_model::{errors::GenSeqErrorModelError, utils::config::RunConfiguration};
use std::path::PathBuf;

pub fn main(config_file: &PathBuf) -> Result<(), GenSeqErrorModelError> {
    let config = RunConfiguration::from(config_file)?;
    utils::runner::runner(&config)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_main_end_to_end() {
        // Smoke test: config → RunConfiguration → runner → serialized model on disk.
        let manifest_dir = env!("CARGO_MANIFEST_DIR");
        let fastq = format!("{}/test_data/H1N1_read1.fq.gz", manifest_dir);
        let out_dir = tempfile::tempdir().unwrap();
        let output = out_dir.path().join("model.json.gz");
        let yaml = format!(
            "fastq_file: {}\noutput_file: {}\noverwrite_output: true\nmax_reads: 50\n",
            fastq,
            output.display()
        );
        let mut tmp = NamedTempFile::new().unwrap();
        write!(tmp, "{}", yaml).unwrap();
        main(&tmp.path().to_path_buf()).unwrap();
        assert!(
            output.exists(),
            "output model file should have been written"
        );
    }
}

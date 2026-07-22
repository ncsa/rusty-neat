use thiserror::Error;

use crate::gen_reads::errors::GenerateReadsError;

/// Errors from the native cancer (tumor/normal) read simulator. Wraps the
/// underlying gen-reads errors (each pass is a `run_neat` call) plus the
/// orchestration- and merge-specific failures.
#[derive(Error, Debug)]
pub enum GenCancerReadsError {
    #[error("Cancer config error: {0}")]
    ConfigError(String),
    #[error("Purity must be strictly in (0, 1); got {0}. For purity 0 or 1 use plain gen-reads.")]
    PurityOutOfRange(f64),
    #[error(
        "Per-pass coverage rounds to zero (normal={normal}x, tumor={tumor}x). \
         Raise total_coverage or pick a less extreme purity."
    )]
    PerPassCoverageZero { normal: usize, tumor: usize },
    #[error("gen-reads pass '{pass}' failed: {source}")]
    GenReadsPass {
        pass: &'static str,
        #[source]
        source: GenerateReadsError,
    },
    #[error("RNG init failed: {0}")]
    Rng(String),
    #[error("Merge failed: {0}")]
    MergeError(String),
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
}

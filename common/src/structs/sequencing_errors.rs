//! This enum covers all error types. It's unclear yet what, if anything, we need this for. 

pub enum SequencingErrorType {
    SNPError,
    InsertionError,
    DeletionError,
}

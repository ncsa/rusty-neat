use serde::{Deserialize, Serialize};
use thiserror::Error;
use log::*;

#[derive(Debug, Error)]
pub enum BedErrors {
    #[error("Error creating bed record: {0}")]
    BedRecordError(String)
}

#[derive(Debug, Serialize, Deserialize)]
pub struct BedRecord {
    contig: String,
    start: usize,
    end: usize,
}

impl BedRecord {
    pub fn new(contig: String, start: usize, end: usize) -> Result<Self, BedErrors> {
        if start >= end {
            return Err(BedErrors::BedRecordError("Region with start pos >= end_pos".to_string()))
        }

        Ok(BedRecord { contig, start, end })
    }

    pub fn len(&self) -> usize {
        self.end - self.start
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn new_bedrecord() {
        let contig = "chr1".to_string();
        let start = 0;
        let end = 900;
        let result = BedRecord::new(contig, start, end).unwrap();
        assert_eq!(result.start, 0);
        assert_eq!(result.end, 900)
    }

    #[test]
    #[should_panic]
    fn test_bad_record() {
        let contig = "chr1".to_string();
        let start = 900;
        let end = 0;
        BedRecord::new(contig, start, end).unwrap();
    }
}
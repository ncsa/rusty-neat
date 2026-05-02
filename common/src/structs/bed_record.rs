use serde::{Deserialize, Serialize};
use thiserror::Error;

#[derive(Debug, Error)]
pub enum BedErrors {
    #[error("Error creating bed record: {0}")]
    BedRecordError(String)
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct BedRecord {
    contig: String,
    pub start: usize,
    pub end: usize,
    other: Vec<String>,
}

impl BedRecord {
    pub fn new(contig: String, start: usize, end: usize, other: Vec<String>) -> Result<Self, BedErrors> {
        if start >= end {
            return Err(BedErrors::BedRecordError("Region with start pos >= end_pos".to_string()))
        }

        Ok(BedRecord { contig, start, end, other })
    }

    pub fn len(&self) -> usize {
        self.end - self.start
    }

    pub fn get_contig(&self) -> String {
        self.contig.to_owned()
    }

    pub fn contains(&self, contig: &str, position: usize) -> bool {
        if self.contig == contig {
            if (position >= self.start) && (position < self.end) {
                true
            } else {
                false
            }
        } else {
            false
        }
    }

    pub fn overlaps(&self, contig: &str, start: usize, end: usize) -> bool {
        // This if says if either start or end are in range, then it is an overlap.
        if self.contains(contig, start) || self.contains(contig, end) {
            return true
        }
        false
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_overlaps() {
        let record = BedRecord::new("chr1".to_string(), 100, 200, vec![]).unwrap();
        assert!(record.overlaps("chr1", 150, 250));  // start inside
        assert!(record.overlaps("chr1", 50, 100));   // end == start of record (contains(50)=false, contains(100)=true)
        assert!(!record.overlaps("chr1", 200, 300)); // fully after (end is exclusive)
        assert!(!record.overlaps("chr1", 0, 50));    // fully before
        assert!(!record.overlaps("chr2", 150, 250)); // wrong contig
    }

    #[test]
    fn test_contains() {
        let record = BedRecord::new("chr1".to_string(), 100, 200, vec![]).unwrap();
        assert!(record.contains("chr1", 100));  // start is inclusive
        assert!(record.contains("chr1", 199));  // last valid position
        assert!(!record.contains("chr1", 200)); // end is exclusive
        assert!(!record.contains("chr1", 99));  // before start
        assert!(!record.contains("chr2", 150)); // wrong contig
    }

    #[test]
    fn new_bedrecord() {
        let contig = "chr1".to_string();
        let start = 0;
        let end = 900;
        let other = vec!["b0001".to_string(), "65".to_string(), "+".to_string()];
        let result = BedRecord::new(contig, start, end, other).unwrap();
        assert_eq!(result.start, 0);
        assert_eq!(result.end, 900)
    }

    #[test]
    #[should_panic]
    fn test_bad_record() {
        let contig = "chr1".to_string();
        let start = 900;
        let end = 0;
        let other = vec!["b0001".to_string(), "65".to_string(), "+".to_string()];
        BedRecord::new(contig, start, end, other).unwrap();
    }
}
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
    pub mut_rate: Option<f64>,
}

impl BedRecord {
    pub fn new(contig: String, start: usize, end: usize, other: &str) -> Result<Self, BedErrors> {
        if start >= end {
            return Err(BedErrors::BedRecordError("Region with start pos >= end_pos".to_string()))
        }
        let mut_rate = Self::parse_other_for_mut(other);
        Ok(BedRecord { contig, start, end, mut_rate})
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

    pub fn parse_other_for_mut(other: &str) -> Option<f64> {
        if let Some(index) = other.find("mut_rate=") {
            let start_index = index + "mut_rate=".len();
            let remainder = &other[start_index..];
            // The value might be followed by a space, tab, or semicolon if there are other fields
            let end_index = remainder.find(|c: char| c.is_whitespace() || c == ';').unwrap_or(remainder.len());
            let value_str = &remainder[..end_index];
            
            value_str.parse::<f64>().ok()
        } else {
            None
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_parse_other_finds_mut_rate() {
        let other = "mut_rate=0.003";
        assert_eq!(BedRecord::parse_other_for_mut(&other), Some(0.003));
        
        let other = "name=gene1;mut_rate=0.005;score=100";
        assert_eq!(BedRecord::parse_other_for_mut(&other), Some(0.005));
        
        let other = "some_other_field mut_rate=0.001\tmore_fields";
        assert_eq!(BedRecord::parse_other_for_mut(&other), Some(0.001));
        
        let other = "no_rate_here";
        assert_eq!(BedRecord::parse_other_for_mut(&other), None);
        
        let other = "mut_rate=not_a_number";
        assert_eq!(BedRecord::parse_other_for_mut(&other), None);
    }

    #[test]
    fn test_overlaps() {
        let record = BedRecord::new(
            "chr1".to_string(),
            100,
            200,
            "").unwrap();
        assert!(record.overlaps("chr1", 150, 250));  // start inside
        assert!(record.overlaps("chr1", 50, 100));   // end == start of record (contains(50)=false, contains(100)=true)
        assert!(!record.overlaps("chr1", 200, 300)); // fully after (end is exclusive)
        assert!(!record.overlaps("chr1", 0, 50));    // fully before
        assert!(!record.overlaps("chr2", 150, 250)); // wrong contig
    }

    #[test]
    fn test_contains() {
        let record = BedRecord::new(
            "chr1".to_string(),
            100,
            200,
            "").unwrap();
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
        let other = "b0001\t65\t+";
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
        let other = "b0001\t65\t+";
        BedRecord::new(contig, start, end, other).unwrap();
    }
}
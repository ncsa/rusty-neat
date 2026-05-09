use thiserror::Error;
use crate::structs::nucleotides::Nucleotide;
use log::error;

#[derive(Error, Debug)]
pub enum SequenceBlockError {
    #[error("SequenceBlock: bad coordinates (start >= end)")]
    BadCoordinatesError,
    #[error("SequenceBlock: requested range out of bounds")]
    OutOfBoundsError,
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct SequenceMap {
    pub region_type: RegionType,
    pub start: usize,
    pub end: usize,
}

impl SequenceMap {
    pub fn from(region_type: RegionType, start: usize, end: usize) -> SequenceMap {
        SequenceMap { region_type, start, end }
    }

    pub fn get_len(&self) -> usize {
        self.end - self.start
    }
}

#[derive(Clone, Default, Debug, PartialEq)]
pub enum RegionType {
    #[default]
    NRegion,
    NonNRegion,
}

#[derive(Debug, Default)]
pub struct SequenceBlock {
    pub contig: String,
    pub ref_start: usize,
    pub ref_end: usize,
    pub sequence: Vec<Nucleotide>,
    pub sequence_map: Vec<SequenceMap>,
}

impl SequenceBlock {
    pub fn get_non_n_regions(&self) -> Vec<&SequenceMap> {
        self.sequence_map
            .iter()
            .filter(|m| m.region_type == RegionType::NonNRegion)
            .collect()
    }

    pub fn get_len(&self) -> usize {
        self.ref_end - self.ref_start
    }

    pub fn get_subseq(&self, request_start: usize, request_end: usize) -> Result<Vec<Nucleotide>, SequenceBlockError> {
        let block_start = request_start + self.ref_start;
        let block_end = request_end + self.ref_start;
        if request_start >= request_end {
            error!("Bad coordinates: start {} >= end {}", request_start, request_end);
            return Err(SequenceBlockError::BadCoordinatesError);
        }
        if (block_end <= self.ref_start)
            || (block_start > self.ref_end)
            || (block_end > self.ref_end)
            || (block_start < self.ref_start)
        {
            error!(
                "Coordinates out of bounds: request ({}, {}), block ({}, {})",
                request_start, request_end, self.ref_start, self.ref_end,
            );
            return Err(SequenceBlockError::OutOfBoundsError);
        }
        Ok(self.sequence[request_start..request_end].to_owned())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structs::nucleotides::Nucleotide::{A, C, G, N, T};

    fn make_block() -> SequenceBlock {
        SequenceBlock {
            contig: "chr1".to_string(),
            ref_start: 0,
            ref_end: 20,
            sequence: vec![A, A, A, A, A, C, C, C, C, C, G, G, G, G, G, T, T, T, T, T],
            sequence_map: vec![],
        }
    }

    #[test]
    fn test_get_subseq_valid() {
        let block = make_block();
        assert_eq!(block.get_subseq(5, 10).unwrap(), vec![C, C, C, C, C]);
    }

    #[test]
    fn test_get_subseq_full() {
        let block = make_block();
        let sub = block.get_subseq(0, 20).unwrap();
        assert_eq!(sub, block.sequence);
    }

    #[test]
    fn test_get_subseq_bad_coords() {
        let block = make_block();
        assert!(matches!(block.get_subseq(10, 5), Err(SequenceBlockError::BadCoordinatesError)));
    }

    #[test]
    fn test_get_subseq_out_of_bounds() {
        let block = make_block();
        assert!(matches!(block.get_subseq(15, 25), Err(SequenceBlockError::OutOfBoundsError)));
    }

    #[test]
    fn test_get_non_n_regions_mixed() {
        let block = SequenceBlock {
            contig: "chr1".to_string(),
            ref_start: 0,
            ref_end: 20,
            sequence: vec![N; 20],
            sequence_map: vec![
                SequenceMap::from(RegionType::NRegion, 0, 5),
                SequenceMap::from(RegionType::NonNRegion, 5, 15),
                SequenceMap::from(RegionType::NRegion, 15, 20),
            ],
        };
        let regions = block.get_non_n_regions();
        assert_eq!(regions.len(), 1);
        assert_eq!(regions[0].start, 5);
        assert_eq!(regions[0].end, 15);
    }

    #[test]
    fn test_get_non_n_regions_all_non_n() {
        let block = SequenceBlock {
            contig: "chr1".to_string(),
            ref_start: 0,
            ref_end: 20,
            sequence: vec![A; 20],
            sequence_map: vec![SequenceMap::from(RegionType::NonNRegion, 0, 20)],
        };
        assert_eq!(block.get_non_n_regions().len(), 1);
    }

    #[test]
    fn test_get_len() {
        let block = make_block();
        assert_eq!(block.get_len(), 20);
    }
}
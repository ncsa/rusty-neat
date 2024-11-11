use std::collections::HashMap;
use std::fs::read_to_string;
use std::str::FromStr;

// The SequenceBlock has coordinates that map back to the original reference. Local manipulations
// of the block can use local coordinates, but any variant locations etc will be modified by the
// ref start point
#[derive(Clone, Debug)]
pub struct SequenceBlock {
    // start is the index for when this feature begins on the reference sequence. Start is always
    // relative to the original sequence
    pub ref_start: usize,
    // end is where this feature ends, relative to the reference.
    pub ref_end: usize,
    pub length: usize,
    pub file: String,
}

impl SequenceBlock {
    pub fn new(ref_start: usize, ref_end:usize, file: String) -> Result<Self, ()> {
        if ref_start >= ref_end {
            panic!(
                "Bad coordinates for sequence block start: {} >= end: {}", ref_start, ref_end
            )
        } else {
            Ok(SequenceBlock {
                ref_start,
                ref_end,
                length: ref_end - ref_start,
                file,
            })
        }
    }

    pub fn len(&self) -> usize {
        self.length
    }

    pub fn retrieve(&mut self, start: usize, end: usize) -> Result<Vec<u8>, &'static str> {
        if start >= end {
            panic!(
                "Bad coordinates for sequence block retrieval: {} >= end: {}", start, end
            );
        } else if end > self.ref_end || start < self.ref_start {
            panic!(
                "Bad coordinates in sequence block retrieval -\
                 start: {}, end: {}, block start: {}, block end: {}",
                start, end, self.ref_start, self.ref_end,
            );
        } else {
            // There should only be one line, and since these are programmatically generated, there
            // shouldn't be an issue with security. If the format is wrong, then the map part should
            // throw an error. Had to use vec::from rather than collect for reasons unknown (can't
            // clone a u8?
            let buffer: Vec<u8> = read_to_string(&self.file)
                .unwrap()
                .lines()
                .map(|x: &str| u8::from_str(x).unwrap())
                .collect();
            // start for this function is expected to be a coordinate from the overall contig,
            // so we modify by the known start and end point of this block to get the final
            // coordinates
            let block_start = start - self.ref_start;
            let block_end = end - self.ref_start;
            Ok(buffer[block_start..block_end].to_owned())
        }

    }
}

#[derive(Debug)]
// Each ContigBlock contains multiple sequence blocks, which map to files containing sequences that
// are 128Kb long.
pub struct ContigBlock {
    pub name: String,
    pub len: usize,
    pub blocks: Vec<SequenceBlock>
}

impl ContigBlock {
    fn new(name: String, len: usize, blocks: Vec<SequenceBlock>) -> Result<Self, &'static str> {
        Ok(ContigBlock {
            name,
            len,
            blocks
        })
    }
}

#[derive(Debug)]
// The fasta map holds all the data for retrieving items from the fasta files
pub struct FastaMap {
    pub map: HashMap<String, ContigBlock>
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_basic_map() {
        // use crate::structs::nucleotides::Nuc::*;
        //
        // let read = vec![
        //     Unk, Unk, Unk, Unk, Unk, Unk, Unk, Unk, Unk, Unk,
        //     Unk, Unk, Ada, Cyt, Thy, Unk, Unk, Unk, Ada, Unk
        // ];
        // let fasta_sequence = Box::new(HashMap::from([
        //     ("chr 1".to_string(), read)
        // ]));
        // let fasta_map = FastaMap::build_from_reference(&fasta_sequence).unwrap();
        //
        // let expected_output = "FastaMap { map: {\"chr 1\": \
        // (20, [\
        // FastaBlock { start: 0, block_type: NBlock }, \
        // FastaBlock { start: 12, block_type: Allowed }, \
        // FastaBlock { start: 15, block_type: NBlock }, \
        // FastaBlock { start: 18, block_type: Allowed }, \
        // FastaBlock { start: 19, block_type: NBlock }\
        // ])} }".to_string();
        //
        // let test_out = format!("{:?}", fasta_map);
        //
        // assert_eq!(test_out, expected_output)
    }
}
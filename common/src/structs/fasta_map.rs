use std::cmp::PartialEq;
use std::collections::HashMap;
use log::debug;
use structs::variants::{Variant, VariantType};
use structs::nucleotides::Nuc;
use structs::nucleotides::Nuc::*;
use structs::variants::VariantType::Deletion;

#[derive(Debug)]
// The fasta map is a collection of blocks.
pub struct FastaMap {
    pub map: Vec<ContigBlock>
}

#[derive(Debug)]
// Each block contains the name, length, and the list of blocks that cover the sequence
pub struct ContigBlock {
    pub name: String,
    pub len: usize,
    pub blocks: Vec<SequenceBlock>
}

// The SequenceBlock coordinates map directly back to the reference. All positions are relative
// to that sequence, or (potentially) subsequence
#[derive(Copy, Clone, Debug)]
pub struct SequenceBlock {
    // start is the index for when this feature begins on the reference sequence. Start is always
    // relative to the original sequence
    pub start: usize,
    // end is where this feature ends, relative to the reference.
    pub end: usize,
    pub block_type: BlockType,
}

#[derive(Copy, Clone, Debug, PartialEq)]
pub enum BlockType {
    // Allowed = either A, C, G, or T.
    Standard,
    // N block is a block of N's or any unknown character
    NBlock,
    // Variant is a segment with several possible values
    Variant,
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

impl SequenceBlock {
    pub fn new(start: usize, end:usize, block_type: BlockType) -> Result<Self, &'static str> {
        Ok(SequenceBlock {
            start,
            end,
            block_type
        })
    }

    pub fn len(&self) -> usize {
        // absolute values protect against sequences that are numbered backwards for valid reasons
        (self.start - self.end).abs()
    }
}

impl FastaMap {
    pub fn new(chrom: &str, start: usize, chrom_size: usize) -> Result<Self, &'static str> {
        // By default, we'll create a map with a single entry, with start and end points at start
        // and end and all bases allowed.
        let default_block = SequenceBlock::new(
            start,
            chrom_size,
            BlockType::Standard
        ).unwrap();

        let default_contig_block = ContigBlock::new(
            chrom.to_string(),
            chrom_size,
            vec![default_block],
        ).unwrap();

        Ok(FastaMap {
            map: vec![default_contig_block]
        })
    }

    pub fn build_from_reference(
        fasta_sequences: &Box<HashMap<String, Vec<Nuc>>>
    ) -> Result<Self, &'static str> {
        // We are specifically calling this function build from "reference" because we are assuming
        // that a reference will contain only (case-insensitive) A, C, G, T, or N. Anything else is
        // considered an "N" in NEAT.
        let mut fasta_map = Vec::with_capacity(fasta_sequences.len());
        // we'll organize things in a vector and then export it as an array.
        // which means any updates will need to overwrite the array, which is fine because this
        // should be a very small object
        let mut sequence_blocks = Vec::new();
        let mut current_start = 0;
        let mut current_end = 1;
        let mut in_allowed_block = false;
        for (contig, sequence) in fasta_sequences.iter() {
            // we'll store n-block size here
            let mut n_block = 0;
            for i in 0..sequence.len() {
                let mut base = sequence.get(i);
                match base {
                    Some(N) => {
                        if in_allowed_block {
                            // This is the case when we just finished reading an Allowed block
                            sequence_blocks.push(SequenceBlock::new(
                                current_start,
                                current_end,
                                BlockType::Standard,
                            ).unwrap());
                            in_allowed_block = false;
                            current_start = current_end;
                            current_end = current_start;
                        }
                        else {
                            n_block += 1;
                            continue
                        }
                    },
                    Some(A) | Some(C) | Some(G) | Some(T) => {
                        if !in_allowed_block {
                            // start of a new allowed_block
                            in_allowed_block = true;
                            // This also means we just finished an N-block
                            let block_to_add = SequenceBlock::new(
                                current_start,
                                current_start + n_block,
                                BlockType::NBlock,
                            ).unwrap();
                            sequence_blocks.push(block_to_add);
                            current_start += n_block;
                            current_end = current_start;
                            n_block = 0
                        }
                        if i == sequence.len() - 1 {
                            // If we're at the end of a sequence, still in the Allowed block, we
                            // need to push or else we'll wind up missing the tail.
                            sequence_blocks.push(SequenceBlock::new(
                                current_start,
                                sequence.len(),
                                BlockType::Standard,
                            ).unwrap());

                        }
                        current_end += 1;
                    },
                    _ => {
                        debug!("Technically, FASTA format allows nucleotides other than A, C, G,\
                        T, and N, but for the purposes of NEAT, they are invalid. Please supply \
                        a FASTA file with only these characters");
                        return Err(panic!("Reference contained invalid nucleotides."));
                    }
                }
            }
            let contig_block = ContigBlock::new(
                contig.to_string(),
                sequence.len(),
                sequence_blocks,
            );
            // reset the vector
            sequence_blocks = Vec::new();
        }

        Ok(FastaMap {
            map: fasta_map,
        })
    }
}

pub fn insert_variant(
    fasta_blocks: &ContigBlock, variant: Variant, variant_location: usize
) -> ContigBlock {
    // Todo: re-write this using ContigBlock
    let mut return_blocks = Vec::with_capacity(fasta_blocks.len() + 2);
    for block in fasta_blocks {
        if (block.start < variant_location) && (variant_location < block.end_point()) {
            if block.block_type != BlockType::Standard {
                panic!("BUG! Trying to insert a variant into an invalid section of the genome.")
            }

            let new_block_before =
                SequenceBlock::new(
                    block.start,
                    variant_location - block.start,
                    block.block_type,
                ).unwrap();

            let variant_block =
                SequenceBlock::new(
                    variant_location,
                    variant.reference.len(),
                    BlockType::Variant,
                ).unwrap();

            let start = match variant.variant_type {
                Deletion => variant_location + variant.reference.len(),
                _ => variant_location + 1,
            };
            let new_block_after = SequenceBlock::new(
                start,
                block.end_point() - start,
                block.block_type,
            ).unwrap();
            // The lengths, after accounting for any indels, should be the same
            let total_length = new_block_before.len() + new_block_after.len() + {
                if variant.variant_type == Deletion  {
                    variant_block.len() - 1
                } else {
                    1
                }
            };
            debug_assert_eq!(total_length, block.len);
            return_blocks.push(new_block_before);
            return_blocks.push(variant_block);
            return_blocks.push(new_block_after)
        } else {
            return_blocks.push(*block);
        }
    }
    return_blocks
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structs::nucleotides::Nuc::*;

    #[test]
    fn test_basic_map() {
        use crate::structs::nucleotides::Nuc::*;

        let read = vec![N, N, N, N, N, N, N, N, N, N, N, N, A, C, T, N, N, N, A, N];
        let fasta_sequence = Box::new(HashMap::from([
            ("chr 1".to_string(), read)
        ]));
        let fasta_map = FastaMap::build_from_reference(&fasta_sequence).unwrap();

        let expected_output = "FastaMap { map: {\"chr 1\": \
        (20, [\
        FastaBlock { start: 0, block_type: NBlock }, \
        FastaBlock { start: 12, block_type: Allowed }, \
        FastaBlock { start: 15, block_type: NBlock }, \
        FastaBlock { start: 18, block_type: Allowed }, \
        FastaBlock { start: 19, block_type: NBlock }\
        ])} }".to_string();

        let test_out = format!("{:?}", fasta_map);

        assert_eq!(test_out, expected_output)
    }
}
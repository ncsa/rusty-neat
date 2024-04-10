use std::cmp::PartialEq;
use std::collections::HashMap;
use log::debug;
use structs::variants::{Variant, VariantType};
use structs::nucleotides::Nuc;
use structs::nucleotides::Nuc::*;
use structs::variants::VariantType::Deletion;

#[derive(Debug)]
pub struct FastaMap {
    pub map: HashMap<String, (usize, Vec<FastaBlock>)>
}

#[derive(Copy, Clone, Debug)]
pub struct FastaBlock {
    // start is the index for when this feature begins on the sequence. End is always the start of
    // the next feature or the length of the sequence.
    pub start: usize,
    // length doesn't prescribe where this block ends. This is intentional, since
    pub length: usize,
    pub block_type: BlockType,
}

#[derive(Copy, Clone, Debug, PartialEq)]
pub enum BlockType {
    // Allowed = either A, C, G, or T.
    Allowed,
    // N block is a block of N's
    NBlock,
    // Variant is a segment with several possible values
    Variant,
}

impl FastaBlock {
    pub fn new(start: usize, length:usize, block_type: BlockType) -> Result<Self, &'static str> {
        Ok(FastaBlock {
            start,
            length,
            block_type
        })
    }

    pub fn end_point(&self) -> usize {
        // Because we store the length, not the end point, we leave the end point to be calculated.
        // The logic here is that some blocks will point to specific places on the reference, but
        // variants the end point will be meaningless, since they are inserted orthogonal to the
        // fasta map. No matter how long the variant, it only maps to one position on the reference.
        todo!()
    }
}

impl FastaMap {
    pub fn new(chrom: &str, start: usize, chrom_size: usize) -> Result<Self, &'static str> {
        // By default, we'll create a map with a single entry, with start and end points at start
        // and end and all bases allowed.
        let default_block = FastaBlock::new(
            start,
            chrom_size,
            BlockType::Allowed
        ).unwrap();

        Ok(FastaMap {
            map: HashMap::from([(chrom.to_string(), (chrom_size, vec![default_block]))])
        })
    }

    pub fn build_from_reference(
        fasta_sequences: &Box<HashMap<String, Vec<Nuc>>>
    ) -> Result<Self, &'static str> {
        // We are specifically calling this function build from "reference" because we are assuming
        // that a reference will contain only A, C, G, T, or N. Technically, FASTA format allows
        // nucleotides outside these choices, but for the sake of NEAT, these are invalid.
        let mut fasta_map = HashMap::new();
        // we'll organize things in a vector and then export it as an array.
        // which means any updates will need to overwrite the array, which is fine because this
        // should be a very small object
        let mut temp_map = Vec::new();
        let mut current_start = 0;
        let mut current_end = 1;
        let mut in_allowed_block = false;
        for (contig, sequence) in fasta_sequences.iter() {
            let mut size = sequence.len();
            for i in 0..sequence.len() {
                let base = sequence.get(i);
                match base {
                    Some(N(num)) => {
                        if in_allowed_block {
                            // This is the case when we just finished reading an Allowed block
                            temp_map.push(FastaBlock::new(
                                current_start,
                                current_end,
                                BlockType::Allowed,
                            ).unwrap());
                            in_allowed_block = false;
                            current_start = current_end;
                            current_end = current_start;
                        }
                        let block_to_add = FastaBlock::new(
                            current_start,
                            current_start + num,
                            BlockType::NBlock,
                        ).unwrap();
                        temp_map.push(block_to_add);
                        current_start += num;
                        // We need to account for N's to determine the full size of the reference
                        size += num - 1;
                        current_end = current_start;
                    },
                    Some(A) | Some(C) | Some(G) | Some(T) => {
                        if !in_allowed_block {
                            // start of a new allowed_block
                            in_allowed_block = true;
                        }
                        if i == sequence.len() - 1 {
                            // If we're at the end of a sequence, still in the Allowed block, we
                            // need to push or else we'll wind up missing the tail.
                            temp_map.push(FastaBlock::new(
                                current_start,
                                sequence.len(),
                                BlockType::Allowed,
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
            // update the map
            fasta_map.insert(contig.clone(), (size, temp_map));
            // reset the vector
            temp_map = Vec::new();
        }

        Ok(FastaMap {
            map: fasta_map,
        })
    }
}

pub fn insert_variant(
    fasta_blocks: &Vec<FastaBlock>, variant: Variant, variant_location: usize
) -> Vec<FastaBlock> {
    let mut return_blocks = Vec::with_capacity(fasta_blocks.len() + 2);
    for block in fasta_blocks {
        if (block.start < variant_location) && (variant_location < block.end_point()) {
            if block.block_type != BlockType::Allowed {
                panic!("BUG! Trying to insert a variant into an invalid section of the genome.")
            }
            // The one hiccup in the idea of "inserting" variant is that a deletion is subtractive.
            // Therefore, deletions require splitting differently than the additive/nondestructive
            // insert or SNP (note that this isn't a comment on their biological function, merely
            // how they interact with the reference from a code perspective).
            //
            // Split the old block in two. The first step is the same regardless of variant type.
            let new_block_before =
                FastaBlock::new(
                    block.start,
                    variant_location - block.start,
                    block.block_type,
                ).unwrap();

            let variant_block =
                FastaBlock::new(
                    variant_location,
                    variant.reference.len(),
                    BlockType::Variant,
                ).unwrap();
            // Because deletions eliminate reference bases, we need to account for that here.
            // for example, let's say our reference sequence is [A, C, T, A, G]
            // The map then looks like this {0, 5, Allowed}. If we add an insertion of 1 base at 2,
            // [A, C, T, T, A, G], gives a ref sequence as [T] and an alt of [T, T], so we have
            // {0, 2, Allowed}, {2, 2, Variant}, {3, 2, Allowed}. You can read this as bases 0 and 1
            // are normal nucleotides, at base 2 on the reference, we have either T or TT in the
            // mutated subject, and then at base 3 we have the last 2 bases from the Allowed set.
            // If we reverse that, where [A, C, T, T, A, G] is the reference sequence and
            // [A, C, T, A, G] is the alt. We started then with {0, 6, Allowed}, and we instead have
            // {0, 2, Allowed}, {2, 2, Variant}, {4, 2, Allowed}.
            // this reads as the mutated sequence matches the first 2 bases of the ref, then there
            // is either a "TT" or a "T" next, followed by the reference from 4-6 (i.e., the rest
            // of it). That way, position 3 is no longer directly referenced by the map and there
            // is no chance of accidentally mutating a deleted base.
            let start = match variant.variant_type {
                Deletion => variant_location + variant.reference.len(),
                _ => variant_location + 1,
            };
            let new_block_after = FastaBlock::new(
                start,
                block.end_point() - start,
                block.block_type,
            ).unwrap();
            // The lengths, after accounting for any indels, should be the same
            let total_length = new_block_before.length + new_block_after.length + {
                if variant.variant_type == Deletion  {
                    variant_block.length - 1
                } else {
                    1
                }
            };
            debug_assert_eq!(total_length, block.length);
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

        let read = vec![N(12), A, C, T, N(3), A, N(1)];
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
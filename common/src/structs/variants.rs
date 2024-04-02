use std::collections::HashMap;
use rand::Rng;
use rand_chacha::ChaChaRng;
use crate::structs::nucleotides::Nuc;

#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash, Ord, PartialOrd)]
pub enum VariantType {
    SNP,
    Indel,
}

#[derive(Debug, PartialEq, Eq, Clone)]
pub struct Variant {
    // This is a basic holder for a variant. There are several aspects of the variant that we
    // don't need to store, such as which ploid the variant appears on and the context sequence,
    // which we only need for the output files, we can generate at the time we need it, such as the
    // insertion sequence.
    //
    // This is the type of variant. Any new types will need to be added to VariantType
    // to be implemented here.
    pub variant_type: VariantType,
    // Position on the reference where the variant starts
    pub position: usize,
    // The reference allele of interest. This is either one base or several bases for a deletion.
    pub reference: Vec<Nuc>,
    // The alternate allele of interest. This is either one base or several bases for an insertion.
    pub alternate: Vec<Nuc>,
    // the genotype will tell us which ploid is mutated in vcf, fastq and bam.
    pub genotype: Vec<u8>,
    // Saved later for convenience
    pub is_homozygous: bool,
}

impl Variant {
    pub fn new(
        variant_type: VariantType,
        position: usize,
        reference: &Vec<Nuc>,
        alternate: &Vec<Nuc>,
        genotype: Vec<u8>,
        is_homozygous: bool,
    ) -> Self {
        Variant {
            variant_type,
            position,
            reference: reference.clone(),
            alternate: alternate.clone(),
            genotype,
            is_homozygous,
        }
    }

    pub fn is_insertion(&self) -> bool { self.reference.len() < self.alternate.len() }

    pub fn apply(&self, sequence: Vec<Nuc>, start: usize) -> Vec<Nuc> {
        let location = self.position - start;
        debug_assert_eq!(self.reference[0], sequence[location]);
        match self.variant_type {
            VariantType::Indel => {
                if self.is_insertion() {
                    // since the alternate includes the base from the reference in the first
                    // position, we don't need to get it from the sequence.
                    [
                        sequence.get(..location).unwrap(),
                        self.alternate.as_slice(),
                        sequence.get(location+1..).unwrap()
                    ].concat()
                } else {
                    // +1 in the first segment so that we grab the reference base,
                    // what we are actually doing in the second is
                    //     (self.reference.len() - 1) // to account for not deleting the first base
                    //   + (location + 1) // to skip the reference base
                    // since -1 and +1 cancel out, we just have length + location.
                    // example:
                    //     sequence = [A, C, C, G, T, C, A, A, A, T, T, G, C, T]
                    //     start = 11243
                    //     variant = { variant_type: Indel, position: 11245
                    //                 reference: [C, G, T], alternate: [C] // Deletion of 2 bases
                    //                 ... } // other info not relevant for now.
                    // our location in this case is 2 (the index where the variant starts)
                    // part 1 gives us [A, C, C] (all indexes <= 2)
                    // part 2 gives us [C, A, A, A, T, T, G, C, T]
                    // concatenating the 2 together gives us the sequence with a deletion of 2 bases
                    [
                        sequence.get(..location+1).unwrap(),
                        sequence.get((self.reference.len() + location)..).unwrap(),
                    ].concat()
                }
            },
            VariantType::SNP => {
                [
                    sequence.get(..location).unwrap(),
                    self.alternate.as_slice(),
                    sequence.get((location + 1)..).unwrap(),
                ].concat()
            },
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use structs::nucleotides::Nuc::*;

    #[test]
    fn test_get_len() {
        let my_variant = Variant::new(
            VariantType::Indel,
            &vec![A, A, A, T, C, G, T, T, T, A],
            &vec![T, G, A, A, T, G],
            vec![0, 1],
            true,
        );
        assert_eq!(my_variant.get_length(), 5)
    }

    #[test]
    fn test_contains() {
        let test_variant = Variant::new(
            VariantType::Indel,
            &vec![A, A, A, T, C, G, T, T, T, A],
            &vec![T, G, A, A, T, G],
            vec![0, 1],
            false,
        );
        let test_pos = 13;
        assert!(test_variant.contains(test_pos));
    }
}

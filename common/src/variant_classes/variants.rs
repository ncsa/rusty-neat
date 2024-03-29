use crate::models;

use models::nucleotides::Nuc;
use std::cmp::Ordering;
use std::borrow::Cow;

#[derive(Debug, Eq, PartialEq, Hash, Ord, PartialOrd)]
pub enum VariantType {
    SNP,
    Indel,
}

#[derive(Debug, PartialEq, Eq)]
pub struct Variant<'a> {
    // This is a basic holder for a variant. There are several aspects of the variant that we
    // don't need to store, such as which ploid the variant appears on and the context sequence,
    // which we only need for the output files, we can generate at the time we need it.
    //
    // chromosome: uniquely identifies the contig from the fasta file. Stored as a "clone on write"
    //      object for convenience of writing out later. The `'a` I think means it lives as long as
    //      the variant struct lives, which is what we want.
    // position: the position within the reference (0-indexed) where the first ref base is
    //      located. This is the base before the event, in the case of indels, or of the event,
    //      in the case of SNPs.
    // reference: the vector of nucelotides on chromosome starting at position which are altered
    //      for this variant
    // alternate: the vector of nucleotides that take the place of reference.
    chromosome: Cow<'a, str>,
    position: usize,
    length: usize,
    reference: Vec<Nuc>,
    alternate: Vec<Nuc>,
    variant_type: VariantType,
}

impl Variant<'_> {
    fn new(
        chromosome: &str,
        position: usize,
        reference: Vec<Nuc>,
        alternate: Vec<Nuc>,
        variant_type: VariantType,
    ) -> Self {
        let length = Self::get_length(&reference, &alternate);
        Variant {
            chromosome: Cow::from(chromosome),
            position,
            length,
            reference,
            alternate,
            variant_type,
        }
    }
    fn get_length(reference: &Vec<Nuc>, alternate: &Vec<Nuc>) -> usize {
        // Returns the length of the variant. Since all we have at the moment are SNPs and indels,
        // this is pretty basic, but we may need to expand this out a bit more for different variant
        // types.
        //
        // Different lengths indicates an indel of some variety, or an SV.
        // Minus one because the first base of an indel/SV is unchanged
        if reference.len() == alternate.len() {
            reference.len()
        } else {
            if reference.len() > alternate.len() {
                reference.len() - 1
            } else {
                alternate.len() - 1
            }
        }
    }

    fn contains(&self, position: usize) -> bool {
        // Checks whether the given position is within the variant. If the position is in between
        // start and end, that ought to do it.
        self.position <= position && (self.position + self.get_length() > position)
    }
}

impl PartialOrd for Variant<'_> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        if self.chromosome != other.chromosome {
            return None;
        }
        Some(self.position.cmp(&other.position))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use lib::nucleotides::Nuc::*;

    #[test]
    fn test_get_len() {
        let my_variant = Variant::new(
            "chr1",
            10,
            vec![A, A, A, G, G, T],
            vec![A],
            VariantType::Indel,
        );

        assert_eq!(my_variant.length, 5)
    }
}
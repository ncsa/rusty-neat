//! These enums and structs help us keep track of variants added to the reference
//! It stores all the necessary data to write a variant out.

#[derive(Debug)]
pub enum VariantError {
    MalformedIndel,
    MalformedSnp,
    MalformedRef,
    MalformedAlt,
}

#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash, Ord, PartialOrd)]
pub enum VariantType {
    SNP,
    Insertion,
    Deletion,
}

#[derive(Debug, PartialEq, Eq, Clone, Hash)]
pub struct Variant {
    // This is a basic holder for a variant. There are several aspects of the variant that we
    // don't need to store, such as which ploid the variant appears on and the context sequence,
    // which we only need for the output files, we can generate at the time we need it, such as the
    // insertion sequence.
    //
    // This is the type of variant. Any new types will need to be added to VariantType
    // to be implemented here.
    pub variant_type: VariantType,
    // the location, in reference coordinates, of the start of this variant
    pub location: usize,
    // The reference allele of interest. This is either one base or several bases for a deletion.
    pub reference: Vec<u8>,
    // The alternate allele of interest. This is either one base or several bases for an insertion.
    pub alternate: Vec<u8>,
    // the genotype will tell us if it is heterozygous or homozygous.
    pub genotype: Vec<u8>,
}

impl Variant {
    pub fn new(
        variant_type: VariantType,
        location: usize,
        reference: &Vec<u8>,
        alternate: &Vec<u8>,
        genotype: &Vec<u8>,
    ) -> Result<Self, VariantError> {
        // This will generate a variant, storing the reference, alternate, start point relative to the contig, and a genotype
        // but first a quick sanity check
        match variant_type {
            VariantType::Insertion | 
            VariantType::Deletion => {
                if (reference.len() == alternate.len()) ||
                    reference[0] != alternate[0] {
                    return Err(VariantError::MalformedIndel)
                }
            },
            VariantType::SNP => { 
                if reference.len() != alternate.len() {
                    return Err(VariantError::MalformedSnp)
                }
            },
        }

        // Enforcing a policy that we assume elsewhere in the code
        if reference.len() < 1 {
            return Err(VariantError::MalformedRef)
        }
        if reference.len() < 1 {
            return Err(VariantError::MalformedAlt)
        }

        Ok(Variant {
            variant_type,
            location,
            reference: reference.clone(),
            alternate: alternate.clone(),
            genotype: genotype.clone(),
        })
    }

    pub fn is_homozygous(&self) -> Result<bool, VariantError> {
        // We won't have any 0/0[/0/0...] genotypes in NEAT, so genotype will either be 1/0 or 1/1 (or 1/1/0 v 1/1/1)
        // if there is a 0 in the genotype, it must be heterozygous (false)
        if self.genotype.contains(&0) { Ok(false) } else { Ok(true) }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structs::variants::VariantType::{SNP, Insertion, Deletion};

    #[test]
    fn test_variant_creation() {
        let variant = Variant {
            // Not actually a valid variant, but just testing the constructor
            variant_type: SNP,
            location: 222,
            reference: vec![0, 1, 3, 2],
            alternate: vec![0],
            genotype: vec![1, 1, 1],
        };
        assert_eq!(variant.variant_type, SNP);
        assert_eq!(variant.variant_type, SNP);
        assert!(variant.is_homozygous().unwrap());
    }

    #[test]
    fn test_variant_new() {
        let variant = Variant::new(
            Insertion,
            10,
            &vec![0],
            &vec![0, 1, 3, 2],
            &vec![0, 1],
        ).unwrap();
        assert!(!variant.is_homozygous().unwrap())
    }

    #[test]
    fn test_deletion_new() {
        let variant = Variant::new(
            Deletion,
            22,
            &vec![0, 1, 3, 2],
            &vec![0],
            &vec![0, 1],
        ).unwrap();
        assert!(!variant.is_homozygous().unwrap())
    }

    #[test]
    #[should_panic]
    fn test_bad_variant_creation() {
        Variant::new(
            // with new it should catch this error
            SNP,
            22,
            &vec![0, 1, 3, 2],
            &vec![0],
            &vec![1, 1, 1],
        ).unwrap();
    }
}

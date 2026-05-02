//! These enums and structs help us keep track of variants added to the reference
//! It stores all the necessary data to write a variant out.
use thiserror::Error;
use crate::structs::nucleotides::Nucleotide;

#[derive(Error, Debug)]
pub enum VariantError {
    #[error("ERROR: Malformed indel variant!")]
    MalformedIndel,
    #[error("ERROR: Malformed SNP variant!")]
    MalformedSnp,
    #[error("ERROR: Malformed variant ref!")]
    MalformedRef,
    #[error("ERROR: Malformed variant alt!")]
    MalformedAlt,
    #[error("Error Generating the genotype string!")]
    GenoStringError,
}

#[derive(Debug, Clone, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub enum Genotypes {
    Heterozygous,
    Homozygous,
}

#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash, Ord, PartialOrd)]
pub enum VariantType {
    SNP,
    Insertion,
    Deletion,
}

impl From<usize> for VariantType {
    fn from(i: usize) -> Self {
        match i {
            0 => Self::SNP,
            1 => Self::Insertion,
            2 => Self::Deletion,
            _ => panic!("Index out of range!")
        }
    }
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
    pub reference: Vec<Nucleotide>,
    // The alternate allele of interest. This is either one base or several bases for an insertion.
    pub alternate: Vec<Nucleotide>,
    // the genotype string is for writing in the vcf. 
    pub genotype_str: String,
    // The genotype is either heterozygous (not on all alleles) or homozygous (on all alleles)
    pub genotype: Genotypes,
}

impl Variant {

    pub fn new(
        variant_type: VariantType,
        location: usize,
        reference: &Vec<Nucleotide>,
        alternate: &Vec<Nucleotide>,
        genotype: &mut Vec<usize>,
    ) -> Result<Self, VariantError> {
        // Check for empty vecs before any index access
        if reference.is_empty() {
            return Err(VariantError::MalformedRef)
        }
        if alternate.is_empty() {
            return Err(VariantError::MalformedAlt)
        }

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

        let mut genotype_label = Genotypes::Homozygous;
        // Homozygous would mean all alleles have the variant, so all are 1's
        // the presence of a single 0 indicates heterozygous.
        // This is a little shakey biologically, but we're going with this.
        if genotype.contains(&0) {
            genotype_label = Genotypes::Heterozygous;
        }
        
        Ok(Variant {
            variant_type,
            location,
            reference: reference.clone(),
            alternate: alternate.clone(),
            genotype_str: genotype_to_string(genotype)?,
            genotype: genotype_label,
        })
    }

    pub fn get_loc(&self) -> Result<usize, VariantError> {
        Ok(self.location.clone())
    }
}

fn genotype_to_string(genotype: &mut Vec<usize>) -> Result<String, VariantError> {
    // Converts a vector of 0s and 1s representing genotype to a standard
    // vcf genotype string.
    //
    // In order to make the vcf look a little more realistic, we shuffle the genotypes, that way
    // it doesn't look like one ploid got all the mutations. I actually don't know if this is
    // realistic or not, but that's what we're doing.
    let geno_str = genotype.iter().map(|x| x.to_string() + "/").collect::<String>();
    match geno_str.strip_suffix("/") {
        Some(str) => Ok(str.to_string()),
        None => Err(VariantError::GenoStringError)
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
            reference: vec![Nucleotide::A, Nucleotide::C, Nucleotide::T, Nucleotide::G],
            alternate: vec![Nucleotide::A],
            genotype_str: "1/1/1".to_string(),
            genotype: Genotypes::Homozygous,
        };
        assert_eq!(variant.variant_type, SNP);
        assert_eq!(variant.reference, vec![Nucleotide::A, Nucleotide::C, Nucleotide::T, Nucleotide::G]);
        assert_eq!(variant.genotype, Genotypes::Homozygous);
    }

    #[test]
    fn test_variant_new() {
        let variant = Variant::new(
            Insertion,
            10,
            &vec![Nucleotide::A],
            &vec![Nucleotide::A, Nucleotide::C, Nucleotide::T, Nucleotide::G],
            &mut vec![0, 1],
        ).unwrap();
        assert_eq!(variant.genotype, Genotypes::Heterozygous)
    }

    #[test]
    fn test_deletion_new() {
        let variant = Variant::new(
            Deletion,
            22,
            &vec![Nucleotide::A, Nucleotide::C, Nucleotide::T, Nucleotide::G],
            &vec![Nucleotide::A],
            &mut vec![1, 1],
        ).unwrap();
        assert_eq!(variant.genotype, Genotypes::Homozygous)
    }

    #[test]
    #[should_panic(expected = "MalformedSnp")]
    fn test_bad_variant_creation() {
        Variant::new(
            SNP,
            22,
            &vec![Nucleotide::A, Nucleotide::C, Nucleotide::T, Nucleotide::G],
            &vec![Nucleotide::A],
            &mut vec![1, 1, 1],
        ).unwrap();
    }

    #[test]
    fn test_genotype_to_string() {
        let mut genotype = vec![0, 1, 0];
        assert_eq!(String::from("0/1/0"), genotype_to_string(&mut genotype).unwrap());
    }

    #[test]
    fn test_empty_reference_returns_malformed_ref() {
        let result = Variant::new(
            SNP,
            0,
            &vec![],
            &vec![Nucleotide::G],
            &mut vec![1, 1],
        );
        assert!(matches!(result, Err(VariantError::MalformedRef)));
    }

    #[test]
    fn test_empty_alternate_returns_malformed_alt() {
        let result = Variant::new(
            SNP,
            0,
            &vec![Nucleotide::A],
            &vec![],
            &mut vec![1, 1],
        );
        assert!(matches!(result, Err(VariantError::MalformedAlt)));
    }
}

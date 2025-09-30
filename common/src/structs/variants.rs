//! These enums and structs help us keep track of variants added to the reference
//! It stores all the necessary data to write a variant out.
use serde::{Deserialize, Serialize};
use thiserror::Error;
use crate::structs::nucleotides::Nucleotide;
use log::*;

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
    #[error("Invalid Vcf Input: {0}")]
    InvalidVcf(String),
}

#[derive(Debug, Clone, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub enum Genotype {
    Heterozygous,
    Homozygous,
}

fn gt_from_str(input: &str) -> Genotype {
    let mut result: Vec<&str> = input.split("/").collect();
    if result.len() == 1 {
        // Split may have failed, trying a different delimiter
        result = input.split("|").collect();
    }
    let mut found_zero = false;
    for item in result {
        if item.parse::<u64>().unwrap() == 0 {
            found_zero = true;
            break
        }
    }
    if found_zero {
        Genotype::Heterozygous
    } else {
        Genotype::Homozygous
    }
}

#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash, Ord, PartialOrd, Serialize, Deserialize)]
pub enum VariantType {
    SNP,
    Insertion,
    Deletion,
    Complex,
}

impl From<usize> for VariantType {
    fn from(i: usize) -> Self {
        match i {
            0 => Self::SNP,
            1 => Self::Insertion,
            2 => Self::Deletion,
            3 => Self::Complex,
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
    // the genotype string is feor writing in the vcf. 
    pub genotype_str: String,
    // The genotype is either heterozygous (not on all alleles) or homozygous (on all alleles)
    pub genotype: Genotype,

    // VCF files specific. These are variables specific to VCF files read as input,
    // either for variant inputs or for generating models.

    // Included id, if present
    pub id: Option<String>,
    // Included quality score
    pub quality_score: Option<usize>,
    // Included filter, if present
    pub filter: Option<String>,
    // Included Info, if present
    pub info: Option<String>,
    // Other included information, if present
    pub format: Vec<String>,
    // Sample, if present, must be same length as format
    pub sample: Vec<String>
}

impl Variant {

    pub fn from_file_params(
        location: usize,
        id: &str,
        filter: &str,
        info_field: &str,
        vcf_reference: &str,
        vcf_alternate: &str,
        quality_score: usize,
        format: Vec<String>,
        sample: Vec<String>,
    ) -> Result<Self, VariantError> {

        let genotype_str = {
            if format.contains(&String::from("GT")) {
                info!("Found genotype in vcf file");
                let gt_pos = format
                    .clone()
                    .into_iter()
                    .position(|x| x == String::from("GT"))
                    .unwrap();
                sample[gt_pos].to_string()
            } else {
                return Err(
                    VariantError::InvalidVcf(
                        "Missing GT field. Cannot process input VCF without GT in FORMAT field".to_string()
                    )
                )
            }
        };
    
        let mut variant_type = VariantType::SNP;
        if vcf_reference.len() > 1 {
            if vcf_alternate.len() > 1 {
                variant_type = VariantType::Complex;
            } else if vcf_alternate.len() == 1 {
                variant_type = VariantType::Deletion
            } else {
                return Err(
                    VariantError::InvalidVcf(
                        "ALT field in vcf with empty alt. Invalid VCF.".to_string()
                    )
                )
            }
        } else if vcf_reference.len() == 1 {
            if vcf_alternate.len() > 1 {
                variant_type = VariantType::Insertion;
            } else {
                // snp
            }
        } else {
            return Err(
                VariantError::InvalidVcf(
                    "REF field in vcf with empty ref. Invalid VCF.".to_string()
                )
            )
        }
        let mut reference = Vec::new();
        for char in vcf_reference.chars() {
            reference.push(Nucleotide::from(char))
        }
        let mut alternate = Vec::new();
        for char in vcf_alternate.chars() {
            alternate.push(Nucleotide::from(char))
        }
        let genotype = gt_from_str(&genotype_str);
        Ok(Variant{
            variant_type,
            location,
            reference,
            alternate,
            genotype_str,
            genotype,
            id: Some(id.to_string()),
            quality_score: Some(quality_score),
            filter: Some(filter.to_string()),
            info: Some(info_field.to_string()),
            format,
            sample
        })

    }

    pub fn new(
        variant_type: VariantType,
        location: usize,
        reference: &Vec<Nucleotide>,
        alternate: &Vec<Nucleotide>,
        genotype: &mut Vec<usize>,
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
            _ => {
                panic!("Complex variants not yet supported.")
            }
        }

        // Enforcing a policy that we assume elsewhere in the code
        if reference.len() < 1 {
            return Err(VariantError::MalformedRef)
        }
        if reference.len() < 1 {
            return Err(VariantError::MalformedAlt)
        }

        let mut genotype_label = Genotype::Homozygous;
        // Homozygous would mean all alleles have the variant, so all are 1's
        // the presence of a single 0 indicates heterozygous.
        // This is a little shakey biologically, but we're going with this.
        if genotype.contains(&0) {
            genotype_label = Genotype::Heterozygous;
        }
        
        Ok(Variant {
            variant_type,
            location,
            reference: reference.clone(),
            alternate: alternate.clone(),
            genotype_str: genotype_to_string(genotype)?,
            genotype: genotype_label,
            id: None,
            quality_score: None,
            filter: None,
            info: None,
            format: Vec::new(),
            sample: Vec::new()
        })
    }

    pub fn get_loc(&self) -> Result<usize, VariantError> {
        Ok(self.location)
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
            genotype: Genotype::Homozygous,
            id: None,
            quality_score: None,
            filter: None,
            info: None,
            format: Vec::new(),
            sample: Vec::new(),
        };
        assert_eq!(variant.variant_type, SNP);
        assert_eq!(variant.reference, vec![Nucleotide::A, Nucleotide::C, Nucleotide::T, Nucleotide::G]);
        assert_eq!(variant.genotype, Genotype::Homozygous);
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
        assert_eq!(variant.genotype, Genotype::Heterozygous)
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
        assert_eq!(variant.genotype, Genotype::Homozygous)
    }

    #[test]
    #[should_panic]
    fn test_bad_variant_creation() {
        Variant::new(
            // with new it should catch this error
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
}

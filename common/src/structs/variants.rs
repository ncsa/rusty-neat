//! These enums and structs help us keep track of variants added to the reference
//! It stores all the necessary data to write a variant out.
use crate::structs::nucleotides::{Nucleotide, is_acgtn};
use log::*;
use serde::{Deserialize, Serialize};
use thiserror::Error;

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
        if item == "." {
            continue;
        }
        if item.parse::<u64>().unwrap_or(1) == 0 {
            found_zero = true;
            break;
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
            _ => panic!("Index out of range!"),
        }
    }
}

#[derive(Debug, PartialEq, Eq, Clone, Hash)]
pub struct SvData {
    pub alternate: String
}

#[derive(Debug, PartialEq, Eq, Clone, Hash)]
pub enum AlternateType {
    Literal(Vec<Nucleotide>),
    Symbolic(SvData),
}

impl AlternateType {
    pub fn get_vec(&self) -> Option<Vec<Nucleotide>> {
        match self {
            AlternateType::Literal(vector) => Some(vector.clone()),
            _ => None,
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
<<<<<<< Updated upstream
    // The alternate allele of interest. This is either one base or several bases for an insertion.
    // For complex variants, we may need to make this an Option. Nones indicate special cases.
    // or add another enum for alternate types
    pub alternate: Vec<Nucleotide>,
=======
    // The alternate allele of interest. This is either one base or several bases for an indel, or
    // it is a key string from the vcf spec indicating a complex variant.
    pub alternate: AlternateType,
>>>>>>> Stashed changes
    // the genotype string is for writing in the vcf.
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
    pub sample: Vec<String>,
}

impl Variant {
    pub fn from_file(
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
                debug!("Found genotype in vcf file");
                let gt_pos = format.clone().into_iter().position(|x| x == "GT").unwrap();
                sample[gt_pos].to_string()
            } else {
                return Err(VariantError::InvalidVcf(
                    "Missing GT field. Cannot process input VCF without GT in FORMAT field"
                        .to_string(),
                ));
            }
        };

        let variant_type = parse_alternate(vcf_reference, vcf_alternate)?;
        let mut reference = Vec::new();
        for char in vcf_reference.chars() {
            reference.push(Nucleotide::from(char))
        }
        let mut alternate = Vec::new();
        for char in vcf_alternate.chars() {
            alternate.push(Nucleotide::from(char))
        }
        let genotype = gt_from_str(&genotype_str);
        Ok(Variant {
            variant_type,
            location,
            reference,
            alternate: AlternateType::Literal(alternate),
            genotype_str,
            genotype,
            id: Some(id.to_string()),
            quality_score: Some(quality_score),
            filter: Some(filter.to_string()),
            info: Some(info_field.to_string()),
            format,
            sample,
        })
    }

    pub fn new(
        variant_type: VariantType,
        location: usize,
        reference: &Vec<Nucleotide>,
        alternate: &Vec<Nucleotide>,
        genotype: &mut Vec<usize>,
    ) -> Result<Self, VariantError> {
        // Check for empty vecs before any index access
        if reference.is_empty() {
            return Err(VariantError::MalformedRef);
        }
        if alternate.is_empty() {
            return Err(VariantError::MalformedAlt);
        }

        match variant_type {
            VariantType::Insertion | VariantType::Deletion => {
                if (reference.len() == alternate.len()) || reference[0] != alternate[0] {
                    return Err(VariantError::MalformedIndel);
                }
            }
            VariantType::SNP => {
                if reference.len() != alternate.len() {
                    return Err(VariantError::MalformedSnp);
                }
            }
            _ => {
                panic!("Complex variants not yet supported.")
            }
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
            alternate: AlternateType::Literal(alternate.clone()),
            genotype_str: genotype_to_string(genotype)?,
            genotype: genotype_label,
            id: None,
            quality_score: None,
            filter: None,
            info: None,
            format: Vec::new(),
            sample: Vec::new(),
        })
    }

    pub fn get_loc(&self) -> Result<usize, VariantError> {
        Ok(self.location)
    }
}

fn parse_alternate(reference: &str, alt: &str) -> Result<VariantType, VariantError> {
    if reference.is_empty() {
        return Err(VariantError::InvalidVcf(
            "REF field in vcf with empty ref. Invalid VCF.".to_string(),
        ));
    }
    if alt.is_empty() {
        return Err(VariantError::InvalidVcf(
            "ALT field in vcf with empty alt. Invalid VCF.".to_string(),
        ));
    }
    if !is_acgtn(reference) {
        return Err(VariantError::InvalidVcf(format!(
            "REF must contain only A/C/G/T/N; got {reference}"
        )));
    }
    // Symbolic / structural ALTs (e.g. <DUP:TANDEM>, <DEL>, <INS:ME:ALU>, breakends
    // like `G]17:198982]`) are valid per VCF 4.2 but not yet representable by
    // VariantType. Surface them as an InvalidVcf rather than silently coercing.
    if !is_acgtn(alt) {
        return Err(VariantError::InvalidVcf(format!(
            "symbolic / structural ALT not yet supported: {alt}"
        )));
    }
    match (reference.len(), alt.len()) {
        (1, 1) => Ok(VariantType::SNP),
        (1, _) => Ok(VariantType::Insertion),
        (_, 1) => Ok(VariantType::Deletion),
        (_, _) => Ok(VariantType::Complex),
    }
}

fn genotype_to_string(genotype: &mut Vec<usize>) -> Result<String, VariantError> {
    // Converts a vector of 0s and 1s representing genotype to a standard
    // vcf genotype string.
    //
    // In order to make the vcf look a little more realistic, we shuffle the genotypes, that way
    // it doesn't look like one ploid got all the mutations. I actually don't know if this is
    // realistic or not, but that's what we're doing.
    let geno_str = genotype
        .iter()
        .map(|x| x.to_string() + "/")
        .collect::<String>();
    match geno_str.strip_suffix("/") {
        Some(str) => Ok(str.to_string()),
        None => Err(VariantError::GenoStringError),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structs::variants::VariantType::{Complex, Deletion, Insertion, SNP};

    #[test]
    fn test_variant_creation() {
        let variant = Variant {
            // Not actually a valid variant, but just testing the constructor
            variant_type: SNP,
            location: 222,
            reference: vec![Nucleotide::A, Nucleotide::C, Nucleotide::T, Nucleotide::G],
            alternate: AlternateType::Literal(vec![Nucleotide::A]),
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
        assert_eq!(
            variant.reference,
            vec![Nucleotide::A, Nucleotide::C, Nucleotide::T, Nucleotide::G]
        );
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
        )
        .unwrap();
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
        )
        .unwrap();
        assert_eq!(variant.genotype, Genotype::Homozygous)
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
        )
        .unwrap();
    }

    #[test]
    fn test_genotype_to_string() {
        let mut genotype = vec![0, 1, 0];
        assert_eq!(
            String::from("0/1/0"),
            genotype_to_string(&mut genotype).unwrap()
        );
    }

    #[test]
    fn test_empty_reference_returns_malformed_ref() {
        let result = Variant::new(SNP, 0, &vec![], &vec![Nucleotide::G], &mut vec![1, 1]);
        assert!(matches!(result, Err(VariantError::MalformedRef)));
    }

    #[test]
    fn test_empty_alternate_returns_malformed_alt() {
        let result = Variant::new(SNP, 0, &vec![Nucleotide::A], &vec![], &mut vec![1, 1]);
        assert!(matches!(result, Err(VariantError::MalformedAlt)));
    }

    #[test]
    fn parse_alternate_snp() {
        assert_eq!(parse_alternate("A", "C").unwrap(), SNP);
    }

    #[test]
    fn parse_alternate_insertion() {
        assert_eq!(parse_alternate("A", "ACGT").unwrap(), Insertion);
    }

    #[test]
    fn parse_alternate_deletion() {
        assert_eq!(parse_alternate("ACGT", "A").unwrap(), Deletion);
    }

    #[test]
    fn parse_alternate_complex_multibase_substitution() {
        assert_eq!(parse_alternate("AC", "GT").unwrap(), Complex);
    }

    #[test]
    fn parse_alternate_complex_unequal_length() {
        // Both > 1 but different lengths is still Complex (e.g., 2bp ref → 3bp alt).
        assert_eq!(parse_alternate("AC", "GTA").unwrap(), Complex);
    }

    #[test]
    fn parse_alternate_accepts_lowercase_and_n() {
        // Soft-masked bases and N should still classify as literal-base variants.
        assert_eq!(parse_alternate("a", "n").unwrap(), SNP);
        assert_eq!(parse_alternate("A", "acgN").unwrap(), Insertion);
    }

    #[test]
    fn parse_alternate_empty_reference_is_invalid_vcf() {
        let err = parse_alternate("", "A").unwrap_err();
        match err {
            VariantError::InvalidVcf(msg) => assert!(msg.contains("REF")),
            other => panic!("expected InvalidVcf, got {other:?}"),
        }
    }

    #[test]
    fn parse_alternate_empty_alt_is_invalid_vcf() {
        let err = parse_alternate("A", "").unwrap_err();
        match err {
            VariantError::InvalidVcf(msg) => assert!(msg.contains("ALT")),
            other => panic!("expected InvalidVcf, got {other:?}"),
        }
    }

    #[test]
    fn parse_alternate_symbolic_alt_is_invalid_vcf() {
        // VCF 4.2 symbolic ALTs aren't representable as a Vec<Nucleotide> yet.
        for alt in ["<DUP:TANDEM>", "<DEL>", "<INS:ME:ALU>", "<CNV>"] {
            let err = parse_alternate("A", alt).unwrap_err();
            match err {
                VariantError::InvalidVcf(msg) => assert!(
                    msg.contains("symbolic") || msg.contains("structural"),
                    "expected symbolic/structural message for {alt}, got: {msg}"
                ),
                other => panic!("expected InvalidVcf for {alt}, got {other:?}"),
            }
        }
    }

    #[test]
    fn parse_alternate_breakend_alt_is_invalid_vcf() {
        // VCF 4.2 breakend notation like `G]17:198982]` or `]13:123456]T`.
        let err = parse_alternate("G", "G]17:198982]").unwrap_err();
        assert!(matches!(err, VariantError::InvalidVcf(_)));
    }

    #[test]
    fn parse_alternate_non_iupac_ref_is_invalid_vcf() {
        // REF must be literal bases — a stray bracket here means a malformed VCF row.
        let err = parse_alternate("A<", "C").unwrap_err();
        match err {
            VariantError::InvalidVcf(msg) => assert!(msg.contains("REF")),
            other => panic!("expected InvalidVcf, got {other:?}"),
        }
    }
}

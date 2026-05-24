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
    #[error("Complex variants are not yet supported.")]
    ComplexNotSupported,
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

/// Classifies a symbolic / structural ALT (VCF 4.2 §1.4.1 — §1.4.5).
///
/// `Unknown` is the catch-all for symbolic ALTs we recognize as symbolic but
/// can't classify — e.g. a tag like `<NOVEL_THING>` that isn't in the spec. We
/// still round-trip the raw ALT string verbatim from `SvData::raw_alt`.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash, Serialize, Deserialize)]
pub enum SvType {
    Del,
    Dup,
    Cnv,
    Ins,
    Inv,
    Bnd,
    Unknown,
}

impl SvType {
    /// Parse a VCF symbolic / breakend ALT string into a coarse type tag.
    /// Returns `None` for literal-base ALTs (callers should not invoke this
    /// on `AlternateType::Literal` payloads).
    ///
    /// Recognizes:
    ///   - `<DEL>`, `<DEL:...>` → Del
    ///   - `<DUP>`, `<DUP:TANDEM>`, `<DUP:...>` → Dup
    ///   - `<CNV>`, `<CNV:...>` → Cnv
    ///   - `<INS>`, `<INS:ME:ALU>`, `<INS:...>` → Ins
    ///   - `<INV>` → Inv
    ///   - Breakend notation (`G]17:198982]`, `]13:123456]T`, etc.) → Bnd
    ///   - Other `<TAG>` forms → Unknown
    pub fn from_alt_string(alt: &str) -> Option<Self> {
        if alt.is_empty() {
            return None;
        }
        if alt.starts_with('<') && alt.ends_with('>') {
            let body = &alt[1..alt.len() - 1];
            let head = body.split(':').next().unwrap_or(body);
            return Some(match head {
                "DEL" => SvType::Del,
                "DUP" => SvType::Dup,
                "CNV" => SvType::Cnv,
                "INS" => SvType::Ins,
                "INV" => SvType::Inv,
                _ => SvType::Unknown,
            });
        }
        // Breakend notation contains `[` or `]` and a `:`-separated mate locus.
        if (alt.contains('[') || alt.contains(']')) && alt.contains(':') {
            return Some(SvType::Bnd);
        }
        None
    }
}

/// Parsed payload for a symbolic / structural ALT.
///
/// Stores the original ALT string for verbatim round-trip, plus optional
/// INFO-derived fields (`END=`, `SVLEN=`, `CN=`). The Phase 2 reader will
/// populate the optional fields from the VCF INFO column; today, anything
/// constructing an `SvData` must fill them in itself.
#[derive(Debug, PartialEq, Eq, Clone, Hash)]
pub struct SvData {
    /// The original ALT field, verbatim — used to round-trip back to the
    /// output VCF without re-synthesizing the angle-bracket notation.
    pub raw_alt: String,
    pub sv_type: SvType,
    /// 1-based inclusive end coordinate from `INFO/END`, if present.
    pub end: Option<usize>,
    /// `INFO/SVLEN` if present. Signed per VCF 4.2: negative for deletions,
    /// positive for insertions/duplications.
    pub svlen: Option<i64>,
    /// `INFO/CN` if present. Total copy number for `<CNV>` / `<DUP>` records.
    pub copy_number: Option<u32>,
}

impl SvData {
    /// Construct an `SvData` with only the raw ALT and type tag populated.
    /// Optional INFO fields default to `None` — set them on the returned
    /// value if available.
    pub fn new(raw_alt: impl Into<String>, sv_type: SvType) -> Self {
        SvData {
            raw_alt: raw_alt.into(),
            sv_type,
            end: None,
            svlen: None,
            copy_number: None,
        }
    }

    /// Compute the reference span (in bases) of this SV, given the variant's
    /// 1-based start `location`. Prefers `END` (span = end - location + 1);
    /// falls back to `|SVLEN|`. Returns `None` if neither is set.
    pub fn span(&self, location: usize) -> Option<usize> {
        if let Some(end) = self.end {
            return Some(end.saturating_sub(location).saturating_add(1));
        }
        if let Some(svlen) = self.svlen {
            return Some(svlen.unsigned_abs() as usize);
        }
        None
    }
}

#[derive(Debug, PartialEq, Eq, Clone, Hash)]
pub enum AlternateType {
    Literal(Vec<Nucleotide>),
    Symbolic(SvData),
}

impl AlternateType {
    /// Returns the literal base sequence if this is a `Literal` variant,
    /// or `None` for symbolic / structural ALTs.
    pub fn as_literal(&self) -> Option<&[Nucleotide]> {
        match self {
            AlternateType::Literal(vector) => Some(vector),
            AlternateType::Symbolic(_) => None,
        }
    }

    /// Returns the parsed symbolic payload if this is a `Symbolic` variant,
    /// or `None` for literal-base ALTs.
    pub fn as_symbolic(&self) -> Option<&SvData> {
        match self {
            AlternateType::Symbolic(sv) => Some(sv),
            AlternateType::Literal(_) => None,
        }
    }

    pub fn is_literal(&self) -> bool {
        matches!(self, AlternateType::Literal(_))
    }

    pub fn is_symbolic(&self) -> bool {
        matches!(self, AlternateType::Symbolic(_))
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
    // The alternate allele of interest. This is either one base or several bases for an indel, or
    // it is a key string from the vcf spec indicating a complex variant.
    pub alternate: AlternateType,
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
                return Err(VariantError::ComplexNotSupported);
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
    // like `G]17:198982]`) are valid per VCF 4.2 but not yet wired through
    // AlternateType::Symbolic — downstream code assumes literal bases. Surface
    // them as an InvalidVcf rather than silently coercing.
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
    fn alternate_type_literal_as_literal_returns_slice() {
        let alt = AlternateType::Literal(vec![Nucleotide::A, Nucleotide::C]);
        assert_eq!(alt.as_literal(), Some(&[Nucleotide::A, Nucleotide::C][..]));
    }

    #[test]
    fn alternate_type_symbolic_as_literal_returns_none() {
        let alt = AlternateType::Symbolic(SvData::new("<DUP:TANDEM>", SvType::Dup));
        assert!(alt.as_literal().is_none());
        // The dispatcher round-trip: symbolic ALTs go through as_symbolic().
        assert!(alt.is_symbolic());
        assert_eq!(alt.as_symbolic().unwrap().raw_alt, "<DUP:TANDEM>");
        assert_eq!(alt.as_symbolic().unwrap().sv_type, SvType::Dup);
    }

    #[test]
    fn alternate_type_literal_helpers() {
        let alt = AlternateType::Literal(vec![Nucleotide::A]);
        assert!(alt.is_literal());
        assert!(!alt.is_symbolic());
        assert!(alt.as_symbolic().is_none());
    }

    #[test]
    fn sv_type_from_alt_string_recognizes_standard_tags() {
        assert_eq!(SvType::from_alt_string("<DEL>"), Some(SvType::Del));
        assert_eq!(SvType::from_alt_string("<DEL:ME>"), Some(SvType::Del));
        assert_eq!(SvType::from_alt_string("<DUP>"), Some(SvType::Dup));
        assert_eq!(SvType::from_alt_string("<DUP:TANDEM>"), Some(SvType::Dup));
        assert_eq!(SvType::from_alt_string("<CNV>"), Some(SvType::Cnv));
        assert_eq!(SvType::from_alt_string("<INS>"), Some(SvType::Ins));
        assert_eq!(SvType::from_alt_string("<INS:ME:ALU>"), Some(SvType::Ins));
        assert_eq!(SvType::from_alt_string("<INV>"), Some(SvType::Inv));
    }

    #[test]
    fn sv_type_from_alt_string_breakend() {
        // VCF 4.2 breakend forms — any of these should be tagged Bnd.
        assert_eq!(
            SvType::from_alt_string("G]17:198982]"),
            Some(SvType::Bnd)
        );
        assert_eq!(
            SvType::from_alt_string("]13:123456]T"),
            Some(SvType::Bnd)
        );
        assert_eq!(
            SvType::from_alt_string("[2:321682[A"),
            Some(SvType::Bnd)
        );
    }

    #[test]
    fn sv_type_from_alt_string_unknown_tag() {
        // Angle-bracket form we don't have a dedicated variant for —
        // still classified as symbolic (Unknown), not literal.
        assert_eq!(
            SvType::from_alt_string("<NOVEL_THING>"),
            Some(SvType::Unknown)
        );
    }

    #[test]
    fn sv_type_from_alt_string_literal_returns_none() {
        // Literal-base ALTs shouldn't pass through this parser at all;
        // callers gate it on AlternateType, but defend the contract.
        assert_eq!(SvType::from_alt_string("A"), None);
        assert_eq!(SvType::from_alt_string("ACGT"), None);
        assert_eq!(SvType::from_alt_string(""), None);
    }

    #[test]
    fn sv_data_span_from_end_field() {
        // END is 1-based inclusive — for a deletion at POS=100, END=109
        // covers 10 bases (100..=109).
        let mut sv = SvData::new("<DEL>", SvType::Del);
        sv.end = Some(109);
        assert_eq!(sv.span(100), Some(10));
    }

    #[test]
    fn sv_data_span_falls_back_to_svlen() {
        // SVLEN is signed per VCF 4.2: negative for deletions, positive
        // for insertions/duplications. span() takes the absolute value.
        let mut sv = SvData::new("<DEL>", SvType::Del);
        sv.svlen = Some(-50);
        assert_eq!(sv.span(100), Some(50));

        let mut sv = SvData::new("<DUP>", SvType::Dup);
        sv.svlen = Some(200);
        assert_eq!(sv.span(100), Some(200));
    }

    #[test]
    fn sv_data_span_prefers_end_over_svlen() {
        // If both are populated and disagree, END wins (it's the more
        // authoritative coordinate in VCF 4.2).
        let mut sv = SvData::new("<DEL>", SvType::Del);
        sv.end = Some(109);
        sv.svlen = Some(-999);
        assert_eq!(sv.span(100), Some(10));
    }

    #[test]
    fn sv_data_span_none_when_neither_set() {
        let sv = SvData::new("<DEL>", SvType::Del);
        assert!(sv.span(100).is_none());
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

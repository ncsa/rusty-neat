//! This higher level struct will combine variant data and struct data, and will keep track
//! of which sequence blocks are affected
use crate::structs::fasta_map::{Contig, FastaMapError};
use crate::structs::variants::{Variant, VariantError};
use std::collections::HashMap;
use log::debug;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum ContigVariantError {
    #[error("ContigVariant reported a FastaMap error: {0}")]
    ContigError(FastaMapError),
    #[error("ContigVariant reported a Variant error: {0}")]
    VariantError(VariantError),
    #[error("ContiVariant reported a mapping error")]
    MappingError,
}

#[derive(Debug)]
pub struct ContigVariants {
    // This struct will keep track of all variants on a contig.
    //
    // The contig containing the variants
    contig: Contig,
    // The list of variants on the contig
    variants_list: Vec<Variant>,
    block_map: Option<HashMap<Variant, Vec<String>>>,
}

impl ContigVariants {
    pub fn new_unmapped(
        contig: Contig,
        variants_list: Vec<Variant>,
    ) -> Result<Self, ContigVariantError> {
        // This creates a basic contig variants box without the mapping of variants to sequence blocks
        Ok(ContigVariants { contig, variants_list, block_map: None })
    }

    pub fn new_mapped(
        contig: Contig,
        variants_list: Vec<Variant>,
    ) -> Result<Self, ContigVariantError> {
        // This will create a HashMap mapping variants to all the sequence blocks they affect
        let mut variant_map: HashMap<Variant, Vec<String>> = HashMap::new();
        for variant in variants_list.clone().iter() {
            let request_start = variant.location;
            // Guaranteed to work because a reference vector must have at least 1 element
            let request_end = request_start + variant.reference.len();
            let result = contig.get_sequence_block(request_start, request_end);
            match result {
                Ok(result) => { variant_map.insert(variant.clone(), result.clone()); },
                // This is a little weird, but we'll debug later
                Err(e) => { debug!("{}", e) }
            }
        }
        if variant_map.is_empty() {
            return Err(ContigVariantError::MappingError)
        }
        Ok(ContigVariants { 
            contig: contig,
            variants_list: variants_list.to_owned(),
            block_map: Some(variant_map).to_owned(),
        })
    }
}

#[cfg(test)]
mod tests {
    use crate::structs::contig_variants::ContigVariants;
    use crate::structs::fasta_map::{Contig, RegionType};
    use crate::structs::variants::{Variant, VariantType};

    #[test]
    fn test_new_unmapped() {
        // create contig
        let contig_name = "chr1".to_string();
        let len: usize = 1000;
        let blocks: Vec<String> = vec!["chr1_001000_002000.json".to_string()];
        let map: Vec<(usize, usize, RegionType)> = Vec::new();
        let contig = Contig::new(contig_name.to_owned(), len, blocks.to_owned(), map.to_owned()).unwrap();
        // create a variant
        let variant_type = VariantType::SNP;
        let location = 55;
        let reference: Vec<u8> = vec![1];
        let alternate: Vec<u8> = vec![3];
        let mut genotype: Vec<u8> = vec![1,0];
        let variant = Variant::new(variant_type, location, &reference, &alternate, &mut genotype).unwrap();
        let variants_list = vec![variant];
        // Now create the contig variant
        let cv = ContigVariants::new_unmapped(contig, variants_list.to_owned()).unwrap();
        assert_eq!(cv.contig.name, contig_name);
        assert_eq!(cv.variants_list, variants_list);
        assert!(cv.block_map.is_none());
    }
}
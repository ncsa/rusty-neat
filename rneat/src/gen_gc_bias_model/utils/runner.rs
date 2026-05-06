use std::{collections::HashMap, path::PathBuf};
use itertools::Itertools;
use log::*;

use common::{
    models::{mutation_model::MutationModel, snp_trinuc_model::{TrinucFrame}},
    structs::{
        bed_record::BedRecord,
        fasta_map::{FastaMap},
        nucleotides::{Nucleotide, ALLOWED_NUCS},
        transition_matrix::TransitionMatrix,
        variants::{Genotype, Variant, VariantType}
    }
};
use crate::gen_gc_bias_model::errors::GenGcBiasModelError;
use crate::gen_mut_model::errors::GenMutationModelError;


pub fn runner(
    path: &PathBuf
) -> Result<(), GenGcBiasModelError> {
    Ok(())
}
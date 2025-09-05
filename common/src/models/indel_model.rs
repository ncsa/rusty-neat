//! The following section are the models for each type of variant. In order to create the variant,
//! we need to model its statistical property. NEAT included two types of variants: SNPs and Indels.
//! Note that indels are actually two types, insertions and deletions, but they are usually classed
//! together in the literature and are usually short. Insertions are often slips that cause sections
//! to be duplicated, but NEAT made no attempt to distinguish between the types of insertions or
//! deletions, since they are treated similarly by variant calling software.
//!
//! We have the variants separated out in the variants struct, into Insertion and Deletion, but they
//! will share this model to avoid duplicating code, since code-wise they are closely linked.
use serde::{Deserialize, Serialize};
use simple_rng::{NeatRng, NeatRngError};
use thiserror::Error;
use std::{io, path::PathBuf};
use flate2::read::GzDecoder;
use crate::{models::lib::{model_reader, model_writer}, structs::{
    distributions::{DiscreteDistribution, DistributionErrors}, 
    nucleotides::{allowed_vec, Nucleotide}
}};

#[derive(Debug, Error)]
pub enum IndelModelError {
    #[error("Indel Model returned an Rng Error: {0}")]
    RngError(#[from] NeatRngError),
    #[error("Indel model returnd an error while initializng distributions.")]
    DistributionInitializationError(#[from] DistributionErrors),
    #[error("Error during input/output: {0}")]
    IoError(#[from] io::Error),
    #[error("Serde error building default model: {0}")]
    SerdeError(#[from] serde_json::Error),
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IndelModel {
    // Based what was in the original NEAT
    pub (crate) insertion_probability: f64,
    ins_dist: DiscreteDistribution,
    del_dist: DiscreteDistribution,
}

static DATA_FILE: &'static [u8] = include_bytes!("model_data/default_indel_model.json.gz");

impl IndelModel {
    pub fn default() -> Result<Self, IndelModelError> {
        // Default Indel model from the original NEAT, scaled by the lowest value and rounded.
        let reader = GzDecoder::new(DATA_FILE);
        let data: IndelModel = serde_json::from_reader(reader)
            .map_err(IndelModelError::SerdeError)?;
        Ok(data)
    }

    pub fn from(filename: &PathBuf) -> Result<Self, IndelModelError> {
        let data: IndelModel = model_reader(filename)?;
        Ok(data)
    }

    pub fn write_to_file(&self, filename: &PathBuf) -> Result<(), IndelModelError> {
        model_writer(self, filename)?;
        Ok(())
    }

    pub fn is_insertion(&self, rand: f64) -> Result<bool, IndelModelError> {
        Ok(rand < self.insertion_probability)
    }

    pub fn new_insert_length(&self, rand: f64) -> Result<usize, IndelModelError> {
        // This function uses the insertion weights to choose an insertion length from the list.
        Ok(self.ins_dist.sample(rand)?)
    }

    pub fn new_delete_length(&self, rand: f64) -> Result<usize, IndelModelError> {
        // This function is the same as above, but uses the deletion lengths instead.
        Ok(self.del_dist.sample(rand)?)
    }
    
    pub fn generate_random_insertion(
        &self, 
        length: usize, 
        rng: &mut NeatRng
    ) -> Result<Vec<Nucleotide>, IndelModelError> {
        // We could refine this with a nucleotide bias matrix. Maybe it would make a difference,
        // but probably not, since the presence of the insertion is more important than it's content,
        // for this use. If there was a call for it, maybe.
        let mut insertion_vec = Vec::new();
        let allowed_nucs = allowed_vec();
        for _ in 0..length {
            // since our range is restricted, we are covered
            let nuc = rng.choose(&allowed_nucs).unwrap();
            insertion_vec.push(nuc)
        }
        Ok(insertion_vec)
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::structs::nucleotides::Nucleotide::*;

    #[test]
    fn test_indel_model() {
        let mut rng = NeatRng::new_from_seed(&vec![
            "Hello".to_string(),
            "Cruel".to_string(),
            "World".to_string(),
        ]).unwrap();

        let indel_model = IndelModel {
            insertion_probability: 0.75,
            ins_dist: DiscreteDistribution::new(&vec![100.0, 1.0], &vec![1, 2]).unwrap(),
            del_dist: DiscreteDistribution::new(&vec![1.0, 1000.0], &vec![3, 4]).unwrap(),
        };

        let length = indel_model.new_insert_length(rng.random().unwrap()).unwrap();
        if length > 0 {
            assert!((length == 1) || (length == 2));
        }
    }

    #[test]
    fn test_generate_insertion () {
        let mut rng = NeatRng::new_from_seed(&vec![
            "Hello".to_string(),
            "Cruel".to_string(),
            "World".to_string(),
        ]).unwrap();

        let indel_model = IndelModel {
            insertion_probability: 0.75,
            ins_dist: DiscreteDistribution::new(&vec![100.0, 1.0], &vec![1, 2]).unwrap(),
            del_dist: DiscreteDistribution::new(&vec![1.0, 1000.0], &vec![3, 4]).unwrap(),
        };

        let insertion = indel_model.generate_random_insertion(10, &mut rng).unwrap();
        assert_eq!(insertion.len(), 10);
        assert!(!insertion.contains(&Nucleotide::A));

        let insertion2 = indel_model.generate_random_insertion(12, &mut rng).unwrap();
        assert_eq!(insertion2.len(), 12);
        assert_eq!(insertion2, Vec::from([G, T, T, G, G, G, T, C, T, C, T, T]));
    }
}

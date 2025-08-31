//! In DNA sequencing, a fragment is a bit of DNA, roughly uniform in length that is sequenced
//! by the machine. Sometimes these fragments have special molecules attached to the end for ID
//! purposes. How this is done is a process called Chemistry Magic. For our purposes, we exect
//! the data to be uniform enough that a mean and standard deviation will describe the set.

use simple_rng::NeatRngError;
use thiserror::Error;
use std::io;
use std::path::Path;
use serde_json;
use serde::{Deserialize, Serialize};
use crate::structs::distributions::{DiscreteDistribution, DistributionErrors, NormalDistribution};
use crate::models::lib::{model_reader, model_writer};

#[derive(Error, Debug)]
pub enum FragmentModelError {
    #[error("Fragment model returned an RNG error: {0}")]
    RngError(NeatRngError),
    #[error("Fragment model returned an IO error: {0}")]
    IoError(io::Error),
    #[error("Fragment Model attempted to load a file that it could not find: {0}")]
    FileNotFound(String),
    #[error("Fragment model reported a distribution initiation error: {0}")]
    DistributionInitError(DistributionErrors)
}

impl From<DistributionErrors> for FragmentModelError {
    fn from(error: DistributionErrors) -> Self {
        FragmentModelError::DistributionInitError(error)
    }
}

impl From<NeatRngError> for FragmentModelError {
    fn from(error: NeatRngError) -> Self {
        FragmentModelError::RngError(error)
    }
}

impl From<io::Error> for FragmentModelError {
    fn from(error: io::Error) -> Self {
        FragmentModelError::IoError(error)
    }
}

#[derive(Debug)]
pub enum FragmentLengthModel {
    Discrete(DiscreteFragmentLengthModel),
    Normal(NormalFragmentLengthModel),
}

#[derive(Debug, Serialize, Deserialize)]
pub struct DiscreteFragmentLengthModel {
    // The mean length of a fragment for this simulation
    // fragment_mean: f64, // maybe we don't want to store this
    // The standard deviation of a fragment for this simulation
    // fragment_std: f64, // maybe we don't want to store this
    // The Normal Distribution for this model
    distrbution: DiscreteDistribution,
}

impl DiscreteFragmentLengthModel {
    pub fn new(lengths: Vec<usize>, weights: Vec<f64>) -> Result<Self, FragmentModelError> {
        // These were numbers routinely used for testing in NEAT genReads
        let fragment_dist = DiscreteDistribution::new_from_values(&weights, &lengths)?;
        Ok(DiscreteFragmentLengthModel { 
            distrbution: fragment_dist
        })
    }

    pub fn discrete_from_file(filename: &str) -> Result<Self, FragmentModelError> {
        // The baseline model is really just a mathematical equation and can be reconstructed by the two input parameters
        // But reading from data, it may be better to use the discrete distribution to maintain outliers. This will load such
        // a model from a json file. The file must be of the format:
        // {
        //   "fragment_lengths": [ ... ],
        //   "fragment_weights": [ ... ]
        // }
        // Where fragment lengths are of type usize and fragment weights are of (or can be cast as) type f64.
        let path = Path::new(&filename);
        match path.exists() {
            false => return Err(FragmentModelError::FileNotFound(filename.to_string())),
            _ => {},
        }

        let data: DiscreteFragmentLengthModel = model_reader(&filename).unwrap();

        Ok(data)
    }

    pub fn write_file(&self, filename: &str) -> Result<(), FragmentModelError> {
        // serialize a model with serde and write it to file
        let data = serde_json::to_string(self).unwrap();
        model_writer(data, filename)?;
        Ok(())
    }

    pub fn generate_fragment(&self, rand: f64) -> Result<usize, FragmentModelError> {
        // This function generates a fragment length based on mean and standard deviation.
        // This is basically a random number generator.
        Ok(self.distrbution.sample_index(rand)? as usize)
    }
}

#[derive(Debug)]
pub struct NormalFragmentLengthModel {
    // The mean length of a fragment for this simulation
    // fragment_mean: f64, // maybe we don't want to store this
    // The standard deviation of a fragment for this simulation
    // fragment_std: f64, // maybe we don't want to store this
    // The Normal Distribution for this model
    distrbution: NormalDistribution,
}

impl NormalFragmentLengthModel {
    pub fn default() -> Result<Self, FragmentModelError> {
        // These were numbers routinely used for testing in NEAT genReads
        let fragment_dist = NormalDistribution::new(300.0, 30.0)?.into();
        Ok(NormalFragmentLengthModel { 
            distrbution: fragment_dist,
        })
    }

    pub fn new_from_mean(fragment_mean: f64, fragment_std: f64) -> Result<Self, FragmentModelError> {
        let fragment_dist = NormalDistribution::new(
            fragment_mean,
            fragment_std,
        )?.into();
        Ok(NormalFragmentLengthModel {
            distrbution: fragment_dist,
        })
    }

    pub fn params(&self) -> Result<(f64, f64), FragmentModelError> {
        // This returns the parameters used to initiate the model
        Ok(self.distrbution.params()?)
    }

    pub fn generate_fragment(&self, rand: f64) -> Result<usize, FragmentModelError> {
        // This function generates a fragment length based on mean and standard deviation.
        // This is basically a random number generator.
        Ok(self.distrbution.sample(rand)?.trunc() as usize)
    }
}

#[cfg(test)]
mod tests{
    use super::*;

    #[test]
    fn test_default() {
        let model = NormalFragmentLengthModel::default().unwrap();
        assert_eq!(model.params().unwrap().0, 300.0);
        assert_eq!(model.params().unwrap().1, 30.0);
    }

    #[test]
    fn test_new_from_mean() {
        let mean = 34.33;
        let std_dev = 1.232;
        let model = NormalFragmentLengthModel::new_from_mean(
            mean, 
            std_dev
        ).unwrap();
        assert_eq!(model.distrbution.params().unwrap().0, 34.33);
        assert_eq!(model.distrbution.params().unwrap().1, 1.232);
    }

    #[test]
    fn test_new_discrete() {
        let l_vec = vec![1, 8, 9, 10];
        let w_vec = vec![1.0, 3.0, 2.0, 1.2];
        let model = DiscreteFragmentLengthModel::new(
            l_vec.clone(),
            w_vec.clone(),
        ).unwrap();
        assert_eq!(model.distrbution.values().unwrap(), l_vec);
        assert_eq!(model.distrbution.weights().unwrap(), w_vec);
    }

    #[test]
    fn test_generate_fragment() {
        todo!()
    }

    #[test]
    fn test_from_file() {
        todo!()
    }

    #[test]
    fn test_write_file() {
        todo!()
    }
}
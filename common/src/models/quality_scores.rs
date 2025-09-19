//! Walking through Zac Stephens original algorithm, to try to make sure I replicate it correctly.
//!   * For position 1, there is a vector of weights for each score, extracted from data.
//!   * For each position in the read length after that
//!         * For each possible quality score, a distribution is constructed with weights and
//!            scores, as determined by a matrix of weights
//!   * For read length N and # possible quality scores Q, this creates a vector with length N
//!         * first element is a 1-D vector of weights with length Q
//!         * each subsequent element is a vector of length Q,
//!           each element of which is a vector of length N.
//! To generate quality scores, they follow the following procedure:
//!   * Sample the first element (1D vector) for initial quality score.
//!   * For next position, the previous Q score determines which N-length set of weights to use to
//!     determine the next quality score
//! Advantages of this approach:
//!   * Does a fairly effective job of modeling shapes of the quality scores for a set read length
//! Disadvantages of this approach:
//!   * The fact that we're working with a matrix with a different first element is
//!     extremely confusing
//!   * In addition to the difficulty keeping track of indexes, I'm not sure how well this will
//!     translate to Rust. May need a custom data structure. Like seed + subsequent.
//!   * Assumes a fixed read length, meaning you have to extrapolate for longer read lengths.
//!   * In Python, at least, this was slow, although in retrospect it didn't eat up much memory.
use serde_json;
use std::{io, path::PathBuf};
use thiserror::Error;
use flate2::read::GzDecoder;
use serde::{Deserialize, Serialize};
use std::{fmt::{Display, Formatter}};
use simple_rng::{NeatRng, NeatRngError};

use crate::models::lib::{model_reader, model_writer};
use crate::models::sequencing_error_model::SeqModelError;
use crate::structs::distributions::{DiscreteDistribution, DistributionErrors};

pub const QUALITY_OFFSET: usize = 33;

#[derive(Debug, Error)]
pub enum QualityModelError {
    #[error("Quality model initiation returned a distribution error: {0}")]
    DistributionError(DistributionErrors),
    #[error("Quality score creation reported an RNG Error: {0}")]
    RngError(NeatRngError),
    #[error("Quality model return an IO error: {0}")]
    IoError(io::Error),
    #[error("Serde error building default model: {0}")]
    SerdeError(serde_json::Error),
}

impl From<serde_json::Error> for QualityModelError {
    fn from(error: serde_json::Error) -> Self {
        QualityModelError::SerdeError(error)
    }
}

impl From<DistributionErrors> for QualityModelError {
    fn from(error: DistributionErrors) -> Self {
        QualityModelError::DistributionError(error)
    }
}

impl From<NeatRngError> for QualityModelError {
    fn from(error: NeatRngError) -> Self {
        QualityModelError::RngError(error)
    }
}

impl From<std::io::Error> for QualityModelError {
    fn from(error: std::io::Error) -> Self {
        QualityModelError::IoError(error)
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct QualityScoreModel {
    // This is the vector of the quality scores possible in this dataset. This could be a list
    // of numbers from 1-42, for example, or bins of scores, [2, 13, 27, 33] or whatever the
    // dataset uses. This list is expected to be sorted.
    pub quality_score_options: Vec<usize>,
    // True for binned scores, false for continuous
    pub binned_scores: bool,
    // The assumed read length of this dataset. The model will assume this read length and adjust
    // on a per-run basis in a deterministic way (doubling positional weight arrays)
    pub assumed_read_length: usize,
    // Weights for the first position in the read length.
    pub seed_dist: DiscreteDistribution,
    // A matrix for each subsequent position along the read length after the first. Each row is a
    // discrete distribution, keyed the previous score. For example, for possible scores 0-41,
    // inclusive, there would be 42 vectors (one for each possible previous score), each giving the
    // distribution for the current position (one weight for each of 42 scores). This is based on
    // the original design in NEAT. Previous attempts to simplify this have not been able to
    // successfully reproduce quality scores.
    pub distros_from_one: Vec<Vec<DiscreteDistribution>>,
}

impl Display for QualityScoreModel {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        // Just a basic display showing the possible quality scores and read length of the model.
        // This won't be used by generate_reads, but could be used by a generate quality score model
        write!(
            f,
            "QualityScoreModel: (rl: {})\n\
            \tscores: {:?}\n\
            \tbinned? {:?}\n",
            self.assumed_read_length, self.quality_score_options, self.binned_scores,
        )
    }
}

static DATA_FILE: &'static [u8] = include_bytes!("model_data/default_quality_score_model.json.gz");

impl QualityScoreModel {
    // methods for QualityScoreModel objects

    pub fn default() -> Result<Self, QualityModelError> {
        // This generates the default quality score model based on the original NEAT default. 
        // the parameters are sort of outdated now:
        //  - quality_score_options: 0 - 41
        //  - binned_scores: false
        //  - assumed_read_length: 101
        // I won't list all the NEAT calculated weights. zcat the model file if you are interested, but 
        // it's a lot of data.
        
        // this map_err thing I had help on. I feel like this is something rust-analyzer should have known.
        let reader = GzDecoder::new(DATA_FILE);
        let data: QualityScoreModel = serde_json::from_reader(reader)
            .map_err(QualityModelError::SerdeError)?;
        Ok(data)
    }

    /// we will write subutilities that use these features, eventually
    pub fn from_file(filename: &PathBuf) -> Result<Self, QualityModelError> {
        // Uses the serde_json crate to read a quality model from file
        let data: QualityScoreModel = model_reader(filename).unwrap();
        Ok(data)
    }

    #[allow(dead_code)]
    // need a test
    pub fn display(&self) -> String {
        // Some detailed formats. I thing these will be useful for quality model generation debugging.
        format!(
            "QualityScoreModel: (rl: {})\n\
            \tscores: {:?}\n\
            \tbinned? {:?}\n\
            \tseed distribution: {:?}\n\
            \tfirst weight array: {:?}",
            self.assumed_read_length,
            self.quality_score_options,
            self.binned_scores,
            self.seed_dist,
            self.distros_from_one[1][0],
        )
    }

    #[allow(dead_code)]
    // need a test
    pub fn display_it_all(&self) -> String {
        format!(
            "QualityScoreModel: (rl: {})\n\
            \tscores: {:?}\n\
            \tbinned? {:?}\n\
            \tseed distribution: {:?}\n\
            \tscore weight array: {:?}",
            self.assumed_read_length,
            self.quality_score_options,
            self.binned_scores,
            self.seed_dist,
            self.distros_from_one,
        )
    }

    pub fn generate_quality_scores(
        &self,
        length: usize,
        rng: &mut NeatRng,
    ) -> Result<Vec<usize>, SeqModelError> {
        // Generates a list of quality scores of length run_read_length using the model. If the
        // input read length differs, we do some index magic to extrapolate the model
        // run_read_length: The desired read length for the model to generate.

        // This will be the list of scores generated
        let mut score_list: Vec<usize> = Vec::with_capacity(length);
        // sample the scores list with the seed weights applied to generate the first score.
        // Samples an index based on the weights, which then selects the quality score.
        let mut seed_score = self.seed_dist.sample(rng.random()?)?;
        // We don't want the first score to be 31. That throws off the parser, because it makes an @ symbol
        if seed_score == 31 {
            seed_score -= 1
        }
        // Adding the seed score to the list. This is safe so long as the qual score max is <= 255 (currently 40)
        score_list.push(seed_score);
        // To map from one length to another, we use the algorithm found in the original NEAT 2.0,
        // adapted to rust. See function for implementation details.
        let indexes: Vec<usize> = self.quality_index_remap(length);
        // Sort of annoying, but to account for the remap, in order to get "previous score" from
        // the score list, we need to know the current index we are filling, absolutely, in cases
        // of mismatches between model read length and run read length
        // We can skip the first one, since we already generated it above. On loop 1, we will look
        // at that seed score to get our first set of weights.
        let mut current_index = 1;
        for i in indexes {
            // The weight at this index is the score weights for position i, given a
            // previous score of score_list[i-1]

            // First get the index of the previous score from the original scores list.
            // This will match the index in the score weights table that corresponds to that score.
            let score_position = self
                .quality_score_options
                .iter()
                .position(|&x| x == score_list[current_index - 1])
                .unwrap();
            // Now we have an index (in the default case 0..<4) of a vector for the position, based
            // on the previous score. We have to subtract one from the index because the first matrix
            // (position 0) corresponds to the second quality score (position 1).
            // We need to figure out a better way if we're going to do discrete.
            let index = self.distros_from_one[i-1][score_position]
                .sample(rng.random()?)?;
            let score = self.quality_score_options[index];
            score_list.push(score); 
            current_index += 1;
        }
        Ok(score_list)
    }
    
    fn quality_index_remap(&self, length: usize) -> Vec<usize> {
        // Basically, this function does integer division (truncation) to fill positions
        // in a vector the length of the desired read length.
        // for example. You are mapping from read length 6 to read length 8,
        // A: [0, 1, 2, 3, 4, 5] -> B: [0, 1, 2, 3, 4, 5, 6, 7]
        // The question is, since our model (A) only has 6 columns of values, and for each of the 8
        // items in the desire output (B), we need to know which of the 5 columns from A to use.
        // So we create a vector of length 8, which tells us, at each position, which column from
        // A to use. For example, at position 1, we take (6 * 1) // 8 (where `//` denotes integer
        // division), resulting in 6//8 = 0, so for the first quality score, we select based on
        // The 0 column from the A model. At position 2, (6 * 2) // 8 = 12//8 = 1. Repeating, this
        // for each i from 0 to read length, we end up with: C: [0, 0, 1, 2, 3, 3, 4, 5]
        // Similarly, we can remap this way from larger to smaller:
        // A: [0, 1, 2, 3, 4, 5, 6, 7] -> B: [0, 1, 2, 3, 4, 5]
        // C = [(8*0)//6 = 0, (8*1)//6 = 1, (8*2)//6 = 2, (8*3)//6 = 4, (8*4)//6 = 5]
        //   = [0, 1, 2, 4, 5]
        // Advantages: should be pretty quick. Easy calculations.
        // Disadvantages: Tends to lose info from the back of the read when downsizing. Might need
        //                to check that.
        if length == self.assumed_read_length as usize {
            (1..length).collect()
        } else {
            let mut indexes: Vec<usize> = Vec::new();
            for i in 1..length {
                let index: usize = (self.assumed_read_length as usize * i) / length;
                // This first value(s) will always be zero when run_read_length is longer than
                // assumed read length.
                if index < 1 {
                    indexes.push(1);
                } else if index >= (self.assumed_read_length - 1) {
                    indexes.push(self.assumed_read_length - 1)
                } else {
                    indexes.push(index);
                }
            }
            indexes
        }
    }

    #[allow(unused)]
    /// we will write subutilities that use these features, eventually
    fn write_to_file(&self, filename: &PathBuf) -> std::io::Result<()> {
        // Uses the serde_json crate to write out the json form of the model. This will help us
        // create base datasets from old neat data, and give us a way to write out models that are
        // generated from user data.
        model_writer(self, filename)?;
        Ok(())
    }

}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_model_write_read() {
        // rewrite this test to read default model, so we aren't keeping long.json around
        let temp_dir = tempfile::tempdir().unwrap();
        let mut temp_file = PathBuf::from(temp_dir.path());
        temp_file.push("test.json.gz");
        let model: QualityScoreModel = QualityScoreModel::default().unwrap();
        assert_eq!(model.assumed_read_length, 101);
        let result = model.write_to_file(&temp_file);
        assert_eq!(result.unwrap(), ());
        temp_dir.close().unwrap();
    }

    #[test]
    fn test_model_default() {
        let model = QualityScoreModel::default().unwrap();
        assert_eq!(model.assumed_read_length, 101)
    }

    #[test]
    fn test_quality_scores_short() {
        let run_read_length = 100;
        let mut rng = NeatRng::new_from_seed(&vec![
            "Hello".to_string(),
            "Cruel".to_string(),
            "World".to_string(),
        ]).unwrap();
        let model = QualityScoreModel::default().unwrap();
        let scores = model.generate_quality_scores(run_read_length, &mut rng).unwrap();
        assert!(!scores.is_empty());
        assert_eq!(scores.len(), 100);
        scores
            .iter()
            .map(|x| assert!(model.quality_score_options.contains(x)))
            .collect()
    }

    #[test]
    fn test_quality_scores_same() {
        let run_read_length = 150;
        let mut rng = NeatRng::new_from_seed(&vec![
            "Hello".to_string(),
            "Cruel".to_string(),
            "World".to_string(),
        ]).unwrap();
        let model = QualityScoreModel::default().unwrap();
        let scores = model.generate_quality_scores(run_read_length, &mut rng).unwrap();
        assert!(!scores.is_empty());
        assert_eq!(scores.len(), 150);
        scores
            .iter()
            .map(|x| assert!(model.quality_score_options.contains(x)))
            .collect()
    }

    #[test]
    fn test_quality_scores_long() {
        let run_read_length = 200;
        let mut rng = NeatRng::new_from_seed(&vec![
            "Hello".to_string(),
            "Cruel".to_string(),
            "World".to_string(),
        ]).unwrap();
        let model = QualityScoreModel::default().unwrap();
        let scores = model.generate_quality_scores(run_read_length, &mut rng).unwrap();
        assert!(!scores.is_empty());
        assert_eq!(scores.len(), 200);
        scores
            .iter()
            .map(|x| assert!(model.quality_score_options.contains(x)))
            .collect()
    }

    #[test]
    fn test_quality_scores_vast_difference() {
        let run_read_length = 2000;
        let mut rng = NeatRng::new_from_seed(&vec![
            "Hello".to_string(),
            "Cruel".to_string(),
            "World".to_string(),
        ]).unwrap();
        let model = QualityScoreModel::default().unwrap();
        let scores = model.generate_quality_scores(run_read_length, &mut rng).unwrap();
        assert!(!scores.is_empty());
        assert_eq!(scores.len(), 2000);
        scores
            .iter()
            .map(|x| assert!(model.quality_score_options.contains(x)))
            .collect()
    }
}

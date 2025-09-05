//! In DNA sequencing, a fragment is a bit of DNA, roughly uniform in length that is sequenced
//! by the machine. Sometimes these fragments have special molecules attached to the end for ID
//! purposes. How this is done is a process called Chemistry Magic. For our purposes, we exect
//! the data to be uniform enough that a mean and standard deviation will describe the set.

use flate2::read::GzDecoder;
use simple_rng::NeatRngError;
use thiserror::Error;
use std::io;
use std::path::{Path, PathBuf};
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
    DistributionInitError(DistributionErrors),
    #[error("Error building default model!")]
    SerdeError(serde_json::Error)
}

impl From<serde_json::Error> for FragmentModelError {
    fn from(error: serde_json::Error) -> Self {
        FragmentModelError::SerdeError(error)
    }
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

impl From<DiscreteFragmentLengthModel> for FragmentLengthModel {
    fn from(model: DiscreteFragmentLengthModel) -> Self {
        FragmentLengthModel::Discrete(model)
    }
}

impl From<NormalFragmentLengthModel> for FragmentLengthModel {
    fn from(model: NormalFragmentLengthModel) -> Self {
        FragmentLengthModel::Normal(model)
    }
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

static DATA_FILE: &'static [u8] = include_bytes!("model_data/default_fragment_length_model.json.gz");

impl DiscreteFragmentLengthModel {
    pub fn new(lengths: Vec<usize>, weights: Vec<f64>) -> Result<Self, FragmentModelError> {
        // These were numbers routinely used for testing in NEAT genReads
        let fragment_dist = DiscreteDistribution::new(&weights, &lengths)?;
        Ok(DiscreteFragmentLengthModel { 
            distrbution: fragment_dist
        })
    }

    pub fn default() -> Result<Self, FragmentModelError> {
        // The parameters of the default model from the original neat 
        // The lengths range from 1 to 799, though it skips from 1 to 32 before counting up.
        // The weights are just numbers between 0 and 1.
        // These are data gathered from publicly availble human data, but should reflect
        // whatever chemistry was used at the time
        let reader = GzDecoder::new(DATA_FILE);
        let data: DiscreteFragmentLengthModel = serde_json::from_reader(reader)
            .map_err(FragmentModelError::SerdeError)?;
        Ok(data)
    }

    pub fn discrete_from_file(filename: &PathBuf) -> Result<Self, FragmentModelError> {
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
            false => return Err(FragmentModelError::FileNotFound(filename.display().to_string())),
            _ => {},
        }

        let data: DiscreteFragmentLengthModel = model_reader(&filename).unwrap();

        Ok(data)
    }

    pub fn write_file(&self, filename: &PathBuf) -> Result<(), FragmentModelError> {
        // serialize a model with serde and write it to file
        model_writer(self, filename)?;
        Ok(())
    }

    pub fn generate_fragment(&self, rand: f64) -> Result<usize, FragmentModelError> {
        // This function generates a fragment length based on mean and standard deviation.
        // This is basically a random number generator.
        Ok(self.distrbution.sample(rand)? as usize)
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
    use std::fs;

    #[test]
    fn test_reader_writer() {
        let input_file = PathBuf::from("/home/joshfactorial/code/rusty-neat/common/src/models/model_data/default_fragment_length_model.json.gz");
        let output_file = PathBuf::from("/home/joshfactorial/code/rusty-neat/common/src/models/model_data/test.json.gz");
        let model: DiscreteFragmentLengthModel = model_reader(&input_file).unwrap();
        let l_vec = vec![1,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,255,256,257,258,259,260,261,262,263,264,265,266,267,268,269,270,271,272,273,274,275,276,277,278,279,280,281,282,283,284,285,286,287,288,289,290,291,292,293,294,295,296,297,298,299,300,301,302,303,304,305,306,307,308,309,310,311,312,313,314,315,316,317,318,319,320,321,322,323,324,325,326,327,328,329,330,331,332,333,334,335,336,337,338,339,340,341,342,343,344,345,346,347,348,349,350,351,352,353,354,355,356,357,358,359,360,361,362,363,364,365,366,367,368,369,370,371,372,373,374,375,376,377,378,379,380,381,382,383,384,385,386,387,388,389,390,391,392,393,394,395,396,397,398,399,400,401,402,403,404,405,406,407,408,409,410,411,412,413,414,415,416,417,418,419,420,421,422,423,424,425,426,427,428,429,430,431,432,433,434,435,436,437,438,439,440,441,442,443,444,445,446,447,448,449,450,451,452,453,454,455,456,457,458,459,460,461,462,463,464,465,466,467,468,469,470,471,472,473,474,475,476,477,478,479,480,481,482,483,484,485,486,487,488,489,490,491,492,493,494,495,496,497,498,499,500,501,502,503,504,505,506,507,508,509,510,511,512,513,514,515,516,517,518,519,520,521,522,523,524,525,526,527,528,529,530,531,532,533,534,535,536,537,538,539,540,541,542,543,544,545,546,547,548,549,550,551,552,553,554,555,556,557,558,559,560,561,562,563,564,565,566,567,568,569,570,571,572,573,574,575,576,577,578,579,580,581,582,583,584,585,586,587,588,589,590,591,592,593,594,595,596,597,598,599,600,601,602,603,604,605,606,607,608,609,610,611,612,613,614,615,616,617,618,619,620,621,622,623,624,625,626,627,628,629,630,631,632,633,634,635,636,637,638,639,640,641,642,643,644,645,646,647,648,649,650,651,652,653,654,655,656,657,658,659,660,661,662,663,664,665,666,667,668,669,670,671,672,673,674,675,676,677,678,679,680,681,682,683,684,685,686,687,688,689,690,691,692,693,694,695,696,697,698,699,700,701,702,703,704,705,706,707,708,709,710,711,712,713,714,715,716,717,718,719,720,721,722,723,724,725,726,727,728,729,730,731,732,733,734,735,736,737,738,739,740,741,742,743,744,745,746,747,748,749,750,751,752,753,754,755,756,757,758,759,760,761,762,763,764,765,766,767,768,769,770,771,772,773,774,775,776,777,778,779,780,781,782,783,784,785,786,787,788,789,790,791,793,795,796,797,799];
        assert_eq!(model.distrbution.values().unwrap(), l_vec);

        let result = model.write_file(&output_file);
        assert_eq!(result.unwrap(), ());
        fs::remove_file(output_file).unwrap();
    }

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
        assert_eq!(model.distrbution.weights().unwrap(), [0.1388888888888889, 0.5555555555555556, 0.8333333333333334, 1.0]);
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
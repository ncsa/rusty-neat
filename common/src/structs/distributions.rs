//! This DiscreteDistribution is an implementation of Zach Stephen's original neat-genReads code
//! from the py/probability.py file in tag 2.1 of github.com/ncsa/neat
//! (see also github.com/zstephens/neat-genreads). We may try the statrs Categorical distribution
//! as well, as I think it does the same thing.
use simple_rng::NeatRngError;
use log::debug;
use thiserror::Error;
use serde::{
    Deserialize, 
    Serialize
};
use statrs::{
    distribution::{ContinuousCDF, Normal},
    statistics::Distribution
};

#[derive(Debug, Error)]
pub enum DistributionErrors {
    #[error("A distribution returned a sampling error: {0}.")]
    SamplingError(&'static str),
    #[error("A distribution returned an Rng Error: {0}")]
    RngError(#[from] NeatRngError),
    #[error("Requested a value vector from an index-only model")]
    VectorNotFound,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DiscreteDistribution {
    pub(crate) values: Vec<usize>,
    pub(crate) weights: Vec<f64>,
}

impl DiscreteDistribution {
    pub fn new(w_vec: &Vec<f64>, v_vec: &Vec<usize>) -> Result<Self, DistributionErrors> {
        let cumulative_probability = {
            let sum_weights: f64 = w_vec.iter().sum();
            if sum_weights > 0.0 {
                let mut normalized_weights = Vec::with_capacity(w_vec.len());
                // we no longer need the w_vec_64 after this, so we consume it
                for weight in w_vec{
                    normalized_weights.push(weight / sum_weights);
                }
                cumulative_sum(&mut normalized_weights).unwrap()
            } else {
                // Edge case. Really, this shouldn't be allowed in the data. We'll fix it when 
                // we rebuild the models.
                vec![1.0_f64; w_vec.len()]
            }
        };

        Ok(DiscreteDistribution {
            values: v_vec.to_owned(),
            weights: cumulative_probability,
        })
    }

    pub fn values(&self) -> Result<Vec<usize>, DistributionErrors> {
        // fetches the value matrix
        Ok(self.values.clone())
    }

    pub fn weights(&self) -> Result<Vec<f64>, DistributionErrors> {
        // fetches the value matrix, if it exists
        Ok(self.weights.to_vec())
    }

    pub fn sample(&self, rand: f64) -> Result<usize, DistributionErrors> {
        // returns a random value for the distribution, based on cumulative probability
        // This is basically an icdf for a discrete distribution
        // We take a random number as an input to avoid copying the RNG around the program.
        let mut lo: usize = 0;
        // There was no minus 1 in the original, but I ran into an issue with this, for example,
        // lo = 47, hi = 48
        // mid = (lo + hi) / 2 = 47
        // lo = mid + 1 = 48
        // now we're out of bounds
        let mut hi: usize = self.weights.len();
        // bisect left
        while lo < hi {
            let mid = (lo + hi) / 2;
            if self.weights[mid] < rand {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        if lo == self.values.len() {
            // All right, the edge case in question is when the quaity array has no values, because the original quality array wasn't 100% filled
            // out. In the future, when we fix the quality arrays so that this doesn't happen, we will be able to remove this safety check, I think
            debug!("Errored on {:?}", self);
            return Ok(self.values[0].clone())
        }
        Ok(self.values[lo].clone())
    }
}

/// TODO: Figure out how to implement serialize and deserialize for normal. I think it should be easy 
/// enough to do, but I'm tired
#[derive(Debug, Clone)]
pub struct NormalDistribution {
    distribution: Normal,
}

impl NormalDistribution {
    // This is little more than a wrapper for the statrs normal distribution. Other distributions
    // may be more complicated, as needed. I would copy in the code for custom tailoring, but the
    // inverse function requires tons of hardcoded coefficient tables
    // and no one has time for all that
    pub fn new(mean: f64, std_dev: f64) -> Result<Self, DistributionErrors> {
        Ok(NormalDistribution {
            distribution: Normal::new(mean, std_dev).unwrap()
        })
    }

    pub fn params(&self) -> Result<(f64, f64), DistributionErrors> {
        // returns the initiating parameters of the underlying distribution
        Ok((self.distribution.mean().unwrap(), self.distribution.std_dev().unwrap()))
    }

    pub fn inverse_cdf(&self, x: f64) -> Result<f64, DistributionErrors> {
        Ok(self.distribution.inverse_cdf(x))
    }

    pub fn sample(&self, rand: f64) -> Result<f64, DistributionErrors> {
        // Takes a statrs NormalDistribution object, uses the ICDF to start with a probability (our
        // RNG, which generates numbers between 0 and 1 with approximately normal distribution) and
        // generate the corresponding Y value from the PDF. Very handy!
        Ok(self.distribution.inverse_cdf(rand))
    }
}

fn cumulative_sum(a: &mut Vec<f64>) -> Result<Vec<f64>, DistributionErrors> {
    // This calculation allows us to sample the dataset with a simple bisect,
    // which is very fast an convenient.
    let mut acc = 0.0;
    let mut cumvec = Vec::with_capacity(a.len());
    for x in a {
        acc += *x;
        cumvec.push(acc);
    };
    Ok(cumvec)
}

#[cfg(test)]
mod test {
    use simple_rng::NeatRng;
    use super::*;

    #[test]
    fn new_from_values() {
        let l_vec = vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
        let w_vec = vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        assert_eq!(l_vec.len(), w_vec.len());
        let distribution = DiscreteDistribution::new(&w_vec, &l_vec).unwrap();
        assert_eq!(distribution.values().unwrap(), l_vec);
        assert_eq!(distribution.weights().unwrap(), vec![1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]);
    }

    #[test]
    fn test_discrete_distribution() {
        let values = vec![1, 2, 3, 4, 5, 6];
        let weights: Vec<f64> = vec![1.1, 2.0, 1.0, 8.0, 0.2, 2.0];
        let d = DiscreteDistribution::new(&weights, &values).unwrap();
        let mut rng = NeatRng::new_from_seed(&vec![
            "Hello".to_string(),
            "Cruel".to_string(),
            "World".to_string(),
        ]).unwrap();
        let rand = rng.random().unwrap();
        let x = d.sample(rand).unwrap();
        assert_eq!(x, 6);
        let rand = rng.random().unwrap();
        let x = d.sample(rand).unwrap();
        assert_eq!(x, 4);
        let rand = rng.random().unwrap();
        let x = d.sample(rand).unwrap();
        assert_eq!(x, 4);
        let rand = rng.random().unwrap();
        let x = d.sample(rand).unwrap();
        assert_eq!(x, 2);
    }

    #[test]
    fn test_discrete_distribution_with_values() {
        let weights: Vec<f64> = vec![1.1, 2.0, 1.0, 8.0, 0.2, 2.0];
        let values: Vec<usize> = vec![2, 3, 7, 8, 9, 11];
        let d = DiscreteDistribution::new(&weights, &values).unwrap();
        let mut rng = NeatRng::new_from_seed(&vec![
            "Hello".to_string(),
            "Cruel".to_string(),
            "World".to_string(),
        ]).unwrap();
        let rand = rng.random().unwrap();
        let x = d.sample(rand).unwrap();
        assert_eq!(x, 11);
        let rand = rng.random().unwrap();
        let x = d.sample(rand).unwrap();
        assert_eq!(x, 8);
        let rand = rng.random().unwrap();
        let x = d.sample(rand).unwrap();
        assert_eq!(x, 8);
        let rand = rng.random().unwrap();
        let x = d.sample(rand).unwrap();
        assert_eq!(x, 3);
    }

    #[test]
    fn test_normal_distribution_params() {
        let nd = NormalDistribution::new(100.0, 15.0).unwrap();
        let (mean, std_dev) = nd.params().unwrap();
        assert_eq!(mean, 100.0);
        assert_eq!(std_dev, 15.0);
    }

    #[test]
    fn test_normal_distribution_sample_median() {
        // At p=0.5 the inverse CDF returns the mean exactly
        let nd = NormalDistribution::new(200.0, 30.0).unwrap();
        let median = nd.sample(0.5).unwrap();
        assert!((median - 200.0).abs() < 1e-9);
    }

    #[test]
    fn test_normal_distribution_sample_range() {
        // Samples between 0 and 1 should all land within ±4 standard deviations
        let nd = NormalDistribution::new(0.0, 1.0).unwrap();
        for p in [0.001, 0.01, 0.1, 0.5, 0.9, 0.99, 0.999] {
            let val = nd.sample(p).unwrap();
            assert!(val.abs() < 4.0, "sample({}) = {} is outside ±4 sigma", p, val);
        }
    }

    #[test]
    fn test_normal_distribution_inverse_cdf_symmetry() {
        // inverse_cdf(p) and inverse_cdf(1-p) should be symmetric around the mean
        let nd = NormalDistribution::new(50.0, 10.0).unwrap();
        let lo = nd.inverse_cdf(0.1).unwrap();
        let hi = nd.inverse_cdf(0.9).unwrap();
        assert!((lo + hi - 100.0).abs() < 1e-9);
    }
}
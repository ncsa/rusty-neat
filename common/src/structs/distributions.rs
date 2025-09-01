//! This DiscreteDistribution is an implementation of Zach Stephen's original neat-genReads code
//! from the py/probability.py file in tag 2.1 of github.com/ncsa/neat
//! (see also github.com/zstephens/neat-genreads). We may try the statrs Categorical distribution
//! as well, as I think it does the same thing.
use simple_rng::NeatRngError;
use thiserror::Error;
use serde::{Deserialize, Serialize};
use statrs::{distribution::{ContinuousCDF, Normal}, statistics::Distribution};

#[derive(Debug, Error)]
pub enum DistributionErrors {
    #[error("A distribution returned a sampling error: {0}.")]
    SamplingError(&'static str),
    #[error("A distribution returned an Rng Error: {0}")]
    RngError(NeatRngError),
    #[error("Requested a value vector from an index-only model")]
    VectorNotFound,
}

impl From<NeatRngError> for DistributionErrors {
    fn from(error: NeatRngError) -> Self {
        DistributionErrors::RngError(error)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DiscreteDistribution {
    values: Option<Vec<usize>>,
    weights: Vec<f64>,
}

impl DiscreteDistribution {
    pub fn new_index_only(w_vec: &Vec<f64>) -> Result<Self, DistributionErrors> {
        let cumulative_probability = {
            let sum_weights: f64 = w_vec.iter().sum();
            let mut normalized_weights = Vec::with_capacity(w_vec.len());
            // we no longer need the w_vec_64 after this, so we consume it
            for weight in w_vec{
                normalized_weights.push(weight / sum_weights);
            }
            cumulative_sum(&mut normalized_weights)
        }?;

        Ok(DiscreteDistribution {
            values: None,
            weights: cumulative_probability,
        })
    }

    pub fn new_from_values(w_vec: &Vec<f64>, l_vec: &Vec<usize>) -> Result<Self, DistributionErrors> {
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
                vec![0.0_f64; w_vec.len()]
            }
        };

        Ok(DiscreteDistribution {
            values: Some(l_vec.clone()),
            weights: cumulative_probability,
        })
    }

    pub fn values(&self) -> Result<Vec<usize>, DistributionErrors> {
        // fetches the value matrix, if it exists
        match &self.values {
            Some(vector) => Ok(vector.to_vec()),
            None => return Err(DistributionErrors::VectorNotFound),
        }
    }

    pub fn weights(&self) -> Result<Vec<f64>, DistributionErrors> {
        // fetches the value matrix, if it exists
        Ok(self.weights.to_vec())
    }

    pub fn sample_index(&self, rand: f64) -> Result<usize, DistributionErrors> {
        // returns a random index for the distribution, based on cumulative probability
        // This is basically an icdf for a discrete distribution
        // We take a random number as an input to avoid copying the RNG around the program.
        let mut lo: usize = 0;
        let mut hi: usize = self.weights.len();
        // bisect left
        while lo < hi {
            let mid = (lo + hi) / 2;
            if rand > self.weights[mid] {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        Ok(lo)
    }

    pub fn sample_values(&self, rand: f64) -> Result<usize, DistributionErrors> {
        // returns a random value for the distribution, based on cumulative probability
        // This is basically an icdf for a discrete distribution
        // We take a random number as an input to avoid copying the RNG around the program.
        let mut lo: usize = 0;
        let mut hi: usize = self.weights.len();
        // bisect left
        while lo < hi {
            let mid = (lo + hi) / 2;
            if rand > self.weights[mid] {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        match &self.values {
            Some(v) => Ok(v[lo]),
            None => Err(DistributionErrors::SamplingError(
                "DiscreteDistribution tried to sample a value from an instance with no values. Use the sample method instead."
            ))
        }
    }
}

fn cumulative_sum(a: &mut Vec<f64>) -> Result<Vec<f64>, DistributionErrors> {
    let mut acc = 0.0;
    let mut cumvec = Vec::with_capacity(a.len());
    for x in a {
        acc += *x;
        cumvec.push(acc);
    };
    Ok(cumvec)
}

#[derive(Debug, Clone, Copy)]
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


#[cfg(test)]
mod test {
    use simple_rng::NeatRng;
    use super::*;

    #[test]
    fn new_from_values() {
        let l_vec = vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
        let w_vec = vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        assert_eq!(l_vec.len(), w_vec.len());
        let distribution = DiscreteDistribution::new_from_values(&w_vec, &l_vec);

    }

    #[test]
    fn test_discrete_distribution() {
        let weights: Vec<f64> = vec![1.1, 2.0, 1.0, 8.0, 0.2, 2.0];
        let d = DiscreteDistribution::new_index_only(&weights).unwrap();
        let mut rng = NeatRng::new_from_seed(&vec![
            "Hello".to_string(),
            "Cruel".to_string(),
            "World".to_string(),
        ]).unwrap();
        let rand = rng.random().unwrap();
        let x = d.sample_index(rand).unwrap();
        assert_eq!(x, 5);
        let rand = rng.random().unwrap();
        let x = d.sample_index(rand).unwrap();
        assert_eq!(x, 3);
        let rand = rng.random().unwrap();
        let x = d.sample_index(rand).unwrap();
        assert_eq!(x, 3);
        let rand = rng.random().unwrap();
        let x = d.sample_index(rand).unwrap();
        assert_eq!(x, 1);
    }

    #[test]
    fn test_discrete_distribution_with_values() {
        let weights: Vec<f64> = vec![1.1, 2.0, 1.0, 8.0, 0.2, 2.0];
        let values: Vec<usize> = vec![2, 3, 7, 8, 9, 11];
        let d = DiscreteDistribution::new_from_values(&weights, &values).unwrap();
        let mut rng = NeatRng::new_from_seed(&vec![
            "Hello".to_string(),
            "Cruel".to_string(),
            "World".to_string(),
        ]).unwrap();
        let rand = rng.random().unwrap();
        let x = d.sample_index(rand).unwrap();
        assert_eq!(x, 5);
        let rand = rng.random().unwrap();
        let x = d.sample_index(rand).unwrap();
        assert_eq!(x, 3);
        let rand = rng.random().unwrap();
        let x = d.sample_index(rand).unwrap();
        assert_eq!(x, 3);
        let rand = rng.random().unwrap();
        let x = d.sample_index(rand).unwrap();
        assert_eq!(x, 1);
    }
}
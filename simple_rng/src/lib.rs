//! This RNG is based on Alea (https://github.com/nquinlan/better-random-numbers-for-javascript-mirror/
//! by Johannes Baagøe. This is a rust implementation. Some of the behavior of javascript I had
//! to guess at a little bit. I got the mash function to reproduce the results of the
//! javascript one (as tested in a browser runner online) and I used the TypeScript version
//! (https://rampantmonkey.com/writing/ts-prng/) to help me figure out some of the typing. Because
//! of overflow, I kept the numbers to within the range of u32, even though they are u64 and f64
//! outputs. I'm hoping the simplicity of this overall makes it very fast for NEAT
mod mash;
use thiserror::Error;
use mash::Mash;

#[derive(Error, Debug)]
pub enum NeatRngError {
    #[error("Neat RNG suffered an input error")]
    InputError,
    #[error("Neat RNG sufferd an output error")]
    OutputError,
    #[error("Neat RNG suffered a sampling error {0}")]
    SamplingError(String),
    #[error("Neat RNG suffered an invalid range error")]
    InvalidRangeError(String),
    #[error("Rand reported an error!")]
    RandRetrievalError,
}

#[derive(Debug, Clone, Copy)]
pub struct NeatRng {
    // This will be a simple seeded random number generator and some associated
    // functions. The other RNGs I have tried were too complicated and designed
    // with crypto in mind, which is much more complicated than we need for genomic
    // simulation. s0 and s1 are the seeds of our simulation. The original text mentioned
    // 2 calculations, and these are the intermediate factors. s2 is the current random number,
    // and it becomes s1 for the next iteration. c is a placeholder for a u32 that is preserved
    // over each iteration.
    s0: f64,
    s1: f64,
    s2: f64,
    c: u32,
}


impl NeatRng {
    pub fn new_from_seed(seed_list: &Vec<String>) -> Result<NeatRng, NeatRngError> {
        // The seed list is assumed to be a vector of strings. We'll have
        // to figure out a clever way to construct seed strings from random seeds
        // For the default string, we present the abstract of the initial paper up to the first

        // Initialize seeds.
        let mut masher = Mash::new();
        let mut s0 = masher.mash(&vec![' ']);
        let mut s1 = masher.mash(&vec![' ']);
        let mut s2 = masher.mash(&vec![' ']);
        let c = 1;

        // Update seeds.
        let mut seed_vec: Vec<Vec<char>> = Vec::new();
        for seed in seed_list {
            // vectorize the seed and add it to the master list
            let vector_seed = seed.chars().collect::<Vec<char>>();
            seed_vec.push(vector_seed.clone());
            // update the parameters
            s0 -= masher.mash(&vector_seed);
            if s0 < 0f64 {
                s0 += 1f64;
            }
            s1 -= masher.mash(&vector_seed);
            if s1 < 0f64 {
                s1 += 1f64;
            }
            s2 -= masher.mash(&vector_seed);
            if s2 < 0f64 {
                s2 += 1f64;
            }
        }

        Ok(NeatRng {
            s0,
            s1,
            s2,
            c,
        })
    }

<<<<<<< HEAD
    pub fn random(&mut self) -> Result<f64, NeatRngError> {
=======
    /// Generates a random value from the range [0, 1].
    pub fn random(&mut self) -> f64 {
>>>>>>> main
        // Not sure where 2091639 comes from, but let's roll with it
        // The other factor is the reciprocal of the max of u32
        let t = (2091639_f64 * self.s0) + (self.c as f64 * (1.0 / (u32::MAX as f64)));
        self.s0 = self.s1;
        self.s1 = self.s2;
        self.c = t.floor() as u32;
        self.s2 = t - self.c as f64;
        Ok(self.s2)
    }

<<<<<<< HEAD
    pub fn gen_bool(&mut self, frac: f64) -> Result<bool, NeatRngError> {
        // Checks if a random number is less than the requested frac, which is taken as the
        // percent chance of true.
        Ok(self.random().unwrap() < frac)
    }

    pub fn shuffle_in_place<T: Clone>(&mut self, a: &mut Vec<T>) -> Result<(), NeatRngError> {
        // reverse iterator
        for i in (0..=(a.len() - 1)).rev() {
            let j = (self.random().unwrap() * i as f64).floor() as usize;
=======
    /// Returns a bool with a probability `p` of being true.
    ///
    /// # Panics
    ///
    /// If `p` does *not* fall in range [0, 1].

    // This method is lifted from rand 0.85, but with restrictions edited so we can use our
    // custom RNG
    pub fn gen_bool(&mut self, p: f64) -> bool {
        if !(0.0..1.0).contains(&p) {
            if p == 1.0 {
                return true;
            }
            panic!("p={:?} is outside range [0.0, 1.0]", p);
        }

        // This is just `2.0.powi(64)`, but written this way because it is not available
        // in `no_std` mode. (from rand 0.8.5 docs). We used u64 max + 1, the equivalent.
        let scale = u64::MAX as f64 + 1.0;

        let p_int = (p * scale) as u64;
        let x = self.rand_u64();
        x < p_int
    }

    /// Shuffles elements in a sequence in place.
    pub fn shuffle<T: Clone>(&mut self, a: &mut Vec<T>) {
        for i in (1..a.len()).rev() {
            let j = (self.random() * i as f64).floor() as usize;
>>>>>>> main
            let temp = a[i].clone();
            a[i] = a[j].clone();
            a[j] = temp.clone();
        }
        Ok(())
    }

<<<<<<< HEAD
    pub fn range_i64(&mut self, min:i64, max:i64) -> Result<i64, String> {
        if max > i64::MAX {
            return Err("Max value out of range for i64".to_string())
        }
        Ok(((min.clone() as f64) + (self.random().unwrap() * ((max - min) as f64).round())) as i64)
    }

    pub fn rand_int(&mut self) -> Result<u64, NeatRngError> {
        let ret_num: f64 = self.random().unwrap() * (u64::MAX as f64);
        Ok(ret_num.round() as u64)
    }

    pub fn rand_u32(&mut self) -> Result<u32, NeatRngError> {
        let ret_num = self.rand_int().unwrap();
        // Modulo the number to u32 max to guarantee a valid u32
        Ok((ret_num % (u32::MAX as u64)) as u32)
    }

    pub fn choose<T: Clone>(&mut self, a: &Vec<T>) -> Result<T, NeatRngError> {
        // Randomly select based on which calculation comes up as 0 first
        // since i = 0 will force j = 0, this will always return an element
        for i in (0..=(a.len() - 1)).rev() {
            let j = (self.random().unwrap() * i as f64).floor() as usize;
            if j == 0 {
                return Ok(a[i].clone())
            }
        }
        Ok(a[0].clone())
    }
}

#[cfg(test)]
mod tests {
    use super::*;    

    #[test]
    fn test_gen_bool() {
        let mut rng = NeatRng::new_from_seed(&vec![
            "Hello".to_string(),
            "Cruel".to_string(),
            "World".to_string(),
        ]).unwrap();
        assert_eq!(rng.gen_bool(0.5).unwrap(), false);
        assert_eq!(rng.gen_bool(0.5).unwrap(), true);
        assert_eq!(rng.gen_bool(0.5).unwrap(), false);
        assert_eq!(rng.gen_bool(0.5).unwrap(), true);

=======
    /// Returns a random `i64` from a closed interval [min, max].
    pub fn range_i64(&mut self, min: i64, max: i64) -> i64 {
        ((min.clone() as f64) + (self.random() * ((max - min) as f64))) as i64
    }

    /// Returns a random `u64`.
    pub fn rand_u64(&mut self) -> u64 {
        let temp = self.random();
        // Keeping it within range of u32 to prevent overflow
        let ret_num: f64 = temp * (u64::MAX as f64);
        ret_num.trunc() as u64
    }

    /// Returns a random `u32`.
    pub fn rand_u32(&mut self) -> u32 {
        let ret_num = self.rand_u64();
        (ret_num % (u32::MAX as u64)) as u32
    }

    /// Picks an element from a sequence at random.
    pub fn choose<T: Clone>(&mut self, a: &[T]) -> T {
        let i = (self.random() * (a.len() - 1) as f64).floor() as usize;
        a[i].clone()
    }
}

/// Implements the normal distribution.
pub struct NormalDistribution {
    distribution: Normal,
}

impl NormalDistribution {
    /// Constructs a new normal distribution with mean of `mean` and standard deviation of
    /// `std_dev`.
    ///
    /// # Panics
    ///
    /// If `mean` or `std_dev` are `NaN` or `std_dev <= 0`.
    pub fn new(mean: f64, std_dev: f64) -> Self {
        NormalDistribution {
            distribution: Normal::new(mean, std_dev).unwrap(),
        }
    }

    /// Calculates the inverse cumulative distribution function for the normal distribution at `x`.
    ///
    /// # Panics
    ///
    /// If `x < 0` or `x > 1`.
    pub fn inverse_cdf(&self, x: f64) -> f64 {
        self.distribution.inverse_cdf(x)
    }

    /// Generates a random value using the given source of randomness.
    pub fn sample(&self, rng: &mut Rng) -> f64 {
        let x = rng.random();
        self.distribution.inverse_cdf(x)
    }
}

/// Implements the discrete distribution.
///
/// This distribution is an implementation of Zach Stephen's original code from
/// the `py/probability.py` file in version 2.1 of [neat-genreads] (see also [Neat]).
///
/// [neat-genreads]: https://github.com/zstephens/neat-genreads
/// [Neat]: https://github.com/ncsa/neat

// We may try the statrs Categorical distribution as well, as I think it does the same thing.
pub struct DiscreteDistribution {
    degenerate: bool,
    cumulative_probability: Vec<f64>,
}

impl DiscreteDistribution {
    /// Constructs a new discrete distribution with the probability masses defined by `prob_mass`.
    pub fn new<T>(prob_mass: &[T], degenerate: bool) -> Self
    where
        f64: From<T>,
        T: Copy,
    {
        let cumulative_probability = if !degenerate {
            // Convert probability masses to floats.
            let mut prob_mass: Vec<f64> = prob_mass.iter().map(|&w| f64::from(w)).collect();

            // Normalize them.
            let sum: f64 = prob_mass.iter().sum();
            prob_mass = prob_mass.iter().map(|&x| x / sum).collect();

            // Calculate the their cumulative sum.
            prob_mass
                .iter()
                .scan(0.0, |acc, &x| {
                    *acc += x;
                    Some(*acc)
                })
                .collect::<Vec<_>>()
        } else {
            vec![1.0]
        };

        DiscreteDistribution {
            degenerate,
            cumulative_probability,
        }
    }

    /// Generates a random value using the give source of randomness.
    pub fn sample(&self, rng: &mut Rng) -> usize {
        // returns a random index for the distribution, based on cumulative probability
        // This is basically an icdf for a discrete distribution
        let mut lo: usize = 0;
        let mut hi: usize = self.cumulative_probability.len();
        // bisect left
        if !self.degenerate {
            let r = rng.random();
            while lo < hi {
                let mid = (lo + hi) / 2;
                if r > self.cumulative_probability[mid] {
                    lo = mid + 1;
                } else {
                    hi = mid;
                }
            }
        }
        lo
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_discrete_distribution_sample() {
        let mut rng = Rng::from_seed(vec![
            "Hello".to_string(),
            "Cruel".to_string(),
            "World".to_string(),
        ]);
        let weights: Vec<f64> = vec![1.1, 2.0, 1.0, 8.0, 0.2, 2.0];
        let d = DiscreteDistribution::new(&weights, false);
        assert_eq!(d.sample(&mut rng), 5);
        assert_eq!(d.sample(&mut rng), 3);
        assert_eq!(d.sample(&mut rng), 3);
    }

    #[test]
    #[should_panic]
    fn test_normal_distribution_panic() {
        let _n = NormalDistribution::new(0.0, 0.0);
    }

    #[test]
    fn test_normal_distribution_inverse_cdf() {
        let n = NormalDistribution::new(0.0, 1.0);
        assert_eq!(n.inverse_cdf(0.5), 0.0);
    }

    #[test]
    fn test_normal_distribution_sample() {
        let mut rng = Rng::from_seed(vec![
            "Hello".to_string(),
            "Cruel".to_string(),
            "World".to_string(),
        ]);
        let n = NormalDistribution::new(0.0, 1.0);
        assert_eq!(n.sample(&mut rng), 1.181326790364034);
        assert_eq!(n.sample(&mut rng), -0.005968986117806898);
        assert_eq!(n.sample(&mut rng), 0.11272207687563508);
>>>>>>> main
    }

    #[test]
    fn test_random() {
<<<<<<< HEAD
        let mut rng = NeatRng::new_from_seed(&vec![
            "hello".to_string(),
            "cruel".to_string(),
            "world".to_string(),
        ]).unwrap();
        let test = rng.random().unwrap();
        assert_eq!(test, 0.8797469853889197);
        let test2 = rng.random().unwrap();
        assert_eq!(test2, 0.5001547893043607);
        let test3 = rng.random().unwrap();
        assert_eq!(test3, 0.6195652585010976);
    }

    #[test]
    fn test_shuffle_in_place() {
        let mut my_vec = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0];
        let mut rng = NeatRng::new_from_seed(&vec![
            "hello".to_string(),
            "cruel".to_string(),
            "world".to_string(),
        ]).unwrap();
        rng.shuffle_in_place(&mut my_vec).unwrap();
        assert_eq!(my_vec, vec![5.0, 7.0, 6.0, 3.0, 2.0, 1.0, 9.0, 4.0, 8.0]);
        rng.shuffle_in_place(&mut my_vec).unwrap();
        assert_eq!(my_vec, vec![4.0, 9.0, 1.0, 8.0, 3.0, 7.0, 2.0, 6.0, 5.0]);
    }

    #[test]
    fn test_range() {
        let min = 0;
        let max = 10;
        let mut rng = NeatRng::new_from_seed(&vec![
            "hello".to_string(),
            "cruel".to_string(),
            "world".to_string(),
        ]).unwrap();
        let num = rng.range_i64(min, max).unwrap();
        assert_eq!(num, 8);
        let num2 = rng.range_i64(min, max).unwrap();
        assert_eq!(num2, 5);
        assert_eq!(rng.range_i64(-3, 5), Ok(1));
    }

    #[test]
    fn test_rand_int() {
        let mut rng = NeatRng::new_from_seed(&vec![
            "hello".to_string(),
            "cruel".to_string(),
            "world".to_string(),
        ]).unwrap();
        assert_eq!(rng.rand_int().unwrap(), 16228467489086898176);
        assert_eq!(rng.rand_int().unwrap(), 9226227395537666048);
        assert_eq!(rng.rand_int().unwrap(), 11428961760531447808);
=======
        let mut rng = Rng::from_seed(vec![
            "hello".to_string(),
            "cruel".to_string(),
            "world".to_string(),
        ]);
        assert_eq!(rng.random(), 0.8797469853889197);
        assert_eq!(rng.random(), 0.5001547893043607);
        assert_eq!(rng.random(), 0.6195652585010976);
    }

    #[test]
    fn test_gen_bool() {
        let mut rng = Rng::from_seed(vec![
            "Hello".to_string(),
            "Cruel".to_string(),
            "World".to_string(),
        ]);
        assert_eq!(rng.gen_bool(0.5), false);
        assert_eq!(rng.gen_bool(0.5), true);
        assert_eq!(rng.gen_bool(0.5), false);
        assert_eq!(rng.gen_bool(0.5), true);
    }

    #[test]
    fn test_shuffle() {
        let mut vec = vec![1, 2, 3, 4, 5, 6, 7, 8, 9];
        let mut rng = Rng::from_seed(vec![
            "hello".to_string(),
            "cruel".to_string(),
            "world".to_string(),
        ]);
        rng.shuffle(&mut vec);
        assert_eq!(vec, vec![5, 7, 6, 3, 2, 1, 9, 4, 8]);
        rng.shuffle(&mut vec);
        assert_eq!(vec, vec![8, 1, 4, 9, 7, 3, 6, 5, 2]);
    }

    #[test]
    fn test_range_i64() {
        let mut rng = Rng::from_seed(vec![
            "hello".to_string(),
            "cruel".to_string(),
            "world".to_string(),
        ]);
        assert_eq!(rng.range_i64(0, 10), 8);
        assert_eq!(rng.range_i64(0, 10), 5);
        assert_eq!(rng.range_i64(-3, 5), 1);
    }

    #[test]
    fn test_rand_u64() {
        let mut rng = Rng::from_seed(vec![
            "hello".to_string(),
            "cruel".to_string(),
            "world".to_string(),
        ]);
        assert_eq!(rng.rand_u64(), 16228467489086898176);
        assert_eq!(rng.rand_u64(), 9226227395537666048);
        assert_eq!(rng.rand_u64(), 11428961760531447808);
    }

    #[test]
    fn test_rand_u32() {
        let mut rng = Rng::from_seed(vec![
            "hello".to_string(),
            "cruel".to_string(),
            "world".to_string(),
        ]);
        assert_eq!(rng.rand_u32(), 3778484531);
        assert_eq!(rng.rand_u32(), 2148148463);
        assert_eq!(rng.rand_u32(), 2661012523);
    }

    #[test]
    fn test_choose() {
        let mut rng = Rng::from_seed(vec![
            "hello".to_string(),
            "cruel".to_string(),
            "world".to_string(),
        ]);
        let vec: Vec<u32> = (0..10).collect();
        assert_eq!(rng.choose(&vec), 7);
        assert_eq!(rng.choose(&vec), 4);
        assert_eq!(rng.choose(&vec), 5);
>>>>>>> main
    }
}

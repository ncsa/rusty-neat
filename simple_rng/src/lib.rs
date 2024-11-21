/// This RNG is based on Alea (https://github.com/nquinlan/better-random-numbers-for-javascript-mirror/
/// by Johannes Baagøe. This is a rust implementation. Some of the behavior of javascript I had
/// to guess at a little bit. I got the mash function to reproduce the results of the
/// javascript one (as tested in a browser runner online) and I used the TypeScript version
/// (https://rampantmonkey.com/writing/ts-prng/) to help me figure out some of the typing. Because
/// of overflow, I kept the numbers to within the range of u32, even though they are u64 and f64
/// outputs. I'm hoping the simplicity of this overall makes it very fast for NEAT

mod mash;

use mash::Mash;
use statrs::distribution::{ContinuousCDF, Normal};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone)]
pub struct NeatRng {
    // This will be a simple seeded random number generator and some associated
    // functions. The other RNGs I have tried were too complicated and designed
    // with crypto in mind, which is much more complicated than we need for genomic
    // simulation. s0 and s1 are the seeds of our simulation. The original text mentioned
    // 2 calculations, and these are the intermediate factors. s2 is the current random number,
    // and it becomes s1 for the next iteration. c is a placeholder for a u32 that is preserved
    // over each iteration.
    pub seed_vec: Vec<Vec<char>>,
    s0: f64,
    s1: f64,
    s2: f64,
    c: u32,
}


impl NeatRng {
    pub fn new_from_seed(seed_list: Vec<String>) -> NeatRng {
        // The seed list is assumed to be a vector of strings. We'll have
        // to figure out a clever way to construct seed strings from random seeds
        // For the default string, we present the abstract of the initial paper up to the first

        // initialize seeds
        let mut masher = Mash::new();
        let mut s0 = masher.mash(&vec![' ']);
        let mut s1 = masher.mash(&vec![' ']);
        let mut s2 = masher.mash(&vec![' ']);
        let c = 1;
        let mut seed_vec: Vec::<Vec<char>> = Vec::new();
        // update seeds
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

        NeatRng {
            seed_vec,
            s0,
            s1,
            s2,
            c,
        }
    }

    pub fn random(&mut self) -> f64 {
        // Not sure where 2091639 comes from, but let's roll with it
        // The other factor is the reciprocal of the max of u32
        let t = (2091639_f64 * self.s0) + (self.c as f64 * (1.0/(u32::MAX as f64)));
        self.s0 = self.s1;
        self.s1 = self.s2;
        self.c = t.floor() as u32;
        self.s2 = t - self.c as f64;
        self.s2
    }

    pub fn gen_bool(&mut self, frac: f64) -> bool {
        // Checks if a random number is less than the requested frac, which is taken as the
        // percent chance of true.
        if !(0.0..1.0).contains(&frac) {
            if frac == 1.0 {
                return true
            }
            panic!("Invalid proportion for gen_bool {} (must be in [0.0, 1.0))", frac)
        }

        self.random() < frac
    }

    pub fn shuffle_in_place<T: Clone>(&mut self, a: &mut Vec<T>) {
        // reverse iterator
        for i in (0..=(a.len() - 1)).rev() {
            let j = (self.random() * i as f64).floor() as usize;
            let temp = a[i].clone();
            a[i] = a[j].clone();
            a[j] = temp.clone();
        }
    }

    pub fn range_i64(&mut self, min:i64, max:i64) -> Result<i64, String> {
        if max > i64::MAX {
            return Err("Max value out of range for i64".to_string())
        }
        Ok(((min.clone() as f64) + (self.random() * ((max - min) as f64).round())) as i64)
    }

    pub fn rand_int(&mut self) -> u64 {
        let temp = self.random();
        let ret_num: f64 = temp * (u64::MAX as f64);
        ret_num.round() as u64
    }

    pub fn rand_u32(&mut self) -> u32 {
        let ret_num = self.rand_int();
        // Modulo the number to u32 max to guarantee a valid u32
        (ret_num % (u32::MAX as u64)) as u32
    }

    pub fn choose<T: Clone>(&mut self, a: &Vec<T>) -> T {
        // Randomly select based on which calculation comes up as 0 first
        // since i = 0 will force j = 0, this will always return an element
        for i in (0..=(a.len() - 1)).rev() {
            let j = (self.random() * i as f64).floor() as usize;
            if j == 0 {
                return a[i].clone()
            }
        }
        a[0].clone()
    }
}

#[derive(Debug, Clone)]
pub struct NormalDistribution {
    distribution: Normal,
}

impl NormalDistribution {
    // This is little more than a wrapper for the statrs normal distribution. Other distributions
    // may be more complicated, as needed. I would copy in the code for custom tailoring, but the
    // inverse function requires tons of hardcoded coefficient tables
    // and no one has time for all that
    pub fn new(mean: f64, std_dev: f64) -> Self {
        NormalDistribution {
            distribution: Normal::new(mean, std_dev).unwrap()
        }
    }

    pub fn inverse_cdf(&self, x: f64) -> f64 {
        self.distribution.inverse_cdf(x)
    }

    pub fn sample(&self, rng: &mut NeatRng) -> f64 {
        // Takes a statrs NormalDistribution object, uses the ICDF to start with a probability (our
        // RNG, which generates numbers between 0 and 1 with approximately normal distribution) and
        // generate the corresponding Y value from the PDF. Very handy!
        let x = rng.random();
        self.distribution.inverse_cdf(x)
    }
}


/// This DiscreteDistribution is an implementation of Zach Stephen's original neat-genReads code
/// from the py/probability.py file in tag 2.1 of github.com/ncsa/neat
/// (see also github.com/zstephens/neat-genreads). We may try the statrs Categorical distribution
/// as well, as I think it does the same thing.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DiscreteDistribution {
    cumulative_probability: Vec<f64>,
}

impl DiscreteDistribution {
    pub fn new(w_vec: &Vec<f64>) -> Self
    {
        let cumulative_probability = {
            let sum_weights: f64 = w_vec.iter().sum();
            let mut normalized_weights = Vec::with_capacity(w_vec.len());
            // we no longer need the w_vec_64 after this, so we consume it
            for weight in w_vec{
                normalized_weights.push(weight / sum_weights);
            }
            cumulative_sum(&mut normalized_weights)
        };

        DiscreteDistribution {
            cumulative_probability,
        }
    }

    pub fn sample(&self, rng: &mut NeatRng) -> usize {
        // returns a random index for the distribution, based on cumulative probability
        // This is basically an icdf for a discrete distribution
        let mut lo: usize = 0;
        let mut hi: usize = self.cumulative_probability.len();
        // bisect left
        let r = rng.random();
        while lo < hi {
            let mid = (lo + hi) / 2;
            if r > self.cumulative_probability[mid] {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        lo
    }
}

fn cumulative_sum(a: &mut Vec<f64>) -> Vec<f64> {
    let mut acc = 0.0;
    let mut cumvec = Vec::with_capacity(a.len());
    for x in a {
        acc += *x;
        cumvec.push(acc);
    };
    cumvec
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_discrete_distribution() {
        let weights: Vec<f64> = vec![1.1, 2.0, 1.0, 8.0, 0.2, 2.0];
        let d = DiscreteDistribution::new(&weights);
        let mut rng = NeatRng::new_from_seed(vec![
            "Hello".to_string(),
            "Cruel".to_string(),
            "World".to_string(),
        ]);
        let x = d.sample(&mut rng);
        assert_eq!(x, 5);
        let x = d.sample(&mut rng);
        assert_eq!(x, 3);
        let x = d.sample(&mut rng);
        assert_eq!(x, 3);
        let x = d.sample(&mut rng);
        assert_eq!(x, 1);
    }

    #[test]
    fn test_gen_bool() {
        let mut rng = NeatRng::new_from_seed(vec![
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
    fn test_random() {
        let mut rng = NeatRng::new_from_seed(vec![
            "hello".to_string(),
            "cruel".to_string(),
            "world".to_string(),
        ]);
        let test = rng.random();
        assert_eq!(test, 0.8797469853889197);
        let test2 = rng.random();
        assert_eq!(test2, 0.5001547893043607);
        let test3 = rng.random();
        assert_eq!(test3, 0.6195652585010976);
    }

    #[test]
    fn test_shuffle_in_place() {
        let mut my_vec = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0];
        let mut rng = NeatRng::new_from_seed(vec![
            "hello".to_string(),
            "cruel".to_string(),
            "world".to_string(),
        ]);
        rng.shuffle_in_place(&mut my_vec);
        assert_eq!(my_vec, vec![5.0, 7.0, 6.0, 3.0, 2.0, 1.0, 9.0, 4.0, 8.0]);
        rng.shuffle_in_place(&mut my_vec);
        assert_eq!(my_vec, vec![4.0, 9.0, 1.0, 8.0, 3.0, 7.0, 2.0, 6.0, 5.0]);
    }

    #[test]
    fn test_range() {
        let min = 0;
        let max = 10;
        let mut rng = NeatRng::new_from_seed(vec![
            "hello".to_string(),
            "cruel".to_string(),
            "world".to_string(),
        ]);
        let num = rng.range_i64(min, max);
        assert_eq!(num, Ok(8));
        let num2 = rng.range_i64(min, max);
        assert_eq!(num2, Ok(5));
        assert_eq!(rng.range_i64(-3, 5), Ok(1));
    }

    #[test]
    fn test_rand_int() {
        let mut rng = NeatRng::new_from_seed(vec![
            "hello".to_string(),
            "cruel".to_string(),
            "world".to_string(),
        ]);
        assert_eq!(rng.rand_int(), 16228467489086898176);
        assert_eq!(rng.rand_int(), 9226227395537666048);
        assert_eq!(rng.rand_int(), 11428961760531447808);
    }
}

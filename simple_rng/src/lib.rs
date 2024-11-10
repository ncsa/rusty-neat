mod mash;

use mash::Mash;

#[derive(Debug)]
pub struct Rng {
    // This will be a simple seeded random number generator and some associated
    // functions. The other RNGs I have tried were too complicated and designed
    // with crypto in mind, which is much more complicated than we need for genomic
    // simulation
    pub seed_vec: Vec<Vec<char>>,
    s0: f64,
    s1: f64,
    s2: f64,
    c: u32,
}

impl Rng {
    pub fn from_string_list(seed_list: Vec<String>) -> Rng {
        // initialize seeds
        let mut masher = Mash::new();
        let mut s0 = masher.mash(&vec![' ']);
        let mut s1 = masher.mash(&vec![' ']);
        let mut s2 = masher.mash(&vec![' ']);
        let c = 1;
        let mut seed_vec: Vec::<Vec<char>> = Vec::new();
        // update seeds
        for seed in &seed_list {
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

        Rng {
            seed_vec,
            s0,
            s1,
            s2,
            c,
        }
    }

    pub fn random(&mut self) -> f64 {
        let t = 2091639f64 * self.s0 + (self.c as f64) * 2.3283064365386963e-10f64;
        self.s0 = self.s1;
        self.s1 = self.s2;
        self.c = ((t as u64) | 0) as u32;
        self.s2 = t - self.c as f64;
        self.s2
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

    pub fn range_i64(&mut self, min:i64, max:i64) -> i64 {
        ((min.clone() as f64) + (self.random() * ((max - min) as f64))) as i64
    }

    pub fn rand_int(&mut self) -> u64 {
        let temp = self.random();
        let ret_num: f64 = temp * (u64::MAX as f64);
        ret_num.trunc() as u64
    }

    pub fn choose<T: Clone>(&mut self, a: &mut Vec<T>) -> T {
        a[0].clone()
    }
}



#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_random() {
        let mut rng = Rng::from_string_list(vec![
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
        let mut my_vec = vec![1, 2, 3, 4, 5, 6, 7, 8, 9];
        let mut rng = Rng::from_string_list(vec![
            "hello".to_string(),
            "cruel".to_string(),
            "world".to_string(),
        ]);
        rng.shuffle_in_place(&mut my_vec);
        assert_eq!(my_vec, vec![5, 7, 6, 3, 2, 1, 9, 4, 8]);
        rng.shuffle_in_place(&mut my_vec);
        assert_eq!(my_vec, vec![4, 9, 1, 8, 3, 7, 2, 6, 5]);
    }

    #[test]
    fn test_range() {
        let min = 0;
        let max = 10;
        let mut rng = Rng::from_string_list(vec![
            "hello".to_string(),
            "cruel".to_string(),
            "world".to_string(),
        ]);
        let num = rng.range_i64(min, max);
        assert_eq!(num, 8);
        let num2 = rng.range_i64(min, max);
        assert_eq!(num2, 5);
        assert_eq!(rng.range_i64(-3, 5), 1);
    }

    #[test]
    fn test_rand_int() {
        let mut rng = Rng::from_string_list(vec![
            "hello".to_string(),
            "cruel".to_string(),
            "world".to_string(),
        ]);
        assert_eq!(rng.rand_int(), 16228467489086898176);
        assert_eq!(rng.rand_int(), 9226227395537666048);
        assert_eq!(rng.rand_int(), 11428961760531447808);
    }
}

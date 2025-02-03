//! A masher for generating a seed for Alea pseudo-random number generator.
//!
//! A masher is a hash function with an initialized internal state that is used to hash
//! the PRNG's seed (strings). Each element of the internal state is altered by a hash of each
//! argument.
//!
//! Current implementation was adapted from a TypeScript version found [here].
//!
//! [here]: https://rampantmonkey.com/writing/ts-prng/

pub struct Mash {
    n: u64,
}

impl Mash {
    /// Construct a new masher.
    pub fn new() -> Mash {
        // Unclear where this number comes from.
        Mash { n: 0xefc8249d }
    }

    /// Create a seed for a PRNG.

    // In case the site we're referring to will go offline, here's a verbatim copy from the source:
    //
    // ```typescript
    // const Mash = (): ((data: string) => number) => {
    //   var n = 0xefc8249d;
    //
    //   return ((data: string): number => {
    //     for (var i = 0; i < data.length; i++) {
    //       n += data.charCodeAt(i);
    //       var h = 0.02519603282416938 * n;
    //       n = h >>> 0;
    //       h -= n;
    //       h *= n;
    //       n = h >>> 0;
    //       h -= n;
    //       n += h * 0x100000000; // 2^32
    //     }
    //     return (n >>> 0) * 2.3283064365386963e-10; // 2^-32
    //  });
    // };
    // ```
    pub fn mash(&mut self, input_data: &Vec<char>) -> f64 {
        // It seemed like n was doing all the heavy lifting, so I created n_copy to bear the weight,
        // so we aren't touching the object as much. No idea if this is a good thought.
        let mut n_copy: u64 = self.n.clone();
        for char in input_data {
            n_copy += *char as u64;
            // Unclear where this number comes from, other than it being in the original
            let mut h: f64 = 0.02519603282416938 * (n_copy as f64);
            n_copy = h.floor() as u64;
            h -= n_copy as f64;
            h *= n_copy as f64;
            n_copy = h.floor() as u64;
            h -= n_copy as f64;
            // maybe just a difference between JS and Rust, but Rust does not allow 0x100000000,
            // since technically that is not a valid u32. Trying to keep the same constraints but
            // in a rust-friendly way. It does produce results without the +1, but they don't match
            // the javascript.
            n_copy += (h * (u32::MAX as f64 + 1.0)).floor() as u64;
        }
        // now update n
        self.n = n_copy.clone() as u64;
        self.n as f64 * (1.0 / (u32::MAX as f64 + 1.0))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mash() {
        let mut masher = Mash::new();
        let test1 = masher.mash(&"hello".chars().collect());
        assert_eq!(test1, 0.7957609200384468);

        let test2 = masher.mash(&"cruel".chars().collect());
        assert_eq!(test2, 0.8173183863982558);

        let test3 = masher.mash(&"world".chars().collect());
        assert_eq!(test3, 0.2441756660118699);
    }
}

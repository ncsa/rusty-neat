// Copyright 2018 Developers of the Rand project.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

// This is basically a direct lift from rand::rngs::StdRng. Why did I do it this way? I want
// flexibility to change the rng behind the scenes, and was having trouble implementing RNGCore
// myself. I could have just stolen the techniques, but decided to give credit.
// I have borrowed this code for Rust because I am not 100% how to borrow it the "right" way.
// And we may improve on this at some point, so this gives us more flexibility.

//! The standard RNG

use rand_core::{CryptoRng, Error, RngCore, SeedableRng};
use rand_chacha::ChaCha12Rng;

/// The standard RNG. The PRNG algorithm in `NeatRng` is chosen to be efficient
/// on the current platform, to be statistically strong and unpredictable
/// (meaning a cryptographically secure PRNG).
///
/// The current algorithm used is the ChaCha block cipher with 12 rounds. Please
/// see this relevant [rand issue] for the discussion. This may change as new
/// evidence of cipher security and performance becomes available.
///
/// The algorithm is deterministic but should not be considered reproducible
/// due to dependence on configuration and possible replacement in future
/// library versions. For a secure reproducible generator, we recommend use of
/// the [rand_chacha] crate directly.
///
/// [rand_chacha]: https://crates.io/crates/rand_chacha
/// [rand issue]: https://github.com/rust-random/rand/issues/932
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct NeatRng(ChaCha12Rng);

impl RngCore for NeatRng {
    fn next_u32(&mut self) -> u32 { self.0.next_u32() }

    fn next_u64(&mut self) -> u64 { self.0.next_u64() }

    fn fill_bytes(&mut self, dest: &mut [u8]) { self.0.fill_bytes(dest); }

    fn try_fill_bytes(&mut self, dest: &mut [u8]) -> Result<(), Error> {
        self.0.try_fill_bytes(dest)
    }
}

impl SeedableRng for NeatRng {
    type Seed = <ChaCha12Rng as SeedableRng>::Seed;


    fn from_seed(seed: Self::Seed) -> Self { NeatRng(ChaCha12Rng::from_seed(seed)) }


    fn from_rng<R: RngCore>(rng: R) -> Result<Self, Error> {
        ChaCha12Rng::from_rng(rng).map(NeatRng)
    }
}

impl CryptoRng for NeatRng {}

#[cfg(test)]
mod test {
    use super::*;
    use rand::{RngCore, SeedableRng};

    #[test]
    fn test_neat_rng_construction() {
        // Test value-stability of NeatRng. This is expected to break any time
        // the algorithm is changed.
        #[rustfmt::skip]
            let seed = [
            1,0,0,0, 23,0,0,0, 200,1,0,0, 210,30,0,0,
            0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0
        ];

        let target = [10719222850664546238, 14064965282130556830];

        let mut rng0 = NeatRng::from_seed(seed);
        let x0 = rng0.next_u64();

        let mut rng1 = NeatRng::from_rng(rng0).unwrap();
        let x1 = rng1.next_u64();

        assert_eq!([x0, x1], target);
    }
}

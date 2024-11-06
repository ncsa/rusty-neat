// Copyright 2017 10x Genomics

/// This code is borrowed from https://github.com/10XGenomics/rust-debruijn, where the contributors
/// used it to accelerate a debruijn graph construction. I thought about just using the crate, but
/// I needed to genericize some of the interactions so I think it is easier to copy the code here.
/// Their project is under an MIT license:
///
/// The MIT License (MIT)
///
/// Copyright (c) 2014-2018 10x Genomics, Inc.
///
/// Permission is hereby granted, free of charge, to any person obtaining a copy
/// of this software and associated documentation files (the "Software"), to deal
/// in the Software without restriction, including without limitation the rights
/// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
/// copies of the Software, and to permit persons to whom the Software is
/// furnished to do so, subject to the following conditions:
///
/// The above copyright notice and this permission notice shall be included in
/// all copies or substantial portions of the Software.
///
/// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
/// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
/// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
/// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
/// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
/// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
/// THE SOFTWARE.
///
/// All credit to that team for sorting out the details of this procedure, which I only kind of
/// understand, but I think I get it enough to start trying to use it.

/// Generate random genomes (with lots of re-used substrings), reassemble them, and check sanity

use crate::complement;
use crate::Kmer;
use crate::Vmer;

use rand::distributions::{Distribution, Gamma, Range};
use rand::{self, Rng, RngCore};
use std::cmp::{max, min};

/// Generate a uniformly random base
pub fn random_base(r: &mut impl Rng) -> u8 {
    (r.next_u64() % 4) as u8
}

/// Generate uniformly random DNA sequences
pub fn random_dna(len: usize) -> Vec<u8> {
    let mut r = rand::thread_rng();
    (0..len)
        .into_iter()
        .map(|_| (r.next_u64() % 4) as u8)
        .collect()
}

/// Randomly mutate each base with probability `p`
pub fn edit_dna<R: Rng>(seq: &mut [u8], p: f64, r: &mut R) {
    for b in seq.iter_mut() {
        if r.gen_range(0.0,1.0) < p {
            *b = random_base(r);
        }
    }
}

pub fn random_kmer<K: Kmer>() -> K {
    let mut r = rand::thread_rng();
    let mut kmer = K::empty();
    for pos in 0..K::k() {
        let b = (r.next_u64() % 4) as u8;
        kmer.set_mut(pos, b);
    }
    kmer
}

pub fn random_vmer<K: Kmer, V: Vmer>() -> V {
    let mut r = rand::thread_rng();
    let len = r.gen_range(K::k(), min(200, V::max_len()));
    let mut lmer = V::new(len);

    for pos in 0..len {
        let b = (r.next_u64() % 4) as u8;
        lmer.set_mut(pos, b);
    }
    lmer
}

pub fn simple_random_contigs() -> Vec<Vec<u8>> {
    let p1 = random_dna(40);
    let p2 = random_dna(30);

    let pc = random_dna(100);

    let p3 = random_dna(30);
    let p4 = random_dna(40);

    // Simulate contigs
    let mut c1 = Vec::new();
    c1.extend(p1);
    c1.extend(pc.clone());
    c1.extend(p3);

    let mut c2 = Vec::new();
    c2.extend(p2);
    c2.extend(pc);
    c2.extend(p4);

    let mut c3 = Vec::new();
    c3.extend(random_dna(30));

    // Stick a palindrome in the middle
    let palindrome1 = random_dna(33);
    let mut palindrome2 = palindrome1.clone();
    palindrome2.reverse();
    for v in &mut palindrome2 {
        *v = complement(*v);
    }

    c3.extend(palindrome1);
    c3.extend(palindrome2);
    c3.extend(random_dna(50));

    let contigs = vec![c1, c2, c3];
    contigs
}

// Generate random contigs with complicated repeats
pub fn random_contigs() -> Vec<Vec<u8>> {
    // Generate a bunch of sequence chunks
    let mut rng = rand::thread_rng();

    let gamma_dist = Gamma::new(0.6, 25.0);

    let nchunks = max(5, gamma_dist.sample(&mut rng) as u32);
    let chunk_sample = Range::new(0, nchunks);

    let length_dist = Gamma::new(1.5, 200.0);

    let mut chunks: Vec<Vec<u8>> = Vec::with_capacity(nchunks as usize);
    for _ in 0..nchunks {
        let len = max(10, length_dist.sample(&mut rng) as usize);
        let seq = random_dna(len);
        chunks.push(seq);
    }

    // Now make a bunch of chromosomes by pasting together chunks
    let nchrom = max(4, gamma_dist.sample(&mut rng) as u32);

    let mut chroms = Vec::with_capacity(nchrom as usize);
    for _ in 0..nchrom {
        let chrom_chunks = max(4, gamma_dist.sample(&mut rng) as u32);

        let mut chrom_seq = Vec::new();
        for _ in 0..chrom_chunks {
            let chunk_idx = chunk_sample.sample(&mut rng) as usize;
            chrom_seq.extend(chunks[chunk_idx].clone());
        }
        chroms.push(chrom_seq);
    }

    chroms
}

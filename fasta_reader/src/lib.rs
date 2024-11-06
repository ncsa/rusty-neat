use std::collections::{HashMap, VecDeque};
/// Most of this code is borrowed from https://github.com/10XGenomics/rust-debruijn, where the contributors
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

use serde_derive::{Deserialize, Serialize};
use std::{fmt, io};
use std::hash::Hash;
use crate::dna_string::DnaString;
use crate::NucleotideBlock::{NonN, NucBlock};
use common::file_tools::read_lines;

pub mod dna_string;

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
mod bitops_avx2;
#[cfg(test)]
pub mod test;

/// Enum to distinguish between nucleotide regions and unknown regions
#[derive(Debug, PartialEq)]
pub enum NucleotideBlock {
    // This is a region of unreadable bases with the length of the segment
    NonN(u64),
    // This is a block of actual DnaString data
    NucBlock(DnaString),
}


/// Convert a 2-bit representation of a base to a char
#[inline]
pub fn bits_to_ascii(c: u8) -> u8 {
    match c {
        0u8 => b'A',
        1u8 => b'C',
        2u8 => b'G',
        3u8 => b'T',
        _ => b'X',
    }
}

/// Convert an ASCII-encoded DNA base to a 2-bit representation
#[inline]
pub fn base_to_bits(c: u8) -> u8 {
    match c {
        b'A' | b'a' => 0u8,
        b'C' | b'c' => 1u8,
        b'G' | b'g' => 2u8,
        b'T' | b't' => 3u8,
        _ => 4u8,
    }
}


/// Convert an ASCII-encoded DNA base to a 2-bit representation
#[inline]
pub fn is_valid_base(c: u8) -> bool {
    matches!(c, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't')
}

/// Convert a 2-bit representation of a base to a char
#[inline]
pub fn bits_to_base(c: u8) -> char {
    match c {
        0u8 => 'A',
        1u8 => 'C',
        2u8 => 'G',
        3u8 => 'T',
        _ => 'X',
    }
}

/// The complement of a 2-bit encoded base
#[inline(always)]
pub fn complement(base: u8) -> u8 {
    (!base) & 0x3u8
}

/// Trait for interacting with DNA sequences
pub trait Sequence: Sized + fmt::Debug {
    /// Length of DNA sequence
    fn len(&self) -> usize;

    /// True if the sequence is empty.
    fn is_empty(&self) -> bool;

    /// Get 2-bit encoded base at position `pos`
    fn get(&self, pos: usize) -> u8;

    /// Set base at `pos` to 2-bit encoded base `val`
    fn set_mut(&mut self, pos: usize, val: u8);

    /// Set `nbases` positions in the sequence, starting at `pos`.
    /// Values must  be packed into the upper-most bits of `value`.
    fn set_slice_mut(&mut self, pos: usize, nbases: usize, value: u64);

    /// Return a new object containing the reverse complement of the sequence
    fn rc(&self) -> Self;

    /// Iterate over the bases in the sequence
    fn iter(&self) -> SeqIter<Self> {
        SeqIter {
            sequence: self,
            i: 0,
        }
    }

    /// Count the number of A/T bases in the kmer
    fn at_count(&self) -> u32 {
        let mut count = 0;
        for i in 0..self.len() {
            let base = self.get(i);
            if base == 0 || base == 3 {
                count += 1;
            }
        }
        count
    }

    /// Count the number of G/C bases in the kmer
    fn gc_count(&self) -> u32 {
        let mut count = 0;
        for i in 0..self.len() {
            let base = self.get(i);
            if base == 1 || base == 2 {
                count += 1;
            }
        }
        count
    }
}

/// Iterator over bases of a DNA sequence (bases will be unpacked into bytes).
pub struct SeqIter<'a, M: 'a + Sequence> {
    sequence: &'a M,
    i: usize,
}

impl<'a, M: 'a + Sequence> Iterator for SeqIter<'a, M> {
    type Item = u8;

    fn next(&mut self) -> Option<u8> {
        if self.i < self.sequence.len() {
            let value = self.sequence.get(self.i);
            self.i += 1;
            Some(value)
        } else {
            None
        }
    }
}

/// Encapsulates a Kmer sequence with statically known K.
pub trait Kmer: Sequence + Sized + Copy + PartialEq + PartialOrd + Eq + Ord + Hash {
    /// Create a Kmer initialized to all A's
    fn empty() -> Self;

    /// K value for this concrete type.
    fn k() -> usize;

    /// Return the rank of this kmer in an lexicographic ordering of all kmers
    /// E.g. 'AAAA' -> 0, 'AAAT' -> 1, etc. This will panic if K > 32.
    fn to_u64(&self) -> u64;

    /// Construct a kmer from the given lexicographic rank of the kmer.
    /// If K > 32, the leads bases will be A's.
    fn from_u64(value: u64) -> Self;

    // Compute Hamming distance between self and other
    fn hamming_dist(&self, other: Self) -> u32;

    /// Add the base `v` to the left side of the sequence, and remove the rightmost base
    fn extend_left(&self, v: u8) -> Self;

    /// Add the base `v` to the right side of the sequence, and remove the leftmost base
    fn extend_right(&self, v: u8) -> Self;

    /// Add the base `v` to the side of sequence given by `dir`, and remove a base at the opposite side
    fn extend(&self, v: u8, dir: Dir) -> Self {
        match dir {
            Dir::Left => self.extend_left(v),
            Dir::Right => self.extend_right(v),
        }
    }

    /// Generate all the extension of this sequence given by `exts` in direction `Dir`
    fn get_extensions(&self, exts: Exts, dir: Dir) -> Vec<Self> {
        let ext_bases = exts.get(dir);
        ext_bases.iter().map(|b| self.extend(*b, dir)).collect()
    }

    /// Return the minimum of the kmer and it's reverse complement, and a flag indicating if sequence was flipped
    fn min_rc_flip(&self) -> (Self, bool) {
        let rc = self.rc();
        if *self < rc {
            (*self, false)
        } else {
            (rc, true)
        }
    }

    /// Return the minimum of the kmer and it's reverse complement
    fn min_rc(&self) -> Self {
        let rc = self.rc();
        if *self < rc {
            *self
        } else {
            rc
        }
    }

    /// Test if this Kmer and it's reverse complement are the same
    fn is_palindrome(&self) -> bool {
        self.len() % 2 == 0 && *self == self.rc()
    }

    /// Create a Kmer from the first K bytes of `bytes`, which must be encoded as the integers 0-4.
    fn from_bytes(bytes: &[u8]) -> Self {
        if bytes.len() < Self::k() {
            panic!("bytes not long enough to form kmer")
        }

        let mut k0 = Self::empty();

        for (i, b) in bytes.iter().take(Self::k()).enumerate() {
            k0.set_mut(i, *b)
        }

        k0
    }

    /// Create a Kmer from the first K bytes of `bytes`, which must be encoded as ASCII letters A,C,G, or T.
    fn from_ascii(bytes: &[u8]) -> Self {
        if bytes.len() < Self::k() {
            panic!("bytes not long enough to form kmer")
        }

        let mut k0 = Self::empty();

        for (i, b) in bytes.iter().take(Self::k()).enumerate() {
            k0.set_mut(i, base_to_bits(*b))
        }

        k0
    }

    /// Return String containing Kmer sequence
    fn to_string(&self) -> String {
        let mut s = String::with_capacity(self.len());
        for pos in 0..self.len() {
            s.push(bits_to_base(self.get(pos)))
        }
        s
    }

    /// Generate vector of all kmers contained in `str` encoded as 0-4.
    fn kmers_from_bytes(str: &[u8]) -> Vec<Self> {
        if str.len() < Self::k() {
            return Vec::default();
        }
        let mut k0 = Self::empty();
        for (i, v) in str.iter().take(Self::k()).enumerate() {
            k0.set_mut(i, *v);
        }

        let mut r = Vec::with_capacity(str.len() - Self::k() + 1);
        r.push(k0);

        for v in str.iter().skip(Self::k()) {
            k0 = k0.extend_right(*v);
            r.push(k0);
        }

        r
    }

    /// Generate vector of all kmers contained in `str`, encoded as ASCII ACGT.
    fn kmers_from_ascii(str: &[u8]) -> Vec<Self> {
        if str.len() < Self::k() {
            return Vec::default();
        }
        let mut k0 = Self::empty();
        for (i, b) in str.iter().take(Self::k()).enumerate() {
            k0.set_mut(i, base_to_bits(*b));
        }

        let mut r = Vec::with_capacity(str.len() - Self::k() + 1);
        r.push(k0);

        for v in str.iter().skip(Self::k()) {
            k0 = k0.extend_right(base_to_bits(*v));
            r.push(k0);
        }

        r
    }
}

/// An immutable interface to a Mer sequence.
pub trait MerImmut: Sequence + Clone {
    fn set(&self, pos: usize, val: u8) -> Self {
        let mut new = self.clone();
        new.set_mut(pos, val);
        new
    }

    fn set_slice(&self, pos: usize, nbases: usize, bits: u64) -> Self {
        let mut new = self.clone();
        new.set_slice_mut(pos, nbases, bits);
        new
    }
}

impl<T> MerImmut for T where T: Sequence + Clone {}

/// A DNA sequence with run-time variable length, up to a statically known maximum length
pub trait Vmer: Sequence + PartialEq + Eq {
    /// Create a new sequence with length `len`, initialized to all A's
    fn new(len: usize) -> Self;

    /// Maximum sequence length that can be stored in this type
    fn max_len() -> usize;

    /// Create a Vmer from a sequence of bytes
    fn from_slice(seq: &[u8]) -> Self {
        let mut vmer = Self::new(seq.len());
        for (i, v) in seq.iter().enumerate() {
            vmer.set_mut(i, *v);
        }

        vmer
    }

    /// Efficiently extract a Kmer from the sequence
    fn get_kmer<K: Kmer>(&self, pos: usize) -> K;

    /// Get the first Kmer from the sequence
    fn first_kmer<K: Kmer>(&self) -> K {
        self.get_kmer(0)
    }

    /// Get the last kmer in the sequence
    fn last_kmer<K: Kmer>(&self) -> K {
        self.get_kmer(self.len() - K::k())
    }

    /// Get the terminal kmer of the sequence, on the both side of the sequence
    fn both_term_kmer<K: Kmer>(&self) -> (K, K) {
        (self.first_kmer(), self.last_kmer())
    }

    /// Get the terminal kmer of the sequence, on the side of the sequence given by dir
    fn term_kmer<K: Kmer>(&self, dir: Dir) -> K {
        match dir {
            Dir::Left => self.first_kmer(),
            Dir::Right => self.last_kmer(),
        }
    }

    /// Iterate over the kmers in the sequence
    fn iter_kmers<K: Kmer>(&self) -> KmerIter<'_, K, Self> {
        let kmer = if self.len() >= K::k() {
            self.first_kmer()
        } else {
            // Default kmer, will not be used
            K::empty()
        };

        KmerIter {
            bases: self,
            kmer,
            pos: K::k(),
        }
    }

    /// Iterate over the kmers and their extensions, given the extensions of the whole sequence
    fn iter_kmer_exts<K: Kmer>(&self, seq_exts: Exts) -> KmerExtsIter<'_, K, Self> {
        let kmer = if self.len() >= K::k() {
            self.first_kmer()
        } else {
            // Default kmer, will not be used
            K::empty()
        };

        KmerExtsIter {
            bases: self,
            exts: seq_exts,
            kmer,
            pos: K::k(),
        }
    }
}

/// A newtype wrapper around a `Vec<u8>` with implementations
/// of the `Mer` and `Vmer` traits.
#[derive(Debug, Clone, Eq, PartialEq, Ord, PartialOrd)]
pub struct DnaBytes(pub Vec<u8>);

impl Sequence for DnaBytes {
    fn len(&self) -> usize {
        self.0.len()
    }

    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    fn get(&self, pos: usize) -> u8 {
        self.0[pos]
    }

    /// Set base at `pos` to 2-bit encoded base `val`
    fn set_mut(&mut self, pos: usize, val: u8) {
        self.0[pos] = val
    }

    /// Set `nbases` positions in the sequence, starting at `pos`.
    /// Values must  be packed into the upper-most bits of `value`.
    fn set_slice_mut(&mut self, _pos: usize, _nbases: usize, _value: u64) {
        unimplemented!();
        //for i in pos .. (pos + nbases) {
        //
        //}
    }

    /// Return a new object containing the reverse complement of the sequence
    fn rc(&self) -> Self {
        unimplemented!();
    }
}

impl Vmer for DnaBytes {
    /// Create a new sequence with length `len`, initialized to all A's
    fn new(len: usize) -> Self {
        DnaBytes(vec![0u8; len])
    }

    /// Maximum sequence length that can be stored in this type
    fn max_len() -> usize {
        1 << 48
    }

    /// Efficiently extract a Kmer from the sequence
    fn get_kmer<K: Kmer>(&self, pos: usize) -> K {
        K::from_bytes(&self.0[pos..pos + K::k()])
    }
}

/// A newtype wrapper around a `&[u8]` with implementations
/// of the `Mer` and `Vmer` traits.
#[derive(Debug, Eq, PartialEq, Ord, PartialOrd)]
pub struct DnaSlice<'a>(pub &'a [u8]);

impl<'a> Sequence for DnaSlice<'a> {
    fn len(&self) -> usize {
        self.0.len()
    }

    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    fn get(&self, pos: usize) -> u8 {
        self.0[pos]
    }

    /// Set base at `pos` to 2-bit encoded base `val`
    fn set_mut(&mut self, _pos: usize, _val: u8) {
        unimplemented!()
    }

    /// Set `nbases` positions in the sequence, starting at `pos`.
    /// Values must  be packed into the upper-most bits of `value`.
    fn set_slice_mut(&mut self, _pos: usize, _nbases: usize, _value: u64) {
        unimplemented!();
        //for i in pos .. (pos + nbases) {
        //
        //}
    }

    /// Return a new object containing the reverse complement of the sequence
    fn rc(&self) -> Self {
        unimplemented!();
    }
}

impl<'a> Vmer for DnaSlice<'a> {
    /// Create a new sequence with length `len`, initialized to all A's
    fn new(_len: usize) -> Self {
        unimplemented!();
    }

    /// Maximum sequence length that can be stored in this type
    fn max_len() -> usize {
        1 << 48
    }

    /// Efficiently extract a Kmer from the sequence
    fn get_kmer<K: Kmer>(&self, pos: usize) -> K {
        K::from_bytes(&self.0[pos..pos + K::k()])
    }
}

/// Direction of motion in a DeBruijn graph
#[derive(Copy, Clone, Debug, Serialize, Deserialize)]
pub enum Dir {
    Left,
    Right,
}

impl Dir {
    /// Return a fresh Dir with the opposite direction
    pub fn flip(&self) -> Dir {
        match *self {
            Dir::Left => Dir::Right,
            Dir::Right => Dir::Left,
        }
    }

    /// Return a fresh Dir opposite direction if do_flip == True
    pub fn cond_flip(&self, do_flip: bool) -> Dir {
        if do_flip {
            self.flip()
        } else {
            *self
        }
    }

    /// Pick between two alternatives, depending on the direction
    pub fn pick<T>(&self, if_left: T, if_right: T) -> T {
        match self {
            Dir::Left => if_left,
            Dir::Right => if_right,
        }
    }
}

/// Store single-base extensions for a DNA Debruijn graph.
///
/// 8 bits, 4 higher order ones represent extensions to the right, 4 lower order ones
/// represent extensions to the left. For each direction the bits (from lower order
/// to higher order) represent whether there exists an extension with each of the
/// letters A, C, G, T. So overall the bits are:
///  right   left
/// T G C A T G C A
#[derive(Eq, PartialEq, Copy, Clone, Ord, PartialOrd, Hash, Serialize, Deserialize)]
pub struct Exts {
    pub val: u8,
}

impl Exts {
    pub fn new(val: u8) -> Self {
        Exts { val }
    }

    pub fn empty() -> Exts {
        Exts { val: 0u8 }
    }

    pub fn from_single_dirs(left: Exts, right: Exts) -> Exts {
        Exts {
            val: (right.val << 4) | (left.val & 0xf),
        }
    }

    pub fn merge(left: Exts, right: Exts) -> Exts {
        Exts {
            val: left.val & 0x0f | right.val & 0xf0,
        }
    }

    pub fn add(&self, v: Exts) -> Exts {
        Exts {
            val: self.val | v.val,
        }
    }

    pub fn set(&self, dir: Dir, pos: u8) -> Exts {
        let shift = pos
            + match dir {
            Dir::Right => 4,
            Dir::Left => 0,
        };

        let new_val = self.val | (1u8 << shift);
        Exts { val: new_val }
    }

    #[inline]
    fn dir_bits(&self, dir: Dir) -> u8 {
        match dir {
            Dir::Right => self.val >> 4,
            Dir::Left => self.val & 0xf,
        }
    }

    pub fn get(&self, dir: Dir) -> Vec<u8> {
        let bits = self.dir_bits(dir);
        let mut v = Vec::with_capacity(4);
        for i in 0..4 {
            if bits & (1 << i) > 0 {
                v.push(i);
            }
        }

        v
    }

    pub fn has_ext(&self, dir: Dir, base: u8) -> bool {
        let bits = self.dir_bits(dir);
        (bits & (1 << base)) > 0
    }

    pub fn from_slice_bounds(src: &[u8], start: usize, length: usize) -> Exts {
        let l_extend = if start > 0 {
            1u8 << (src[start - 1])
        } else {
            0u8
        };
        let r_extend = if start + length < src.len() {
            1u8 << src[start + length]
        } else {
            0u8
        };

        Exts {
            val: (r_extend << 4) | l_extend,
        }
    }

    pub fn from_dna_string(src: &dna_string::DnaString, start: usize, length: usize) -> Exts {
        let l_extend = if start > 0 {
            1u8 << (src.get(start - 1))
        } else {
            0u8
        };
        let r_extend = if start + length < src.len() {
            1u8 << src.get(start + length)
        } else {
            0u8
        };

        Exts {
            val: (r_extend << 4) | l_extend,
        }
    }

    pub fn num_exts_l(&self) -> u8 {
        self.num_ext_dir(Dir::Left)
    }

    pub fn num_exts_r(&self) -> u8 {
        self.num_ext_dir(Dir::Right)
    }

    pub fn num_ext_dir(&self, dir: Dir) -> u8 {
        let e = self.dir_bits(dir);
        (e & 1u8) + ((e & 2u8) >> 1) + ((e & 4u8) >> 2) + ((e & 8u8) >> 3)
    }

    pub fn mk_left(base: u8) -> Exts {
        Exts::empty().set(Dir::Left, base)
    }

    pub fn mk_right(base: u8) -> Exts {
        Exts::empty().set(Dir::Right, base)
    }

    pub fn mk(left_base: u8, right_base: u8) -> Exts {
        Exts::merge(Exts::mk_left(left_base), Exts::mk_right(right_base))
    }

    pub fn get_unique_extension(&self, dir: Dir) -> Option<u8> {
        if self.num_ext_dir(dir) != 1 {
            None
        } else {
            let e = self.dir_bits(dir);
            for i in 0..4 {
                if (e & (1 << i)) > 0 {
                    return Some(i);
                }
            }

            None
        }
    }

    pub fn single_dir(&self, dir: Dir) -> Exts {
        match dir {
            Dir::Right => Exts { val: self.val >> 4 },
            Dir::Left => Exts {
                val: self.val & 0xfu8,
            },
        }
    }

    /// Complement the extension bases for each direction
    pub fn complement(&self) -> Exts {
        let v = self.val;

        // swap bits
        let mut r = (v & 0x55u8) << 1 | ((v >> 1) & 0x55u8);

        // swap pairs
        r = (r & 0x33u8) << 2 | ((r >> 2) & 0x33u8);
        Exts { val: r }
    }

    pub fn reverse(&self) -> Exts {
        let v = self.val;
        let r = (v & 0xf) << 4 | (v >> 4);
        Exts { val: r }
    }

    pub fn rc(&self) -> Exts {
        self.reverse().complement()
    }
}

impl fmt::Debug for Exts {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut s = String::new();

        for b in self.get(Dir::Left) {
            s.push(bits_to_base(b));
        }
        s.push('|');

        for b in self.get(Dir::Right) {
            s.push(bits_to_base(b));
        }

        write!(f, "{}", s)
    }
}

/// Iterate over the `Kmer`s of a DNA sequence efficiently
pub struct KmerIter<'a, K: Kmer, D>
where
    D: 'a,
{
    bases: &'a D,
    kmer: K,
    pos: usize,
}

impl<'a, K: Kmer, D: Sequence> Iterator for KmerIter<'a, K, D> {
    type Item = K;

    #[inline]
    fn next(&mut self) -> Option<K> {
        if self.pos <= self.bases.len() {
            let retval = self.kmer;

            if self.pos < self.bases.len() {
                self.kmer = self.kmer.extend_right(self.bases.get(self.pos));
            }

            self.pos += 1;
            Some(retval)
        } else {
            None
        }
    }
}

/// Iterate over the `(Kmer, Exts)` tuples of a sequence and it's extensions efficiently
pub struct KmerExtsIter<'a, K: Kmer, D>
where
    D: 'a,
{
    bases: &'a D,
    exts: Exts,
    kmer: K,
    pos: usize,
}

impl<'a, K: Kmer, D: Sequence> Iterator for KmerExtsIter<'a, K, D> {
    type Item = (K, Exts);

    fn next(&mut self) -> Option<(K, Exts)> {
        if self.pos <= self.bases.len() {
            let next_base = if self.pos < self.bases.len() {
                self.bases.get(self.pos)
            } else {
                0u8
            };

            let cur_left = if self.pos == K::k() {
                self.exts
            } else {
                Exts::mk_left(self.bases.get(self.pos - K::k() - 1))
            };

            let cur_right = if self.pos < self.bases.len() {
                Exts::mk_right(next_base)
            } else {
                self.exts
            };

            let cur_exts = Exts::merge(cur_left, cur_right);

            let retval = self.kmer;
            self.kmer = self.kmer.extend_right(next_base);
            self.pos += 1;
            Some((retval, cur_exts))
        } else {
            None
        }
    }
}

pub fn read_fasta(
    fasta_path: &str,
) -> Result<(Box<HashMap<String, HashMap<u64, NucleotideBlock>>>, VecDeque<String>), io::Error> {
    // Reads a fasta file and turns it into a HashMap and puts it in the heap
    let mut fasta_map: HashMap<String, HashMap<u64, NucleotideBlock>> = HashMap::new();
    let mut fasta_order: VecDeque<String> = VecDeque::new();
    let mut current_key = String::new();
    let mut block_map: HashMap<u64, NucleotideBlock>;
    let lines = read_lines(fasta_path)?;
    let mut temp_seq = String::new();
    // Index will be used to create hashmaps of readable areas and non-readable areas
    let mut index: u64 = 0;
    // For the purposes of NEAT, everything that isn't ACTG or actg is treated like an unknown
    // base (N). NEAT does not yet have the capability of sorting out some of the alternate bases
    // like R indicating a purine (A or G).

    // If we encounter N's or unknown bases, this will count them up for us.
    let mut n_run: u64 = 0;
    let mut within_n_block = false;
    lines.for_each(|line| match line {
        Ok(l) => {
            if l.starts_with('>') {
                if !current_key.is_empty() {
                    if !temp_seq.is_empty() {
                        // if temp seq is not empty, then we need to append a nucleotide block
                        block_map
                            .entry(index.clone())
                            .or_insert(NucBlock(DnaString::from_dna_string(&temp_seq)));
                    } else {
                        // Otherwise we append a non-nucleotide block
                        block_map
                            .entry(index.clone())
                            .or_insert(NonN(n_run));
                        n_run = 0
                        // We don't update the index here, since we are starting a new contig
                    }
                    fasta_map
                        .entry(current_key.clone())
                        .or_insert(block_map.clone());
                    block_map = HashMap::new();
                    // Reset index to zero because it is a new contig
                    index = 0;
                    // For the same reason, reset within_n_block
                    within_n_block = false

                }
                // Set up the info for the contig we have found
                current_key = String::from(l.strip_prefix('>').unwrap());
                fasta_order.push_back(current_key.clone());
                temp_seq = String::new();
            } else {
                for char in l.chars() {
                    match char {
                        'A' | 'C' | 'G' | 'T' => {
                            // If we just exited an N block we need to clean up that
                            if within_n_block {
                                // Since we found Ns within the sequence, we need to append those
                                // blocks to the fasta map and now that we have passed the N's,
                                // we need to start counting up the bases
                                within_n_block = false;
                                // This gives the start location and length of the N-run
                                block_map
                                    .entry(index)
                                    .or_insert(NonN(n_run));
                                index += n_run;
                                n_run = 0;
                            }
                            // Either way, we push the new nucleotide
                            temp_seq.push(char)
                        }
                        _ => {
                            if !within_n_block {
                                // This is the start of a new N-block, so we reset the counter,
                                // set the flag to true, and proceed.
                                within_n_block = true;
                                n_run += 1;
                                // adding any previous nucleotide block
                                if !temp_seq.is_empty() {
                                    block_map
                                        .entry(index)
                                        .or_insert(NucBlock(DnaString::from_dna_string(&temp_seq)));
                                    index += temp_seq.len() as u64;
                                    temp_seq = String::new();
                                }

                            } else {
                                // Each new N encountered simply adds 1 to the total.
                                n_run += 1;
                            }
                        }
                    };
                }
            }
        }
        Err(error) => std::panic!("Problem reading fasta file: {}", error),
    });
    // Need to pick up the last block.
    if within_n_block {
        block_map
            .entry(index)
            .or_insert(NonN(n_run));
    } else {
        if !temp_seq.is_empty() {
            block_map
                .entry(index)
                .or_insert(NucBlock(DnaString::from_dna_string(&temp_seq)));
        }
    }
    Ok((Box::new(fasta_map), fasta_order))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_fasta() {
        let test_fasta = "test_data/references/small.fa";
        let (test_map, map_order) = read_fasta(test_fasta).unwrap();
        assert_eq!(map_order[0], "H1N1_HA".to_string());
        let expected_map = HashMap::from(
            ("H1N1_HA", HashMap::from([
                (0, NonN(8)),
                (8, NucBlock(DnaString::from_dna_string("TTCT")))
            ]))
        );
        assert_eq!(*test_map, expected_map)
    }
}

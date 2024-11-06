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
///
/// Copyright 2014 Johannes Köster and 10x Genomics
/// Licensed under the MIT license (http://opensource.org/licenses/MIT)
/// This file may not be copied, modified, or distributed
/// except according to those terms.

/// A 2-bit encoding of arbitrary length DNA sequences.
///
/// Store arbitrary-length DNA strings in a packed 2-bit encoding. Individual base values are encoded
/// as the integers 0,1,2,3 corresponding to A,C,G,T.


use serde_derive::{Deserialize, Serialize};
use std::borrow::Borrow;
use std::cmp::min;
use std::collections::hash_map::DefaultHasher;
use std::fmt;
use std::hash::{Hash, Hasher};

use crate::base_to_bits;
use crate::bits_to_ascii;
use crate::bits_to_base;

const BLOCK_BITS: usize = 64;
const WIDTH: usize = 2;

const MASK: u64 = 0x3;

/// A container for sequence of DNA bases.
#[derive(Ord, PartialOrd, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct DnaString {
    storage: Vec<u64>,
    len: usize,
}

impl DnaString {
    /// Create an empty DNA string
    pub fn new() -> DnaString {
        DnaString {
            storage: Vec::new(),
            len: 0,
        }
    }

    /// Length of the sequence
    pub fn len(&self) -> usize {
        self.len
    }

    /// Create a new instance with a given capacity.
    pub fn with_capacity(n: usize) -> Self {
        let blocks = ((n * WIDTH) >> 6) + (if (n * WIDTH) & 0x3F > 0 { 1 } else { 0 });
        let storage = Vec::with_capacity(blocks);

        DnaString { storage, len: 0 }
    }

    /// Create a DnaString of length n initialized to all A's
    pub fn blank(n: usize) -> Self {
        let blocks = ((n * WIDTH) >> 6) + (if (n * WIDTH) & 0x3F > 0 { 1 } else { 0 });
        let storage = vec![0; blocks];

        DnaString { storage, len: n }
    }

    /// Create a DnaString corresponding to an ACGT-encoded str.
    pub fn from_dna_string(dna: &str) -> DnaString {
        let mut dna_string = DnaString {
            storage: Vec::new(),
            len: 0,
        };

        dna_string.extend(dna.chars().map(|c| base_to_bits(c as u8)));
        dna_string
    }

    /// Create a DnaString from an ASCII ACGT-encoded byte slice.
    /// Non ACGT positions will be converted to 'A'
    pub fn from_acgt_bytes(bytes: &[u8]) -> DnaString {
        let mut dna_string = DnaString::with_capacity(bytes.len());

        // Accelerated avx2 mode. Should run on most machines made since 2013.
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        {
            if is_x86_feature_detected!("avx2") {
                for chunk in bytes.chunks(32) {
                    if chunk.len() == 32 {
                        let (conv_chunk, _) = unsafe { crate::bitops_avx2::convert_bases(chunk) };
                        let packed = unsafe { crate::bitops_avx2::pack_32_bases(conv_chunk) };
                        dna_string.storage.push(packed);
                    } else {
                        let b = chunk.iter().map(|c| base_to_bits(*c));
                        dna_string.extend(b);
                    }
                }

                dna_string.len = bytes.len();
                return dna_string;
            }
        }

        let b = bytes.iter().map(|c| base_to_bits(*c));
        dna_string.extend(b);
        dna_string
    }

    /// Create a DnaString from an ACGT-encoded byte slice,
    /// Non ACGT positions will be converted to repeatable random base determined
    /// by a hash of the read name and the position within the string.
    pub fn from_acgt_bytes_hashn(bytes: &[u8], read_name: &[u8]) -> DnaString {
        let mut hasher = DefaultHasher::new();
        read_name.hash(&mut hasher);

        let mut dna_string = DnaString::with_capacity(bytes.len());

        for (pos, c) in bytes.iter().enumerate() {
            let v = match c {
                b'A' | b'a' => 0u8,
                b'C' | b'c' => 1u8,
                b'G' | b'g' => 2u8,
                b'T' | b't' => 3u8,
                _ => {
                    let mut hasher_clone = hasher.clone();
                    pos.hash(&mut hasher_clone);
                    (hasher_clone.finish() % 4) as u8
                }
            };

            dna_string.push(v);
        }

        dna_string
    }

    /// Create a DnaString from a 0-4 encoded byte slice
    pub fn from_bytes(bytes: &[u8]) -> DnaString {
        let mut dna_string = DnaString {
            storage: Vec::new(),
            len: 0,
        };

        dna_string.extend(bytes.iter().cloned());
        dna_string
    }

    /// Convert sequence to a Vector of 0-4 encoded bytes
    pub fn to_bytes(&self) -> Vec<u8> {
        self.iter().collect()
    }

    /// Convert sequence to a Vector of ascii-encoded bytes
    pub fn to_ascii_vec(&self) -> Vec<u8> {
        self.iter().map(bits_to_ascii).collect()
    }

    /// Append a 0-4 encoded base.
    #[inline]
    pub fn push(&mut self, value: u8) {
        let (block, bit) = self.addr(self.len);
        if bit == 0 && block >= self.storage.len() {
            self.storage.push(0);
        }
        self.set_by_addr(block, bit, value);
        self.len += 1;
    }

    pub fn extend(&mut self, mut bytes: impl Iterator<Item = u8>) {
        // fill the last incomplete u64 block
        while self.len % 32 != 0 {
            match bytes.next() {
                Some(b) => self.push(b),
                None => return,
            }
        }

        let mut bytes = bytes.peekable();

        // chunk the remaining items into groups of at most 32 and handle them together
        while bytes.peek().is_some() {
            let mut val: u64 = 0;
            let mut offset = 62;
            let mut n_added = 0;

            for _ in 0..32 {
                if let Some(b) = bytes.next() {
                    assert!(b < 4);
                    val |= (b as u64) << offset;
                    offset -= 2;
                    n_added += 1;
                } else {
                    break;
                }
            }

            self.storage.push(val);
            self.len += n_added;
        }
    }

    /// Push 0-4 encoded bases from a byte array.
    ///
    /// # Arguments
    /// `bytes`: byte array to read values from
    /// `seq_length`: how many values to read from the byte array. Note that this
    /// is number of values not number of elements of the byte array.
    pub fn push_bytes(&mut self, bytes: &[u8], seq_length: usize) {
        assert!(
            seq_length <= bytes.len() * 8 / WIDTH,
            "Number of elements to push exceeds array length"
        );

        for i in 0..seq_length {
            let byte_index = (i * WIDTH) / 8;
            let byte_slot = (i * WIDTH) % 8;

            let v = bytes[byte_index];
            let bits = (v >> byte_slot) & (MASK as u8);

            self.push(bits);
        }
    }

    /// Iterate over stored values (values will be unpacked into bytes).
    pub fn iter(&self) -> DnaStringIter<'_> {
        DnaStringIter {
            dna_string: self,
            i: 0,
        }
    }

    /// Clear the sequence.
    pub fn clear(&mut self) {
        self.storage.clear();
        self.len = 0;
    }

    // Get the reverse complement of the string
    pub fn rc(&self) -> DnaString {
        let mut dna_string = DnaString::new();
        let rc = (0..self.len()).rev().map(|i| 3 - self.get(i));

        dna_string.extend(rc);
        dna_string
    }

    // Get the value at position 'i',
    pub fn get(&self, i: usize) -> u8 {
        let (block, bit) = self.addr(i);
        self.get_by_addr(block, bit)
    }

    #[inline(always)]
    fn get_by_addr(&self, block: usize, bit: usize) -> u8 {
        ((self.storage[block] >> (62 - bit)) & MASK) as u8
    }

    #[inline(always)]
    fn set_by_addr(&mut self, block: usize, bit: usize, value: u8) {
        let mask = MASK << (62 - bit);
        self.storage[block] |= mask;
        self.storage[block] ^= mask;
        self.storage[block] |= (value as u64 & MASK) << (62 - bit);
    }

    #[inline(always)]
    fn addr(&self, i: usize) -> (usize, usize) {
        let k = i * WIDTH;
        (k / BLOCK_BITS, k % BLOCK_BITS)
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Get the length `k` prefix of the DnaString
    pub fn prefix(&self, k: usize) -> DnaStringSlice<'_> {
        assert!(k <= self.len, "Prefix size exceeds number of elements.");
        DnaStringSlice {
            dna_string: self,
            start: 0,
            length: k,
            is_rc: false,
        }
    }

    /// Get the length `k` suffix of the DnaString
    pub fn suffix(&self, k: usize) -> DnaStringSlice<'_> {
        assert!(k <= self.len, "Suffix size exceeds number of elements.");

        DnaStringSlice {
            dna_string: self,
            start: self.len() - k,
            length: k,
            is_rc: false,
        }
    }

    /// Get slice containing the interval [`start`, `end`) of `self`
    pub fn slice(&self, start: usize, end: usize) -> DnaStringSlice<'_> {
        assert!(start <= self.len, "coordinate exceeds number of elements.");
        assert!(end <= self.len, "coordinate exceeds number of elements.");

        DnaStringSlice {
            dna_string: self,
            start,
            length: end - start,
            is_rc: false,
        }
    }

    /// Create a fresh DnaString containing the reverse of `self`
    pub fn reverse(&self) -> DnaString {
        let values: Vec<u8> = self.iter().collect();
        let mut dna_string = DnaString::new();
        for v in values.iter().rev() {
            dna_string.push(*v);
        }
        dna_string
    }

    /// Compute Hamming distance between this DnaString and another DnaString. The two strings must have the same length.
    pub fn hamming_distance(&self, other: &DnaString) -> usize {
        ndiffs(self, other)
    }

    // pub fn complement(&self) -> DnaString {
    //    assert!(self.width == 2, "Complement only supported for 2bit encodings.");
    //    let values: Vec<u32> = Vec::with_capacity(self.len());
    //    for i, v in self.storage.iter() {
    //        values[i] = v;
    //    }
    //    values[values.len() - 1] =
    // }
}

impl fmt::Display for DnaString {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for v in self.iter() {
            write!(f, "{}", bits_to_base(v))?;
        }
        Ok(())
    }
}

impl fmt::Debug for DnaString {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut s = String::new();
        for pos in 0..self.len() {
            s.push(bits_to_base(self.get(pos)))
        }

        write!(f, "{}", s)
    }
}

impl Default for DnaString {
    fn default() -> Self {
        Self::new()
    }
}

/// Iterator over values of a DnaStringoded sequence (values will be unpacked into bytes).
pub struct DnaStringIter<'a> {
    dna_string: &'a DnaString,
    i: usize,
}

impl<'a> Iterator for DnaStringIter<'a> {
    type Item = u8;

    fn next(&mut self) -> Option<u8> {
        if self.i < self.dna_string.len() {
            let value = self.dna_string.get(self.i);
            self.i += 1;
            Some(value)
        } else {
            None
        }
    }
}

impl<'a> IntoIterator for &'a DnaString {
    type Item = u8;
    type IntoIter = DnaStringIter<'a>;

    fn into_iter(self) -> DnaStringIter<'a> {
        self.iter()
    }
}

/// count Hamming distance between 2 2-bit DNA packed u64s
#[inline]
fn count_diff_2_bit_packed(a: u64, b: u64) -> u32 {
    let bit_diffs = a ^ b;
    let two_bit_diffs = (bit_diffs | bit_diffs >> 1) & 0x5555555555555555;
    two_bit_diffs.count_ones()
}

/// Compute the number of base positions at which two DnaStrings differ, assuming
/// that they have the same length.
pub fn ndiffs(b1: &DnaString, b2: &DnaString) -> usize {
    assert_eq!(b1.len(), b2.len());
    let mut diffs = 0;
    let (s1, s2) = (&b1.storage, &b2.storage);
    for i in 0..s1.len() {
        diffs += count_diff_2_bit_packed(s1[i], s2[i])
    }
    diffs as usize
}

/// An immutable slice into a DnaString
#[derive(Clone)]
pub struct DnaStringSlice<'a> {
    pub dna_string: &'a DnaString,
    pub start: usize,
    pub length: usize,
    pub is_rc: bool,
}

impl<'a> PartialEq for DnaStringSlice<'a> {
    fn eq(&self, other: &DnaStringSlice) -> bool {
        if other.length != self.length {
            return false;
        }
        for i in 0..self.length {
            if self.get(i) != other.get(i) {
                return false;
            }
        }
        true
    }
}
impl<'a> Eq for DnaStringSlice<'a> {}

impl<'a> DnaStringSlice<'a> {
    pub fn is_palindrome(&self) -> bool {
        unimplemented!();
    }

    /// Get the base at position `i`.
    #[inline(always)]
    pub fn get(&self, i: usize) -> u8 {
        if !self.is_rc {
            self.dna_string.get(i + self.start)
        } else {
            crate::complement(self.dna_string.get(self.start + self.length - 1 - i))
        }
    }

    pub fn rc(&self) -> DnaStringSlice<'a> {
        DnaStringSlice {
            dna_string: self.dna_string,
            start: self.start,
            length: self.length,
            is_rc: !self.is_rc,
        }
    }

    pub fn bytes(&self) -> Vec<u8> {
        let mut v = Vec::with_capacity(self.length);
        for pos in 0..self.length {
            v.push(self.get(pos));
        }
        v
    }

    pub fn ascii(&self) -> Vec<u8> {
        let mut v = Vec::with_capacity(self.length);
        for pos in 0..self.length {
            v.push(bits_to_ascii(self.get(pos)));
        }
        v
    }

    pub fn to_dna_string(&self) -> String {
        let mut dna: String = String::with_capacity(self.length);
        for pos in 0..self.length {
            dna.push(bits_to_base(self.get(pos)));
        }
        dna
    }

    pub fn to_owned(&self) -> DnaString {
        // FIXME make this faster
        let mut be = DnaString::with_capacity(self.length);
        for pos in 0..self.length {
            be.push(self.get(pos));
        }

        be
    }
    /// Get slice containing the interval [`start`, `end`) of `self`
    pub fn slice(&self, start: usize, end: usize) -> DnaStringSlice<'_> {
        assert!(
            start <= self.length,
            "coordinate exceeds number of elements."
        );
        assert!(end <= self.length, "coordinate exceeds number of elements.");
        assert!(end >= start, "invalid interval");

        if !self.is_rc {
            DnaStringSlice {
                dna_string: self.dna_string,
                start: self.start + start,
                length: end - start,
                is_rc: self.is_rc,
            }
        } else {
            // remap coords for RC
            let new_start = self.start + self.length - end;
            let new_length = end - start;

            DnaStringSlice {
                dna_string: self.dna_string,
                start: new_start,
                length: new_length,
                is_rc: self.is_rc,
            }
        }
    }
}

impl<'a> fmt::Display for DnaStringSlice<'a> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for pos in 0..self.length {
            write!(f, "{}", bits_to_base(self.get(pos)))?;
        }
        Ok(())
    }
}

impl<'a> fmt::Debug for DnaStringSlice<'a> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut s = String::new();
        if self.length < 256 {
            for pos in self.start..(self.start + self.length) {
                s.push(bits_to_base(self.dna_string.get(pos)))
            }
            write!(f, "{}", s)
        } else {
            write!(
                f,
                "start: {}, len: {}, is_rc: {}",
                self.start, self.length, self.is_rc
            )
        }
    }
}

/// Container for many distinct sequences, concatenated into a single DnaString.  Each
/// sequence is accessible by index as a DnaStringSlice.
#[derive(Debug, Default, Clone, Serialize, Deserialize)]
pub struct PackedDnaStringSet {
    pub sequence: DnaString,
    pub start: Vec<usize>,
    pub length: Vec<u32>,
}

impl PackedDnaStringSet {
    /// Create an empty `PackedDnaStringSet`
    pub fn new() -> Self {
        PackedDnaStringSet {
            sequence: DnaString::new(),
            start: Vec::new(),
            length: Vec::new(),
        }
    }

    /// Get a `DnaStringSlice` containing `i`th sequence in the set
    pub fn get(&self, i: usize) -> DnaStringSlice {
        DnaStringSlice {
            dna_string: &self.sequence,
            start: self.start[i],
            length: self.length[i] as usize,
            is_rc: false,
        }
    }

    /// Get a `DnaStringSlice` containing `i`th sequence in the set
    pub fn slice(&self, i: usize, start: usize, end: usize) -> DnaStringSlice {
        assert!(start <= self.length[i] as usize);
        assert!(end <= self.length[i] as usize);

        DnaStringSlice {
            dna_string: &self.sequence,
            start: self.start[i] + start,
            length: end - start,
            is_rc: false,
        }
    }

    /// Number of sequences in the set
    pub fn len(&self) -> usize {
        self.start.len()
    }

    pub fn is_empty(&self) -> bool {
        self.start.is_empty()
    }

    pub fn add<R: Borrow<u8>, S: IntoIterator<Item = R>>(&mut self, sequence: S) {
        let start = self.sequence.len();
        self.start.push(start);

        let mut length = 0;
        for b in sequence {
            self.sequence.push(*b.borrow());
            length += 1;
        }
        self.length.push(length as u32);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test;
    use rand::{self, Rng};

    fn hamming_dist_slow(s1: &DnaStringSlice, s2: &DnaStringSlice) -> u32 {
        assert_eq!(s1.len(), s2.len());

        let mut ndiff = 0;
        for pos in 0..s1.len() {
            if s1.get(pos) != s2.get(pos) {
                ndiff += 1;
            }
        }

        ndiff
    }

    /// Randomly mutate each base with probability `p`
    pub fn edit_dna_string(string: &mut DnaString, p: f64, r: &mut impl Rng) {
        for pos in 0..string.len() {
            if r.gen_range(0.0, 1.0) < p {
                let base = test::random_base(r);
                string.set_mut(pos, base)
            }
        }
    }

    fn random_dna_string(len: usize) -> DnaString {
        let bytes = test::random_dna(len);
        DnaString::from_bytes(&bytes)
    }


    fn random_slice(len: usize) -> (usize, usize) {
        if len == 0 {
            return (0, 0);
        }
        let mut r = rand::thread_rng();
        let a = (rand::RngCore::next_u64(&mut r) as usize) % len;
        let b = (rand::RngCore::next_u64(&mut r) as usize) % len;
        (std::cmp::min(a, b), std::cmp::max(a, b))
    }

    #[test]
    fn test_push_bytes() {
        let in_values: Vec<u8> = vec![2, 20];

        let mut dna_string = DnaString::new();
        dna_string.push_bytes(&in_values, 8);
        // Contents should be [01000000, 00101000]
        let values: Vec<u8> = dna_string.iter().collect();
        assert_eq!(values, [2, 0, 0, 0, 0, 1, 1, 0]);

        let mut dna_string = DnaString::new();
        dna_string.push_bytes(&in_values, 2);
        // Contents should be 01000000
        let values: Vec<u8> = dna_string.iter().collect();
        assert_eq!(values, [2, 0]);
    }

    // #[test]
    // fn test_from_dna_string() {
    //     dna_string_test("");
    //     dna_string_test("A");
    //     dna_string_test("C");
    //     dna_string_test("G");
    //     dna_string_test("T");
    //
    //     dna_string_test("GC");
    //     dna_string_test("ATA");
    //
    //     dna_string_test("ACGTACGT");
    //     dna_string_test("ACGTAAAAAAAAAATTATATAACGT");
    //     dna_string_test("AACGTAAAAAAAAAATTATATAACGT");
    // }

    // fn dna_string_test(dna: &str) {
    //     let dna_string_a = DnaString::from_dna_string(dna);
    //     let dna_string = DnaString::from_acgt_bytes(dna.as_bytes());
    //     assert_eq!(dna_string_a, dna_string);
    //
    //     let rc = dna_string_a.rc();
    //     let rc2 = rc.rc();
    //     assert_eq!(dna_string_a, rc2);
    //
    //     assert_eq!(dna_string.iter().count(), dna.len());
    //
    //     assert_eq!(dna_string.len, dna.len());
    //
    //     let dna_cp = dna_string.to_string();
    //     assert_eq!(dna, dna_cp);
    // }

    #[test]
    fn test_prefix() {
        let in_values: Vec<u8> = vec![2, 20];
        let mut dna_string = DnaString::new();
        dna_string.push_bytes(&in_values, 8);
        // Contents should be [01000000, 00101000]

        let pref_dna_string = dna_string.prefix(0).to_owned();
        assert_eq!(pref_dna_string.len(), 0);

        let pref_dna_string = dna_string.prefix(8).to_owned();
        assert_eq!(pref_dna_string, dna_string);

        let pref_dna_string = dna_string.prefix(4).to_owned();
        let values: Vec<u8> = pref_dna_string.iter().collect();
        assert_eq!(values, [2, 0, 0, 0]);

        let pref_dna_string = dna_string.prefix(6).to_owned();
        let values: Vec<u8> = pref_dna_string.iter().collect();
        assert_eq!(values, [2, 0, 0, 0, 0, 1]);

        dna_string.push_bytes(&in_values, 8);
        dna_string.push_bytes(&in_values, 8);

        let pref_dna_string = dna_string.prefix(17).to_owned();
        let values: Vec<u8> = pref_dna_string.iter().collect();
        assert_eq!(values, [2, 0, 0, 0, 0, 1, 1, 0, 2, 0, 0, 0, 0, 1, 1, 0, 2]);
    }

    #[test]
    fn test_suffix() {
        let in_values: Vec<u8> = vec![2, 20];
        let mut dna_string = DnaString::new();
        dna_string.push_bytes(&in_values, 8);
        // Contents should be [01000000, 00101000]

        let suf_dna_string = dna_string.suffix(0).to_owned();
        assert_eq!(suf_dna_string.len(), 0);

        let suf_dna_string = dna_string.suffix(8).to_owned();
        assert_eq!(suf_dna_string, dna_string);

        let suf_dna_string = dna_string.suffix(4).to_owned();
        let values: Vec<u8> = suf_dna_string.iter().collect();
        assert_eq!(values, [0, 1, 1, 0]);

        // 000101000000 64+256
        let suf_dna_string = dna_string.suffix(6).to_owned();
        let values: Vec<u8> = suf_dna_string.iter().collect();
        assert_eq!(values, [0, 0, 0, 1, 1, 0]);

        dna_string.push_bytes(&in_values, 8);
        dna_string.push_bytes(&in_values, 8);

        let suf_dna_string = dna_string.suffix(17).to_owned();
        let values: Vec<u8> = suf_dna_string.iter().collect();
        assert_eq!(values, [0, 2, 0, 0, 0, 0, 1, 1, 0, 2, 0, 0, 0, 0, 1, 1, 0]);
    }

    #[test]
    fn test_reverse() {
        let in_values: Vec<u8> = vec![2, 20];
        let mut dna_string = DnaString::new();
        let rev_dna_string = dna_string.reverse();
        assert_eq!(dna_string, rev_dna_string);

        dna_string.push_bytes(&in_values, 8);
        // Contents should be 00010100 00000010

        let rev_dna_string = dna_string.reverse();
        let values: Vec<u8> = rev_dna_string.iter().collect();
        assert_eq!(values, [0, 1, 1, 0, 0, 0, 0, 2]);
    }


    #[test]
    fn test_ndiffs() {
        let x1 = DnaString::from_dna_string("TGCATTAGAAAACTCCTTGCCTGTCTAGAAACTCATTAATCCACACATTGA");
        let x2 = DnaString::from_dna_string("TGCATTAGTAAACTCCTTCGCTGTCTAGAAAATCATTAAGCCACACATTGA");
        assert_eq!(ndiffs(&x1, &x2), 5);

        let x1 = DnaString::from_dna_string("TGCATT");
        let x2 = DnaString::from_dna_string("TGCATT");
        assert_eq!(ndiffs(&x1, &x2), 0);

        let x1 = DnaString::from_dna_string("");
        let x2 = DnaString::from_dna_string("");
        assert_eq!(ndiffs(&x1, &x2), 0);

        let x1 = DnaString::from_dna_string("TGCATTAGAAAACTCCTTGCCTGTCTAGAAACTCATTAATCCACACATTGA\
            TGCATTAGAAAACTCCTTGCCTGTCTAGAAACTCATTAATCCACACATTGATGCATTAGAAAACTCCTTGCCTGTCTAGAAACTCATTAATCCACACATTGA");
        let x2 = DnaString::from_dna_string("TGCATTAGTAAACTCCTTCGCTGTCTAGAAAATCATTAAGCCACACATTGA\
            TGCATTAGTAAACTCCTTCGCTGTCTAGAAAATCATTAAGCCACACATTGATGCATTAGTAAACTCCTTCGCTGTCTAGAAAATCATTAAGCCACACATTGA");
        assert_eq!(ndiffs(&x1, &x2), 15);
    }
}

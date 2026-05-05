use std::{io, path::PathBuf};

use noodles::bam;
use noodles::sam::alignment::record::{
    cigar::op::Kind as CigarKind,
    data::field::{Tag, Value},
    Flags,
};
use thiserror::Error;

use crate::structs::nucleotides::Nucleotide;

#[derive(Debug, Error)]
pub enum BamReaderError {
    #[error("I/O error reading BAM file: {0}")]
    IoError(#[from] io::Error),
}

const SKIP_FLAGS: Flags = Flags::UNMAPPED
    .union(Flags::SECONDARY)
    .union(Flags::SUPPLEMENTARY);

/// MD tag token used during mismatch extraction.
enum MdToken {
    Matches(usize),
    Mismatch(u8),
    Deletion,
}

fn parse_md(bytes: &[u8]) -> Vec<MdToken> {
    let mut tokens = Vec::new();
    let mut i = 0;
    while i < bytes.len() {
        match bytes[i] {
            b'0'..=b'9' => {
                let start = i;
                while i < bytes.len() && bytes[i].is_ascii_digit() {
                    i += 1;
                }
                let n: usize = std::str::from_utf8(&bytes[start..i])
                    .ok()
                    .and_then(|s| s.parse().ok())
                    .unwrap_or(0);
                tokens.push(MdToken::Matches(n));
            }
            b'^' => {
                i += 1;
                while i < bytes.len() && bytes[i].is_ascii_alphabetic() {
                    i += 1;
                }
                tokens.push(MdToken::Deletion);
            }
            b if b.is_ascii_alphabetic() => {
                tokens.push(MdToken::Mismatch(bytes[i]));
                i += 1;
            }
            _ => {
                i += 1;
            }
        }
    }
    tokens
}

/// Walks the parsed MD token list in sync with CIGAR M/=/X operations.
struct MdWalker {
    tokens: Vec<MdToken>,
    idx: usize,
    match_remaining: usize,
}

impl MdWalker {
    fn new(tokens: Vec<MdToken>) -> Self {
        Self {
            tokens,
            idx: 0,
            match_remaining: 0,
        }
    }

    /// Advance one M/=/X base. Returns the reference base byte for a mismatch, None for a match.
    fn next_alignment_base(&mut self) -> Option<u8> {
        loop {
            if self.match_remaining > 0 {
                self.match_remaining -= 1;
                return None;
            }
            match self.tokens.get(self.idx) {
                None => return None,
                Some(MdToken::Matches(0)) => {
                    self.idx += 1;
                }
                Some(MdToken::Matches(n)) => {
                    self.match_remaining = n - 1;
                    self.idx += 1;
                    return None;
                }
                Some(MdToken::Mismatch(b)) => {
                    let b = *b;
                    self.idx += 1;
                    return Some(b);
                }
                Some(MdToken::Deletion) => {
                    self.idx += 1;
                }
            }
        }
    }

    /// Advance past the next Deletion token (called once per D/N CIGAR operation).
    fn skip_deletion(&mut self) {
        while self.idx < self.tokens.len() {
            match &self.tokens[self.idx] {
                MdToken::Matches(0) => {
                    self.idx += 1;
                }
                MdToken::Deletion => {
                    self.idx += 1;
                    return;
                }
                _ => return,
            }
        }
    }
}

/// Minimum mapping quality for a read to contribute to the fragment length model.
const FRAG_FILTER_MAPQUAL: u8 = 10;

/// Reads a BAM or SAM file and returns the absolute template lengths (TLEN) for
/// paired, first-in-pair reads that are confidently mapped to the same reference
/// as their mate and have mapping quality > FRAG_FILTER_MAPQUAL.
///
/// Files with a `.sam` extension are read as plain-text SAM; all others are
/// treated as BGZF-compressed BAM.
pub fn read_fragment_lengths(path: &PathBuf) -> Result<Vec<usize>, BamReaderError> {
    let is_sam = path
        .extension()
        .and_then(|e| e.to_str())
        .map(|s| s.eq_ignore_ascii_case("sam"))
        .unwrap_or(false);

    if is_sam {
        read_fragment_lengths_sam(path)
    } else {
        read_fragment_lengths_bam(path)
    }
}

fn read_fragment_lengths_bam(path: &PathBuf) -> Result<Vec<usize>, BamReaderError> {
    let mut tlens = Vec::new();
    let mut reader = std::fs::File::open(path).map(bam::io::Reader::new)?;
    reader.read_header()?;
    for result in reader.records() {
        let record = result?;
        let flags = record.flags();
        // Must be paired and first in pair
        if !flags.is_segmented() || !flags.is_first_segment() {
            continue;
        }
        // Skip unmapped, secondary, supplementary, or mate-unmapped reads
        if flags.intersects(SKIP_FLAGS) || flags.is_mate_unmapped() {
            continue;
        }
        let mq: u8 = match record.mapping_quality() {
            Some(mq) => u8::from(mq),
            None => continue,
        };
        if mq <= FRAG_FILTER_MAPQUAL {
            continue;
        }
        let ref_id = record.reference_sequence_id().transpose()?;
        let mate_ref_id = record.mate_reference_sequence_id().transpose()?;
        if ref_id != mate_ref_id {
            continue;
        }
        let tlen = record.template_length().unsigned_abs() as usize;
        if tlen > 0 {
            tlens.push(tlen);
        }
    }
    Ok(tlens)
}

fn read_fragment_lengths_sam(path: &PathBuf) -> Result<Vec<usize>, BamReaderError> {
    use noodles::sam;
    use std::io::BufReader;

    let file = std::fs::File::open(path)?;
    let mut reader = sam::io::Reader::new(BufReader::new(file));
    let header = reader.read_header()?;
    let mut tlens = Vec::new();
    for result in reader.records() {
        let record = result?;
        let flags = record.flags()?;
        if !flags.is_segmented() || !flags.is_first_segment() {
            continue;
        }
        if flags.intersects(SKIP_FLAGS) || flags.is_mate_unmapped() {
            continue;
        }
        let mq: u8 = match record.mapping_quality() {
            Some(Ok(mq)) => u8::from(mq),
            _ => continue,
        };
        if mq <= FRAG_FILTER_MAPQUAL {
            continue;
        }
        let ref_id = record.reference_sequence_id(&header).transpose()?;
        let mate_ref_id = record.mate_reference_sequence_id(&header).transpose()?;
        if ref_id != mate_ref_id {
            continue;
        }
        let tlen = record.template_length()?.unsigned_abs() as usize;
        if tlen > 0 {
            tlens.push(tlen);
        }
    }
    Ok(tlens)
}

/// Reads a BAM file and accumulates a raw 4×4 SNP mismatch count matrix.
///
/// `counts[ref_base][read_base]` is incremented for each position where a read
/// base differs from the reference base. Both axes follow ALLOWED_NUCS order
/// (A=0, C=1, G=2, T=3). Bases that are not ACGT (N, IUPAC ambiguity codes,
/// etc.) are silently ignored.
///
/// Unmapped, secondary, and supplementary records are skipped. Records that
/// lack an MD tag are silently skipped; if no records with MD tags are found
/// the returned matrix will be all-zeros.
pub fn read_bam_transitions(path: &PathBuf) -> Result<[[usize; 4]; 4], BamReaderError> {
    let mut counts = [[0usize; 4]; 4];

    let mut reader = std::fs::File::open(path).map(bam::io::Reader::new)?;
    reader.read_header()?;

    for result in reader.records() {
        let record = result?;

        if record.flags().intersects(SKIP_FLAGS) {
            continue;
        }

        let md_bytes: Vec<u8> = match record.data().get(&Tag::MISMATCHED_POSITIONS) {
            Some(Ok(Value::String(s))) => s.iter().copied().collect(),
            _ => continue,
        };

        let tokens = parse_md(&md_bytes);
        let sequence = record.sequence();
        let mut walker = MdWalker::new(tokens);
        let mut read_pos = 0usize;

        for op_result in record.cigar().iter() {
            let op = op_result?;
            let len = op.len();

            match op.kind() {
                CigarKind::Match
                | CigarKind::SequenceMatch
                | CigarKind::SequenceMismatch => {
                    for _ in 0..len {
                        if let Some(ref_b) = walker.next_alignment_base() {
                            if let Some(read_b) = sequence.get(read_pos) {
                                let ri: usize = Nucleotide::from(ref_b as char).into();
                                let wi: usize = Nucleotide::from(read_b as char).into();
                                if ri < 4 && wi < 4 {
                                    counts[ri][wi] += 1;
                                }
                            }
                        }
                        read_pos += 1;
                    }
                }
                CigarKind::Insertion | CigarKind::SoftClip => {
                    read_pos += len;
                }
                CigarKind::Deletion | CigarKind::Skip => {
                    walker.skip_deletion();
                }
                CigarKind::HardClip | CigarKind::Pad => {}
            }
        }
    }

    Ok(counts)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_md_simple() {
        // "10A5G3": 10 matches, mismatch ref=A, 5 matches, mismatch ref=G, 3 matches
        let tokens = parse_md(b"10A5G3");
        let mut walker = MdWalker::new(tokens);
        // 10 matches
        for _ in 0..10 {
            assert!(walker.next_alignment_base().is_none());
        }
        // mismatch ref=A
        assert_eq!(walker.next_alignment_base(), Some(b'A'));
        // 5 matches
        for _ in 0..5 {
            assert!(walker.next_alignment_base().is_none());
        }
        // mismatch ref=G
        assert_eq!(walker.next_alignment_base(), Some(b'G'));
        // 3 matches
        for _ in 0..3 {
            assert!(walker.next_alignment_base().is_none());
        }
        // exhausted
        assert!(walker.next_alignment_base().is_none());
    }

    #[test]
    fn test_parse_md_with_deletion() {
        // "5^AT3": 5 matches, deletion AT, 3 matches
        let tokens = parse_md(b"5^AT3");
        let mut walker = MdWalker::new(tokens);
        for _ in 0..5 {
            assert!(walker.next_alignment_base().is_none());
        }
        walker.skip_deletion();
        for _ in 0..3 {
            assert!(walker.next_alignment_base().is_none());
        }
        assert!(walker.next_alignment_base().is_none());
    }

    #[test]
    fn test_parse_md_leading_zero_mismatch() {
        // "0T3": immediate mismatch at position 0, ref=T, 3 matches
        let tokens = parse_md(b"0T3");
        let mut walker = MdWalker::new(tokens);
        assert_eq!(walker.next_alignment_base(), Some(b'T'));
        for _ in 0..3 {
            assert!(walker.next_alignment_base().is_none());
        }
    }
}
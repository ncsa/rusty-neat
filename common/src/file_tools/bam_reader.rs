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
/// A read is kept only if its MAPQ is strictly greater than this value.
const FRAG_FILTER_MAPQUAL: u8 = 10;

// ── Walker abstraction ───────────────────────────────────────────────────────

/// Per-record filter applied by `walk_bam` before any observer sees a record.
///
/// A record is kept iff every enabled predicate passes. `min_mapq` is treated
/// as a strict lower bound (records with `mq <= min_mapq` are dropped);
/// `min_mapq = 0` disables the MAPQ check entirely, including the implicit
/// "drop records with no MAPQ" behavior.
pub struct BamWalkFilter {
    pub min_mapq: u8,
    pub skip_flags: Flags,
    pub require_paired: bool,
    pub require_first_in_pair: bool,
    pub require_mate_mapped: bool,
    pub require_same_ref_as_mate: bool,
}

impl BamWalkFilter {
    /// Matches the legacy `read_fragment_lengths_bam` filter: paired, first-in-pair,
    /// mate mapped to the same reference, MAPQ > FRAG_FILTER_MAPQUAL.
    pub fn for_frag_length() -> Self {
        Self {
            min_mapq: FRAG_FILTER_MAPQUAL,
            skip_flags: SKIP_FLAGS,
            require_paired: true,
            require_first_in_pair: true,
            require_mate_mapped: true,
            require_same_ref_as_mate: true,
        }
    }

    /// Matches the legacy `read_bam_transitions` filter: skip unmapped/secondary/
    /// supplementary only — no MAPQ filter, no pairing requirements.
    pub fn for_transitions() -> Self {
        Self {
            min_mapq: 0,
            skip_flags: SKIP_FLAGS,
            require_paired: false,
            require_first_in_pair: false,
            require_mate_mapped: false,
            require_same_ref_as_mate: false,
        }
    }
}

/// Visitor invoked once per record that passes the `BamWalkFilter`.
pub trait RecordObserver {
    fn observe(&mut self, record: &bam::Record) -> Result<(), BamReaderError>;
}

/// Counts returned by `walk_bam`.
#[derive(Debug, Default, Clone, Copy)]
pub struct WalkStats {
    pub records_seen: u64,
    pub records_kept: u64,
}

/// Single-pass walk over a BGZF BAM file. Each record that survives `filter`
/// is dispatched to every observer in order. Errors from any observer abort
/// the walk.
pub fn walk_bam(
    path: &PathBuf,
    filter: &BamWalkFilter,
    observers: &mut [&mut dyn RecordObserver],
) -> Result<WalkStats, BamReaderError> {
    let mut reader = std::fs::File::open(path).map(bam::io::Reader::new)?;
    reader.read_header()?;

    let mut stats = WalkStats::default();

    for result in reader.records() {
        let record = result?;
        stats.records_seen += 1;

        let flags = record.flags();
        if flags.intersects(filter.skip_flags) {
            continue;
        }
        if filter.require_paired && !flags.is_segmented() {
            continue;
        }
        if filter.require_first_in_pair && !flags.is_first_segment() {
            continue;
        }
        if filter.require_mate_mapped && flags.is_mate_unmapped() {
            continue;
        }
        if filter.min_mapq > 0 {
            let mq = match record.mapping_quality() {
                Some(mq) => u8::from(mq),
                None => continue,
            };
            if mq <= filter.min_mapq {
                continue;
            }
        }
        if filter.require_same_ref_as_mate {
            let ref_id = record.reference_sequence_id().transpose()?;
            let mate_ref_id = record.mate_reference_sequence_id().transpose()?;
            if ref_id != mate_ref_id {
                continue;
            }
        }

        stats.records_kept += 1;
        for obs in observers.iter_mut() {
            obs.observe(&record)?;
        }
    }

    Ok(stats)
}

// ── Observers ────────────────────────────────────────────────────────────────

/// Collects absolute template lengths (TLEN > 0). Pair with
/// `BamWalkFilter::for_frag_length()`.
#[derive(Debug, Default)]
pub struct FragLengthObserver {
    pub tlens: Vec<usize>,
}

impl RecordObserver for FragLengthObserver {
    fn observe(&mut self, record: &bam::Record) -> Result<(), BamReaderError> {
        let tlen = record.template_length().unsigned_abs() as usize;
        if tlen > 0 {
            self.tlens.push(tlen);
        }
        Ok(())
    }
}

/// Accumulates a 4×4 read-vs-reference mismatch count matrix from MD tags.
/// `counts[ref_base][read_base]` follows ALLOWED_NUCS order (A=0, C=1, G=2, T=3).
/// Pair with `BamWalkFilter::for_transitions()`.
#[derive(Debug, Default)]
pub struct TransitionObserver {
    pub counts: [[usize; 4]; 4],
}

impl RecordObserver for TransitionObserver {
    fn observe(&mut self, record: &bam::Record) -> Result<(), BamReaderError> {
        let md_bytes: Vec<u8> = match record.data().get(&Tag::MISMATCHED_POSITIONS) {
            Some(Ok(Value::String(s))) => s.iter().copied().collect(),
            _ => return Ok(()),
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
                                    self.counts[ri][wi] += 1;
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
        Ok(())
    }
}

// ── Public API (wrappers over walk_bam) ──────────────────────────────────────

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
        let mut obs = FragLengthObserver::default();
        walk_bam(path, &BamWalkFilter::for_frag_length(), &mut [&mut obs])?;
        Ok(obs.tlens)
    }
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
    let mut obs = TransitionObserver::default();
    walk_bam(path, &BamWalkFilter::for_transitions(), &mut [&mut obs])?;
    Ok(obs.counts)
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

    /// Builds a minimal BGZF BAM at `path` with one record per `records` entry.
    /// Fields: (flags, mapq, tlen, reference_sequence_id, mate_reference_sequence_id).
    /// Two reference sequences (chr1, chr2) are declared so mate-on-different-ref
    /// scenarios are expressible.
    fn write_test_bam(
        path: &std::path::PathBuf,
        records: &[(Flags, Option<u8>, i32, Option<usize>, Option<usize>)],
    ) {
        use noodles::sam::{
            self as sam,
            alignment::{
                RecordBuf,
                io::Write as _,
                record::{cigar::{op::Kind, Op}, MappingQuality},
                record_buf::{Cigar, Sequence},
            },
            header::record::value::{Map, map::ReferenceSequence},
        };
        let header = sam::Header::builder()
            .add_reference_sequence(
                b"chr1".to_vec(),
                Map::<ReferenceSequence>::new(std::num::NonZero::<usize>::new(1_000_000).unwrap()),
            )
            .add_reference_sequence(
                b"chr2".to_vec(),
                Map::<ReferenceSequence>::new(std::num::NonZero::<usize>::new(1_000_000).unwrap()),
            )
            .build();
        let file = std::fs::File::create(path).unwrap();
        let mut writer = bam::io::Writer::new(file);
        writer.write_header(&header).unwrap();
        for &(flags, mapq, tlen, ref_id, mate_ref_id) in records {
            let cigar: Cigar = [Op::new(Kind::Match, 4)].into_iter().collect();
            let mut record = RecordBuf::default();
            *record.flags_mut() = flags;
            *record.cigar_mut() = cigar;
            *record.sequence_mut() = Sequence::from(b"ACGT".as_ref());
            if let Some(mq) = mapq {
                *record.mapping_quality_mut() = Some(MappingQuality::try_from(mq).unwrap());
            }
            *record.reference_sequence_id_mut() = ref_id;
            *record.mate_reference_sequence_id_mut() = mate_ref_id;
            *record.template_length_mut() = tlen;
            writer.write_alignment_record(&header, &record).unwrap();
        }
    }

    #[test]
    fn test_walk_bam_dispatches_to_multiple_observers() {
        // One kept record reaches every observer in the slice.
        let temp = tempfile::tempdir().unwrap();
        let path = temp.path().join("t.bam");
        write_test_bam(
            &path,
            &[(
                Flags::SEGMENTED | Flags::FIRST_SEGMENT,
                Some(30),
                200,
                Some(0),
                Some(0),
            )],
        );
        let mut frag = FragLengthObserver::default();
        let mut trans = TransitionObserver::default();
        let stats = walk_bam(
            &path,
            &BamWalkFilter::for_frag_length(),
            &mut [&mut frag, &mut trans],
        )
        .unwrap();
        assert_eq!(stats.records_seen, 1);
        assert_eq!(stats.records_kept, 1);
        assert_eq!(frag.tlens, vec![200]);
        // No MD tag → no mismatches counted, but the observer was reached.
        assert_eq!(trans.counts, [[0usize; 4]; 4]);
    }

    #[test]
    fn test_walk_bam_min_mapq_is_strict_lower_bound() {
        // mq == min_mapq is dropped; mq > min_mapq is kept. Pins the `mq <= min_mapq`
        // semantics inherited from the original FRAG_FILTER_MAPQUAL check.
        let temp = tempfile::tempdir().unwrap();
        let path = temp.path().join("t.bam");
        let flags = Flags::SEGMENTED | Flags::FIRST_SEGMENT;
        write_test_bam(
            &path,
            &[
                (flags, Some(10), 100, Some(0), Some(0)),
                (flags, Some(11), 200, Some(0), Some(0)),
            ],
        );
        let mut frag = FragLengthObserver::default();
        let stats =
            walk_bam(&path, &BamWalkFilter::for_frag_length(), &mut [&mut frag]).unwrap();
        assert_eq!(stats.records_seen, 2);
        assert_eq!(stats.records_kept, 1);
        assert_eq!(frag.tlens, vec![200]);
    }

    #[test]
    fn test_walk_bam_min_mapq_zero_keeps_records_with_no_mapq() {
        // The transitions filter sets min_mapq = 0, which must disable the MAPQ
        // check entirely — including the "no MAPQ → drop" behavior that applies
        // when min_mapq > 0.
        let temp = tempfile::tempdir().unwrap();
        let path = temp.path().join("t.bam");
        write_test_bam(
            &path,
            &[(Flags::empty(), None, 0, Some(0), Some(0))],
        );
        let mut trans = TransitionObserver::default();
        let stats =
            walk_bam(&path, &BamWalkFilter::for_transitions(), &mut [&mut trans]).unwrap();
        assert_eq!(stats.records_seen, 1);
        assert_eq!(stats.records_kept, 1);
    }

    #[test]
    fn test_walk_bam_skip_flags_drops_secondary() {
        let temp = tempfile::tempdir().unwrap();
        let path = temp.path().join("t.bam");
        write_test_bam(
            &path,
            &[(Flags::SECONDARY, Some(30), 100, Some(0), Some(0))],
        );
        let mut trans = TransitionObserver::default();
        let stats =
            walk_bam(&path, &BamWalkFilter::for_transitions(), &mut [&mut trans]).unwrap();
        assert_eq!(stats.records_seen, 1);
        assert_eq!(stats.records_kept, 0);
    }

    #[test]
    fn test_walk_bam_require_same_ref_as_mate() {
        // Paired record with mate on a different reference is dropped by the
        // frag-length filter.
        let temp = tempfile::tempdir().unwrap();
        let path = temp.path().join("t.bam");
        write_test_bam(
            &path,
            &[(
                Flags::SEGMENTED | Flags::FIRST_SEGMENT,
                Some(30),
                100,
                Some(0),
                Some(1),
            )],
        );
        let mut frag = FragLengthObserver::default();
        let stats =
            walk_bam(&path, &BamWalkFilter::for_frag_length(), &mut [&mut frag]).unwrap();
        assert_eq!(stats.records_kept, 0);
    }
}
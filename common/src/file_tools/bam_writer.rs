use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::num::NonZero;
use std::path::PathBuf;
use std::sync::Arc;

use noodles::bam;
use noodles::bgzf;
use noodles::core::Position;
use noodles::sam::{
    self as sam,
    alignment::{
        io::Write as AlignmentWrite,
        record::{
            cigar::op::{Kind as CigarKind, Op},
Flags,
            MappingQuality,
        },
        record_buf::{Cigar, QualityScores, RecordBuf, Sequence},
    },
    header::record::value::{map::ReferenceSequence, Map},
};
use thiserror::Error;

use crate::structs::read_record::ReadRecord;

#[derive(Debug, Error)]
pub enum BamWriterError {
    #[error(transparent)]
    IoError(#[from] std::io::Error),
    #[error("Unknown contig: {0}")]
    UnknownContig(String),
}

pub struct BamWriter {
    writer: bam::io::Writer<bgzf::io::Writer<std::fs::File>>,
    header: sam::Header,
    pub contig_index: HashMap<String, usize>,
    carry: Vec<RecordBuf>,
}

impl BamWriter {
    /// Creates a new BAM file with a header derived from `contigs` (ordered name + length pairs).
    pub fn new(path: &PathBuf, contigs: &[(String, usize)]) -> Result<Self, BamWriterError> {
        let mut builder = sam::Header::builder();
        let mut contig_index = HashMap::new();
        for (i, (name, length)) in contigs.iter().enumerate() {
            let len = NonZero::<usize>::new(*length)
                .unwrap_or_else(|| NonZero::<usize>::new(1).unwrap());
            builder = builder.add_reference_sequence(
                name.as_bytes().to_vec(),
                Map::<ReferenceSequence>::new(len),
            );
            contig_index.insert(name.clone(), i);
        }
        let header = builder.build();
        let file = std::fs::File::create(path)?;
        let mut writer = bam::io::Writer::new(file);
        writer.write_header(&header)?;
        Ok(Self {
            writer,
            header,
            contig_index,
            carry: Vec::new(),
        })
    }

    pub fn write_record(&mut self, record: &RecordBuf) -> Result<(), BamWriterError> {
        self.writer.write_alignment_record(&self.header, record)?;
        Ok(())
    }

    /// Returns the 0-based index of `name` in the BAM header @SQ table.
    pub fn contig_id(&self, name: &str) -> Result<usize, BamWriterError> {
        self.contig_index
            .get(name)
            .copied()
            .ok_or_else(|| BamWriterError::UnknownContig(name.to_string()))
    }

    pub fn write_read_record(&mut self, record: &ReadRecord) -> Result<(), BamWriterError> {
        let ref_id = self.contig_id(&record.contig)?;
        let mate_ref_id = self.contig_id(&record.mate_contig)?;
        let cigar = rle_to_cigar(&record.cigar_ops);
        let flags = read_flags(record.is_paired, record.is_reverse);
        let bam_record = build_bam_record(
            &record.name,
            ref_id,
            record.position,
            flags,
            60,
            cigar,
            mate_ref_id,
            record.mate_position,
            record.template_length,
            &record.sequence,
            &record.quality_scores,
        );
        self.write_record(&bam_record)
    }

    /// Buffer a read for deferred coordinate-sorted output.
    /// Call `flush_up_to` after each block and `flush_all` after each contig.
    pub fn stage_read_record(&mut self, record: &ReadRecord) -> Result<(), BamWriterError> {
        let ref_id = self.contig_id(&record.contig)?;
        let mate_ref_id = self.contig_id(&record.mate_contig)?;
        let cigar = rle_to_cigar(&record.cigar_ops);
        let flags = read_flags(record.is_paired, record.is_reverse);
        let bam_record = build_bam_record(
            &record.name,
            ref_id,
            record.position,
            flags,
            60,
            cigar,
            mate_ref_id,
            record.mate_position,
            record.template_length,
            &record.sequence,
            &record.quality_scores,
        );
        self.carry.push(bam_record);
        Ok(())
    }

    /// Sort the carry buffer and write all records whose 0-based alignment start
    /// is less than `flush_pos`. Records at or beyond `flush_pos` are retained
    /// for the next flush (they may be joined by records from the next block).
    pub fn flush_up_to(&mut self, flush_pos: usize) -> Result<(), BamWriterError> {
        self.carry.sort_unstable_by_key(|r| {
            r.alignment_start().map(|p| p.get()).unwrap_or(0)
        });
        // alignment_start is 1-based; 0-based pos = p.get() - 1 < flush_pos ↔ p.get() <= flush_pos
        let split = self.carry.partition_point(|r| {
            r.alignment_start().map(|p| p.get()).unwrap_or(0) <= flush_pos
        });
        for record in self.carry.drain(..split) {
            self.writer.write_alignment_record(&self.header, &record)?;
        }
        Ok(())
    }

    /// Sort and write all remaining buffered records, then clear the carry.
    /// Call once after all blocks of a contig are processed.
    pub fn flush_all(&mut self) -> Result<(), BamWriterError> {
        self.carry.sort_unstable_by_key(|r| {
            r.alignment_start().map(|p| p.get()).unwrap_or(0)
        });
        for record in self.carry.drain(..) {
            self.writer.write_alignment_record(&self.header, &record)?;
        }
        Ok(())
    }

    pub fn carry_len(&self) -> usize {
        self.carry.len()
    }
}

const BGZF_EOF_LEN: usize = 28;
const BGZF_EOF_BLOCK: [u8; BGZF_EOF_LEN] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00,
    0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00,
];

/// Shared BAM header context for parallel per-contig writers.
/// Build once per run and wrap in `Arc` to share across threads.
pub struct BamContext {
    pub(crate) header: sam::Header,
    pub(crate) contig_index: HashMap<String, usize>,
}

impl BamContext {
    pub fn new(contigs: &[(String, usize)]) -> Self {
        let mut builder = sam::Header::builder();
        let mut contig_index = HashMap::new();
        for (i, (name, length)) in contigs.iter().enumerate() {
            let len = NonZero::<usize>::new(*length)
                .unwrap_or_else(|| NonZero::<usize>::new(1).unwrap());
            builder = builder.add_reference_sequence(
                name.as_bytes().to_vec(),
                Map::<ReferenceSequence>::new(len),
            );
            contig_index.insert(name.clone(), i);
        }
        Self { header: builder.build(), contig_index }
    }
}

/// Abstraction over BAM record staging for both single-file and per-contig writers.
pub trait BamRecordStager {
    fn stage_read_record(&mut self, record: &ReadRecord) -> Result<(), BamWriterError>;
    fn flush_up_to(&mut self, flush_pos: usize) -> Result<(), BamWriterError>;
}

impl BamRecordStager for BamWriter {
    fn stage_read_record(&mut self, record: &ReadRecord) -> Result<(), BamWriterError> {
        let ref_id = self.contig_id(&record.contig)?;
        let mate_ref_id = self.contig_id(&record.mate_contig)?;
        let cigar = rle_to_cigar(&record.cigar_ops);
        let flags = read_flags(record.is_paired, record.is_reverse);
        let bam_record = build_bam_record(
            &record.name, ref_id, record.position, flags, 60, cigar,
            mate_ref_id, record.mate_position, record.template_length,
            &record.sequence, &record.quality_scores,
        );
        self.carry.push(bam_record);
        Ok(())
    }

    fn flush_up_to(&mut self, flush_pos: usize) -> Result<(), BamWriterError> {
        self.carry.sort_unstable_by_key(|r| {
            r.alignment_start().map(|p| p.get()).unwrap_or(0)
        });
        let split = self.carry.partition_point(|r| {
            r.alignment_start().map(|p| p.get()).unwrap_or(0) <= flush_pos
        });
        for record in self.carry.drain(..split) {
            self.writer.write_alignment_record(&self.header, &record)?;
        }
        Ok(())
    }
}

/// Per-contig BAM writer that writes alignment records only (no header).
/// Intended to be concatenated with a header via `concat_temp_bams`.
pub struct BamBodyWriter {
    writer: bam::io::Writer<bgzf::io::Writer<std::fs::File>>,
    context: Arc<BamContext>,
    carry: Vec<RecordBuf>,
    pub path: PathBuf,
}

impl BamBodyWriter {
    pub fn new(path: PathBuf, context: Arc<BamContext>) -> Result<Self, BamWriterError> {
        let file = std::fs::File::create(&path)?;
        // Do not write a header — this file holds only alignment record BGZF blocks.
        let writer = bam::io::Writer::new(file);
        Ok(Self { writer, context, carry: Vec::new(), path })
    }

    fn contig_id(&self, name: &str) -> Result<usize, BamWriterError> {
        self.context.contig_index
            .get(name)
            .copied()
            .ok_or_else(|| BamWriterError::UnknownContig(name.to_string()))
    }

    pub fn flush_all(&mut self) -> Result<(), BamWriterError> {
        self.carry.sort_unstable_by_key(|r| {
            r.alignment_start().map(|p| p.get()).unwrap_or(0)
        });
        for record in self.carry.drain(..) {
            self.writer.write_alignment_record(&self.context.header, &record)?;
        }
        Ok(())
    }
}

impl BamRecordStager for BamBodyWriter {
    fn stage_read_record(&mut self, record: &ReadRecord) -> Result<(), BamWriterError> {
        let ref_id = self.contig_id(&record.contig)?;
        let mate_ref_id = self.contig_id(&record.mate_contig)?;
        let cigar = rle_to_cigar(&record.cigar_ops);
        let flags = read_flags(record.is_paired, record.is_reverse);
        let bam_record = build_bam_record(
            &record.name, ref_id, record.position, flags, 60, cigar,
            mate_ref_id, record.mate_position, record.template_length,
            &record.sequence, &record.quality_scores,
        );
        self.carry.push(bam_record);
        Ok(())
    }

    fn flush_up_to(&mut self, flush_pos: usize) -> Result<(), BamWriterError> {
        self.carry.sort_unstable_by_key(|r| {
            r.alignment_start().map(|p| p.get()).unwrap_or(0)
        });
        let split = self.carry.partition_point(|r| {
            r.alignment_start().map(|p| p.get()).unwrap_or(0) <= flush_pos
        });
        for record in self.carry.drain(..split) {
            self.writer.write_alignment_record(&self.context.header, &record)?;
        }
        Ok(())
    }
}

/// Concatenate per-contig temp BAM body files into a final coordinate-sorted BAM.
///
/// Writes the header once, then appends each body file (stripping its BGZF EOF block),
/// and terminates the output with a single BGZF EOF block.
pub fn concat_temp_bams(
    context: &BamContext,
    body_files: &[PathBuf],
    output_path: &PathBuf,
) -> Result<(), BamWriterError> {
    // Serialize the header to bytes, finishing with a BGZF EOF block.
    let mut header_bytes: Vec<u8> = Vec::new();
    {
        let mut writer = bam::io::Writer::new(&mut header_bytes);
        writer.write_header(&context.header)?;
        let _ = writer.into_inner().finish()?;
    }
    // Strip the trailing BGZF EOF so body files can be appended directly.
    let header_end = header_bytes.len().saturating_sub(BGZF_EOF_LEN);

    let output_file = std::fs::File::create(output_path)?;
    let mut out = BufWriter::new(output_file);
    out.write_all(&header_bytes[..header_end])?;

    if body_files.is_empty() {
        out.write_all(&BGZF_EOF_BLOCK)?;
    } else {
        for (i, body_path) in body_files.iter().enumerate() {
            let f = File::open(body_path)?;
            let file_len = f.metadata()?.len();
            let is_last = i == body_files.len() - 1;
            let bytes_to_copy = if is_last {
                file_len
            } else {
                file_len.saturating_sub(BGZF_EOF_LEN as u64)
            };
            let mut reader = BufReader::new(f).take(bytes_to_copy);
            std::io::copy(&mut reader, &mut out)?;
        }
    }
    Ok(())
}

fn rle_to_cigar(char_ops: &[char]) -> Cigar {
    if char_ops.is_empty() {
        return [Op::new(CigarKind::Match, 1)].into_iter().collect();
    }
    let mut ops: Vec<Op> = Vec::new();
    let mut current = char_ops[0];
    let mut count = 1usize;
    for &ch in &char_ops[1..] {
        if ch == current {
            count += 1;
        } else {
            ops.push(Op::new(char_to_cigar_kind(current), count));
            current = ch;
            count = 1;
        }
    }
    ops.push(Op::new(char_to_cigar_kind(current), count));
    ops.into_iter().collect()
}

fn char_to_cigar_kind(c: char) -> CigarKind {
    match c {
        'I' => CigarKind::Insertion,
        'D' => CigarKind::Deletion,
        _ => CigarKind::Match,
    }
}

/// Builds SAM flags for a simulated read.
///
/// For single-ended runs, returns `Flags::empty()`. For paired-end:
/// - R1 (forward): SEGMENTED | PROPERLY_SEGMENTED | MATE_REVERSE_COMPLEMENTED | FIRST_SEGMENT
/// - R2 (reverse): SEGMENTED | PROPERLY_SEGMENTED | REVERSE_COMPLEMENTED | LAST_SEGMENT
pub fn read_flags(is_paired_run: bool, is_reverse: bool) -> Flags {
    if !is_paired_run {
        return Flags::empty();
    }
    if is_reverse {
        Flags::SEGMENTED
            | Flags::PROPERLY_SEGMENTED
            | Flags::REVERSE_COMPLEMENTED
            | Flags::LAST_SEGMENT
    } else {
        Flags::SEGMENTED
            | Flags::PROPERLY_SEGMENTED
            | Flags::MATE_REVERSE_COMPLEMENTED
            | Flags::FIRST_SEGMENT
    }
}

/// Assembles a `RecordBuf` from alignment parameters.
///
/// `pos` and `mate_pos` are 0-based reference coordinates.
/// `qual_scores` are raw Phred values (not ASCII-offset).
pub fn build_bam_record(
    name: &str,
    ref_id: usize,
    pos: usize,
    flags: Flags,
    mapq: u8,
    cigar: Cigar,
    mate_ref_id: usize,
    mate_pos: usize,
    template_length: i32,
    sequence: &str,
    qual_scores: &[usize],
) -> RecordBuf {
    let mut record = RecordBuf::default();

    *record.name_mut() = Some(name.as_bytes().to_vec().into());
    *record.flags_mut() = flags;
    *record.reference_sequence_id_mut() = Some(ref_id);
    *record.alignment_start_mut() = Position::new(pos + 1);
    *record.mapping_quality_mut() = MappingQuality::try_from(mapq).ok();
    *record.cigar_mut() = cigar;
    *record.mate_reference_sequence_id_mut() = Some(mate_ref_id);
    *record.mate_alignment_start_mut() = Position::new(mate_pos + 1);
    *record.template_length_mut() = template_length;
    *record.sequence_mut() = Sequence::from(sequence.as_bytes());

    let qual_u8: Vec<u8> = qual_scores.iter().map(|&q| q as u8).collect();
    *record.quality_scores_mut() = QualityScores::from(qual_u8);

    record
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;
    use std::sync::Arc;
    use crate::structs::read_record::ReadRecord;
    use tempfile::tempdir;

    fn make_read(name: &str, pos: usize) -> ReadRecord {
        ReadRecord {
            name: name.to_string(),
            sequence: "ACGT".to_string(),
            quality_scores: vec![30, 30, 30, 30],
            cigar_ops: vec!['M', 'M', 'M', 'M'],
            is_paired: false,
            is_reverse: false,
            contig: "chr1".to_string(),
            position: pos,
            mate_contig: "chr1".to_string(),
            mate_position: 0,
            template_length: 0,
        }
    }

    fn read_bam_positions(path: &PathBuf) -> Vec<usize> {
        let file = std::fs::File::open(path).unwrap();
        let mut reader = bam::io::Reader::new(file);
        reader.read_header().unwrap();
        reader.records()
            .map(|r| {
                r.unwrap()
                    .alignment_start()
                    .unwrap()  // Option<Result<Position>> → Result<Position>
                    .unwrap()  // Result<Position>         → Position
                    .get() - 1
            })
            .collect()
    }

    #[test]
    fn test_flush_all_sorts_and_writes() {
        let temp = tempdir().unwrap();
        let bam_path = temp.path().join("sorted.bam");
        let contigs = vec![("chr1".to_string(), 100_000usize)];
        {
            let mut writer = BamWriter::new(&bam_path, &contigs).unwrap();
            writer.stage_read_record(&make_read("r3", 300)).unwrap();
            writer.stage_read_record(&make_read("r1", 100)).unwrap();
            writer.stage_read_record(&make_read("r2", 200)).unwrap();
            assert_eq!(writer.carry_len(), 3);
            writer.flush_all().unwrap();
            assert_eq!(writer.carry_len(), 0);
        } // drop writes BGZF EOF block
        assert_eq!(read_bam_positions(&bam_path), vec![100, 200, 300]);
    }

    #[test]
    fn test_flush_up_to_partitions_carry() {
        let temp = tempdir().unwrap();
        let bam_path = temp.path().join("windowed.bam");
        let contigs = vec![("chr1".to_string(), 100_000usize)];
        {
            let mut writer = BamWriter::new(&bam_path, &contigs).unwrap();
            writer.stage_read_record(&make_read("r1", 100)).unwrap();
            writer.stage_read_record(&make_read("r3", 900)).unwrap();
            writer.stage_read_record(&make_read("r2", 500)).unwrap();
            // flush_pos=600: 0-based pos < 600 → writes 100 and 500, carries 900
            writer.flush_up_to(600).unwrap();
            assert_eq!(writer.carry_len(), 1);
            writer.flush_all().unwrap();
            assert_eq!(writer.carry_len(), 0);
        }
        assert_eq!(read_bam_positions(&bam_path), vec![100, 500, 900]);
    }

    #[test]
    fn test_flush_up_to_empty_carry_is_noop() {
        let temp = tempdir().unwrap();
        let bam_path = temp.path().join("empty.bam");
        let contigs = vec![("chr1".to_string(), 1_000usize)];
        let mut writer = BamWriter::new(&bam_path, &contigs).unwrap();
        writer.flush_up_to(500).unwrap();
        writer.flush_all().unwrap();
        assert_eq!(writer.carry_len(), 0);
    }

    #[test]
    fn test_read_flags_single_ended() {
        assert_eq!(read_flags(false, false), Flags::empty());
        assert_eq!(read_flags(false, true), Flags::empty());
    }

    #[test]
    fn test_read_flags_paired_r1() {
        let flags = read_flags(true, false);
        assert!(flags.is_segmented());
        assert!(flags.contains(Flags::FIRST_SEGMENT));
        assert!(!flags.contains(Flags::REVERSE_COMPLEMENTED));
    }

    #[test]
    fn test_read_flags_paired_r2() {
        let flags = read_flags(true, true);
        assert!(flags.is_segmented());
        assert!(flags.contains(Flags::LAST_SEGMENT));
        assert!(flags.contains(Flags::REVERSE_COMPLEMENTED));
    }

    #[test]
    fn test_bam_writer_creates_file() {
        let temp = tempdir().unwrap();
        let bam_path: PathBuf = temp.path().join("test.bam");
        let contigs = vec![("chr1".to_string(), 100_000usize)];
        let mut writer = BamWriter::new(&bam_path, &contigs).unwrap();
        assert_eq!(writer.contig_id("chr1").unwrap(), 0);
        assert!(writer.contig_id("chrX").is_err());

        let record = ReadRecord {
            name: "read1/1".to_string(),
            sequence: "ACGT".to_string(),
            quality_scores: vec![30, 30, 30, 30],
            cigar_ops: vec!['M', 'M', 'M', 'M'],
            is_paired: false,
            is_reverse: false,
            contig: "chr1".to_string(),
            position: 100,
            mate_contig: "chr1".to_string(),
            mate_position: 0,
            template_length: 0,
        };
        writer.write_read_record(&record).unwrap();
        assert!(bam_path.exists());
    }

    #[test]
    fn test_concat_temp_bams_single_body_file() {
        // Exercises the branch where is_last=true fires on the first (and only) iteration,
        // i.e. the EOF block is kept intact rather than stripped.
        let temp = tempdir().unwrap();
        let contigs = vec![("chr1".to_string(), 100_000usize)];
        let context = Arc::new(BamContext::new(&contigs));

        let body = temp.path().join("body_only.bam");
        {
            let mut bw = BamBodyWriter::new(body.clone(), Arc::clone(&context)).unwrap();
            // Stage out-of-order to confirm flush_all sorts.
            bw.stage_read_record(&make_read("r2", 200)).unwrap();
            bw.stage_read_record(&make_read("r1", 100)).unwrap();
            bw.flush_all().unwrap();
        }

        let out = temp.path().join("single.bam");
        concat_temp_bams(&context, &[body], &out).unwrap();

        let positions = read_bam_positions(&out);
        assert_eq!(positions, vec![100, 200],
            "single-body concat: expected [100, 200], got: {:?}", positions);
    }

    #[test]
    fn test_concat_temp_bams_produces_readable_bam() {
        let temp = tempdir().unwrap();
        let contigs = vec![
            ("chr1".to_string(), 100_000usize),
            ("chr2".to_string(), 50_000usize),
        ];
        let context = Arc::new(BamContext::new(&contigs));

        // Write two body files
        let body1 = temp.path().join("body_chr1.bam");
        {
            let mut bw = BamBodyWriter::new(body1.clone(), Arc::clone(&context)).unwrap();
            bw.stage_read_record(&make_read("r1", 100)).unwrap();
            bw.stage_read_record(&make_read("r2", 200)).unwrap();
            bw.flush_all().unwrap();
        } // EOF written on drop

        let body2 = temp.path().join("body_chr2.bam");
        {
            let mut bw = BamBodyWriter::new(body2.clone(), Arc::clone(&context)).unwrap();
            bw.stage_read_record(&make_read("r3", 300)).unwrap();
            bw.flush_all().unwrap();
        }

        let out = temp.path().join("final.bam");
        concat_temp_bams(&context, &[body1, body2], &out).unwrap();

        // Verify the result is a readable BAM with 3 records
        let positions = read_bam_positions(&out);
        assert_eq!(positions, vec![100, 200, 300],
            "expected [100, 200, 300], got: {:?}", positions);
    }

    #[test]
    fn test_concat_temp_bams_empty_body_files() {
        let temp = tempdir().unwrap();
        let contigs = vec![("chr1".to_string(), 1_000usize)];
        let context = Arc::new(BamContext::new(&contigs));
        let out = temp.path().join("empty.bam");
        concat_temp_bams(&context, &[], &out).unwrap();

        let file = std::fs::File::open(&out).unwrap();
        let mut reader = bam::io::Reader::new(file);
        reader.read_header().unwrap();
        let count = reader.records().count();
        assert_eq!(count, 0);
    }
}
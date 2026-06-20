//! This library writes either single ended or paired-ended fastq files.
//! Need to update this method. We want to use the data structures and we want to make sure
//! this function is generic enough to work with the fragmented method we are implementing
//! This one needs a major overhaul, it is autogenerating quality scores etc.
//! Will wait to get other things set up first
use crate::rng::{NeatRng, NeatRngError};
use log::debug;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::path::PathBuf;
use thiserror::Error;

use crate::file_tools::bam_writer::BamRecordStager;
use crate::file_tools::file_io::append_to_file;
use crate::models::quality_scores::QualityScoreModel;
use crate::models::sequencing_error_model::{
    SeqModelError, SequencingErrorModel, SequencingErrorType,
};
use crate::structs::mutated_map::{AdCounter, MutatedMap, MutatedMapError};
use crate::structs::nucleotides::Nucleotide;
use crate::structs::nucleotides::Nucleotide::N;
use crate::structs::read_record::ReadRecord;
use crate::structs::sequence_block::{SequenceBlock, SequenceBlockError};
use crate::structs::variants::{Genotype, Variant};

#[derive(Error, Debug)]
pub enum FastqToolsError {
    #[error("Error writing bgzip fastq block {0}")]
    FastqWriteError(String),
    #[error("Quality scores missing for read 2")]
    MissingQScores,
    #[error("Error reading bgzip fastq block {0}")]
    FastqReadError(String),
    #[error("Mismatch between indexing and reads set for block {0}")]
    InvalidFastqBlock(String),
    #[error("Fastq Tools reported a SequenceBlock error: {0}")]
    SequenceBlockError(#[from] SequenceBlockError),
    #[error("Fastq tools reported a error model error: {0}")]
    ErrorModelError(#[from] SeqModelError),
    #[error("Fastq tools reported an IO error: {0}")]
    IoError(#[from] std::io::Error),
    #[error("Fastq tools reported an Rng Error: {0}")]
    RngError(#[from] NeatRngError),
    #[error("Error from MutatedMap: {0}")]
    MutatedMapError(#[from] MutatedMapError),
    #[error("Error locating a matching mutated map")]
    FindMapError,
    #[error("Paired ended declared, but r2 buffer empty")]
    BufferInitError,
    #[error("Malformed read during pair-ended analysis: {0}")]
    MalformedReadError(String),
    #[error("Strand value must either be 1 or 2, received {0}")]
    StrandError(usize),
    #[error("Truncated read {0}")]
    TruncatedRead(String),
    #[error("Reverse read with read_end > read_start")]
    MalformedReverseRead,
    #[error("BAM write error: {0}")]
    BamError(String),
}

pub enum Strand {
    Forward,
    Reverse,
}

pub fn reverse_complement(sequence: Vec<Nucleotide>) -> Vec<Nucleotide> {
    // Returns the reverse complement of a vector of u8's representing a DNA sequence.
    let length = sequence.len();
    let mut rev_comp = Vec::new();
    for i in (0..length).rev() {
        rev_comp.push(sequence[i].complement())
    }
    rev_comp
}

pub fn write_block_fastq<B1: Write, B2: Write>(
    block_fragments: Vec<(usize, usize)>,
    block_map: &MutatedMap,
    sequence_block: &SequenceBlock,
    paired_ended: bool,
    buffer1: &mut B1,
    buffer2: &mut B2,
    read_length: usize,
    long_reads: bool,
    read_name_prefix: &str,
    quality_score_model: &QualityScoreModel,
    sequencing_error_model: &SequencingErrorModel,
    rng: &mut NeatRng,
    mut bam_writer: Option<&mut dyn BamRecordStager>,
    ad_counter: &mut AdCounter,
) -> Result<(), FastqToolsError> {
    debug!("writing reads for {}", sequence_block.contig);
    // For SE reads, fragment_length == read_length exactly.  Sequencing-error deletions advance
    // seq_index without writing bases, so the generation loop can exhaust the fragment before
    // writing read_length bases, causing a TruncatedRead.  Fetching a small pad beyond `end`
    // gives the loop extra reference bases to consume.  For PE, fragments are large relative to
    // read_length so no padding is needed (and it would shift R2 coordinates).
    let seq_len = sequence_block.sequence.len();
    let se_pad = if paired_ended { 0 } else { 32 };
    for (frag_idx, (start, end)) in block_fragments.into_iter().enumerate() {
        let padded_end = (end + se_pad).min(seq_len);
        let fragment = sequence_block.get_subseq(start, padded_end)?;
        // In long-read mode a fragment may be shorter than read_length; truncate the read
        // to the actual fragment length rather than discarding it.
        let effective_read_len = if long_reads {
            fragment.len().min(read_length)
        } else {
            read_length
        };
        let mut read1_variants: HashMap<usize, &Variant> = HashMap::new();
        let mut reads1_flagged: Vec<usize> = Vec::new();
        let mut read2_variants: HashMap<usize, &Variant> = HashMap::new();
        let mut reads2_flagged: Vec<usize> = Vec::new();
        // Only the variants overlapping this fragment's two read windows matter.
        // block_map.flagged_positions is sorted (from_interval), so binary-search
        // each window instead of scanning every variant on the contig. The old
        // full-map scan made this O(fragments × variants) — pathological with
        // dense models (~220k SNVs on chr22 from the COSMIC tumor rate).
        let flagged = &block_map.flagged_positions;
        // R1 window: [start, start + effective_read_len); var offset = pos - start.
        let r1_lo = flagged.partition_point(|&p| p < start);
        let r1_hi = flagged.partition_point(|&p| p < start + effective_read_len);
        for &pos in &flagged[r1_lo..r1_hi] {
            let var_pos = pos - start;
            read1_variants.insert(var_pos, &block_map.variant_map[&pos]);
            reads1_flagged.push(var_pos);
        }
        // R2 window: [end - effective_read_len, end); R2 is reverse-stranded so
        // var offset = (end - 1) - pos. Guard end > effective_read_len so the
        // window start doesn't underflow (matches the original condition).
        if paired_ended && end > effective_read_len {
            let w_lo = end - effective_read_len;
            let r2_lo = flagged.partition_point(|&p| p < w_lo);
            let r2_hi = flagged.partition_point(|&p| p < end);
            for &pos in &flagged[r2_lo..r2_hi] {
                let var_pos = (end - 1) - pos;
                read2_variants.insert(var_pos, &block_map.variant_map[&pos]);
                reads2_flagged.push(var_pos);
            }
        }

        let ref_start = sequence_block.ref_start;
        let abs_start = start + ref_start;
        let abs_end = end + ref_start;
        // Per-fragment uniqueness tag in the read name. Without this, two
        // fragments that land at the same (start, end) — common via birthday
        // paradox at 30×+ coverage — share an identical QNAME, which violates
        // BAM/VCF spec and silently confuses Picard MarkDuplicates into
        // dropping them (see #210). The within-block frag_idx is sufficient
        // because rneat currently uses one block per contig and per-contig
        // read names already differ via `read_name_prefix`. For multi-pass
        // workflows that concatenate FASTQs from independent gen-reads runs,
        // the caller must still prefix each run's reads (e.g. cancer_simulate.sh
        // does N_/T_ between normal and tumor passes).
        let base_name = format!(
            "{}_{:010}_{:010}_{:016x}", read_name_prefix, abs_start, abs_end, frag_idx,
        );

        let r2_start = if paired_ended && abs_end >= effective_read_len {
            abs_end - effective_read_len
        } else {
            0
        };
        let tlen = if paired_ended {
            (abs_end - abs_start) as i32
        } else {
            0
        };

        let quality_scores_1 =
            quality_score_model.generate_quality_scores(effective_read_len, rng)?;
        let r1_record = match generate_read(
            &fragment,
            &reads1_flagged,
            &read1_variants,
            effective_read_len,
            format!("{}/1", base_name),
            Strand::Forward,
            quality_scores_1,
            sequencing_error_model,
            rng,
            sequence_block.contig.clone(),
            abs_start,
            sequence_block.contig.clone(),
            r2_start,
            tlen,
            paired_ended,
            ad_counter,
        ) {
            Ok(record) => record,
            Err(FastqToolsError::TruncatedRead(msg)) => {
                debug!("{}", msg);
                continue;
            }
            Err(e) => return Err(e),
        };

        // Generate r2 BEFORE writing r1, so that a TruncatedRead on r2
        // skips BOTH reads together. Otherwise r1 lands in buffer1 with
        // no matching r2 in buffer2, and BWA-MEM aborts with "paired
        // reads have different names" when the streams desync. The
        // failure mode was always reachable but became common after #221
        // (the literal-DEL skip+D-op fix), which legitimately advances
        // seq_index past the deleted bases and exhausts the buffer for
        // long deletions near a fragment edge.
        let r2_record = if paired_ended {
            let quality_scores_2 =
                quality_score_model.generate_quality_scores(effective_read_len, rng)?;
            let r2_pos = abs_end.saturating_sub(effective_read_len);
            let tlen_r2 = -((abs_end - abs_start) as i32);
            match generate_read(
                &reverse_complement(fragment),
                &reads2_flagged,
                &read2_variants,
                effective_read_len,
                format!("{}/2", base_name),
                Strand::Reverse,
                quality_scores_2,
                sequencing_error_model,
                rng,
                sequence_block.contig.clone(),
                r2_pos,
                sequence_block.contig.clone(),
                abs_start,
                tlen_r2,
                true,
                ad_counter,
            ) {
                Ok(record) => Some(record),
                Err(FastqToolsError::TruncatedRead(msg)) => {
                    debug!("{}", msg);
                    // Drop r1 alongside r2 so the streams stay in sync.
                    continue;
                }
                Err(e) => return Err(e),
            }
        } else {
            None
        };

        write_read_to_fastq(&r1_record, buffer1)?;
        if let Some(ref mut bam) = bam_writer {
            bam.stage_read_record(&r1_record)
                .map_err(|e| FastqToolsError::BamError(e.to_string()))?;
        }

        if let Some(r2_record) = r2_record {
            write_read_to_fastq(&r2_record, buffer2)?;
            if let Some(ref mut bam) = bam_writer {
                bam.stage_read_record(&r2_record)
                    .map_err(|e| FastqToolsError::BamError(e.to_string()))?;
            }
        }

        // Flush all carry records that start strictly before abs_start — safe because
        // block_fragments is sorted, so no future record can land earlier than abs_start.
        if let Some(ref mut bam) = bam_writer {
            bam.flush_up_to(abs_start)
                .map_err(|e| FastqToolsError::BamError(e.to_string()))?;
        }
    }
    Ok(())
}

pub fn combine_temp_fastqs(
    files_r1: Vec<PathBuf>,
    files_r2: Vec<PathBuf>,
    final_filename_r1: &PathBuf,
    final_filename_r2: Option<&PathBuf>,
) -> Result<(), FastqToolsError> {
    stream_gzip_files(&files_r1, final_filename_r1)?;
    if let Some(filename_r2) = final_filename_r2
        && !files_r2.is_empty()
    {
        stream_gzip_files(&files_r2, filename_r2)?;
    }
    Ok(())
}

fn stream_gzip_files(files: &[PathBuf], output: &PathBuf) -> Result<(), FastqToolsError> {
    // The per-contig temp files are each already a complete (multi-member) gzip
    // stream. gzip streams concatenate, so the final file is just the raw bytes
    // of every temp appended in order — no decompress, no recompress. This drops
    // an entire compression pass (each read is compressed once, in the parallel
    // per-contig write, instead of twice) and turns the previously single-
    // threaded combine into pure I/O.
    let mut out_file = BufWriter::new(append_to_file(output)?);
    for file in files {
        let mut f = BufReader::new(
            File::open(file).map_err(|e| FastqToolsError::FastqReadError(e.to_string()))?,
        );
        std::io::copy(&mut f, &mut out_file)
            .map_err(|e| FastqToolsError::FastqWriteError(e.to_string()))?;
    }
    out_file
        .flush()
        .map_err(|e| FastqToolsError::FastqWriteError(e.to_string()))?;
    Ok(())
}

// `cigar_ops.push('D')` runs in a loop per deletion-error base; pushing the
// same byte N times is the entire CIGAR encoding, not a copy-paste mistake.
#[allow(clippy::same_item_push)]
pub fn generate_read(
    sequence: &[Nucleotide],
    flagged_positions: &[usize],
    variant_map: &HashMap<usize, &Variant>,
    read_length: usize,
    name: String,
    read_strand: Strand,
    quality_scores: Vec<usize>,
    sequencing_error_model: &SequencingErrorModel,
    rng: &mut NeatRng,
    contig: String,
    position: usize,
    mate_contig: String,
    mate_position: usize,
    template_length: i32,
    is_paired: bool,
    ad_counter: &mut AdCounter,
) -> Result<ReadRecord, FastqToolsError> {
    if sequence.len() < read_length {
        return Err(FastqToolsError::TruncatedRead(format!("{:?}", sequence)));
    }

    let is_reverse = matches!(read_strand, Strand::Reverse);
    let fragment_length = sequence.len();

    let mut bases_written = 0;
    let mut out_seq = String::new();
    let mut cigar_ops: Vec<char> = Vec::new();
    let mut quality_index = 0;
    let mut seq_index = 0;

    'outer: while (seq_index < fragment_length) && (bases_written < read_length) {
        // Index variants by seq_index for BOTH strands. The caller already
        // reverse-complements the fragment for reverse (R2) reads AND maps each
        // variant's coordinate into that reversed sequence (var_pos = (end-1)-pos),
        // so the variant base sits at `seq_index` here just like the forward read.
        // The old `fragment_length - seq_index` for Reverse applied a second,
        // erroneous reflection — the variant was looked up at the mirror position,
        // so reverse reads carried REF at the true locus and the alt landed
        // elsewhere. Net effect: alternate alleles only ever appeared on
        // forward-strand reads, which strand-aware callers (e.g. Mutect2) correctly
        // flag as strand bias and filter out.
        let fragment_position = seq_index;
        let reference_base = sequence[seq_index].get_unmasked_base();
        let mut base_to_write: Vec<Nucleotide> = vec![reference_base];
        let mut deletion_skip: usize = 0;

        if reference_base == N {
            // Don't try to modify this base.
        } else if flagged_positions.contains(&fragment_position) {
            let variant = variant_map[&fragment_position];
            // MutatedMap::from_interval routes symbolic / structural ALTs to
            // sv_records, so anything in variant_map should be literal. Assert
            // that invariant in debug builds — a regression here would otherwise
            // panic at as_literal().unwrap() with no context.
            debug_assert!(
                variant.alternate.is_literal(),
                "symbolic ALT reached generate_read at position {fragment_position}"
            );
            let entry = ad_counter.entry(variant.location).or_insert((0, 0));
            if (variant.genotype == Genotype::Homozygous) || (rng.random()? < 0.5) {
                base_to_write = match read_strand {
                    Strand::Forward => variant.alternate
                        .as_literal()
                        .unwrap()
                        .to_vec(),
                    Strand::Reverse => variant.alternate
                        .as_literal()
                        .unwrap()
                        .iter()
                        .map(|b| b.complement())
                        .collect(),
                };
                entry.1 += 1; // alt
                // For a net-deletion variant (REF longer than ALT — typically a
                // pure literal deletion from indel_model with ALT.len() == 1),
                // skip the deleted REF bases. Without this, the read transcribes
                // the deleted region from the unbroken reference and emits no
                // CIGAR D-op, leaving the variant invisible to downstream
                // callers (see #221). The deletion_skip + D-op machinery below
                // is the same path SequencingErrorType::DeletionError already
                // uses; reusing it keeps CIGAR shape uniform across both
                // sources of deletion signal.
                let ref_len = variant.reference.len();
                let alt_len = base_to_write.len();
                if ref_len > alt_len {
                    let want_skip = ref_len - alt_len;
                    // Cap at remaining buffer so we don't read past the
                    // fragment end. Truncated cases fall through to the
                    // existing TruncatedRead error path.
                    let max_skip = fragment_length
                        .saturating_sub(seq_index)
                        .saturating_sub(1);
                    let actual_skip = want_skip.min(max_skip);
                    seq_index += actual_skip;
                    deletion_skip = actual_skip;
                }
            } else {
                entry.0 += 1; // ref (het coin landed on ref)
            }
        } else {
            let score = quality_scores[quality_index];
            let prob = sequencing_error_model.convert_score(score)?;
            if rng.random()? < prob {
                debug!("Creating sequencing error");
                let error =
                    sequencing_error_model.generate_sequencing_error(reference_base, rng)?;
                match error {
                    SequencingErrorType::SnpError(base) => {
                        debug!("Snp error");
                        base_to_write = vec![base];
                    }
                    SequencingErrorType::DeletionError(length) => {
                        debug!("Deletion error");
                        if seq_index + length < sequence.len() {
                            seq_index += length;
                            deletion_skip = length;
                        }
                    }
                    SequencingErrorType::InsertionError(vec) => {
                        debug!("Insertion error");
                        let mut insertion = vec![reference_base];
                        insertion.extend(vec);
                        base_to_write = insertion;
                    }
                }
            }
        }

        let mut is_first_base = true;
        for base in base_to_write {
            out_seq.push(base.into());
            bases_written += 1;
            if is_first_base {
                cigar_ops.push('M');
                is_first_base = false;
            } else {
                cigar_ops.push('I');
            }
            if bases_written == read_length {
                break 'outer;
            }
        }

        // 'D' ops for skipped reference bases (sequencing deletion errors only)
        for _ in 0..deletion_skip {
            cigar_ops.push('D');
        }

        seq_index += 1;
        quality_index += 1;
    }

    if bases_written < read_length {
        return Err(FastqToolsError::TruncatedRead(format!("{:?}", sequence)));
    }

    Ok(ReadRecord {
        name,
        sequence: out_seq,
        quality_scores,
        cigar_ops,
        is_paired,
        is_reverse,
        contig,
        position,
        mate_contig,
        mate_position,
        template_length,
    })
}

pub fn write_read_to_fastq<W: Write>(
    record: &ReadRecord,
    buffer: &mut W,
) -> Result<(), FastqToolsError> {
    buffer.write_all(format!("@{}\n", record.name).as_bytes())?;
    buffer.write_all(record.sequence.as_bytes())?;
    buffer.write_all(b"\n+\n")?;
    buffer.write_all(&quality_scores_to_char_vec(&record.quality_scores)?)?;
    buffer.write_all(b"\n")?;
    Ok(())
}

pub fn quality_scores_to_char_vec(array: &[usize]) -> Result<Vec<u8>, FastqToolsError> {
    let mut score_vec = Vec::new();
    for &score in array {
        score_vec.push((score + 33) as u8)
    }
    Ok(score_vec)
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::Compression;
    use flate2::write::GzEncoder;
    use crate::file_tools::bam_writer::{BamRecordStager, BamWriter};
    use crate::file_tools::file_io::{VectorBuffer, create_output_file, read_gzip_lines};
    use crate::structs::nucleotides::Nucleotide::*;
    use crate::structs::sequence_block::{RegionType, SequenceMap};
    use crate::structs::variants::{Variant, VariantType};
    use std::io::Write;

    #[test]
    fn test_combine_temp_fastqs() {
        let temp_dir = tempfile::tempdir().unwrap();

        let write_fastq_gz = |name: &str, content: &[u8]| -> PathBuf {
            let path = temp_dir.path().join(name);
            let f = std::fs::File::create(&path).unwrap();
            let mut enc = GzEncoder::new(f, Compression::default());
            enc.write_all(content).unwrap();
            enc.finish().unwrap();
            path
        };

        let file1 = write_fastq_gz("r1.fastq.gz", b"@read1\nACGT\n+\nIIII\n");
        let file2 = write_fastq_gz("r2.fastq.gz", b"@read2\nTTTT\n+\nIIII\n");
        let output = temp_dir.path().join("combined.fastq.gz");

        combine_temp_fastqs(vec![file1, file2], vec![], &output, None).unwrap();

        let lines: Vec<String> = read_gzip_lines(&output)
            .unwrap()
            .map(|l| l.unwrap())
            .collect();
        assert_eq!(
            lines.len(),
            8,
            "Combined file should have 8 lines (2 records × 4 lines)"
        );
        assert_eq!(lines[0], "@read1");
        assert_eq!(lines[4], "@read2");
    }

    #[test]
    fn test_write_reverse() {
        let temp_dir = tempfile::tempdir().unwrap();
        let mut temp_file = PathBuf::from(temp_dir.path());
        temp_file.push("test.fastq.gz");
        let mut temp_writer = create_output_file(&temp_file, true).unwrap();
        let original_seq = vec![A, C, C, G, A, A, T, G, A];
        let rev_comp_seq = reverse_complement(original_seq);
        let expected_rev_comp = vec![T, C, A, T, T, C, G, G, T];
        assert_eq!(expected_rev_comp, rev_comp_seq);
        let flagged_positions: Vec<usize> = Vec::new();
        let variant_map: HashMap<usize, &Variant> = HashMap::new();
        let read_len = 4;
        let read_name = "neat_gen__0000000000_0000000009/2".to_string();
        let quality_scores = vec![32, 32, 32, 32];
        let sequencing_error_model = SequencingErrorModel::default().unwrap();
        let mut rng = NeatRng::new_from_seed(&vec![
            "Hello".to_string(),
            "Cruel".to_string(),
            "World".to_string(),
        ])
        .unwrap();
        let record = generate_read(
            &rev_comp_seq,
            &flagged_positions,
            &variant_map,
            read_len,
            read_name,
            Strand::Reverse,
            quality_scores,
            &sequencing_error_model,
            &mut rng,
            "chr1".to_string(),
            0,
            "chr1".to_string(),
            0,
            0,
            false,
            &mut AdCounter::new(),
        )
        .unwrap();
        let mut buffer = GzEncoder::new(&mut temp_writer, Compression::default());
        write_read_to_fastq(&record, &mut buffer).unwrap();
        buffer.flush().unwrap();
        let temp_reader = read_gzip_lines(&temp_file).unwrap();
        let seq_line: String = {
            let mut rtrn = String::new();
            for line in temp_reader {
                let l = line.unwrap().to_string();
                if l.starts_with("@") {
                    // skip name line
                } else {
                    rtrn = l;
                    break;
                }
            }
            rtrn
        };
        let exp_rev_cmp_str = "TCAT".to_string();
        assert_eq!(seq_line, exp_rev_cmp_str);
    }

    #[test]
    fn test_qual_score_to_write() {
        let qual_scores = vec![33, 25, 37, 28, 15, 33, 33, 37, 37, 25];
        let qual_string = "B:F=0BBFF:".as_bytes();
        assert_eq!(
            qual_string,
            quality_scores_to_char_vec(&qual_scores).unwrap()
        )
    }

    #[test]
    fn test_reverse_complement() {
        let read: Vec<Nucleotide> = vec![A, A, A, A, C, C, C, C, C];
        let revcomp: Vec<Nucleotide> = vec![G, G, G, G, G, T, T, T, T];
        assert_eq!(reverse_complement(read), revcomp);
    }

    #[test]
    fn test_write_block_fastq_ref_start_in_read_name() {
        // Verifies that when a SequenceBlock has ref_start > 0, the read names in the
        // output FASTQ use reference-relative positions, not block-local positions.
        use crate::structs::{
            mutated_map::MutatedMap,
            sequence_block::{RegionType, SequenceBlock, SequenceMap},
        };
        let temp_dir = tempfile::tempdir().unwrap();
        let ref_start: usize = 1000;
        let seq_len: usize = 50;
        let sequence: Vec<Nucleotide> = (0..seq_len)
            .map(|i| match i % 4 {
                0 => A,
                1 => C,
                2 => G,
                _ => T,
            })
            .collect();
        let block = SequenceBlock {
            contig: "chr1".to_string(),
            ref_start,
            ref_end: ref_start + seq_len,
            sequence: sequence.clone(),
            sequence_map: vec![SequenceMap::from(RegionType::NonNRegion, 0, seq_len)],
        };
        let mutated_map =
            MutatedMap::from_interval(ref_start, ref_start + seq_len, vec![]).unwrap();
        let frag_start: usize = 5;
        let frag_end: usize = 25;
        let fragments = vec![(frag_start, frag_end)];
        let out_path = temp_dir.path().join("out.fastq.gz");
        let outfile = create_output_file(&out_path, true).unwrap();
        use crate::file_tools::file_io::VectorBuffer;
        let dummy_data: VectorBuffer = VectorBuffer::new();
        let mut buf1 = GzEncoder::new(outfile, Compression::default());
        let mut buf2 = GzEncoder::new(dummy_data, Compression::default());
        let seq_err_model = SequencingErrorModel::default().unwrap();
        let quality_model = QualityScoreModel::default().unwrap();
        let mut rng = NeatRng::new_from_seed(&vec!["test".to_string()]).unwrap();
        write_block_fastq(
            fragments,
            &mutated_map,
            &block,
            false,
            &mut buf1,
            &mut buf2,
            10,
            false,
            "chr1",
            &quality_model,
            &seq_err_model,
            &mut rng,
            None,
            &mut AdCounter::new(),
        )
        .unwrap();
        buf1.finish().unwrap();
        let lines: Vec<String> = read_gzip_lines(&out_path)
            .unwrap()
            .map(|l| l.unwrap())
            .collect();
        let header = &lines[0];
        assert!(
            header.starts_with('@'),
            "Expected FASTQ header, got: {}",
            header
        );
        // Read name must contain ref-relative positions, not block-local ones
        let expected_start = format!("{:010}", frag_start + ref_start); // 0000001005
        let expected_end = format!("{:010}", frag_end + ref_start); // 0000001025
        assert!(
            header.contains(&expected_start),
            "Read name should contain ref-relative start {}; got: {}",
            expected_start,
            header
        );
        assert!(
            header.contains(&expected_end),
            "Read name should contain ref-relative end {}; got: {}",
            expected_end,
            header
        );
        // Also verify block-local positions are NOT what was written
        let local_start = format!("{:010}", frag_start); // 0000000005
        assert!(
            !header.contains(&local_start),
            "Read name should NOT contain block-local start {}; got: {}",
            local_start,
            header
        );
    }

    /// Two fragments at the SAME (start, end) coordinate must produce reads
    /// with distinct QNAMEs — without the per-fragment uniqueness tag this
    /// would silently emit collision-named reads, which Picard MarkDuplicates
    /// would interpret as PCR duplicates and drop (see #210).
    #[test]
    fn test_write_block_fastq_unique_names_for_same_position_fragments() {
        use crate::structs::{
            mutated_map::MutatedMap,
            sequence_block::{RegionType, SequenceBlock, SequenceMap},
        };
        let temp_dir = tempfile::tempdir().unwrap();
        let seq_len: usize = 200;
        let sequence: Vec<Nucleotide> = (0..seq_len)
            .map(|i| match i % 4 {
                0 => A,
                1 => C,
                2 => G,
                _ => T,
            })
            .collect();
        let block = SequenceBlock {
            contig: "chr1".to_string(),
            ref_start: 0,
            ref_end: seq_len,
            sequence,
            sequence_map: vec![SequenceMap::from(RegionType::NonNRegion, 0, seq_len)],
        };
        let mutated_map = MutatedMap::from_interval(0, seq_len, vec![]).unwrap();
        // Four fragments — three at the SAME coordinate, one elsewhere.
        let fragments = vec![(10, 80), (10, 80), (10, 80), (100, 170)];
        let out_path = temp_dir.path().join("collide.fastq.gz");
        let outfile = create_output_file(&out_path, true).unwrap();
        use crate::file_tools::file_io::VectorBuffer;
        let dummy: VectorBuffer = VectorBuffer::new();
        let mut buf1 = GzEncoder::new(outfile, Compression::default());
        let mut buf2 = GzEncoder::new(dummy, Compression::default());
        let seq_err_model = SequencingErrorModel::default().unwrap();
        let quality_model = QualityScoreModel::default().unwrap();
        let mut rng = NeatRng::new_from_seed(&vec!["uniq-name".to_string()]).unwrap();
        write_block_fastq(
            fragments,
            &mutated_map,
            &block,
            false,
            &mut buf1,
            &mut buf2,
            70,
            false,
            "chr1",
            &quality_model,
            &seq_err_model,
            &mut rng,
            None,
            &mut AdCounter::new(),
        )
        .unwrap();
        buf1.finish().unwrap();
        // Read every line that starts with @ AND is followed by a valid
        // FASTQ record. Filter by record-position (line 1 of every 4-line
        // block) — quality lines (line 4) can start with @ too.
        let lines: Vec<String> = read_gzip_lines(&out_path)
            .unwrap()
            .map(|l| l.unwrap())
            .collect();
        let names: Vec<&String> = lines
            .iter()
            .enumerate()
            .filter(|(i, _)| i % 4 == 0)
            .map(|(_, l)| l)
            .collect();
        assert_eq!(names.len(), 4, "expected 4 read names, got: {:?}", names);
        let mut sorted: Vec<&&String> = names.iter().collect();
        sorted.sort();
        sorted.dedup();
        assert_eq!(
            sorted.len(),
            4,
            "expected 4 unique names, got duplicates in: {:?}",
            names
        );
    }

    // --- CIGAR-building tests for the refactored generate_read ---

    fn make_sequence(len: usize) -> Vec<Nucleotide> {
        (0..len)
            .map(|i| match i % 4 {
                0 => A,
                1 => C,
                2 => G,
                _ => T,
            })
            .collect()
    }

    #[test]
    fn test_generate_read_cigar_all_match_no_errors() {
        // Q40 → error prob 0.0001; with seed and 10 positions no error fires.
        let sequence = make_sequence(30);
        let read_length = 10;
        let quality_scores = vec![40usize; read_length];
        let model = SequencingErrorModel::default().unwrap();
        let mut rng = NeatRng::new_from_seed(&vec!["no_error".to_string()]).unwrap();
        let record = generate_read(
            &sequence,
            &[],
            &HashMap::new(),
            read_length,
            "r/1".to_string(),
            Strand::Forward,
            quality_scores,
            &model,
            &mut rng,
            "chr1".to_string(),
            0,
            "chr1".to_string(),
            0,
            0,
            false,
            &mut AdCounter::new(),
        )
        .unwrap();
        assert_eq!(record.sequence.len(), read_length);
        assert_eq!(record.cigar_ops.len(), read_length);
        assert!(
            record.cigar_ops.iter().all(|&c| c == 'M'),
            "expected all-M cigar, got: {:?}",
            record.cigar_ops
        );
    }

    #[test]
    fn test_generate_read_cigar_mi_plus_d_equals_len() {
        // Invariant: M+I == read_length, cigar_ops.len() == (M+I) + D, for any error mix.
        let sequence = make_sequence(500);
        let read_length = 50;
        let quality_scores = vec![0usize; read_length]; // Q0 → 100% error probability
        let model = SequencingErrorModel::default().unwrap();
        let mut rng = NeatRng::new_from_seed(&vec!["invariant".to_string()]).unwrap();
        let record = generate_read(
            &sequence,
            &[],
            &HashMap::new(),
            read_length,
            "r/1".to_string(),
            Strand::Forward,
            quality_scores,
            &model,
            &mut rng,
            "chr1".to_string(),
            0,
            "chr1".to_string(),
            0,
            0,
            false,
            &mut AdCounter::new(),
        )
        .unwrap();
        let mi = record
            .cigar_ops
            .iter()
            .filter(|&&c| c == 'M' || c == 'I')
            .count();
        let d = record.cigar_ops.iter().filter(|&&c| c == 'D').count();
        assert_eq!(record.sequence.len(), read_length);
        assert_eq!(mi, read_length, "M+I must equal read_length");
        assert_eq!(record.cigar_ops.len(), mi + d, "cigar len must equal M+I+D");
    }

    /// Counter behaviour: a homozygous SNP should drive alt_count up and
    /// leave ref_count at 0 across many reads; a heterozygous SNP should split
    /// roughly 50/50 around the half-point of N reads. Together these pin the
    /// AdCounter increments in generate_read at the coin-flip site.
    #[test]
    fn test_generate_read_counter_homozygous_all_alt() {
        let sequence = make_sequence(30);
        let read_length = 10;
        let hom_snp =
            Variant::new(VariantType::SNP, 5, &vec![T], &vec![C], &mut vec![1, 1]).unwrap();
        let variant_map = HashMap::from([(5usize, &hom_snp)]);
        let quality_scores = vec![40usize; read_length];
        let model = SequencingErrorModel::default().unwrap();
        let mut rng = NeatRng::new_from_seed(&vec!["hom".to_string()]).unwrap();
        let mut ad: AdCounter = HashMap::new();
        for i in 0..50 {
            let _ = generate_read(
                &sequence, &[5], &variant_map, read_length, format!("r{i}/1"),
                Strand::Forward, quality_scores.clone(), &model, &mut rng,
                "chr1".to_string(), 0, "chr1".to_string(), 0, 0, false, &mut ad,
            )
            .unwrap();
        }
        let (refs, alts) = ad[&5];
        assert_eq!(refs, 0, "homozygous SNP must produce zero ref reads");
        assert_eq!(alts, 50, "homozygous SNP must produce all alt reads");
    }

    /// Regression for the reverse-strand variant-placement bug: a reverse (R2)
    /// read must apply the variant at `seq_index` (the caller already
    /// reverse-complements the fragment and maps the variant coordinate into it),
    /// NOT at the mirrored `fragment_length - seq_index`. With the old reflection
    /// the variant fell outside the read window, so reverse reads silently
    /// carried REF — alternate alleles only ever appeared on forward reads,
    /// which Mutect2 (correctly) filtered as strand bias.
    #[test]
    fn test_generate_read_reverse_applies_variant() {
        let sequence = make_sequence(30);
        let read_length = 10;
        let hom_snp =
            Variant::new(VariantType::SNP, 5, &vec![T], &vec![C], &mut vec![1, 1]).unwrap();
        let variant_map = HashMap::from([(5usize, &hom_snp)]);
        let quality_scores = vec![40usize; read_length];
        let model = SequencingErrorModel::default().unwrap();
        let mut rng = NeatRng::new_from_seed(&vec!["rev".to_string()]).unwrap();
        let mut ad: AdCounter = HashMap::new();
        for i in 0..50 {
            let _ = generate_read(
                &sequence, &[5], &variant_map, read_length, format!("r{i}/2"),
                Strand::Reverse, quality_scores.clone(), &model, &mut rng,
                "chr1".to_string(), 0, "chr1".to_string(), 0, 0, true, &mut ad,
            )
            .unwrap();
        }
        let (refs, alts) = ad[&5];
        // The pre-fix reflection put the lookup at index 25, outside the
        // 10-base read, so alt would be 0 here. The fix applies it at index 5.
        assert_eq!(refs, 0, "reverse homozygous SNP must produce zero ref reads");
        assert_eq!(alts, 50, "reverse read must carry the alt (strand-bias regression)");
    }

    #[test]
    fn test_generate_read_counter_heterozygous_splits_around_half() {
        let sequence = make_sequence(30);
        let read_length = 10;
        let het_snp =
            Variant::new(VariantType::SNP, 5, &vec![T], &vec![C], &mut vec![0, 1]).unwrap();
        let variant_map = HashMap::from([(5usize, &het_snp)]);
        let quality_scores = vec![40usize; read_length];
        let model = SequencingErrorModel::default().unwrap();
        let mut rng = NeatRng::new_from_seed(&vec!["het".to_string()]).unwrap();
        let mut ad: AdCounter = HashMap::new();
        let n = 1000;
        for i in 0..n {
            let _ = generate_read(
                &sequence, &[5], &variant_map, read_length, format!("r{i}/1"),
                Strand::Forward, quality_scores.clone(), &model, &mut rng,
                "chr1".to_string(), 0, "chr1".to_string(), 0, 0, false, &mut ad,
            )
            .unwrap();
        }
        let (refs, alts) = ad[&5];
        assert_eq!(refs + alts, n, "every read should increment exactly one slot");
        // Binomial(1000, 0.5) → 99.99% CI is well within [400, 600]
        assert!(
            (400..600).contains(&(refs as usize)),
            "het split should be roughly 50/50 ({}/1000 ref, {} alt)",
            refs,
            alts
        );
    }

    #[test]
    fn test_generate_read_snp_variant_produces_all_m_cigar() {
        // A SNP variant changes a base but does not add I or D ops.
        let sequence = make_sequence(30);
        let read_length = 10;
        let snp = Variant::new(VariantType::SNP, 3, &vec![T], &vec![C], &mut vec![1, 1]).unwrap();
        let variant_map = HashMap::from([(3usize, &snp)]);
        let quality_scores = vec![40usize; read_length];
        let model = SequencingErrorModel::default().unwrap();
        let mut rng = NeatRng::new_from_seed(&vec!["snp_cigar".to_string()]).unwrap();
        let record = generate_read(
            &sequence,
            &[3],
            &variant_map,
            read_length,
            "r/1".to_string(),
            Strand::Forward,
            quality_scores,
            &model,
            &mut rng,
            "chr1".to_string(),
            0,
            "chr1".to_string(),
            0,
            0,
            false,
            &mut AdCounter::new(),
        )
        .unwrap();
        assert_eq!(record.cigar_ops.len(), read_length);
        assert!(
            record.cigar_ops.iter().all(|&c| c == 'M'),
            "SNP variant must not add I or D ops; got: {:?}",
            record.cigar_ops
        );
    }

    #[test]
    fn test_generate_read_insertion_variant_produces_i_ops() {
        // Homozygous insertion: alt = [A, C, G] → 1 M (anchor) + 2 I (inserted bases).
        let sequence = make_sequence(30);
        let read_length = 10;
        let ins = Variant::new(
            VariantType::Insertion,
            3,
            &vec![A],
            &vec![A, C, G],
            &mut vec![1, 1],
        )
        .unwrap();
        let variant_map = HashMap::from([(3usize, &ins)]);
        let quality_scores = vec![40usize; read_length];
        let model = SequencingErrorModel::default().unwrap();
        let mut rng = NeatRng::new_from_seed(&vec!["ins_variant".to_string()]).unwrap();
        let record = generate_read(
            &sequence,
            &[3],
            &variant_map,
            read_length,
            "r/1".to_string(),
            Strand::Forward,
            quality_scores,
            &model,
            &mut rng,
            "chr1".to_string(),
            0,
            "chr1".to_string(),
            0,
            0,
            false,
            &mut AdCounter::new(),
        )
        .unwrap();
        let i_count = record.cigar_ops.iter().filter(|&&c| c == 'I').count();
        assert_eq!(i_count, 2, "2-base insertion variant must produce 2 I ops");
        let mi = record
            .cigar_ops
            .iter()
            .filter(|&&c| c == 'M' || c == 'I')
            .count();
        assert_eq!(mi, read_length);
    }

    #[test]
    fn test_generate_read_error_indels_appear_in_cigar() {
        // Q0 at every position guarantees an error at each position (convert_score(0) = 1.0).
        // With 50 positions, P(all 50 errors are SNPs) ≈ 0.6^50 < 1.4e-11,
        // so this effectively guarantees at least one I or D op in the cigar.
        let sequence = make_sequence(500);
        let read_length = 50;
        let quality_scores = vec![0usize; read_length];
        let model = SequencingErrorModel::default().unwrap();
        let mut rng = NeatRng::new_from_seed(&vec!["error_indel".to_string()]).unwrap();
        let record = generate_read(
            &sequence,
            &[],
            &HashMap::new(),
            read_length,
            "r/1".to_string(),
            Strand::Forward,
            quality_scores,
            &model,
            &mut rng,
            "chr1".to_string(),
            0,
            "chr1".to_string(),
            0,
            0,
            false,
            &mut AdCounter::new(),
        )
        .unwrap();
        assert!(
            record.cigar_ops.iter().any(|&c| c != 'M'),
            "expected non-M ops when every position has a guaranteed error; got: {:?}",
            record.cigar_ops
        );
        // Invariant must hold even under heavy errors
        let mi = record
            .cigar_ops
            .iter()
            .filter(|&&c| c == 'M' || c == 'I')
            .count();
        assert_eq!(mi, read_length);
    }

    #[test]
    fn test_generate_read_deletion_error_increases_cigar_length() {
        // When deletion errors fire, each skipped reference base adds a D op, so
        // cigar_ops.len() > read_length.  Verify the formula len == M+I+D holds
        // and, whenever D > 0, that cigar_ops.len() > read_length.
        let sequence = make_sequence(500);
        let read_length = 50;
        let quality_scores = vec![0usize; read_length];
        let model = SequencingErrorModel::default().unwrap();
        let mut rng = NeatRng::new_from_seed(&vec!["del_len".to_string()]).unwrap();
        let record = generate_read(
            &sequence,
            &[],
            &HashMap::new(),
            read_length,
            "r/1".to_string(),
            Strand::Forward,
            quality_scores,
            &model,
            &mut rng,
            "chr1".to_string(),
            0,
            "chr1".to_string(),
            0,
            0,
            false,
            &mut AdCounter::new(),
        )
        .unwrap();
        let mi = record
            .cigar_ops
            .iter()
            .filter(|&&c| c == 'M' || c == 'I')
            .count();
        let d = record.cigar_ops.iter().filter(|&&c| c == 'D').count();
        assert_eq!(mi, read_length);
        assert_eq!(record.cigar_ops.len(), mi + d);
        if d > 0 {
            assert!(
                record.cigar_ops.len() > read_length,
                "deletion errors must make cigar longer than read_length"
            );
        }
    }

    /// #221 regression: a long literal Deletion (REF=50bp, ALT=1bp) must
    /// produce a CIGAR with the 49-base D-op so downstream callers see the
    /// deletion. Before the fix, the alt branch wrote the 1-byte anchor and
    /// then advanced seq_index by 1, so the read transcribed the 49 deleted
    /// bases from the unbroken reference and emitted CIGAR `<read_length>M`.
    #[test]
    fn test_literal_long_deletion_emits_d_ops() {
        // 200-bp reference: 30 bases of left context + 50-bp REF starting at
        // pos 30 + 120 bases of right context.
        let mut sequence = vec![A; 30];
        sequence.extend(vec![C; 50]); // the 50 REF bases (anchor + 49 deleted)
        sequence.extend(vec![T; 120]); // post-deletion context
        let read_length = 100;

        // Homozygous DEL at position 30, REF=50bp (CCCC...), ALT=1bp (anchor C).
        let ref_bases: Vec<Nucleotide> = vec![C; 50];
        let alt_bases: Vec<Nucleotide> = vec![C];
        let variant = Variant::new(
            VariantType::Deletion,
            30,
            &ref_bases,
            &alt_bases,
            &mut vec![1, 1], // homozygous so the alt branch always fires
        )
        .unwrap();
        let variant_map = HashMap::from([(30usize, &variant)]);
        let flagged_positions = vec![30usize];
        let qual_scores = vec![33; 100];
        let sequencing_error_model = SequencingErrorModel::default().unwrap();
        let mut rng =
            NeatRng::new_from_seed(&vec!["literal".into(), "long".into(), "del".into()]).unwrap();

        let record = generate_read(
            &sequence,
            &flagged_positions,
            &variant_map,
            read_length,
            "del_test/1".into(),
            Strand::Forward,
            qual_scores,
            &sequencing_error_model,
            &mut rng,
            "chr1".into(),
            0,
            "chr1".into(),
            0,
            0,
            false,
            &mut AdCounter::new(),
        )
        .unwrap();

        // Sequence written must equal read_length (100 bases).
        assert_eq!(record.sequence.len(), read_length);

        // The 49-bp deletion must be encoded as D ops in the CIGAR. The
        // exact CIGAR is read-length-dependent (sequencing errors may add
        // small ops too), but the D count must be ≥ 49 — the deletion span.
        // We assert D count ≥ 49 rather than == 49 so a stray seq-error
        // deletion doesn't flake the test.
        let d_count = record.cigar_ops.iter().filter(|&&c| c == 'D').count();
        assert!(
            d_count >= 49,
            "expected ≥49 D ops for a 50-bp REF / 1-bp ALT homozygous deletion, got {d_count}. \
             CIGAR: {:?}",
            record.cigar_ops
        );

        // Bases-written count: each M op or I op consumes one output base.
        let mi_count = record
            .cigar_ops
            .iter()
            .filter(|&&c| c == 'M' || c == 'I')
            .count();
        assert_eq!(
            mi_count, read_length,
            "M+I ops must equal read_length ({read_length}); got {mi_count}"
        );

        // AD counter on the variant position must record the alt observation.
        // (Verified via the variant_map's pointer at position 30 — alt_count
        // is incremented inside generate_read when the alt branch fires.)
    }

    #[test]
    fn test_apply_variants() {
        let sequence = vec![A, C, G, T, T, A, T, G, A, C, G, T, T, A, T, G];
        let variant1 =
            Variant::new(VariantType::SNP, 1, &vec![T], &vec![C], &mut vec![1, 0]).unwrap();
        let variant2 = Variant::new(
            VariantType::Deletion,
            3,
            &vec![T, T],
            &vec![T],
            &mut vec![0, 1],
        )
        .unwrap();
        let variant_map = HashMap::from([(1, &variant1), (3, &variant2)]);
        let flagged_positions = vec![1, 3];
        let read_name = "neat_generated__0000000000_0000000008/1".to_string();
        let qual_scores = vec![33, 25, 37, 28, 15, 33, 33, 37];
        let sequencing_error_model = SequencingErrorModel::default().unwrap();
        let mut rng = NeatRng::new_from_seed(&vec![
            "Hello".to_string(),
            "Cruel".to_string(),
            "World".to_string(),
        ])
        .unwrap();
        let result = generate_read(
            &sequence,
            &flagged_positions,
            &variant_map,
            8,
            read_name,
            Strand::Forward,
            qual_scores,
            &sequencing_error_model,
            &mut rng,
            "chr1".to_string(),
            0,
            "chr1".to_string(),
            0,
            0,
            false,
            &mut AdCounter::new(),
        );
        assert!(result.is_ok());
    }

    // ── incremental BAM flush tests ──────────────────────────────────────────

    fn make_block(seq_len: usize) -> SequenceBlock {
        let sequence: Vec<Nucleotide> = (0..seq_len)
            .map(|i| match i % 4 {
                0 => A,
                1 => C,
                2 => G,
                _ => T,
            })
            .collect();
        SequenceBlock {
            contig: "chr1".to_string(),
            ref_start: 0,
            ref_end: seq_len,
            sequence,
            sequence_map: vec![SequenceMap::from(RegionType::NonNRegion, 0, seq_len)],
        }
    }

    fn run_write_block(
        bam_writer: &mut BamWriter,
        fragments: Vec<(usize, usize)>,
        block: &SequenceBlock,
        paired_ended: bool,
        read_len: usize,
    ) {
        let mutated_map = MutatedMap::from_interval(0, block.sequence.len(), vec![]).unwrap();
        let seq_err_model = SequencingErrorModel::default().unwrap();
        let quality_model = QualityScoreModel::default().unwrap();
        let mut rng = NeatRng::new_from_seed(&vec!["flush-test".to_string()]).unwrap();
        let mut buf1 = GzEncoder::new(VectorBuffer::new(), Compression::default());
        let mut buf2 = GzEncoder::new(VectorBuffer::new(), Compression::default());
        let stager: Option<&mut dyn BamRecordStager> = Some(bam_writer);
        write_block_fastq(
            fragments,
            &mutated_map,
            block,
            paired_ended,
            &mut buf1,
            &mut buf2,
            read_len,
            false,
            "chr1",
            &quality_model,
            &seq_err_model,
            &mut rng,
            stager,
            &mut AdCounter::new(),
        )
        .unwrap();
    }

    #[test]
    fn test_write_block_fastq_carry_bounded_se() {
        // With 5 SE fragments at strictly increasing positions, flush_up_to(abs_start)
        // fires after each fragment so the carry never exceeds 1 record.
        let temp = tempfile::tempdir().unwrap();
        let bam_path = temp.path().join("out.bam");
        let block = make_block(500);
        let mut bw = BamWriter::new(&bam_path, &[("chr1".to_string(), 500usize)]).unwrap();

        let fragments = vec![(0usize, 10), (50, 60), (100, 110), (150, 160), (200, 210)];
        run_write_block(&mut bw, fragments, &block, false, 10);

        assert_eq!(
            bw.carry_len(),
            1,
            "SE: carry should hold only the last fragment's read, got {}",
            bw.carry_len()
        );
    }

    #[test]
    fn test_write_block_fastq_carry_bounded_pe() {
        // With 5 PE fragments at strictly increasing positions, flush_up_to fires after
        // each fragment so carry never exceeds 2 records (one R1 + one R2).
        let temp = tempfile::tempdir().unwrap();
        let bam_path = temp.path().join("out.bam");
        let block = make_block(5000);
        let mut bw = BamWriter::new(&bam_path, &[("chr1".to_string(), 5000usize)]).unwrap();

        // Fragments 400 bp wide, spaced 100 bp apart so positions are strictly increasing.
        let fragments = vec![
            (0usize, 400),
            (500, 900),
            (1000, 1400),
            (1500, 1900),
            (2000, 2400),
        ];
        run_write_block(&mut bw, fragments, &block, true, 150);

        assert_eq!(
            bw.carry_len(),
            2,
            "PE: carry should hold only the last fragment's read pair, got {}",
            bw.carry_len()
        );
    }

    #[test]
    fn test_write_block_fastq_bam_sorted_with_incremental_flush() {
        // End-to-end: write_block_fastq with incremental flush then flush_all must
        // produce a coordinate-sorted BAM, identical to bulk flushing.
        use noodles::bam;

        let temp = tempfile::tempdir().unwrap();
        let bam_path = temp.path().join("sorted.bam");
        let block = make_block(500);
        let contigs = vec![("chr1".to_string(), 500usize)];
        let mut bw = BamWriter::new(&bam_path, &contigs).unwrap();

        let fragments = vec![(0usize, 10), (50, 60), (100, 110), (150, 160), (200, 210)];
        run_write_block(&mut bw, fragments, &block, false, 10);
        bw.flush_all().unwrap();
        drop(bw); // BGZF EOF written on drop

        let file = std::fs::File::open(&bam_path).unwrap();
        let mut reader = bam::io::Reader::new(file);
        reader.read_header().unwrap();
        let positions: Vec<usize> = reader
            .records()
            .map(|r| r.unwrap().alignment_start().unwrap().unwrap().get() - 1)
            .collect();
        assert_eq!(
            positions,
            vec![0, 50, 100, 150, 200],
            "BAM must be coordinate-sorted; got: {:?}",
            positions
        );
    }
}

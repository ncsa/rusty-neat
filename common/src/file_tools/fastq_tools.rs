//! This library writes either single ended or paired-ended fastq files.
//! Need to update this method. We want to use the data structures and we want to make sure
//! this function is generic enough to work with the fragmented method we are implementing
//! This one needs a major overhaul, it is autogenerating quality scores etc.
//! Will wait to get other things set up first
use std::io::{BufWriter, Write};
use log::debug;
use std::path::PathBuf;
use std::collections::HashMap;
use crate::rng::{NeatRng, NeatRngError};
use flate2::Compression;
use flate2::write::GzEncoder;
use thiserror::Error;

use crate::structs::mutated_map::{MutatedMap, MutatedMapError};
use crate::structs::sequence_block::{SequenceBlock, SequenceBlockError};
use crate::structs::nucleotides::Nucleotide;
use crate::structs::read_record::ReadRecord;
use crate::file_tools::file_io::{append_to_file, read_gzip_lines};
use crate::file_tools::bam_writer::BamRecordStager;
use crate::models::quality_scores::QualityScoreModel;
use crate::models::sequencing_error_model::{SeqModelError, SequencingErrorModel, SequencingErrorType};
use crate::structs::nucleotides::Nucleotide::N;
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

fn reverse_complement(sequence: Vec<Nucleotide>) -> Vec<Nucleotide> {
    // Returns the reverse complement of a vector of u8's representing a DNA sequence.
    let length = sequence.len();
    let mut rev_comp = Vec::new();
    for i in (0..length).rev() {
        rev_comp.push(sequence[i].complement())
    }
    rev_comp
}

pub fn write_block_fastq<T: Write, W: Write> (
    block_fragments: Vec<(usize, usize)>,
    block_map: &MutatedMap,
    sequence_block: &SequenceBlock,
    paired_ended: bool,
    buffer1: &mut GzEncoder<T>,
    buffer2: &mut GzEncoder<W>,
    read_length: usize,
    read_name_prefix: &str,
    quality_score_model: &QualityScoreModel,
    sequencing_error_model: &SequencingErrorModel,
    rng: &mut NeatRng,
    mut bam_writer: Option<&mut dyn BamRecordStager>,
) -> Result<(), FastqToolsError> {
    debug!("writing reads for {}", sequence_block.contig);
    for (start, end) in block_fragments {
        let fragment = sequence_block.get_subseq(start, end)?;
        let mut read1_variants: HashMap<usize, &Variant> = HashMap::new();
        let mut reads1_flagged: Vec<usize> = Vec::new();
        let mut read2_variants: HashMap<usize, &Variant> = HashMap::new();
        let mut reads2_flagged: Vec<usize> = Vec::new();
        for (pos, variant) in &block_map.variant_map {
            if (start..(start+read_length)).contains(&pos) {
                let var_pos = pos - start;
                read1_variants.insert(var_pos, variant);
                reads1_flagged.push(var_pos);
            }
            if end > read_length {
                if paired_ended == true && ((end-read_length)..end).contains(&pos) {
                    let var_pos = (end - 1) - pos;
                    read2_variants.insert(var_pos, variant);
                    reads2_flagged.push(var_pos);
                }
            }
        }

        let ref_start = sequence_block.ref_start;
        let abs_start = start + ref_start;
        let abs_end = end + ref_start;
        let base_name = format!("{}_{:010}_{:010}", read_name_prefix, abs_start, abs_end);

        let r2_start = if paired_ended && abs_end >= read_length {
            abs_end - read_length
        } else {
            0
        };
        let tlen = if paired_ended { (abs_end - abs_start) as i32 } else { 0 };

        let quality_scores_1 = quality_score_model.generate_quality_scores(read_length, rng)?;
        let r1_record = match generate_read(
            &fragment,
            &reads1_flagged,
            &read1_variants,
            read_length,
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
        ) {
            Ok(record) => record,
            Err(FastqToolsError::TruncatedRead(msg)) => { debug!("{}", msg); continue }
            Err(e) => return Err(e),
        };

        write_read_to_fastq(&r1_record, buffer1)?;
        if let Some(ref mut bam) = bam_writer {
            bam.stage_read_record(&r1_record)
                .map_err(|e| FastqToolsError::BamError(e.to_string()))?;
        }

        if paired_ended {
            let quality_scores_2 = quality_score_model.generate_quality_scores(read_length, rng)?;
            let r2_pos = if abs_end >= read_length { abs_end - read_length } else { 0 };
            let tlen_r2 = -((abs_end - abs_start) as i32);

            let r2_record = match generate_read(
                &reverse_complement(fragment),
                &reads2_flagged,
                &read2_variants,
                read_length,
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
            ) {
                Ok(record) => record,
                Err(FastqToolsError::TruncatedRead(msg)) => { debug!("{}", msg); continue }
                Err(e) => return Err(e),
            };

            write_read_to_fastq(&r2_record, buffer2)?;
            if let Some(ref mut bam) = bam_writer {
                bam.stage_read_record(&r2_record)
                    .map_err(|e| FastqToolsError::BamError(e.to_string()))?;
            }
        }
    }
    Ok(())
}

pub fn combine_temp_fastqs(
    files_r1: Vec<PathBuf>,
    files_r2: Vec<PathBuf>,
    final_filename_r1: &PathBuf,
    final_filename_r2: Option<&PathBuf>,
    shuffle: bool,
    rng: &mut NeatRng,
) -> Result<(), FastqToolsError> {
    let mut r1_lines: Vec<String> = Vec::new();
    for file in &files_r1 {
        let reader = read_gzip_lines(file)?;
        for line_result in reader {
            match line_result {
                Ok(l) => r1_lines.push(l),
                Err(e) => return Err(FastqToolsError::FastqReadError(e.to_string())),
            }
        }
    }
    let mut records_r1: Vec<Vec<String>> = r1_lines
        .chunks(4)
        .filter(|c| c.len() == 4)
        .map(|c| c.to_vec())
        .collect();

    let paired = !files_r2.is_empty();
    let mut records_r2: Vec<Vec<String>> = if paired {
        let mut r2_lines: Vec<String> = Vec::new();
        for file in &files_r2 {
            let reader = read_gzip_lines(file)?;
            for line_result in reader {
                match line_result {
                    Ok(l) => r2_lines.push(l),
                    Err(e) => return Err(FastqToolsError::FastqReadError(e.to_string())),
                }
            }
        }
        r2_lines
            .chunks(4)
            .filter(|c| c.len() == 4)
            .map(|c| c.to_vec())
            .collect()
    } else {
        Vec::new()
    };

    if shuffle {
        if paired && !records_r2.is_empty() {
            let mut pairs: Vec<(Vec<String>, Vec<String>)> = records_r1
                .into_iter()
                .zip(records_r2.into_iter())
                .collect();
            rng.shuffle_in_place(&mut pairs)?;
            let (new_r1, new_r2): (Vec<Vec<String>>, Vec<Vec<String>>) =
                pairs.into_iter().unzip();
            records_r1 = new_r1;
            records_r2 = new_r2;
        } else {
            rng.shuffle_in_place(&mut records_r1)?;
        }
    }

    {
        let final_file = append_to_file(final_filename_r1)?;
        let enc = GzEncoder::new(final_file, Compression::default());
        let mut writer = BufWriter::new(enc);
        for record in &records_r1 {
            for line in record {
                writer.write_fmt(format_args!("{}\n", line))?;
            }
        }
    }

    if let Some(filename_r2) = final_filename_r2 {
        let final_file = append_to_file(filename_r2)?;
        let enc = GzEncoder::new(final_file, Compression::default());
        let mut writer = BufWriter::new(enc);
        for record in &records_r2 {
            for line in record {
                writer.write_fmt(format_args!("{}\n", line))?;
            }
        }
    }

    Ok(())
}

fn generate_read(
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
) -> Result<ReadRecord, FastqToolsError> {
    if sequence.len() < read_length {
        return Err(FastqToolsError::TruncatedRead(format!("{:?}", sequence)))
    }

    let is_reverse = matches!(read_strand, Strand::Reverse);
    let fragment_length = sequence.len();

    let mut bases_written = 0;
    let mut out_seq = String::new();
    let mut cigar_ops: Vec<char> = Vec::new();
    let mut quality_index = 0;
    let mut seq_index = 0;

    'outer: while (seq_index < fragment_length) && (bases_written < read_length) {
        let fragment_position = match read_strand {
            Strand::Forward => seq_index,
            Strand::Reverse => fragment_length - seq_index,
        };
        let reference_base = sequence[seq_index].get_unmasked_base();
        let mut base_to_write: Vec<Nucleotide> = vec![reference_base];
        let mut deletion_skip: usize = 0;

        if reference_base == N {
            // Don't try to modify this base.
        } else if flagged_positions.contains(&fragment_position) {
            let variant = variant_map[&fragment_position];
            if (variant.genotype == Genotype::Homozygous) || (rng.random()? < 0.5) {
                base_to_write = match read_strand {
                    Strand::Forward => variant.alternate.clone(),
                    Strand::Reverse => variant.alternate.iter().map(|b| b.complement()).collect(),
                };
            }
        } else {
            let score = quality_scores[quality_index];
            let prob = sequencing_error_model.convert_score(score)?;
            if rng.random()? < prob {
                debug!("Creating sequencing error");
                let error = sequencing_error_model.generate_sequencing_error(reference_base, rng)?;
                match error {
                    SequencingErrorType::SnpError(base) => {
                        debug!("Snp error");
                        base_to_write = vec![base];
                    },
                    SequencingErrorType::DeletionError(length) => {
                        debug!("Deletion error");
                        if seq_index + length < sequence.len() {
                            seq_index += length;
                            deletion_skip = length;
                        }
                    },
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
        return Err(FastqToolsError::TruncatedRead(format!("{:?}", sequence)))
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

fn write_read_to_fastq<T: Write>(
    record: &ReadRecord,
    buffer: &mut GzEncoder<T>,
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
    use crate::file_tools::file_io::{create_output_file, read_gzip_lines};
    use crate::structs::variants::{Variant, VariantType};
    use crate::structs::nucleotides::Nucleotide::*;
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

        let mut rng = NeatRng::new_from_seed(&vec!["test".to_string()]).unwrap();
        combine_temp_fastqs(vec![file1, file2], vec![], &output, None, false, &mut rng).unwrap();

        let lines: Vec<String> = read_gzip_lines(&output)
            .unwrap()
            .map(|l| l.unwrap())
            .collect();
        assert_eq!(lines.len(), 8, "Combined file should have 8 lines (2 records × 4 lines)");
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
        ]).unwrap();
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
        ).unwrap();
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
                    break
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
        assert_eq!(qual_string, quality_scores_to_char_vec(&qual_scores).unwrap())
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
            sequence_block::{RegionType, SequenceBlock, SequenceMap},
            mutated_map::MutatedMap,
        };
        let temp_dir = tempfile::tempdir().unwrap();
        let ref_start: usize = 1000;
        let seq_len: usize = 50;
        let sequence: Vec<Nucleotide> = (0..seq_len)
            .map(|i| match i % 4 { 0 => A, 1 => C, 2 => G, _ => T })
            .collect();
        let block = SequenceBlock {
            contig: "chr1".to_string(),
            ref_start,
            ref_end: ref_start + seq_len,
            sequence: sequence.clone(),
            sequence_map: vec![SequenceMap::from(RegionType::NonNRegion, 0, seq_len)],
        };
        let mutated_map = MutatedMap::from_interval(ref_start, ref_start + seq_len, vec![]).unwrap();
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
            fragments, &mutated_map, &block, false,
            &mut buf1, &mut buf2,
            10, "chr1",
            &quality_model, &seq_err_model, &mut rng,
            None,
        ).unwrap();
        buf1.finish().unwrap();
        let lines: Vec<String> = read_gzip_lines(&out_path)
            .unwrap()
            .map(|l| l.unwrap())
            .collect();
        let header = &lines[0];
        assert!(header.starts_with('@'), "Expected FASTQ header, got: {}", header);
        // Read name must contain ref-relative positions, not block-local ones
        let expected_start = format!("{:010}", frag_start + ref_start); // 0000001005
        let expected_end   = format!("{:010}", frag_end   + ref_start); // 0000001025
        assert!(
            header.contains(&expected_start),
            "Read name should contain ref-relative start {}; got: {}", expected_start, header
        );
        assert!(
            header.contains(&expected_end),
            "Read name should contain ref-relative end {}; got: {}", expected_end, header
        );
        // Also verify block-local positions are NOT what was written
        let local_start = format!("{:010}", frag_start); // 0000000005
        assert!(
            !header.contains(&local_start),
            "Read name should NOT contain block-local start {}; got: {}", local_start, header
        );
    }

    // --- CIGAR-building tests for the refactored generate_read ---

    fn make_sequence(len: usize) -> Vec<Nucleotide> {
        (0..len).map(|i| match i % 4 { 0 => A, 1 => C, 2 => G, _ => T }).collect()
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
            &sequence, &[], &HashMap::new(), read_length,
            "r/1".to_string(), Strand::Forward, quality_scores, &model, &mut rng,
            "chr1".to_string(), 0, "chr1".to_string(), 0, 0, false,
        ).unwrap();
        assert_eq!(record.sequence.len(), read_length);
        assert_eq!(record.cigar_ops.len(), read_length);
        assert!(record.cigar_ops.iter().all(|&c| c == 'M'),
            "expected all-M cigar, got: {:?}", record.cigar_ops);
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
            &sequence, &[], &HashMap::new(), read_length,
            "r/1".to_string(), Strand::Forward, quality_scores, &model, &mut rng,
            "chr1".to_string(), 0, "chr1".to_string(), 0, 0, false,
        ).unwrap();
        let mi = record.cigar_ops.iter().filter(|&&c| c == 'M' || c == 'I').count();
        let d  = record.cigar_ops.iter().filter(|&&c| c == 'D').count();
        assert_eq!(record.sequence.len(), read_length);
        assert_eq!(mi, read_length, "M+I must equal read_length");
        assert_eq!(record.cigar_ops.len(), mi + d, "cigar len must equal M+I+D");
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
            &sequence, &[3], &variant_map, read_length,
            "r/1".to_string(), Strand::Forward, quality_scores, &model, &mut rng,
            "chr1".to_string(), 0, "chr1".to_string(), 0, 0, false,
        ).unwrap();
        assert_eq!(record.cigar_ops.len(), read_length);
        assert!(record.cigar_ops.iter().all(|&c| c == 'M'),
            "SNP variant must not add I or D ops; got: {:?}", record.cigar_ops);
    }

    #[test]
    fn test_generate_read_insertion_variant_produces_i_ops() {
        // Homozygous insertion: alt = [A, C, G] → 1 M (anchor) + 2 I (inserted bases).
        let sequence = make_sequence(30);
        let read_length = 10;
        let ins = Variant::new(
            VariantType::Insertion, 3, &vec![A], &vec![A, C, G], &mut vec![1, 1],
        ).unwrap();
        let variant_map = HashMap::from([(3usize, &ins)]);
        let quality_scores = vec![40usize; read_length];
        let model = SequencingErrorModel::default().unwrap();
        let mut rng = NeatRng::new_from_seed(&vec!["ins_variant".to_string()]).unwrap();
        let record = generate_read(
            &sequence, &[3], &variant_map, read_length,
            "r/1".to_string(), Strand::Forward, quality_scores, &model, &mut rng,
            "chr1".to_string(), 0, "chr1".to_string(), 0, 0, false,
        ).unwrap();
        let i_count = record.cigar_ops.iter().filter(|&&c| c == 'I').count();
        assert_eq!(i_count, 2, "2-base insertion variant must produce 2 I ops");
        let mi = record.cigar_ops.iter().filter(|&&c| c == 'M' || c == 'I').count();
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
            &sequence, &[], &HashMap::new(), read_length,
            "r/1".to_string(), Strand::Forward, quality_scores, &model, &mut rng,
            "chr1".to_string(), 0, "chr1".to_string(), 0, 0, false,
        ).unwrap();
        assert!(
            record.cigar_ops.iter().any(|&c| c != 'M'),
            "expected non-M ops when every position has a guaranteed error; got: {:?}", record.cigar_ops
        );
        // Invariant must hold even under heavy errors
        let mi = record.cigar_ops.iter().filter(|&&c| c == 'M' || c == 'I').count();
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
            &sequence, &[], &HashMap::new(), read_length,
            "r/1".to_string(), Strand::Forward, quality_scores, &model, &mut rng,
            "chr1".to_string(), 0, "chr1".to_string(), 0, 0, false,
        ).unwrap();
        let mi = record.cigar_ops.iter().filter(|&&c| c == 'M' || c == 'I').count();
        let d  = record.cigar_ops.iter().filter(|&&c| c == 'D').count();
        assert_eq!(mi, read_length);
        assert_eq!(record.cigar_ops.len(), mi + d);
        if d > 0 {
            assert!(record.cigar_ops.len() > read_length,
                "deletion errors must make cigar longer than read_length");
        }
    }

    #[test]
    fn test_apply_variants() {
        let sequence = vec![A, C, G, T, T, A, T, G, A, C, G, T, T, A, T, G];
        let variant1 = Variant::new(
                VariantType::SNP,
                1,
                &vec![T],
                &vec![C],
                &mut vec![1, 0],
        ).unwrap();
        let variant2 = Variant::new(
                VariantType::Deletion,
                3,
                &vec![T, T],
                &vec![T],
                &mut vec![0,1],
        ).unwrap();
        let variant_map = HashMap::from([
            (1, &variant1),
            (3, &variant2),
        ]);
        let flagged_positions = vec![1,3];
        let read_name = "neat_generated__0000000000_0000000008/1".to_string();
        let qual_scores = vec![33, 25, 37, 28, 15, 33, 33, 37];
        let sequencing_error_model = SequencingErrorModel::default().unwrap();
        let mut rng = NeatRng::new_from_seed(&vec![
            "Hello".to_string(),
            "Cruel".to_string(),
            "World".to_string(),
        ]).unwrap();
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
        );
        assert!(result.is_ok());
    }

}
use std::io;
use std::path::PathBuf;
use thiserror::Error;
use crate::{
    file_tools::file_io::{is_gzipped_file, read_gzip_lines, read_lines},
    rng::{NeatRng, NeatRngError},
    structs::{
        nucleotides::{Nucleotide, Nucleotide::{N, X, Maskeda, Maskedc, Maskedg, Maskedt}, NucleotideSelector},
        sequence_block::{RegionType, SequenceMap},
    },
};

#[derive(Debug, Error)]
pub enum FastaStreamError {
    #[error("I/O error reading FASTA: {0}")]
    IoError(#[from] io::Error),
}

/// Streaming FASTA reader that yields `(contig_name, sequence)` one contig at a time.
///
/// Only one contig's sequence is held in memory at a time, making this suitable for
/// large reference genomes without the temp-file overhead of `read_fasta`.
pub struct FastaStream {
    lines: Box<dyn Iterator<Item = io::Result<String>> + Send>,
    pending_name: Option<String>,
}

impl FastaStream {
    pub fn open(path: &PathBuf) -> Result<Self, FastaStreamError> {
        let lines: Box<dyn Iterator<Item = io::Result<String>> + Send> =
            if is_gzipped_file(path)? {
                Box::new(read_gzip_lines(path)?)
            } else {
                Box::new(read_lines(path)?)
            };

        let mut stream = FastaStream { lines, pending_name: None };

        // Advance to the first header line.
        loop {
            match stream.lines.next() {
                None => return Ok(stream),
                Some(Err(e)) => return Err(e.into()),
                Some(Ok(line)) if line.starts_with('>') => {
                    stream.pending_name = Some(parse_contig_name(&line));
                    break;
                }
                Some(Ok(_)) => {}
            }
        }

        Ok(stream)
    }
}

impl Iterator for FastaStream {
    type Item = Result<(String, Vec<Nucleotide>), FastaStreamError>;

    fn next(&mut self) -> Option<Self::Item> {
        let name = self.pending_name.take()?;
        let mut sequence: Vec<Nucleotide> = Vec::new();

        loop {
            match self.lines.next() {
                None => break,
                Some(Err(e)) => return Some(Err(e.into())),
                Some(Ok(line)) => {
                    if line.starts_with('>') {
                        self.pending_name = Some(parse_contig_name(&line));
                        break;
                    }
                    sequence.extend(line.chars().map(Nucleotide::from));
                }
            }
        }

        Some(Ok((name, sequence)))
    }
}

fn parse_contig_name(header: &str) -> String {
    header[1..]
        .split_whitespace()
        .next()
        .unwrap_or("")
        .to_string()
}

/// Builds a `Vec<SequenceMap>` that partitions `sequence` into contiguous N/non-N runs.
///
/// The returned map covers the full sequence without gaps or overlaps.
/// Never returns an error; the `Result` wrapper is kept for call-site compatibility.
pub fn map_buffer(sequence: &[Nucleotide]) -> Vec<SequenceMap> {
    if sequence.is_empty() {
        return Vec::new();
    }
    let mut map: Vec<SequenceMap> = Vec::new();
    let mut region_start = 0;
    let mut region_end = 1;
    let mut inside_n_region = matches!(sequence[0], N | X | Maskeda | Maskedc | Maskedg | Maskedt);

    if inside_n_region {
        for base in &sequence[1..] {
            match base {
                N | X | Maskeda | Maskedc | Maskedg | Maskedt => region_end += 1,
                _ => {
                    map.push(SequenceMap::from(RegionType::NRegion, 0, region_end));
                    region_start = region_end;
                    region_end = region_start + 1;
                    inside_n_region = false;
                    break;
                }
            }
        }
        if inside_n_region {
            map.push(SequenceMap::from(RegionType::NRegion, 0, region_end));
            return map;
        }
    }

    for base in &sequence[region_end..] {
        match base {
            N | X | Maskeda | Maskedc | Maskedg | Maskedt => {
                if inside_n_region {
                    region_end += 1;
                } else {
                    inside_n_region = true;
                    map.push(SequenceMap::from(RegionType::NonNRegion, region_start, region_end));
                    region_start = region_end;
                    region_end = region_start + 1;
                }
            }
            _ => {
                if inside_n_region {
                    inside_n_region = false;
                    map.push(SequenceMap::from(RegionType::NRegion, region_start, region_end));
                    region_start = region_end;
                    region_end = region_start + 1;
                } else {
                    region_end += 1;
                }
            }
        }
    }

    let region_type = if inside_n_region { RegionType::NRegion } else { RegionType::NonNRegion };
    map.push(SequenceMap::from(region_type, region_start, region_end));
    map
}

/// Replaces each N base in `sequence` with a randomly sampled masked base in-place.
pub fn apply_n_substitution(
    sequence: &mut Vec<Nucleotide>,
    selector: &NucleotideSelector,
    rng: &mut NeatRng,
) -> Result<(), NeatRngError> {
    for base in sequence.iter_mut() {
        if *base == N {
            *base = selector.sample_bases(rng.random()?).get_masked();
        }
    }
    Ok(())
}

/// Returns the contiguous non-N/non-X regions of a sequence as `(start, end)` pairs.
pub fn non_n_regions(sequence: &[Nucleotide]) -> Vec<(usize, usize)> {
    let mut regions = Vec::new();
    let mut region_start: Option<usize> = None;

    for (i, &nuc) in sequence.iter().enumerate() {
        match (nuc, region_start) {
            (Nucleotide::N | Nucleotide::X, Some(s)) => {
                regions.push((s, i));
                region_start = None;
            }
            (Nucleotide::N | Nucleotide::X, None) => {}
            (_, None) => region_start = Some(i),
            (_, Some(_)) => {}
        }
    }
    if let Some(s) = region_start {
        regions.push((s, sequence.len()));
    }
    regions
}

/// Scans a FASTA file and returns `(contig_name, length)` for each contig
/// without storing sequences. Used to build BAM headers before streaming reads.
pub fn scan_fasta_lengths(path: &PathBuf) -> Result<Vec<(String, usize)>, FastaStreamError> {
    let lines: Box<dyn Iterator<Item = io::Result<String>>> =
        if is_gzipped_file(path)? {
            Box::new(read_gzip_lines(path)?)
        } else {
            Box::new(read_lines(path)?)
        };

    let mut result: Vec<(String, usize)> = Vec::new();
    let mut current_name: Option<String> = None;
    let mut current_len: usize = 0;

    for line in lines {
        let line = line?;
        if line.starts_with('>') {
            if let Some(name) = current_name.take() {
                result.push((name, current_len));
            }
            current_name = Some(parse_contig_name(&line));
            current_len = 0;
        } else {
            current_len += line.chars().filter(|c| !c.is_whitespace()).count();
        }
    }
    if let Some(name) = current_name {
        result.push((name, current_len));
    }

    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;
    use flate2::{Compression, write::GzEncoder};

    fn write_fasta(content: &str) -> NamedTempFile {
        let mut f = NamedTempFile::new().unwrap();
        write!(f, "{}", content).unwrap();
        f
    }

    fn write_gzipped_fasta(content: &str) -> NamedTempFile {
        let mut f = NamedTempFile::new().unwrap();
        let mut enc = GzEncoder::new(&mut f, Compression::default());
        enc.write_all(content.as_bytes()).unwrap();
        enc.finish().unwrap();
        f
    }

    fn make_sequence(bases: &str) -> Vec<Nucleotide> {
        bases.chars().map(Nucleotide::from).collect()
    }

    // FastaStream tests

    #[test]
    fn test_reads_single_contig() {
        let f = write_fasta(">chr1\nACGT\n");
        let contigs: Vec<_> = FastaStream::open(&f.path().to_path_buf())
            .unwrap()
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        assert_eq!(contigs.len(), 1);
        assert_eq!(contigs[0].0, "chr1");
        assert_eq!(contigs[0].1, vec![Nucleotide::A, Nucleotide::C, Nucleotide::G, Nucleotide::T]);
    }

    #[test]
    fn test_reads_multiple_contigs() {
        let f = write_fasta(">chr1\nACGT\n>chr2\nTTTT\n");
        let contigs: Vec<_> = FastaStream::open(&f.path().to_path_buf())
            .unwrap()
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        assert_eq!(contigs.len(), 2);
        assert_eq!(contigs[0].0, "chr1");
        assert_eq!(contigs[1].0, "chr2");
    }

    #[test]
    fn test_parses_multiline_sequence() {
        let f = write_fasta(">chr1\nACGT\nACGT\n");
        let contigs: Vec<_> = FastaStream::open(&f.path().to_path_buf())
            .unwrap()
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        assert_eq!(contigs[0].1.len(), 8);
    }

    #[test]
    fn test_strips_description_from_header() {
        let f = write_fasta(">chr1 some description here\nACGT\n");
        let contigs: Vec<_> = FastaStream::open(&f.path().to_path_buf())
            .unwrap()
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        assert_eq!(contigs[0].0, "chr1");
    }

    #[test]
    fn test_reads_gzipped_fasta() {
        let f = write_gzipped_fasta(">chr1\nACGT\n");
        let contigs: Vec<_> = FastaStream::open(&f.path().to_path_buf())
            .unwrap()
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        assert_eq!(contigs.len(), 1);
        assert_eq!(contigs[0].0, "chr1");
        assert_eq!(contigs[0].1.len(), 4);
    }

    #[test]
    fn test_n_bases_decoded_correctly() {
        let f = write_fasta(">chr1\nACNT\n");
        let contigs: Vec<_> = FastaStream::open(&f.path().to_path_buf())
            .unwrap()
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        assert_eq!(contigs[0].1[2], Nucleotide::N);
    }

    #[test]
    fn test_empty_file_yields_no_contigs() {
        let f = write_fasta("");
        let contigs: Vec<_> = FastaStream::open(&f.path().to_path_buf())
            .unwrap()
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        assert!(contigs.is_empty());
    }

    // non_n_regions tests

    #[test]
    fn test_non_n_regions_no_ns() {
        let seq = make_sequence("ACGTACGT");
        assert_eq!(non_n_regions(&seq), vec![(0, 8)]);
    }

    #[test]
    fn test_non_n_regions_leading_n() {
        let seq = make_sequence("NNACGT");
        assert_eq!(non_n_regions(&seq), vec![(2, 6)]);
    }

    #[test]
    fn test_non_n_regions_interior_n() {
        let seq = make_sequence("ACGTNACGT");
        assert_eq!(non_n_regions(&seq), vec![(0, 4), (5, 9)]);
    }

    #[test]
    fn test_non_n_regions_all_n() {
        let seq = make_sequence("NNNN");
        assert!(non_n_regions(&seq).is_empty());
    }

    #[test]
    fn test_non_n_regions_x_closes_region() {
        let mut seq = make_sequence("ACGT");
        seq.push(Nucleotide::X);
        seq.extend(make_sequence("ACGT"));
        assert_eq!(non_n_regions(&seq), vec![(0, 4), (5, 9)]);
    }

    #[test]
    fn test_non_n_regions_x_at_start_is_skipped() {
        let mut seq = vec![Nucleotide::X, Nucleotide::X];
        seq.extend(make_sequence("ACGT"));
        assert_eq!(non_n_regions(&seq), vec![(2, 6)]);
    }

    // scan_fasta_lengths tests

    #[test]
    fn test_scan_fasta_lengths_single_contig() {
        let f = write_fasta(">chr1\nACGT\n");
        let lengths = scan_fasta_lengths(&f.path().to_path_buf()).unwrap();
        assert_eq!(lengths, vec![("chr1".to_string(), 4)]);
    }

    #[test]
    fn test_scan_fasta_lengths_multiple_contigs() {
        let f = write_fasta(">chr1\nACGT\n>chr2\nGGGGGGGGGG\n");
        let lengths = scan_fasta_lengths(&f.path().to_path_buf()).unwrap();
        assert_eq!(lengths, vec![
            ("chr1".to_string(), 4),
            ("chr2".to_string(), 10),
        ]);
    }

    #[test]
    fn test_scan_fasta_lengths_multiline() {
        let f = write_fasta(">chr1\nACGT\nACGT\n");
        let lengths = scan_fasta_lengths(&f.path().to_path_buf()).unwrap();
        assert_eq!(lengths, vec![("chr1".to_string(), 8)]);
    }

    #[test]
    fn test_scan_fasta_lengths_strips_description() {
        let f = write_fasta(">chr1 some description\nACGT\n");
        let lengths = scan_fasta_lengths(&f.path().to_path_buf()).unwrap();
        assert_eq!(lengths[0].0, "chr1");
    }
}
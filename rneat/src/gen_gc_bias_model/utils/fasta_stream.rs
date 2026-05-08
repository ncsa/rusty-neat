use std::path::PathBuf;
use common::{
    file_tools::file_io::{is_gzipped_file, read_gzip_lines, read_lines},
    structs::nucleotides::Nucleotide,
};
use crate::gen_gc_bias_model::errors::GenGcBiasModelError;

/// Streaming FASTA reader that yields `(contig_name, sequence)` one contig at a time.
///
/// Only one contig's sequence is held in memory at a time, making this suitable for
/// large reference genomes without the temp-file overhead of `read_fasta`.
pub struct FastaStream {
    lines: Box<dyn Iterator<Item = std::io::Result<String>>>,
    /// The header line consumed while reading the previous contig's sequence.
    pending_name: Option<String>,
}

impl FastaStream {
    pub fn open(path: &PathBuf) -> Result<Self, GenGcBiasModelError> {
        let lines: Box<dyn Iterator<Item = std::io::Result<String>>> =
            if is_gzipped_file(path)? {
                Box::new(read_gzip_lines(path)?)
            } else {
                Box::new(read_lines(path)?)
            };

        let mut stream = FastaStream { lines, pending_name: None };

        // Advance to the first header line.
        loop {
            match stream.lines.next() {
                None => return Ok(stream), // empty file
                Some(Err(e)) => return Err(e.into()),
                Some(Ok(line)) if line.starts_with('>') => {
                    stream.pending_name = Some(parse_contig_name(&line));
                    break;
                }
                Some(Ok(_)) => {} // skip any leading non-header lines
            }
        }

        Ok(stream)
    }
}

impl Iterator for FastaStream {
    type Item = Result<(String, Vec<Nucleotide>), GenGcBiasModelError>;

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
    // Take the first whitespace-delimited token after '>'
    header[1..]
        .split_whitespace()
        .next()
        .unwrap_or("")
        .to_string()
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
}
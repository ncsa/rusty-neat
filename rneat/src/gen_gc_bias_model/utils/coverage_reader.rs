use std::{
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader, Read, Seek, SeekFrom},
    path::PathBuf,
};
use common::file_tools::file_io::{is_gzipped_file, read_gzip_lines};
use crate::gen_gc_bias_model::{errors::GenGcBiasModelError, utils::config::CoverageFormat};

// ── Plain-text index (byte-offset seek) ──────────────────────────────────────

/// Byte-offset index into a plain-text coverage file.
///
/// Built in a single forward pass. For each contig, stores the start and end byte
/// offsets so we can seek directly to that contig's data when needed.
struct CoverageIndex {
    offsets: HashMap<String, (u64, u64)>,
}

impl CoverageIndex {
    fn build(path: &PathBuf) -> Result<Self, GenGcBiasModelError> {
        let file = File::open(path)?;
        let mut reader = BufReader::new(file);

        let mut offsets: HashMap<String, (u64, u64)> = HashMap::new();
        let mut current_contig = String::new();
        let mut contig_start: u64 = 0;
        let mut byte_offset: u64 = 0;
        let mut line = String::new();

        loop {
            line.clear();
            let n = reader.read_line(&mut line)?;
            if n == 0 {
                break;
            }

            let contig = line.split('\t').next().unwrap_or("").to_string();
            if contig.is_empty() {
                byte_offset += n as u64;
                continue;
            }

            if contig != current_contig {
                if !current_contig.is_empty() {
                    offsets.insert(current_contig.clone(), (contig_start, byte_offset));
                }
                current_contig = contig;
                contig_start = byte_offset;
            }

            byte_offset += n as u64;
        }

        if !current_contig.is_empty() {
            offsets.insert(current_contig, (contig_start, byte_offset));
        }

        Ok(CoverageIndex { offsets })
    }
}

// ── CoverageData ─────────────────────────────────────────────────────────────

/// Per-contig depth array (0-indexed).
///
/// Positions not present in the coverage file (gaps in bedtools-genomecov-dz output)
/// are stored implicitly as 0 and returned as 0 via `depth_at`.
pub struct CoverageData {
    depths: Vec<u32>,
}

impl CoverageData {
    fn load(
        path: &PathBuf,
        index: &CoverageIndex,
        contig: &str,
        format: CoverageFormat,
    ) -> Result<Option<Self>, GenGcBiasModelError> {
        let &(start_byte, end_byte) = match index.offsets.get(contig) {
            Some(offsets) => offsets,
            None => return Ok(None),
        };

        let mut file = File::open(path)?;
        file.seek(SeekFrom::Start(start_byte))?;
        let limited = file.take(end_byte - start_byte);
        let reader = BufReader::new(limited);

        let mut depths: Vec<u32> = Vec::new();
        for line in reader.lines() {
            let line = line?;
            append_depth_record(&line, format, &mut depths)?;
        }
        Ok(Some(CoverageData { depths }))
    }

    /// Depth at a 0-based position. Returns 0 for positions beyond the stored range.
    #[inline]
    pub fn depth_at(&self, pos: usize) -> u32 {
        self.depths.get(pos).copied().unwrap_or(0)
    }
}

// ── Shared line-parsing helper ────────────────────────────────────────────────

/// Parse one `contig\tpos\tdepth` line and push the depth into `out`.
/// Skips the contig field (caller already knows which contig it's reading).
fn append_depth_record(
    line: &str,
    format: CoverageFormat,
    out: &mut Vec<u32>,
) -> Result<(), GenGcBiasModelError> {
    let mut parts = line.splitn(3, '\t');
    let _chrom = parts.next();
    let pos_str = parts.next().ok_or_else(|| {
        GenGcBiasModelError::CoverageParseError(format!(
            "Missing position field in line: {}",
            line
        ))
    })?;
    let depth_str = parts.next().ok_or_else(|| {
        GenGcBiasModelError::CoverageParseError(format!(
            "Missing depth field in line: {}",
            line
        ))
    })?;

    let pos: usize = pos_str.trim().parse().map_err(|_| {
        GenGcBiasModelError::CoverageParseError(format!(
            "Invalid position '{}' in line: {}",
            pos_str, line
        ))
    })?;
    let depth: u32 = depth_str.trim().parse().map_err(|_| {
        GenGcBiasModelError::CoverageParseError(format!(
            "Invalid depth '{}' in line: {}",
            depth_str, line
        ))
    })?;

    let idx = if format.is_zero_based() { pos } else { pos - 1 };

    if idx >= out.len() {
        out.resize(idx, 0);
        out.push(depth);
    } else {
        out[idx] = depth;
    }
    Ok(())
}

// ── CoverageReader ────────────────────────────────────────────────────────────

/// Unified coverage file reader that handles both plain-text and gzip-compressed files.
///
/// **Plain text**: builds a byte-offset index in one forward pass; each contig is then
/// loaded on demand via `seek` — O(1) memory per contig.
///
/// **Gzip**: streams the entire file in a single forward pass, building all contigs at
/// once. Memory scales with total depth data (≈ 4 bytes × genome positions covered).
/// Contigs are evicted from the in-memory map as they are consumed via `get`, so peak
/// memory declines as processing advances. For whole-genome gzip files, plain text uses
/// less peak memory; gzip is convenient for single-chromosome or targeted files.
pub struct CoverageReader {
    inner: ReaderInner,
}

enum ReaderInner {
    Plain {
        path: PathBuf,
        index: CoverageIndex,
        format: CoverageFormat,
    },
    Gzip {
        data: HashMap<String, CoverageData>,
    },
}

impl CoverageReader {
    /// Open a coverage file. Detects gzip by `.gz` extension.
    pub fn open(path: &PathBuf, format: CoverageFormat) -> Result<Self, GenGcBiasModelError> {
        if is_gzipped_file(path)? {
            let data = load_all_gzip(path, format)?;
            Ok(Self { inner: ReaderInner::Gzip { data } })
        } else {
            let index = CoverageIndex::build(path)?;
            Ok(Self { inner: ReaderInner::Plain { path: path.clone(), index, format } })
        }
    }

    /// Return coverage data for `contig`, or `None` if the contig is absent from the file.
    ///
    /// For gzip readers, the contig is removed from the in-memory map after the first call
    /// so that memory is freed as contigs are consumed in FASTA order.
    pub fn get(&mut self, contig: &str) -> Result<Option<CoverageData>, GenGcBiasModelError> {
        match &mut self.inner {
            ReaderInner::Plain { path, index, format } => {
                CoverageData::load(path, index, contig, *format)
            }
            ReaderInner::Gzip { data } => {
                Ok(data.remove(contig))
            }
        }
    }
}

/// Load all contigs from a gzip-compressed coverage file in one forward pass.
fn load_all_gzip(
    path: &PathBuf,
    format: CoverageFormat,
) -> Result<HashMap<String, CoverageData>, GenGcBiasModelError> {
    let mut raw: HashMap<String, Vec<u32>> = HashMap::new();

    for line_result in read_gzip_lines(path)? {
        let line = line_result?;
        if line.is_empty() { continue; }

        let contig = match line.split('\t').next() {
            Some(c) if !c.is_empty() => c.to_string(),
            _ => continue,
        };

        let depths = raw.entry(contig).or_default();
        append_depth_record(&line, format, depths)?;
    }

    Ok(raw.into_iter()
        .map(|(k, depths)| (k, CoverageData { depths }))
        .collect())
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use flate2::{write::GzEncoder, Compression};
    use tempfile::NamedTempFile;

    fn write_plain(content: &str) -> NamedTempFile {
        let mut f = NamedTempFile::new().unwrap();
        write!(f, "{}", content).unwrap();
        f
    }

    fn write_gzip(content: &str) -> tempfile::TempPath {
        let f = tempfile::Builder::new().suffix(".gz").tempfile().unwrap();
        let mut enc = GzEncoder::new(f.as_file(), Compression::default());
        enc.write_all(content.as_bytes()).unwrap();
        enc.finish().unwrap();
        f.into_temp_path()
    }

    // ── CoverageIndex / CoverageData::load (internal, plain-text) ────────────

    #[test]
    fn test_index_finds_contig_offsets() {
        let f = write_plain("chr1\t1\t10\nchr1\t2\t20\nchr2\t1\t5\n");
        let index = CoverageIndex::build(&f.path().to_path_buf()).unwrap();
        assert!(index.offsets.contains_key("chr1"));
        assert!(index.offsets.contains_key("chr2"));
        assert!(!index.offsets.contains_key("chr3"));
    }

    #[test]
    fn test_load_samtools_depth_1based() {
        let f = write_plain("chr1\t1\t10\nchr1\t2\t20\nchr1\t3\t30\n");
        let path = f.path().to_path_buf();
        let index = CoverageIndex::build(&path).unwrap();
        let cov = CoverageData::load(&path, &index, "chr1", CoverageFormat::SamtoolsDepth)
            .unwrap().unwrap();
        assert_eq!(cov.depth_at(0), 10);
        assert_eq!(cov.depth_at(1), 20);
        assert_eq!(cov.depth_at(2), 30);
        assert_eq!(cov.depth_at(3), 0);
    }

    #[test]
    fn test_load_bedtools_genomecov_dz_0based_sparse() {
        let f = write_plain("chr1\t0\t15\nchr1\t2\t25\n");
        let path = f.path().to_path_buf();
        let index = CoverageIndex::build(&path).unwrap();
        let cov = CoverageData::load(&path, &index, "chr1", CoverageFormat::BedtoolsGenomecovDz)
            .unwrap().unwrap();
        assert_eq!(cov.depth_at(0), 15);
        assert_eq!(cov.depth_at(1), 0);
        assert_eq!(cov.depth_at(2), 25);
    }

    #[test]
    fn test_load_returns_none_for_missing_contig() {
        let f = write_plain("chr1\t1\t10\n");
        let path = f.path().to_path_buf();
        let index = CoverageIndex::build(&path).unwrap();
        let result = CoverageData::load(&path, &index, "chr2", CoverageFormat::SamtoolsDepth).unwrap();
        assert!(result.is_none());
    }

    #[test]
    fn test_load_multiple_contigs_independently() {
        let f = write_plain("chr1\t1\t100\nchr1\t2\t200\nchr2\t1\t50\nchr2\t2\t60\n");
        let path = f.path().to_path_buf();
        let index = CoverageIndex::build(&path).unwrap();
        let cov1 = CoverageData::load(&path, &index, "chr1", CoverageFormat::SamtoolsDepth)
            .unwrap().unwrap();
        let cov2 = CoverageData::load(&path, &index, "chr2", CoverageFormat::SamtoolsDepth)
            .unwrap().unwrap();
        assert_eq!(cov1.depth_at(0), 100);
        assert_eq!(cov1.depth_at(1), 200);
        assert_eq!(cov2.depth_at(0), 50);
        assert_eq!(cov2.depth_at(1), 60);
    }

    #[test]
    fn test_load_rejects_malformed_line() {
        let f = write_plain("chr1\tnot_a_number\t10\n");
        let path = f.path().to_path_buf();
        let index = CoverageIndex::build(&path).unwrap();
        let result = CoverageData::load(&path, &index, "chr1", CoverageFormat::SamtoolsDepth);
        assert!(result.is_err());
    }

    // ── CoverageReader — plain path ───────────────────────────────────────────

    #[test]
    fn test_coverage_reader_plain_returns_correct_depths() {
        let f = write_plain("chr1\t1\t10\nchr1\t2\t20\nchr2\t1\t50\n");
        let path = f.path().to_path_buf();
        let mut reader = CoverageReader::open(&path, CoverageFormat::SamtoolsDepth).unwrap();

        let cov1 = reader.get("chr1").unwrap().unwrap();
        assert_eq!(cov1.depth_at(0), 10);
        assert_eq!(cov1.depth_at(1), 20);

        let cov2 = reader.get("chr2").unwrap().unwrap();
        assert_eq!(cov2.depth_at(0), 50);

        assert!(reader.get("chrX").unwrap().is_none());
    }

    // ── CoverageReader — gzip path ────────────────────────────────────────────

    #[test]
    fn test_coverage_reader_gzip_matches_plain() {
        let content = "chr1\t1\t10\nchr1\t2\t20\nchr1\t3\t30\nchr2\t1\t50\nchr2\t2\t60\n";
        let plain_f = write_plain(content);
        let gz_path = write_gzip(content);

        let mut plain_reader = CoverageReader::open(
            &plain_f.path().to_path_buf(), CoverageFormat::SamtoolsDepth
        ).unwrap();
        let mut gz_reader = CoverageReader::open(
            &gz_path.to_path_buf(), CoverageFormat::SamtoolsDepth
        ).unwrap();

        for contig in ["chr1", "chr2"] {
            let plain_cov = plain_reader.get(contig).unwrap().unwrap();
            let gz_cov = gz_reader.get(contig).unwrap().unwrap();
            for pos in 0..3 {
                assert_eq!(
                    plain_cov.depth_at(pos), gz_cov.depth_at(pos),
                    "depth mismatch at {}:{}", contig, pos
                );
            }
        }
    }

    #[test]
    fn test_coverage_reader_gzip_returns_none_for_missing_contig() {
        let gz_path = write_gzip("chr1\t1\t10\n");
        let mut reader = CoverageReader::open(
            &gz_path.to_path_buf(), CoverageFormat::SamtoolsDepth
        ).unwrap();
        assert!(reader.get("chrX").unwrap().is_none());
    }

    #[test]
    fn test_coverage_reader_gzip_evicts_after_get() {
        let gz_path = write_gzip("chr1\t1\t10\n");
        let mut reader = CoverageReader::open(
            &gz_path.to_path_buf(), CoverageFormat::SamtoolsDepth
        ).unwrap();
        assert!(reader.get("chr1").unwrap().is_some());
        // Second call — contig was removed from the map after the first get
        assert!(reader.get("chr1").unwrap().is_none());
    }

    #[test]
    fn test_coverage_reader_gzip_zero_based_format() {
        let gz_path = write_gzip("chr1\t0\t15\nchr1\t2\t25\n");
        let mut reader = CoverageReader::open(
            &gz_path.to_path_buf(), CoverageFormat::BedtoolsGenomecovDz
        ).unwrap();
        let cov = reader.get("chr1").unwrap().unwrap();
        assert_eq!(cov.depth_at(0), 15);
        assert_eq!(cov.depth_at(1), 0); // gap
        assert_eq!(cov.depth_at(2), 25);
    }
}

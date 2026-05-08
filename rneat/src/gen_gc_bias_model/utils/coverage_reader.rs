use std::{
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader, Read, Seek, SeekFrom},
    path::PathBuf,
};
use crate::gen_gc_bias_model::{errors::GenGcBiasModelError, utils::config::CoverageFormat};

/// Byte-offset index into a plain-text coverage file.
///
/// Built in a single forward pass. For each contig, stores the start and end byte
/// offsets so we can seek directly to that contig's data when needed.
pub struct CoverageIndex {
    offsets: HashMap<String, (u64, u64)>,
}

impl CoverageIndex {
    pub fn build(path: &PathBuf) -> Result<Self, GenGcBiasModelError> {
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

    pub fn contains(&self, contig: &str) -> bool {
        self.offsets.contains_key(contig)
    }
}

/// Per-contig depth array (0-indexed).
///
/// Positions not present in the coverage file (gaps in bedtools-genomecov-dz output)
/// are stored implicitly as 0 and returned as 0 via `depth_at`.
pub struct CoverageData {
    depths: Vec<u32>,
}

impl CoverageData {
    /// Load depth data for a single contig. Returns `None` if the contig is absent
    /// from the index.
    pub fn load(
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

            // Fill gaps with 0 (handles sparse bedtools-genomecov-dz output)
            if idx >= depths.len() {
                depths.resize(idx, 0);
                depths.push(depth);
            } else {
                depths[idx] = depth;
            }
        }

        Ok(Some(CoverageData { depths }))
    }

    /// Depth at a 0-based position. Returns 0 for positions beyond the stored range.
    #[inline]
    pub fn depth_at(&self, pos: usize) -> u32 {
        self.depths.get(pos).copied().unwrap_or(0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn write_coverage(content: &str) -> NamedTempFile {
        let mut f = NamedTempFile::new().unwrap();
        write!(f, "{}", content).unwrap();
        f
    }

    #[test]
    fn test_index_finds_contig_offsets() {
        let f = write_coverage("chr1\t1\t10\nchr1\t2\t20\nchr2\t1\t5\n");
        let index = CoverageIndex::build(&f.path().to_path_buf()).unwrap();
        assert!(index.contains("chr1"));
        assert!(index.contains("chr2"));
        assert!(!index.contains("chr3"));
    }

    #[test]
    fn test_load_samtools_depth_1based() {
        // samtools depth: 1-based positions
        let f = write_coverage("chr1\t1\t10\nchr1\t2\t20\nchr1\t3\t30\n");
        let path = f.path().to_path_buf();
        let index = CoverageIndex::build(&path).unwrap();
        let cov = CoverageData::load(&path, &index, "chr1", CoverageFormat::SamtoolsDepth)
            .unwrap()
            .unwrap();
        assert_eq!(cov.depth_at(0), 10);
        assert_eq!(cov.depth_at(1), 20);
        assert_eq!(cov.depth_at(2), 30);
        assert_eq!(cov.depth_at(3), 0); // beyond stored range
    }

    #[test]
    fn test_load_bedtools_genomecov_dz_0based_sparse() {
        // bedtools-genomecov-dz: 0-based, only nonzero depths listed
        let f = write_coverage("chr1\t0\t15\nchr1\t2\t25\n"); // position 1 is absent (depth 0)
        let path = f.path().to_path_buf();
        let index = CoverageIndex::build(&path).unwrap();
        let cov = CoverageData::load(&path, &index, "chr1", CoverageFormat::BedtoolsGenomecovDz)
            .unwrap()
            .unwrap();
        assert_eq!(cov.depth_at(0), 15);
        assert_eq!(cov.depth_at(1), 0); // gap filled with 0
        assert_eq!(cov.depth_at(2), 25);
    }

    #[test]
    fn test_load_returns_none_for_missing_contig() {
        let f = write_coverage("chr1\t1\t10\n");
        let path = f.path().to_path_buf();
        let index = CoverageIndex::build(&path).unwrap();
        let result = CoverageData::load(&path, &index, "chr2", CoverageFormat::SamtoolsDepth).unwrap();
        assert!(result.is_none());
    }

    #[test]
    fn test_load_multiple_contigs_independently() {
        let f = write_coverage("chr1\t1\t100\nchr1\t2\t200\nchr2\t1\t50\nchr2\t2\t60\n");
        let path = f.path().to_path_buf();
        let index = CoverageIndex::build(&path).unwrap();

        let cov1 = CoverageData::load(&path, &index, "chr1", CoverageFormat::SamtoolsDepth)
            .unwrap()
            .unwrap();
        let cov2 = CoverageData::load(&path, &index, "chr2", CoverageFormat::SamtoolsDepth)
            .unwrap()
            .unwrap();

        assert_eq!(cov1.depth_at(0), 100);
        assert_eq!(cov1.depth_at(1), 200);
        assert_eq!(cov2.depth_at(0), 50);
        assert_eq!(cov2.depth_at(1), 60);
    }

    #[test]
    fn test_load_rejects_malformed_line() {
        let f = write_coverage("chr1\tnot_a_number\t10\n");
        let path = f.path().to_path_buf();
        let index = CoverageIndex::build(&path).unwrap();
        let result = CoverageData::load(&path, &index, "chr1", CoverageFormat::SamtoolsDepth);
        assert!(result.is_err());
    }
}
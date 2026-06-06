//! Tag per-pass read names and concatenate FASTQs (replaces the `awk` step in
//! `tools/cancer_simulate.sh`).
//!
//! gen-reads names reads by genomic position (`RNEAT_generated_<contig>_<start>_
//! <end>/<mate>`), so a normal and a tumor read sampled at the same coordinate
//! get IDENTICAL names. On chr22 at 30x combined we saw ~456k collisions out of
//! ~9.6M reads — which Picard MarkDuplicates silently drops (halving effective
//! coverage at those loci). Prefixing each pass's read names (`N_`/`T_`) fixes it.

use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use flate2::Compression;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;

use crate::gen_cancer_reads::errors::GenCancerReadsError;

/// Stream each `(path, tag)` input, prefix the header line of every 4-line FASTQ
/// record with `<tag>_`, and write the concatenation to `out` (gzip).
pub fn tag_and_concat(inputs: &[(&Path, &str)], out: &Path) -> Result<(), GenCancerReadsError> {
    let f = File::create(out)?;
    let mut w = GzEncoder::new(BufWriter::new(f), Compression::default());
    for (path, tag) in inputs {
        let r = BufReader::new(MultiGzDecoder::new(File::open(path)?));
        for (i, line) in r.lines().enumerate() {
            let line = line?;
            // Record header is line 1 of every 4 (NR % 4 == 1). Quality lines
            // (line 4) can also start with '@', so key on position, not char.
            if i % 4 == 0 {
                match line.strip_prefix('@') {
                    Some(rest) => writeln!(w, "@{tag}_{rest}")?,
                    None => writeln!(w, "{line}")?, // defensive passthrough
                }
            } else {
                writeln!(w, "{line}")?;
            }
        }
    }
    w.finish()?;
    Ok(())
}

//! End-to-end test for the BND/INV breakpoint double-counting fix.
//!
//! The chimeric pass emits junction-spanning reads for every BND/INV junction.
//! Before the fix, the regular per-contig pass *also* covered those junctions
//! from the unbroken reference (BND/INV are coverage-neutral), so a homozygous
//! junction sat at ~2x coverage. The fix drops the broken-allele fraction of
//! regular pairs that cross a junction; for a homozygous junction that's all of
//! them, so afterward NO regular read (`RNEAT_generated_`) should span the
//! junction, while interior / flank positions keep normal coverage and the
//! chimeric junction reads (`RNEAT_chimeric_`) are still emitted.

mod common;

use common::{GenReadsConfig, fresh_workdir, h1n1_reference, rneat};
use flate2::read::MultiGzDecoder;
use std::io::{BufRead, BufReader};
use std::io::Write as _;

const READ_LEN: usize = 151;

/// Run a homozygous INV at H1N1_HA:[start,end] (1-based) and return all
/// FASTQ read-name lines (line 1 of every 4-line record).
fn run_hom_inv(test_name: &str, sv_start: usize, sv_end: usize) -> Vec<String> {
    let (_dir, work) = fresh_workdir();
    let input_vcf = work.join("input_inv.vcf");
    {
        let mut f = std::fs::File::create(&input_vcf).unwrap();
        writeln!(f, "##fileformat=VCFv4.2").unwrap();
        writeln!(f, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">").unwrap();
        writeln!(f, "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"t\">").unwrap();
        writeln!(f, "##INFO=<ID=END,Number=1,Type=Integer,Description=\"e\">").unwrap();
        writeln!(f, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS").unwrap();
        writeln!(
            f,
            "H1N1_HA\t{sv_start}\t.\tA\t<INV>\t60\tPASS\tSVTYPE=INV;END={sv_end}\tGT\t1/1"
        )
        .unwrap();
    }

    let mut config = GenReadsConfig::new(h1n1_reference(), work.clone(), test_name);
    config.coverage = 100;
    config.read_len = READ_LEN;
    config.produce_fastq = true;
    config.input_vcf = Some(input_vcf);
    let yaml = config.write_yaml();
    rneat().args(["gen-reads", "-c"]).arg(yaml.path()).assert().success();

    let out = work.join(format!("{test_name}_r1.fastq.gz"));
    assert!(out.exists(), "FASTQ not produced at {:?}", out);
    let r = BufReader::new(MultiGzDecoder::new(std::fs::File::open(&out).unwrap()));
    r.lines()
        .enumerate()
        .filter(|(i, _)| i % 4 == 0)
        .map(|(_, l)| l.unwrap())
        .collect()
}

/// Count regular (`RNEAT_generated_`) reads on H1N1_HA whose span
/// `[abs_start, abs_end)` contains the 0-based position `pos`. Read names are
/// `RNEAT_generated_H1N1_HA_<abs_start>_<abs_end>_<hex>/1`.
fn regular_reads_covering(qnames: &[String], pos: usize) -> usize {
    qnames
        .iter()
        .filter(|l| l.contains("RNEAT_generated_H1N1_HA_"))
        .filter_map(|l| {
            // strip leading '@' and trailing '/1'
            let name = l.trim_start_matches('@').split('/').next()?;
            let parts: Vec<&str> = name.split('_').collect();
            // ...generated_H1N1_HA_<start>_<end>_<hex>
            let n = parts.len();
            if n < 3 {
                return None;
            }
            let start: usize = parts[n - 3].parse().ok()?;
            let end: usize = parts[n - 2].parse().ok()?;
            Some((start, end))
        })
        .filter(|&(start, end)| start <= pos && pos < end)
        .count()
}

#[test]
fn homozygous_inv_junctions_have_no_regular_crossing_reads() {
    // INV at H1N1_HA 1-based [400, 799] → 0-based junctions at 399 and 798.
    let qnames = run_hom_inv("hom_inv_dc", 400, 799);

    // Sanity: the chimeric junction reads are still emitted (we only removed
    // the redundant regular reference reads, not the junction signal).
    let chimeric = qnames
        .iter()
        .filter(|l| l.contains("RNEAT_chimeric_INV_H1N1_HA_"))
        .count();
    assert!(chimeric > 0, "expected chimeric INV junction reads, got 0");

    // Interior / flank control positions keep normal coverage.
    let interior = regular_reads_covering(&qnames, 600); // inside the inversion
    let flank = regular_reads_covering(&qnames, 150); // left of the inversion
    assert!(interior > 0, "interior should retain regular coverage");
    assert!(flank > 0, "flank should retain regular coverage");

    // Homozygous → every regular pair crossing a junction is dropped (no RNG),
    // so EXACTLY zero regular reads span either breakpoint. The suppression
    // junctions match the chimeric pass: POS-1 (0-based start, 399) and the
    // stored END value (799) — see collect_bnd_inv_junctions.
    let cross_start = regular_reads_covering(&qnames, 399);
    let cross_end = regular_reads_covering(&qnames, 799);
    assert_eq!(
        cross_start, 0,
        "homozygous INV start junction (399) must have no regular crossing reads; \
         interior control = {interior}"
    );
    assert_eq!(
        cross_end, 0,
        "homozygous INV end junction (799) must have no regular crossing reads; \
         interior control = {interior}"
    );
}

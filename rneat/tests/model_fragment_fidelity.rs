//! Fidelity test: a fitted model's NUMBERS must actually reach gen-reads output.
//!
//! Two output-fidelity links are already covered elsewhere:
//!   * `gen_reads::...::test_paired_ended_insert_size_matches_model` proves an
//!     EXPLICIT `fragment_mean` reaches the output BAM's TLEN.
//!   * `gen_reads::...::test_input_vcf_snp_appears_in_bam_reads` proves an input
//!     VCF variant appears in the output reads.
//!
//! This proves the remaining link that the Delta `model_builders.sbatch` round-trip
//! only *smoke*-tests (it checks the reads are valid, not that they reflect the
//! model): a fragment_model FILE built by `gen-frag-length-model` from a real BAM
//! actually drives the simulated fragment sizes. This exercises the model-load +
//! precedence path (config `fragment_model:`), not the explicit-number path.

mod common;

use std::fs;
use std::io::Read as _;
use std::path::Path;

use common::{h1n1_reference, rneat};
use flate2::read::GzDecoder;
use noodles::bam;
use noodles::core::Position;
use noodles::sam::{
    self as sam,
    alignment::{
        RecordBuf,
        io::Write as _,
        record::{
            Flags, MappingQuality,
            cigar::{Op, op::Kind},
        },
        record_buf::{Cigar, Sequence},
    },
    header::record::value::{Map, map::ReferenceSequence},
};

/// Write a paired-end BAM whose fragments have the given `tlens` (one proper pair
/// per entry). Mirrors `model_parity::write_test_bam` but varies TLEN per pair so
/// the fitted Normal has a realistic non-zero st_dev.
fn write_paired_bam_with_tlens(path: &Path, contig: &[u8], contig_len: usize, tlens: &[usize], read_len: usize) {
    let header = sam::Header::builder()
        .add_reference_sequence(
            contig.to_vec(),
            Map::<ReferenceSequence>::new(std::num::NonZero::<usize>::new(contig_len).unwrap()),
        )
        .build();
    let file = fs::File::create(path).unwrap();
    let mut writer = bam::io::Writer::new(file);
    writer.write_header(&header).unwrap();

    let cigar: Cigar = [Op::new(Kind::Match, read_len)].into_iter().collect();
    let seq = vec![b'A'; read_len];
    let paired_first = Flags::SEGMENTED | Flags::FIRST_SEGMENT;
    let paired_second = Flags::SEGMENTED;

    for &tlen in tlens {
        let mate_start = tlen - read_len + 1;
        let mut r1 = RecordBuf::default();
        *r1.flags_mut() = paired_first;
        *r1.cigar_mut() = cigar.clone();
        *r1.reference_sequence_id_mut() = Some(0);
        *r1.alignment_start_mut() = Position::new(1);
        *r1.mate_reference_sequence_id_mut() = Some(0);
        *r1.mate_alignment_start_mut() = Position::new(mate_start);
        *r1.template_length_mut() = tlen as i32;
        *r1.sequence_mut() = Sequence::from(seq.as_slice());
        *r1.mapping_quality_mut() = Some(MappingQuality::try_from(30u8).unwrap());
        writer.write_alignment_record(&header, &r1).unwrap();

        let mut r2 = RecordBuf::default();
        *r2.flags_mut() = paired_second;
        *r2.cigar_mut() = cigar.clone();
        *r2.reference_sequence_id_mut() = Some(0);
        *r2.alignment_start_mut() = Position::new(mate_start);
        *r2.mate_reference_sequence_id_mut() = Some(0);
        *r2.mate_alignment_start_mut() = Position::new(1);
        *r2.template_length_mut() = -(tlen as i32);
        *r2.sequence_mut() = Sequence::from(seq.as_slice());
        *r2.mapping_quality_mut() = Some(MappingQuality::try_from(30u8).unwrap());
        writer.write_alignment_record(&header, &r2).unwrap();
    }
}

fn write_yaml(dir: &Path, body: &str) -> std::path::PathBuf {
    let p = dir.join(format!("cfg_{}.yml", body.len()));
    fs::write(&p, body).unwrap();
    p
}

/// Read the fitted mean out of a gzipped `FragmentLengthModel::Normal` JSON.
/// The enum is externally tagged, so the shape is `{"Normal":{"mean":…,"st_dev":…}}`.
fn fitted_normal_mean(path: &Path) -> f64 {
    let mut raw = String::new();
    GzDecoder::new(fs::File::open(path).unwrap())
        .read_to_string(&mut raw)
        .unwrap();
    let v: serde_json::Value = serde_json::from_str(&raw).unwrap();
    v["Normal"]["mean"]
        .as_f64()
        .expect("fragment model should be a Normal with a numeric mean")
}

/// Mean |TLEN| over R1 records of a paired-end BAM (each template counted once).
fn mean_r1_insert_size(bam_path: &Path) -> f64 {
    let mut reader = bam::io::Reader::new(fs::File::open(bam_path).unwrap());
    reader.read_header().unwrap();
    let mut tlens: Vec<f64> = Vec::new();
    for result in reader.records() {
        let record = result.unwrap();
        let flags = record.flags();
        if !flags.is_segmented() || !flags.is_first_segment() {
            continue;
        }
        let tlen = record.template_length().unsigned_abs() as f64;
        if tlen > 0.0 {
            tlens.push(tlen);
        }
    }
    assert!(!tlens.is_empty(), "no positive R1 TLENs in {}", bam_path.display());
    tlens.iter().sum::<f64>() / tlens.len() as f64
}

/// build a fragment model from a BAM whose fragments average 200 bp, then simulate
/// with THAT MODEL FILE and confirm the simulated insert sizes actually average
/// ~200 — not the built-in default. If gen-reads ignored the loaded model, the mean
/// would land elsewhere; ±20% around 200 is the same tolerance the explicit-mean
/// test uses.
#[test]
fn built_fragment_model_drives_output_insert_size() {
    let tmp = tempfile::tempdir().unwrap();

    // 1. input BAM: fragments centered at 200 bp with a real spread (mean 200, sd ~14).
    let tlens: Vec<usize> = [180usize, 190, 200, 210, 220]
        .iter()
        .copied()
        .cycle()
        .take(500)
        .collect();
    let in_bam = tmp.path().join("input.bam");
    write_paired_bam_with_tlens(&in_bam, b"H1N1_HA", 1700, &tlens, 75);

    // 2. build the fragment model.
    let model = tmp.path().join("frag.json.gz");
    let build_cfg = write_yaml(
        tmp.path(),
        &format!(
            "input_file: {}\noutput_file: {}\noverwrite_output: true\nmin_reads: 2\n",
            in_bam.display(),
            model.display()
        ),
    );
    rneat()
        .args(["gen-frag-length-model", "-c"])
        .arg(&build_cfg)
        .assert()
        .success();

    // builder side: the model captured the ~200 bp mean of the input.
    let fitted = fitted_normal_mean(&model);
    assert!(
        (fitted - 200.0).abs() < 5.0,
        "builder mis-fit: fitted mean {fitted:.1} not ~200 (input mean is exactly 200)"
    );

    // 3. simulate using the MODEL FILE — no explicit fragment_mean/st_dev.
    let sim_cfg = write_yaml(
        tmp.path(),
        &format!(
            "reference: {ref}\nread_len: 75\ncoverage: 30\nploidy: 1\npaired_ended: true\n\
             fragment_model: {model}\nproduce_bam: true\nproduce_fastq: false\nproduce_vcf: false\n\
             overwrite_output: true\noutput_dir: {out}\noutput_filename: rt\n\
             rng_seed: frag fidelity\nnum_threads: 1\n",
            ref = h1n1_reference().display(),
            model = model.display(),
            out = tmp.path().display(),
        ),
    );
    rneat()
        .args(["gen-reads", "-c"])
        .arg(&sim_cfg)
        .assert()
        .success();

    // 4. output side: simulated insert sizes reflect the fitted model, not a default.
    let out_bam = tmp.path().join("rt.bam");
    let out_mean = mean_r1_insert_size(&out_bam);
    eprintln!("[fidelity] fitted model mean = {fitted:.1} bp; simulated output insert mean = {out_mean:.1} bp");
    let tolerance = 0.20;
    assert!(
        (out_mean - fitted).abs() <= fitted * tolerance,
        "output insert mean {out_mean:.1} is not within ±20% of the fitted model mean \
         {fitted:.1} — a built fragment_model is NOT reaching gen-reads output"
    );
}

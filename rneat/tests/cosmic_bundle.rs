//! Regression test for the bundled COSMIC tumor model
//! (`tools/cosmic_v104_pancancer_model.json.gz`).
//!
//! Issue #217: the model must carry a populated `sv_model` component so
//! that `cancer_simulate.sh`'s tumor pass generates somatic symbolic SVs
//! when `--sv-rate-scale > 0`. Before v1.12.1 the field was `null` and
//! the tumor pass silently emitted zero `<DEL>` / `<DUP>` / `<CNV>` /
//! `<INV>` / `<BND>` records, leaving `cancer_sv_benchmark.sh` with an
//! empty somatic SV truth.
//!
//! Parameters are injected by `tools/inject_cancer_sv_model.py` from
//! PCAWG pan-cancer rates (Li et al., Nature 578, 2020) — see that file
//! for citations.

mod common;

use std::path::PathBuf;

use ::common::models::mutation_model::MutationModel;
use ::common::structs::variants::SvType;
use common::{GenReadsConfig, fresh_workdir, h1n1_reference, rneat};

fn cosmic_model_path() -> PathBuf {
    PathBuf::from(format!(
        "{}/../tools/cosmic_v104_pancancer_model.json.gz",
        env!("CARGO_MANIFEST_DIR"),
    ))
}

#[test]
fn bundled_cosmic_model_carries_sv_component() {
    let path = cosmic_model_path();
    assert!(
        path.exists(),
        "bundled COSMIC model missing at {}",
        path.display()
    );
    let model = MutationModel::from_file(&path).expect("deserialize COSMIC model");
    let sv = model.sv_model.as_ref().expect(
        "tools/cosmic_v104_pancancer_model.json.gz must have a non-null \
         sv_model — re-run tools/inject_cancer_sv_model.py (issue #217)",
    );

    // All six SV types must be present; otherwise type sampling silently
    // skips whichever is missing.
    for t in [
        SvType::Del,
        SvType::Dup,
        SvType::Cnv,
        SvType::Bnd,
        SvType::Inv,
        SvType::Ins,
    ] {
        assert!(
            sv.type_probabilities.contains_key(&t),
            "missing type_probabilities entry for {:?}",
            t
        );
        assert!(
            sv.length_log_normal.contains_key(&t),
            "missing length_log_normal entry for {:?}",
            t
        );
    }

    let prob_sum: f64 = sv.type_probabilities.values().sum();
    assert!(
        (prob_sum - 1.0).abs() < 1e-6,
        "type_probabilities sum to {prob_sum}, expected 1.0"
    );

    let cn_sum: f64 = sv.cnv_copy_number_distribution.values().sum();
    assert!(
        (cn_sum - 1.0).abs() < 1e-6,
        "cnv_copy_number_distribution sums to {cn_sum}, expected 1.0"
    );

    // Per-base rate should be cancer-typical (~1e-8), not gnomAD-SV
    // germline (~6e-4). A 10,000x difference is the canary that someone
    // mistakenly copied the germline defaults into the cancer bundle.
    assert!(
        sv.per_base_rate < 1e-6,
        "per_base_rate {} looks like a germline rate; cancer should be ≲ 1e-7",
        sv.per_base_rate
    );

    // BND should be over-represented vs germline (PCAWG ~35% vs gnomAD
    // ~19%) — this distinguishes the cancer-tuned bundle from someone
    // accidentally copying default_sv_model() in here.
    let bnd_prob = sv.type_probabilities[&SvType::Bnd];
    assert!(
        bnd_prob >= 0.20,
        "BND proportion {bnd_prob} too low for cancer; PCAWG ~0.35"
    );
}

/// End-to-end: with the bundled COSMIC model wired into gen-reads via
/// `mutation_model:` and a high `sv_rate_scale` to compensate for H1N1's
/// tiny contigs, at least one symbolic SV record must reach the output
/// VCF. Before issue #217 this silently produced zero symbolic SVs
/// because `sv_model` was null in the bundled JSON.
///
/// Why sv_rate_scale=1e6: H1N1 segments average ~1.7 kb and
/// `per_base_rate=3.8e-8` (cancer-typical), so λ at scale=1 is ≈ 7e-5
/// across the whole reference — essentially zero. Scaling to 1e6 gives
/// λ ≈ 70 per segment, which guarantees signal without changing any
/// other code paths.
#[test]
fn bundled_cosmic_model_drives_symbolic_sv_generation() {
    let (_dir, work) = fresh_workdir();

    let mut config = GenReadsConfig::new(h1n1_reference(), work.clone(), "cosmic_sv_smoke");
    config.coverage = 5;
    config.produce_fastq = false;
    config.produce_vcf = true;
    config.mutation_model = Some(cosmic_model_path());
    // High scale to make SV signal abundant on a tiny reference. The
    // bundled cancer per_base_rate is realistic at chr-scale but vanishes
    // on H1N1's 1.7 kb segments without amplification.
    config.sv_rate_scale = Some(1e6);
    config.rng_seed = "v1.12.1-cosmic-sv-smoke".to_string();

    let yaml = config.write_yaml();
    rneat()
        .args(["gen-reads", "-c"])
        .arg(yaml.path())
        .assert()
        .success();

    let out_vcf = work.join("cosmic_sv_smoke.vcf.gz");
    assert!(out_vcf.exists(), "VCF not produced at {:?}", out_vcf);

    let vcf_lines: Vec<String> = {
        use flate2::read::MultiGzDecoder;
        use std::io::{BufRead, BufReader};
        let r = BufReader::new(MultiGzDecoder::new(std::fs::File::open(&out_vcf).unwrap()));
        r.lines().map(|l| l.unwrap()).collect()
    };
    let body: Vec<&String> = vcf_lines.iter().filter(|l| !l.starts_with('#')).collect();
    assert!(!body.is_empty(), "no records emitted at all");

    // Count symbolic-ALT records (`<DEL>`, `<DUP>`, `<CNV>`, `<INV>`) and
    // BND-notation records (`t[p[`, `t]p]`, etc.). Literal Insertions
    // from sv_model also pass through, but checking those overlaps with
    // the SNP/indel path; symbolic + BND records are unambiguously from
    // the SV path.
    let mut symbolic = 0usize;
    let mut bnd = 0usize;
    for rec in &body {
        let fields: Vec<&str> = rec.split('\t').collect();
        if fields.len() < 5 {
            continue;
        }
        let alt = fields[4];
        if alt.starts_with('<') && alt.ends_with('>') {
            symbolic += 1;
        } else if alt.contains('[') || alt.contains(']') {
            bnd += 1;
        }
    }
    assert!(
        symbolic + bnd > 0,
        "expected ≥1 symbolic or BND SV in VCF with bundled COSMIC model + \
         sv_rate_scale=1e6; got 0. {} total records. \
         Check that tools/cosmic_v104_pancancer_model.json.gz has a non-null sv_model.",
        body.len()
    );
}

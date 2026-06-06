//! End-to-end test for the native `rneat gen-cancer-reads` subcommand (#239).
//!
//! Drives a tumor/normal simulation on H1N1 and checks the orchestration +
//! merges: per-pass golden VCFs + tagged/concatenated FASTQs + an origin-tagged
//! truth VCF (germline | somatic | shared).

mod common;

use common::{fresh_workdir, h1n1_reference, rneat};
use flate2::read::MultiGzDecoder;
use std::io::{BufRead, BufReader};

fn read_gz_lines(path: &std::path::Path) -> Vec<String> {
    let r = BufReader::new(MultiGzDecoder::new(std::fs::File::open(path).unwrap()));
    r.lines().map(|l| l.unwrap()).collect()
}

#[test]
fn gen_cancer_reads_produces_tagged_fastqs_and_origin_truth() {
    let (_dir, work) = fresh_workdir();
    let yaml = work.join("cancer.yml");
    // High per-pass mutation rates so H1N1's tiny (~14 kb) reference yields plenty
    // of germline + somatic SNVs to classify. SVs left off (finicky on H1N1's
    // short segments); the merge logic is type-agnostic.
    std::fs::write(
        &yaml,
        format!(
            "reference: {ref}\n\
             output_dir: {out}\n\
             output_prefix: ctest\n\
             total_coverage: 30\n\
             purity: 0.5\n\
             read_len: 70\n\
             paired_ended: true\n\
             fragment_mean: 250\n\
             fragment_st_dev: 30\n\
             normal_mutation_rate: 0.01\n\
             tumor_mutation_rate: 0.01\n\
             overwrite_output: true\n\
             rng_seed: cancer-e2e\n",
            ref = h1n1_reference().display(),
            out = work.display(),
        ),
    )
    .unwrap();

    rneat()
        .args(["gen-cancer-reads", "-c"])
        .arg(&yaml)
        .assert()
        .success();

    // Per-pass + merged outputs exist.
    for f in [
        "ctest_normal.vcf.gz",
        "ctest_tumor.vcf.gz",
        "ctest_merged_r1.fastq.gz",
        "ctest_merged_r2.fastq.gz",
        "ctest_merged_truth.vcf.gz",
    ] {
        assert!(work.join(f).is_file(), "expected output {f}");
    }

    // Merged FASTQ read names carry both N_ and T_ tags.
    let r1 = read_gz_lines(&work.join("ctest_merged_r1.fastq.gz"));
    let headers: Vec<&String> = r1.iter().step_by(4).collect();
    assert!(headers.iter().any(|h| h.starts_with("@N_")), "no N_-tagged reads");
    assert!(headers.iter().any(|h| h.starts_with("@T_")), "no T_-tagged reads");

    // Truth VCF: NEAT_ORIGIN present, with somatic AND shared (germline carried
    // through the tumor pass) classes represented.
    let truth = read_gz_lines(&work.join("ctest_merged_truth.vcf.gz"));
    assert!(
        truth.iter().any(|l| l.contains("##INFO=<ID=NEAT_ORIGIN")),
        "truth VCF missing NEAT_ORIGIN header declaration"
    );
    let body: Vec<&String> = truth.iter().filter(|l| !l.starts_with('#')).collect();
    let has = |tag: &str| body.iter().any(|l| l.contains(&format!("NEAT_ORIGIN={tag}")));
    assert!(has("somatic"), "no somatic records in truth");
    assert!(has("shared"), "no shared (germline-carried) records in truth");
    // every body record must be origin-tagged
    assert!(
        body.iter().all(|l| l.contains("NEAT_ORIGIN=")),
        "some truth records lack NEAT_ORIGIN"
    );
}

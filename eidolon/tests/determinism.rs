//! Determinism guarantees for the gen-reads pipeline.
//!
//! Reproducibility is a load-bearing promise of any read simulator — without it,
//! variant-calling pipelines can't be debugged or regression-tested. These tests run the
//! real binary twice with the same seed and assert the outputs match.
//!
//! The invariant we test is **same seed → same record multiset**: every record (name,
//! sequence, separator, quality) appears in both runs, with no drops, duplicates, or
//! payload differences. Line-by-line output *order* can vary even with `num_threads=1`
//! because eidolon iterates HashMap keys when assembling per-contig outputs, and Rust's
//! HashMap iteration is randomized per process. The multiset claim is the load-bearing
//! one for downstream consumers (variant callers, mappers) — they consume FASTQ records
//! as an unordered set anyway.
//!
//! We exercise the invariant under both num_threads=1 and the rayon default so a
//! single-thread-only determinism leak doesn't slip past us, and also assert that the
//! seed argument actually changes output (otherwise "determinism" is a coincidence).
//!
//! We compare decompressed FASTQ contents rather than raw gzip bytes: gzip headers can
//! embed mtime/OS metadata that varies across runs even when the payload is identical.

mod common;

use common::{GenReadsConfig, eidolon, fresh_workdir, h1n1_reference, read_gzip_fastq_lines};

fn run_gen_reads(seed: &str, output_dir: &std::path::Path, name: &str, threads: Option<usize>) {
    let mut config = GenReadsConfig::new(h1n1_reference(), output_dir.to_path_buf(), name);
    config.rng_seed = seed.to_string();
    config.num_threads = threads;
    let yaml = config.write_yaml();
    eidolon()
        .args(["gen-reads", "-c"])
        .arg(yaml.path())
        .assert()
        .success();
}

/// Group a flat FASTQ line list into 4-line records, then sort lexicographically. This
/// turns line order into a multiset comparison.
fn sorted_records(lines: Vec<String>) -> Vec<[String; 4]> {
    assert!(
        lines.len().is_multiple_of(4),
        "FASTQ line count must be multiple of 4"
    );
    let mut records: Vec<[String; 4]> = lines
        .chunks(4)
        .map(|c| [c[0].clone(), c[1].clone(), c[2].clone(), c[3].clone()])
        .collect();
    records.sort();
    records
}

#[test]
fn same_seed_single_threaded_produces_same_record_multiset() {
    // num_threads=1 removes inter-contig parallelism. If determinism leaks even here,
    // some shared HashMap or atomic counter is consuming RNG non-deterministically and
    // the issue is *not* just a benign per-thread reordering — it's a real bug.
    let (_a, dir_a) = fresh_workdir();
    let (_b, dir_b) = fresh_workdir();

    let seed = "eidolon single threaded";
    run_gen_reads(seed, &dir_a, "st_a", Some(1));
    run_gen_reads(seed, &dir_b, "st_b", Some(1));

    let a = sorted_records(read_gzip_fastq_lines(&dir_a.join("st_a_r1.fastq.gz")));
    let b = sorted_records(read_gzip_fastq_lines(&dir_b.join("st_b_r1.fastq.gz")));

    assert!(!a.is_empty(), "first run produced no records");
    assert_eq!(
        a.len(),
        b.len(),
        "same seed (single-threaded) produced different record counts: a={}, b={}",
        a.len(),
        b.len(),
    );
    assert_eq!(
        a, b,
        "same seed (single-threaded) produced different record multisets. \
         A record was dropped, duplicated, or its content differs.",
    );
}

#[test]
fn same_seed_multi_threaded_produces_same_record_multiset() {
    // With parallel processing the per-contig output order can vary, but the set of
    // records (including names and quality scores) must be identical. If this fails,
    // the RNG is being shared across threads in a way that races, or reads are being
    // dropped/duplicated under contention.
    let (_a, dir_a) = fresh_workdir();
    let (_b, dir_b) = fresh_workdir();

    let seed = "eidolon multi threaded";
    run_gen_reads(seed, &dir_a, "mt_a", None);
    run_gen_reads(seed, &dir_b, "mt_b", None);

    let a = sorted_records(read_gzip_fastq_lines(&dir_a.join("mt_a_r1.fastq.gz")));
    let b = sorted_records(read_gzip_fastq_lines(&dir_b.join("mt_b_r1.fastq.gz")));

    assert!(!a.is_empty(), "first run produced no records");
    assert_eq!(
        a.len(),
        b.len(),
        "same seed produced different record counts (multi-threaded): a={}, b={}",
        a.len(),
        b.len(),
    );
    assert_eq!(
        a, b,
        "same seed produced different record multisets (multi-threaded). \
         Some record was dropped, duplicated, or differs between runs.",
    );
}

#[test]
fn different_seeds_produce_different_output() {
    // Sanity check that the seed argument is actually load-bearing. If this fails, the
    // RNG isn't actually being seeded from the config string and "determinism" is a
    // coincidence rather than a guarantee.
    let (_a, dir_a) = fresh_workdir();
    let (_b, dir_b) = fresh_workdir();

    run_gen_reads("alpha seed text", &dir_a, "diff_a", Some(1));
    run_gen_reads("bravo seed text", &dir_b, "diff_b", Some(1));

    let a = read_gzip_fastq_lines(&dir_a.join("diff_a_r1.fastq.gz"));
    let b = read_gzip_fastq_lines(&dir_b.join("diff_b_r1.fastq.gz"));

    assert!(!a.is_empty() && !b.is_empty());
    assert_ne!(
        a, b,
        "different seeds produced identical FASTQ — seed argument is not load-bearing",
    );
}

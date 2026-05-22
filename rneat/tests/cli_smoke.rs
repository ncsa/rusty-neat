//! CLI-level smoke tests for the `rneat` binary.
//!
//! These tests exercise argument parsing and error propagation through the real binary,
//! not through `run_neat()` and friends directly. They catch breakages that pure-library
//! tests miss: a subcommand renamed in `main.rs`, a required flag dropped, an exit code
//! that silently flips to zero on error.

mod common;

use common::rneat;
use predicates::prelude::*;

/// Subcommands the binary must expose. If any one is missing from `rneat --help`, the
/// public CLI surface has regressed.
const EXPECTED_SUBCOMMANDS: &[&str] = &[
    "gen-reads",
    "filter-reads",
    "gen-mut-model",
    "gen-seq-error-model",
    "gen-frag-length-model",
    "gen-gc-bias-model",
    "gen-bam-models",
    "compare-vcfs",
];

#[test]
fn top_level_help_lists_every_subcommand() {
    let assert = rneat().arg("--help").assert().success();
    let output = std::str::from_utf8(&assert.get_output().stdout)
        .expect("--help must produce valid UTF-8")
        .to_string();
    for sub in EXPECTED_SUBCOMMANDS {
        assert!(
            output.contains(sub),
            "`rneat --help` did not list subcommand `{sub}`. Full output:\n{output}",
        );
    }
}

#[test]
fn each_subcommand_help_exits_clean() {
    // `<sub> --help` must exit 0 and mention the configuration-yaml flag — that flag is
    // the load-bearing input for every subcommand.
    for sub in EXPECTED_SUBCOMMANDS {
        rneat()
            .args([sub, "--help"])
            .assert()
            .success()
            .stdout(predicate::str::contains("--configuration-yaml"));
    }
}

#[test]
fn missing_config_file_exits_nonzero_with_error_message() {
    // Each subcommand asserts the config file exists before invoking the runner. A
    // missing path must surface as a non-zero exit *and* a stderr message that
    // identifies what went wrong, so users / CI logs aren't left guessing.
    let unique = format!("/nonexistent-rneat-smoke-{}.yml", std::process::id(),);
    rneat()
        .args(["gen-reads", "-c", &unique])
        .assert()
        .failure()
        .stderr(predicate::str::contains("Error"));
}

#[test]
fn no_arguments_prints_help_and_exits_nonzero() {
    // The binary is configured with `.arg_required_else_help(true)` at the outer
    // "rneat" subcommand, so running with no arguments must show help text rather
    // than panic, and exit non-zero (clap's standard behavior for the missing-arg case).
    rneat()
        .assert()
        .failure()
        .stderr(predicate::str::contains("SUB-COMMAND").or(predicate::str::contains("Usage")));
}

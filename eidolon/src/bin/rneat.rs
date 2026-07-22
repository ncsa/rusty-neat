//! Deprecation shim for the old `rneat` command name.
//!
//! `rneat` was renamed to `eidolon` in v2.0.0. This alias prints a one-line
//! deprecation warning to stderr and then forwards all arguments to the sibling
//! `eidolon` binary, so existing scripts keep working during the transition.
//! Remove this shim in a future release once callers have migrated.

use std::process::Command;

fn main() {
    eprintln!(
        "warning: `rneat` has been renamed to `eidolon` (v2.0.0). The `rneat` alias \
         is deprecated and will be removed in a future release — please invoke `eidolon` instead."
    );

    // The real binary is installed alongside this shim; fall back to PATH.
    let eidolon = std::env::current_exe()
        .ok()
        .and_then(|p| p.parent().map(|dir| dir.join("eidolon")))
        .unwrap_or_else(|| std::path::PathBuf::from("eidolon"));

    match Command::new(&eidolon)
        .args(std::env::args_os().skip(1))
        .status()
    {
        Ok(status) => std::process::exit(status.code().unwrap_or(1)),
        Err(e) => {
            eprintln!(
                "error: failed to launch `eidolon` ({}): {e}",
                eidolon.display()
            );
            std::process::exit(127);
        }
    }
}

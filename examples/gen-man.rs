//! Render the man page from the live CLI definition, so the two never disagree
//! about a flag's existence or its help text. Regenerate and commit after any
//! flag change with `cargo run --example gen-man`.
//!
//! Writes `man/whittle.1` at the repository root, or under a directory given as
//! the first argument (the release workflow passes one). An example, not a
//! binary, so `clap_mangen` stays a dev-dependency.

use std::io::Write;
use std::path::PathBuf;

fn main() -> std::io::Result<()> {
    let dir = std::env::args_os()
        .nth(1)
        .map_or_else(|| PathBuf::from("."), PathBuf::from)
        .join("man");
    std::fs::create_dir_all(&dir)?;

    let path = dir.join("whittle.1");
    let mut buf = Vec::new();
    clap_mangen::Man::new(whittle::cli::command()).render(&mut buf)?;
    std::fs::File::create(&path)?.write_all(&buf)?;
    println!("{}", path.display());
    Ok(())
}

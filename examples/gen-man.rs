//! Renders the man page from the CLI definition, so the two agree on every flag
//! and its help text. `cargo run --example gen-man` regenerates it after a flag
//! change; the result is committed.
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

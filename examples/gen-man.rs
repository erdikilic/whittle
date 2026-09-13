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
    // clap_mangen prefixes the version with `v`; releases are bare numbers.
    let page = String::from_utf8(buf)
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?
        .replace("\nv0.", "\n0.")
        .replace("\nv1.", "\n1.");
    let mut file = std::fs::File::create(&path)?;
    for line in page.lines() {
        writeln!(file, "{}", line.trim_end())?;
    }
    println!("{}", path.display());
    Ok(())
}

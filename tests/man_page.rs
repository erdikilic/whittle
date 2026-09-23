//! The committed man page is rendered from the CLI definition and must match
//! it; `cargo run --example gen-man` regenerates it after a flag change.

/// Renders the man page the way `examples/gen-man.rs` does.
fn rendered() -> String {
    let mut buf = Vec::new();
    clap_mangen::Man::new(whittle::cli::command())
        .render(&mut buf)
        .unwrap();
    let version = env!("CARGO_PKG_VERSION");
    let page = String::from_utf8(buf)
        .unwrap()
        .replace(&format!("\nv{version}"), &format!("\n{version}"));
    let mut out = String::new();
    for line in page.lines() {
        out.push_str(line.trim_end());
        out.push('\n');
    }
    out
}

#[test]
fn committed_man_page_matches_the_cli() {
    let committed =
        std::fs::read_to_string(concat!(env!("CARGO_MANIFEST_DIR"), "/man/whittle.1")).unwrap();
    assert!(
        committed == rendered(),
        "man/whittle.1 is stale; run `cargo run --example gen-man`"
    );
}

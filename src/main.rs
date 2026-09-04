//! Binary entry point: parses the CLI, initializes observability, runs the
//! workflow, and maps failure to exit status 1 (2 for a CLI parse error).

fn main() {
    // On x86-64 builds compiled for AVX2 (the crate's target-cpu=x86-64-v3
    // default), the running CPU is checked for AVX2 before any SIMD code
    // executes, so an unsupported CPU exits with a message rather than a SIGILL.
    // The check compiles to a no-op on non-AVX2 and aarch64 builds.
    ensure_simd::ensure_simd();

    let mut cfg = match whittle::cli::parse() {
        Ok(c) => c,
        Err(e) => {
            eprintln!("{e:#}");
            std::process::exit(2);
        },
    };
    let mut obs = whittle::obs::init(&mut cfg);

    // The version and command lines precede every other line (the banner, the
    // clamp, mismatch and no-op warnings, and an early hard-error exit), so they
    // are the first two lines of every log. Line mode only; bar mode emits its
    // own one-line start. `args_os` (not `args`, which panics on non-UTF-8 argv)
    // feeds `command_line`.
    if obs.shows_lines() {
        tracing::info!("whittle {}", env!("CARGO_PKG_VERSION"));
        tracing::info!("Command: {}", whittle::command_line(std::env::args_os()));
    }

    let start = std::time::Instant::now();
    if let Err(e) = whittle::run(cfg, &mut obs) {
        // A downstream reader that stops early (`whittle ... | head`) closes the
        // pipe; the failed write is the reader's decision, not a fault in this
        // run. Exiting 0 without a message follows the shell convention for a
        // producer cut off by its consumer.
        let cut_off = e.chain().any(|cause| {
            cause
                .downcast_ref::<std::io::Error>()
                .is_some_and(|io| io.kind() == std::io::ErrorKind::BrokenPipe)
        });
        if !cut_off {
            tracing::error!(
                "Failed after {}: {e:#}",
                whittle::obs::human_dur(start.elapsed())
            );
        }
        // `Drop` stops the ticker and clears any live bar. `process::exit` runs
        // no destructors, so without this drop a mid-run failure would leave a
        // stale bar frame on the terminal in place of the `Failed after` line.
        drop(obs);
        std::process::exit(if cut_off { 0 } else { 1 });
    }
}

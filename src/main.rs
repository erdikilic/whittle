fn main() {
    // On x86-64 builds compiled with AVX2 (the default via the crate's
    // target-cpu=x86-64-v3 config), verify the running CPU actually supports
    // AVX2 before any SIMD code runs, exiting with a clear message instead of
    // a SIGILL. Compiles to a no-op on non-AVX2 / aarch64 builds.
    ensure_simd::ensure_simd();

    let cfg = match whittle::cli::parse() {
        Ok(c) => c,
        Err(e) => {
            eprintln!("{e:#}");
            std::process::exit(2);
        },
    };
    let mut obs = whittle::obs::init(cfg.verbosity, cfg.quiet);

    // `whittle {version}`/`Command: ...` must print before everything else: the
    // banner, any clamp/mismatch/no-op warning, and even an early hard-error bail,
    // so a reader always finds them at the top. Line mode only; bar mode gets its
    // own one-line start instead. `args_os` (not `args`, which panics on non-UTF-8
    // argv) feeds `command_line`.
    if obs.shows_lines() {
        tracing::info!("whittle {}", env!("CARGO_PKG_VERSION"));
        tracing::info!("{}", whittle::command_line(std::env::args_os()));
    }

    let start = std::time::Instant::now();
    if let Err(e) = whittle::run(cfg, &mut obs) {
        tracing::error!(
            "Failed after {}: {e:#}",
            whittle::obs::human_dur(start.elapsed())
        );
        // Explicit drop (not just letting `obs` fall out of scope at `main`'s
        // end): `Drop` stops the ticker and clears any live bar, and that must
        // happen BEFORE `process::exit` (which terminates immediately and
        // runs no destructors) or a mid-run failure would leave a stale bar
        // frame on the terminal instead of a clean "Failed after ..." line.
        drop(obs);
        std::process::exit(1);
    }
}

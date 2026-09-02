fn main() {
    // On x86-64 builds compiled with AVX2 (the default via the crate's
    // target-cpu=x86-64-v3 config), verify the running CPU supports AVX2
    // before any SIMD code runs, exiting with a clear message instead of a
    // SIGILL. Compiles to a no-op on non-AVX2 / aarch64 builds.
    ensure_simd::ensure_simd();

    let mut cfg = match whittle::cli::parse() {
        Ok(c) => c,
        Err(e) => {
            eprintln!("{e:#}");
            std::process::exit(2);
        },
    };
    let mut obs = whittle::obs::init(&mut cfg);

    // `whittle {version}`/`Command: ...` must print before everything else: the
    // banner, any clamp/mismatch/no-op warning, and even an early hard-error bail,
    // so a reader always finds them at the top. Line mode only; bar mode gets its
    // own one-line start instead. `args_os` (not `args`, which panics on non-UTF-8
    // argv) feeds `command_line`.
    if obs.shows_lines() {
        tracing::info!("whittle {}", env!("CARGO_PKG_VERSION"));
        tracing::info!("Command: {}", whittle::command_line(std::env::args_os()));
    }

    let start = std::time::Instant::now();
    if let Err(e) = whittle::run(cfg, &mut obs) {
        // A downstream reader that stops early (`whittle ... | head`) closes the
        // pipe, and the write that fails is the reader's decision, not a fault
        // in this run. Exiting 0 without a message matches the shell convention
        // for a producer cut off by its consumer.
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
        // `Drop` stops the ticker and clears any live bar, and that must
        // happen before `process::exit`, which terminates immediately and runs
        // no destructors; otherwise a mid-run failure leaves a stale bar frame
        // on the terminal instead of a clean "Failed after ..." line.
        drop(obs);
        std::process::exit(if cut_off { 0 } else { 1 });
    }
}

fn main() {
    let cfg = match whittle::cli::parse() {
        Ok(c) => c,
        Err(e) => {
            eprintln!("{e:#}");
            std::process::exit(2);
        },
    };
    let mut obs = whittle::obs::init(cfg.verbosity, cfg.quiet);

    // `whittle {version}`/`Command: ...` must be the very first thing printed —
    // even before the resolved-config banner and any clamp/mismatch/no-op
    // warning `run` emits, and even before an early hard-error bail — so a
    // reader can always find them at the top of a run's output. Line mode
    // only: bar mode gets exactly its own one-line start (see `run`) so the
    // live bar stays otherwise clean. `args_os` (not `args`, which panics on
    // non-UTF-8 argv) feeds `command_line`, which lossily converts and
    // shell-quotes each argument.
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
        // happen BEFORE `process::exit` — which terminates immediately and
        // runs no destructors — or a mid-run failure would leave a stale bar
        // frame on the terminal instead of a clean "Failed after ..." line.
        drop(obs);
        std::process::exit(1);
    }
}

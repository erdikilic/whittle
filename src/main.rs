fn main() {
    let cfg = match whittle::cli::parse() {
        Ok(c) => c,
        Err(e) => {
            eprintln!("{e:#}");
            std::process::exit(2);
        },
    };
    let mut obs = whittle::obs::init(cfg.verbosity, cfg.quiet);
    let start = std::time::Instant::now();
    if let Err(e) = whittle::run(cfg, &mut obs) {
        tracing::error!(
            "Failed after {}: {e:#}",
            whittle::obs::human_dur(start.elapsed())
        );
        std::process::exit(1);
    }
}

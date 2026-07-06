fn main() -> anyhow::Result<()> {
    let cfg = whittle::cli::parse()?;
    let obs = whittle::obs::init(cfg.verbosity, cfg.quiet);
    whittle::run(cfg, &obs)
}

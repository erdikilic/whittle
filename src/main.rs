fn main() -> anyhow::Result<()> {
    let cfg = whittle::cli::parse()?;
    let mut obs = whittle::obs::init(cfg.verbosity, cfg.quiet);
    whittle::run(cfg, &mut obs)
}

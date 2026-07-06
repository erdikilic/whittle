fn main() -> anyhow::Result<()> {
    let cfg = whittle::cli::parse()?;
    whittle::run(cfg)
}

fn main() -> anyhow::Result<()> {
    let cfg = chopping::cli::parse()?;
    chopping::run(cfg)
}

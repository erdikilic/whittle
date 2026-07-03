use crate::config::Config;

// Minimal compile bridge: the real clap-based CLI arrives in Task 8. `Config`
// no longer has a trivial default (it now carries FilterConfig/TrimPlan), so
// there is nothing sensible to construct here yet.
pub fn parse() -> anyhow::Result<Config> {
    anyhow::bail!("CLI wired in Task 8")
}

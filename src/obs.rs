//! Observability: leveled logging (tracing) and progress reporting (indicatif).

use tracing::level_filters::LevelFilter;

/// Map the CLI verbosity/quiet flags to a tracing level. `WHITTLE_LOG`, when set,
/// is applied separately (in `init`) and takes precedence over this.
pub fn level_from(verbosity: u8, quiet: bool) -> LevelFilter {
    if quiet {
        LevelFilter::WARN
    } else {
        match verbosity {
            0 => LevelFilter::INFO,
            1 => LevelFilter::DEBUG,
            _ => LevelFilter::TRACE,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn level_mapping() {
        assert_eq!(level_from(0, false), LevelFilter::INFO);
        assert_eq!(level_from(1, false), LevelFilter::DEBUG);
        assert_eq!(level_from(2, false), LevelFilter::TRACE);
        assert_eq!(level_from(5, false), LevelFilter::TRACE);
        // quiet wins over verbosity
        assert_eq!(level_from(0, true), LevelFilter::WARN);
        assert_eq!(level_from(3, true), LevelFilter::WARN);
    }
}

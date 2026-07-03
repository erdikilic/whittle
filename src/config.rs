use std::path::PathBuf;

use crate::filter::FilterConfig;
use crate::io::Format;
use crate::trim::TrimPlan;

#[derive(Debug, Clone)]
pub struct IoConfig {
    pub input: Option<PathBuf>,
    pub output: Option<PathBuf>,
    pub in_format: Option<Format>,
    pub out_format: Option<Format>,
}

#[derive(Debug, Clone)]
pub struct Config {
    pub io: IoConfig,
    pub filter: FilterConfig,
    pub trim: TrimPlan,
    pub threads: usize,
}

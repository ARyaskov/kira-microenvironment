use clap::{Args, Parser, Subcommand, ValueEnum};
use std::path::PathBuf;

#[derive(Debug, Parser)]
#[command(name = "kira-microenvironment")]
#[command(about = "Microenvironment LR graph builder")]
#[command(version)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Command,
}

#[derive(Debug, Subcommand)]
pub enum Command {
    Run(RunArgs),
    Validate(RunArgs),
}

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum AggModeArg {
    Median,
    TrimmedMean,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum CapModeArg {
    Fixed,
    P99,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum SpecModeArg {
    On,
    Off,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum AutoGroupsModeArg {
    Flat,
    Hierarchical,
}

#[derive(Debug, Clone, Args)]
pub struct RunArgs {
    #[arg(long)]
    pub expr: PathBuf,
    #[arg(long)]
    pub barcodes: PathBuf,
    #[arg(long)]
    pub groups: Option<PathBuf>,
    #[arg(long)]
    pub auto_groups: Option<PathBuf>,
    #[arg(long)]
    pub auto_groups_coarse: Option<PathBuf>,
    #[arg(long)]
    pub auto_groups_fine: Option<PathBuf>,
    #[arg(long)]
    pub auto_groups_anti: Option<PathBuf>,
    #[arg(long, value_enum, default_value_t = AutoGroupsModeArg::Flat)]
    pub auto_groups_mode: AutoGroupsModeArg,
    #[arg(long, default_value_t = 1e-6)]
    pub auto_groups_eps: f32,
    #[arg(long, default_value_t = 0.1)]
    pub auto_groups_min_delta: f32,
    #[arg(long, default_value = "unknown")]
    pub auto_groups_unknown: String,
    #[arg(long, default_value_t = false)]
    pub auto_groups_emit_scores: bool,
    #[arg(long)]
    pub resources: PathBuf,
    #[arg(long)]
    pub secretion: Option<PathBuf>,
    #[arg(long)]
    pub out: PathBuf,
    #[arg(long, default_value_t = false)]
    pub validate_only: bool,
    #[arg(long, default_value = "full")]
    pub lr_profile: String,

    #[arg(long, value_enum, default_value_t = AggModeArg::Median)]
    pub agg: AggModeArg,
    #[arg(long, default_value_t = 0.05)]
    pub trim: f32,
    #[arg(long, default_value_t = 1e-6)]
    pub eps: f32,
    #[arg(long, value_enum, default_value_t = CapModeArg::P99)]
    pub cap_mode: CapModeArg,
    #[arg(long, default_value_t = 0.99)]
    pub cap_p: f32,
    #[arg(long)]
    pub cap_fixed: Option<f32>,

    #[arg(long, default_value_t = 0.10)]
    pub cov_min: f32,
    #[arg(long, default_value_t = 0.05)]
    pub expr_min: f32,
    #[arg(long, value_enum, default_value_t = SpecModeArg::On)]
    pub spec: SpecModeArg,
    #[arg(long, default_value_t = 10.0)]
    pub spec_cap: f32,
    #[arg(long, default_value_t = 200)]
    pub top_n_per_pair: usize,
    #[arg(long, default_value_t = 200)]
    pub top_n_per_source: usize,

    #[arg(long, default_value_t = 0.5)]
    pub loud_thresh: f64,
    #[arg(long, default_value_t = 0.5)]
    pub silent_thresh: f64,
    #[arg(long, default_value_t = 3)]
    pub min_support: usize,
    #[arg(long)]
    pub regime_map: Option<PathBuf>,
}

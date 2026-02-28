use clap::Parser;
use kira_environment::cli::{Cli, Command};
use kira_environment::error::KiraError;
use kira_environment::logging;
use kira_environment::simd::dispatch::detect_simd_level;
use kira_environment::{RunConfig, run};

fn main() {
    let cli = Cli::parse();

    let simd = detect_simd_level();
    logging::info(format!("SIMD: {}", simd.as_str()));

    let result = match cli.command {
        Command::Run(args) => run(RunConfig {
            expr: args.expr,
            barcodes: args.barcodes,
            groups: args.groups,
            auto_groups: args.auto_groups,
            auto_groups_coarse: args.auto_groups_coarse,
            auto_groups_fine: args.auto_groups_fine,
            auto_groups_anti: args.auto_groups_anti,
            auto_groups_mode: args.auto_groups_mode,
            auto_groups_eps: args.auto_groups_eps,
            auto_groups_min_delta: args.auto_groups_min_delta,
            auto_groups_unknown: args.auto_groups_unknown,
            auto_groups_emit_scores: args.auto_groups_emit_scores,
            resources: args.resources,
            embedded_profile: args.embedded_profile,
            secretion: args.secretion,
            out: args.out,
            validate_only: args.validate_only,
            lr_profile: args.lr_profile,
            simd_level: simd,
            agg: args.agg,
            trim: args.trim,
            eps: args.eps,
            cap_mode: args.cap_mode,
            cap_p: args.cap_p,
            cap_fixed: args.cap_fixed,
            cov_min: args.cov_min,
            expr_min: args.expr_min,
            spec: args.spec,
            spec_cap: args.spec_cap,
            top_n_per_pair: args.top_n_per_pair,
            top_n_per_source: args.top_n_per_source,
            loud_thresh: args.loud_thresh,
            silent_thresh: args.silent_thresh,
            min_support: args.min_support,
            regime_map: args.regime_map,
        }),
        Command::Validate(args) => run(RunConfig {
            expr: args.expr,
            barcodes: args.barcodes,
            groups: args.groups,
            auto_groups: args.auto_groups,
            auto_groups_coarse: args.auto_groups_coarse,
            auto_groups_fine: args.auto_groups_fine,
            auto_groups_anti: args.auto_groups_anti,
            auto_groups_mode: args.auto_groups_mode,
            auto_groups_eps: args.auto_groups_eps,
            auto_groups_min_delta: args.auto_groups_min_delta,
            auto_groups_unknown: args.auto_groups_unknown,
            auto_groups_emit_scores: args.auto_groups_emit_scores,
            resources: args.resources,
            embedded_profile: args.embedded_profile,
            secretion: args.secretion,
            out: args.out,
            validate_only: true,
            lr_profile: args.lr_profile,
            simd_level: simd,
            agg: args.agg,
            trim: args.trim,
            eps: args.eps,
            cap_mode: args.cap_mode,
            cap_p: args.cap_p,
            cap_fixed: args.cap_fixed,
            cov_min: args.cov_min,
            expr_min: args.expr_min,
            spec: args.spec,
            spec_cap: args.spec_cap,
            top_n_per_pair: args.top_n_per_pair,
            top_n_per_source: args.top_n_per_source,
            loud_thresh: args.loud_thresh,
            silent_thresh: args.silent_thresh,
            min_support: args.min_support,
            regime_map: args.regime_map,
        }),
    };

    if let Err(err) = result {
        eprintln!("{err}");
        std::process::exit(exit_code(&err));
    }
}

fn exit_code(err: &KiraError) -> i32 {
    err.kind.exit_code()
}

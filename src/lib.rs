pub mod agg;
pub mod auto_groups;
pub mod cli;
pub mod determinism;
pub mod error;
pub mod expr;
pub mod io;
pub mod logging;
pub mod paths;
pub mod resources;
pub mod select;
pub mod simd;
pub mod stages;
pub mod version;

use auto_groups::{AutoGroupsConfig, AutoGroupsMode, AutoGroupsResult, run_auto_groups};
use cli::{AggModeArg, AutoGroupsModeArg, CapModeArg, SpecModeArg};
use error::{ErrorKind, KiraError, Result};
use serde_json::{Map, Value};
use simd::dispatch::SimdLevel;
use stages::stage0_resolve::{Stage0Config, Stage0Result, run_stage0};
use stages::stage1_agg::{Stage1Config, Stage1Result, run_stage1};
use stages::stage2_score::{Stage2Config, Stage2Result, run_stage2};
use stages::stage3_network::{Stage3Config, Stage3Result, run_stage3};
use stages::stage4_link::{Stage4Config, Stage4Result, run_stage4};
use std::path::PathBuf;
use version::{TOOL_NAME, TOOL_VERSION};

#[derive(Debug, Clone)]
pub struct RunConfig {
    pub expr: PathBuf,
    pub barcodes: PathBuf,
    pub groups: Option<PathBuf>,
    pub auto_groups: Option<PathBuf>,
    pub auto_groups_coarse: Option<PathBuf>,
    pub auto_groups_fine: Option<PathBuf>,
    pub auto_groups_anti: Option<PathBuf>,
    pub auto_groups_mode: AutoGroupsModeArg,
    pub auto_groups_eps: f32,
    pub auto_groups_min_delta: f32,
    pub auto_groups_unknown: String,
    pub auto_groups_emit_scores: bool,
    pub resources: PathBuf,
    pub secretion: Option<PathBuf>,
    pub out: PathBuf,
    pub validate_only: bool,
    pub lr_profile: String,
    pub simd_level: SimdLevel,
    pub agg: AggModeArg,
    pub trim: f32,
    pub eps: f32,
    pub cap_mode: CapModeArg,
    pub cap_p: f32,
    pub cap_fixed: Option<f32>,
    pub cov_min: f32,
    pub expr_min: f32,
    pub spec: SpecModeArg,
    pub spec_cap: f32,
    pub top_n_per_pair: usize,
    pub top_n_per_source: usize,
    pub loud_thresh: f64,
    pub silent_thresh: f64,
    pub min_support: usize,
    pub regime_map: Option<PathBuf>,
}

#[derive(Debug, Clone)]
pub struct PipelineResult {
    pub auto_groups: Option<AutoGroupsResult>,
    pub stage0: Stage0Result,
    pub stage1: Option<Stage1Result>,
    pub stage2: Option<Stage2Result>,
    pub stage3: Option<Stage3Result>,
    pub stage4: Option<Stage4Result>,
}

pub fn run(config: RunConfig) -> Result<PipelineResult> {
    crate::logging::info("Pipeline start: kira-microenvironment");
    let has_auto_flags = config.auto_groups.is_some()
        || config.auto_groups_coarse.is_some()
        || config.auto_groups_fine.is_some()
        || config.auto_groups_anti.is_some();

    if config.groups.is_some() && has_auto_flags {
        return Err(KiraError::new(
            ErrorKind::AutoGroupsConfigInvalid,
            "--groups is mutually exclusive with auto-groups flags",
        ));
    }
    if config.groups.is_none() && !has_auto_flags {
        return Err(KiraError::new(
            ErrorKind::InvalidArgument,
            "must provide one of --groups or auto-groups flags",
        ));
    }

    if has_auto_flags {
        match config.auto_groups_mode {
            AutoGroupsModeArg::Flat => {
                if config.auto_groups.is_none() {
                    return Err(KiraError::new(
                        ErrorKind::AutoGroupsConfigInvalid,
                        "--auto-groups is required for --auto-groups-mode flat",
                    ));
                }
                if config.auto_groups_coarse.is_some() || config.auto_groups_fine.is_some() {
                    return Err(KiraError::new(
                        ErrorKind::AutoGroupsConfigInvalid,
                        "--auto-groups and --auto-groups-coarse/--auto-groups-fine are mutually exclusive",
                    ));
                }
            }
            AutoGroupsModeArg::Hierarchical => {
                if config.auto_groups.is_some() {
                    return Err(KiraError::new(
                        ErrorKind::AutoGroupsConfigInvalid,
                        "--auto-groups and --auto-groups-coarse are mutually exclusive",
                    ));
                }
                if config.auto_groups_coarse.is_none() || config.auto_groups_fine.is_none() {
                    return Err(KiraError::new(
                        ErrorKind::AutoGroupsConfigInvalid,
                        "--auto-groups-coarse and --auto-groups-fine are required for hierarchical mode",
                    ));
                }
            }
        }
    }

    let groups_path = match (&config.groups, has_auto_flags) {
        (Some(g), false) => g.clone(),
        (Some(_), true) => {
            return Err(KiraError::new(
                ErrorKind::AutoGroupsConfigInvalid,
                "--groups is mutually exclusive with auto-groups flags",
            ));
        }
        (None, false) => {
            return Err(KiraError::new(
                ErrorKind::InvalidArgument,
                "must provide one of --groups or auto-groups flags",
            ));
        }
        (None, true) => PathBuf::new(),
    };

    let auto_groups = if has_auto_flags {
        crate::logging::info("Stage AG start: auto-groups");
        Some(run_auto_groups(&AutoGroupsConfig {
            expr: config.expr.clone(),
            barcodes: config.barcodes.clone(),
            mode: match config.auto_groups_mode {
                AutoGroupsModeArg::Flat => AutoGroupsMode::Flat,
                AutoGroupsModeArg::Hierarchical => AutoGroupsMode::Hierarchical,
            },
            markers_flat: config.auto_groups.clone(),
            markers_coarse: config.auto_groups_coarse.clone(),
            markers_fine: config.auto_groups_fine.clone(),
            anti_markers: config.auto_groups_anti.clone(),
            resources_dir: config.resources.clone(),
            out_dir: config.out.clone(),
            eps: config.auto_groups_eps,
            min_delta: config.auto_groups_min_delta,
            unknown_label: config.auto_groups_unknown.clone(),
            emit_scores: config.auto_groups_emit_scores,
        })?)
    } else {
        crate::logging::info("Stage AG skipped: using provided --groups");
        None
    };
    if auto_groups.is_some() {
        crate::logging::info("Stage AG done");
    }

    let groups_for_stage0 = if let Some(ag) = &auto_groups {
        ag.groups_path.clone()
    } else {
        groups_path
    };

    crate::logging::info("Stage0 start: resolve");
    let stage0 = run_stage0(Stage0Config {
        expr: config.expr.clone(),
        barcodes: config.barcodes,
        groups: groups_for_stage0,
        resources_dir: config.resources.clone(),
        secretion: config.secretion.clone(),
        out_dir: config.out.clone(),
        validate_only: config.validate_only,
        lr_profile: config.lr_profile.clone(),
        simd_level: config.simd_level,
    })?;
    crate::logging::info("Stage0 done");

    if config.validate_only {
        crate::logging::info("Pipeline done: validate-only");
        return Ok(PipelineResult {
            auto_groups,
            stage0,
            stage1: None,
            stage2: None,
            stage3: None,
            stage4: None,
        });
    }

    crate::logging::info("Stage1 start: aggregate");
    let stage1 = run_stage1(Stage1Config {
        expr: config.expr,
        out_dir: config.out.clone(),
        agg_mode: config.agg,
        trim: config.trim,
        eps: config.eps,
        cap_mode: config.cap_mode,
        cap_p: config.cap_p,
        cap_fixed: config.cap_fixed,
    })?;
    crate::logging::info("Stage1 done");

    crate::logging::info("Stage2 start: score");
    let stage2 = run_stage2(Stage2Config {
        out_dir: config.out.clone(),
        resources_dir: config.resources.clone(),
        lr_profile: config.lr_profile,
        eps: config.eps,
        cov_min: config.cov_min,
        expr_min: config.expr_min,
        spec_on: matches!(config.spec, SpecModeArg::On),
        spec_cap: config.spec_cap,
        top_n_per_pair: config.top_n_per_pair,
        top_n_per_source: config.top_n_per_source,
    })?;
    crate::logging::info("Stage2 done");

    crate::logging::info("Stage3 start: network");
    let stage3 = run_stage3(Stage3Config {
        out_dir: config.out.clone(),
    })?;
    crate::logging::info("Stage3 done");

    let stage4 = if let Some(secretion_dir) = config.secretion {
        crate::logging::info("Stage4 start: secretion linking");
        Some(run_stage4(Stage4Config {
            out_dir: config.out,
            resources_dir: config.resources,
            secretion_dir,
            regime_map_path: config.regime_map,
            loud_thresh: config.loud_thresh,
            silent_thresh: config.silent_thresh,
            min_support: config.min_support,
            eps: config.eps as f64,
        })?)
    } else {
        crate::logging::info("Stage4 skipped: no --secretion provided");
        None
    };
    if stage4.is_some() {
        crate::logging::info("Stage4 done");
    }

    let result = PipelineResult {
        auto_groups,
        stage0,
        stage1: Some(stage1),
        stage2: Some(stage2),
        stage3: Some(stage3),
        stage4,
    };

    write_unified_root_summary(&result)?;
    crate::logging::info("Pipeline done");

    Ok(result)
}

fn write_unified_root_summary(result: &PipelineResult) -> Result<()> {
    let root = result.stage0.out_stage_dir.parent().ok_or_else(|| {
        KiraError::new(ErrorKind::Path, "failed to resolve root output directory")
    })?;
    let summary_path = root.join("summary.json");

    let existing = if summary_path.exists() {
        let bytes = std::fs::read(&summary_path).map_err(|e| {
            KiraError::new(
                ErrorKind::Path,
                format!("failed to read {}: {e}", summary_path.display()),
            )
        })?;
        serde_json::from_slice::<Value>(&bytes).unwrap_or(Value::Object(Map::new()))
    } else {
        Value::Object(Map::new())
    };

    let mut out = Map::new();
    out.insert(
        "tool".to_string(),
        serde_json::json!({
            "name": TOOL_NAME,
            "version": TOOL_VERSION,
        }),
    );
    out.insert(
        "simd".to_string(),
        serde_json::to_value(result.stage0.resolved.simd.clone()).unwrap_or(Value::Null),
    );
    out.insert(
        "inputs".to_string(),
        serde_json::to_value(result.stage0.resolved.inputs.clone())
            .unwrap_or(Value::Array(Vec::new())),
    );
    out.insert(
        "resources".to_string(),
        serde_json::to_value(result.stage0.resolved.resources.clone())
            .unwrap_or(Value::Array(Vec::new())),
    );
    out.insert(
        "stage0".to_string(),
        serde_json::to_value(result.stage0.resolved.clone()).unwrap_or(Value::Null),
    );
    out.insert(
        "stage1".to_string(),
        serde_json::to_value(result.stage1.as_ref().map(|s| s.summary.clone()))
            .unwrap_or(Value::Null),
    );
    out.insert(
        "stage2".to_string(),
        serde_json::to_value(result.stage2.as_ref().map(|s| s.summary.clone()))
            .unwrap_or(Value::Null),
    );
    out.insert(
        "stage3".to_string(),
        serde_json::json!({
            "out_dir": result
                .stage3
                .as_ref()
                .map(|s| s.out_root_dir.to_string_lossy().to_string())
        }),
    );

    if let Some(ag) = &result.auto_groups {
        out.insert(
            "auto_groups".to_string(),
            serde_json::json!({
                "enabled": true,
                "mode": ag.mode.as_str(),
                "markers_file": ag.markers_file,
                "coarse_markers": ag.coarse_markers,
                "fine_markers": ag.fine_markers,
                "anti_markers": ag.anti_markers,
                "eps": ag.eps,
                "min_delta": ag.min_delta,
                "unknown_label": ag.unknown_label,
                "counts": ag.summary.counts,
            }),
        );
    }

    if let Some(stage4) = &result.stage4 {
        let link_path = stage4.out_stage_dir.join("link_summary.json");
        if link_path.exists() {
            let bytes = std::fs::read(&link_path).map_err(|e| {
                KiraError::new(
                    ErrorKind::Path,
                    format!("failed to read {}: {e}", link_path.display()),
                )
            })?;
            let link = serde_json::from_slice::<Value>(&bytes).unwrap_or(Value::Null);
            out.insert("stage4".to_string(), link);
        } else {
            out.insert("stage4".to_string(), Value::Null);
        }
    }

    if let Value::Object(obj) = existing {
        if let Some(v) = obj.get("network") {
            out.insert("network".to_string(), v.clone());
        }
        if let Some(v) = obj.get("top_edges_by_score") {
            out.insert("top_edges_by_score".to_string(), v.clone());
        }
        if let Some(v) = obj.get("linking") {
            out.insert("linking".to_string(), v.clone());
        }
    }

    let regimes = derive_pipeline_regimes(&out);
    out.insert("regimes".to_string(), regimes);
    out.insert(
        "distributions".to_string(),
        derive_pipeline_distributions(&out),
    );

    crate::io::atomic::write_json_atomic(&summary_path, &Value::Object(out))?;
    write_pipeline_step(root)
}

fn derive_pipeline_regimes(summary: &Map<String, Value>) -> Value {
    let counts = summary
        .get("linking")
        .and_then(|v| v.get("counts"))
        .and_then(Value::as_object);

    let n_loud = counts
        .and_then(|c| c.get("n_groups_loud"))
        .and_then(Value::as_u64)
        .unwrap_or(0);
    let n_mixed = counts
        .and_then(|c| c.get("n_groups_mixed"))
        .and_then(Value::as_u64)
        .unwrap_or(0);
    let n_silent = counts
        .and_then(|c| c.get("n_groups_silent"))
        .and_then(Value::as_u64)
        .unwrap_or(0);
    let n_total = counts
        .and_then(|c| c.get("n_groups_total"))
        .and_then(Value::as_u64)
        .unwrap_or(0);
    let n_known = n_loud + n_mixed + n_silent;
    let n_unclassified = n_total.saturating_sub(n_known);
    let denom = if n_total > 0 { n_total as f64 } else { 1.0 };

    serde_json::json!({
        "counts": {
            "LoudCommunication": n_loud,
            "MixedCommunication": n_mixed,
            "SilentCommunication": n_silent,
            "Unclassified": n_unclassified,
        },
        "fractions": {
            "LoudCommunication": (n_loud as f64) / denom,
            "MixedCommunication": (n_mixed as f64) / denom,
            "SilentCommunication": (n_silent as f64) / denom,
            "Unclassified": (n_unclassified as f64) / denom,
        }
    })
}

fn derive_pipeline_distributions(summary: &Map<String, Value>) -> Value {
    let n_edges = summary
        .get("network")
        .and_then(|v| v.get("n_edges"))
        .and_then(Value::as_f64)
        .unwrap_or(0.0);
    let n_groups = summary
        .get("network")
        .and_then(|v| v.get("n_groups"))
        .and_then(Value::as_f64)
        .unwrap_or(0.0);
    let pairs_scored = summary
        .get("linking")
        .and_then(|v| v.get("counts"))
        .and_then(|v| v.get("n_pairs_scored"))
        .and_then(Value::as_f64)
        .unwrap_or(0.0);
    let edges_after_filter = summary
        .get("stage2")
        .and_then(|v| v.get("counts"))
        .and_then(|v| v.get("n_edges_after_filter"))
        .and_then(Value::as_f64)
        .unwrap_or(0.0);

    serde_json::json!({
        "network_edges": { "median": n_edges, "p90": Value::Null, "p99": Value::Null },
        "network_groups": { "median": n_groups, "p90": Value::Null, "p99": Value::Null },
        "pairs_scored": { "median": pairs_scored, "p90": Value::Null, "p99": Value::Null },
        "edges_after_filter": { "median": edges_after_filter, "p90": Value::Null, "p99": Value::Null },
    })
}

fn write_pipeline_step(root: &std::path::Path) -> Result<()> {
    let pipeline_step = serde_json::json!({
        "tool": { "name": TOOL_NAME, "version": TOOL_VERSION },
        "mode": "pipeline",
        "artifacts": {
            "summary": "summary.json"
        },
        "timestamp": "1970-01-01T00:00:00Z"
    });
    crate::io::atomic::write_json_atomic(&root.join("pipeline_step.json"), &pipeline_step)
}

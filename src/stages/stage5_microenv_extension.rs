use crate::error::{ErrorKind, KiraError, Result};
use crate::expr::reader::ExprReader;
use crate::io::atomic::{write_json_atomic, write_tsv_atomic};
use crate::io::format::fmt_f32;
use crate::io::input::read_tsv_input;
use crate::io::tsv::TsvReader;
use crate::metrics::microenv_extension::aggregate::{
    ClusterStatsRow, GlobalStats, Missingness, build_cluster_stats, build_global_stats,
    build_missingness,
};
use crate::metrics::microenv_extension::panels::MICROENV_EXTENSION_PANEL_V1;
use crate::metrics::microenv_extension::scores::{
    CellMicroenvScores, MicroenvThresholds, compute_microenv_scores, resolve_panels,
};
use serde::Serialize;
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

#[derive(Debug, Clone)]
pub struct Stage5Config {
    pub expr: PathBuf,
    pub out_dir: PathBuf,
}

#[derive(Debug, Clone)]
pub struct Stage5Result {
    pub out_stage_dir: PathBuf,
    pub summary: Stage5Summary,
}

#[derive(Debug, Clone, Serialize)]
pub struct Stage5Summary {
    pub panel_version: String,
    pub thresholds: MicroenvThresholds,
    pub global_stats: GlobalStats,
    pub cluster_stats: Vec<ClusterStatsRow>,
    pub missingness: Missingness,
}

pub fn run_stage5(cfg: Stage5Config) -> Result<Stage5Result> {
    crate::logging::info("Stage5: computing per-cell microenvironment extension metrics");
    let root = cfg.out_dir.join("kira-microenvironment");
    let stage0_dir = root.join("stage0_resolve");
    let stage5_dir = root.join("stage5_microenv_extension");
    let metrics_path = root.join("metrics.tsv");

    let groups_normalized_path = stage0_dir.join("groups_normalized.tsv");
    let resolved_path = stage0_dir.join("resolved.json");

    for p in [&cfg.expr, &groups_normalized_path, &resolved_path] {
        if !p.exists() {
            return Err(KiraError::new(
                ErrorKind::InputMissing,
                format!("missing required input: {}", p.display()),
            ));
        }
    }

    let expr = ExprReader::open(&cfg.expr)?;
    let barcodes =
        parse_stage0_barcodes_path(&resolved_path).and_then(|p| parse_barcodes(Path::new(&p)))?;
    let group_by_cell = parse_groups_normalized(&groups_normalized_path)?;

    let mut clusters = Vec::with_capacity(barcodes.len());
    for cell_id in &barcodes {
        clusters.push(
            group_by_cell
                .get(cell_id)
                .cloned()
                .unwrap_or_else(|| "unassigned".to_string()),
        );
    }

    if barcodes.len() != expr.n_cells() {
        return Err(KiraError::new(
            ErrorKind::GroupsFormat,
            format!(
                "barcodes count ({}) does not match expr cells ({})",
                barcodes.len(),
                expr.n_cells()
            ),
        ));
    }

    let resolved = resolve_panels(&expr);
    let scores = compute_microenv_scores(&expr, &resolved);
    write_metrics_tsv(&metrics_path, &barcodes, &clusters, &scores)?;

    let cluster_stats = build_cluster_stats(&scores, &clusters);
    let global_stats = build_global_stats(&scores, &cluster_stats);
    let missingness = build_missingness(&scores, resolved.coverage);
    let summary = Stage5Summary {
        panel_version: MICROENV_EXTENSION_PANEL_V1.to_string(),
        thresholds: MicroenvThresholds::default(),
        global_stats,
        cluster_stats,
        missingness,
    };

    write_json_atomic(&stage5_dir.join("stage5_summary.json"), &summary)?;
    crate::logging::info("Stage5: microenvironment metrics written");

    Ok(Stage5Result {
        out_stage_dir: stage5_dir,
        summary,
    })
}

fn write_metrics_tsv(
    path: &Path,
    cell_ids: &[String],
    clusters: &[String],
    scores: &CellMicroenvScores,
) -> Result<()> {
    let mut lines = Vec::with_capacity(cell_ids.len() + 1);
    lines.push(
        "cell_id\tgroup\thyp_core\tnfkb_core\tifn_core\tcheckpoint_core\tadenosine_core\tstromal_core\tHSI\tIAS\tISS\tMIO\tSII\tMSM\thypoxia_high\tinflammatory_high\timmune_suppression_high\tmetabolic_suppression_high\tstromal_high\tmicroenv_stress_mode"
            .to_string(),
    );

    for i in 0..cell_ids.len() {
        lines.push(format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            cell_ids[i],
            clusters[i],
            fmt_f32(scores.hyp_core[i]),
            fmt_f32(scores.nfkb_core[i]),
            fmt_f32(scores.ifn_core[i]),
            fmt_f32(scores.checkpoint_core[i]),
            fmt_f32(scores.adenosine_core[i]),
            fmt_f32(scores.stromal_core[i]),
            fmt_f32(scores.hsi[i]),
            fmt_f32(scores.ias[i]),
            fmt_f32(scores.iss[i]),
            fmt_f32(scores.mio[i]),
            fmt_f32(scores.sii[i]),
            fmt_f32(scores.msm[i]),
            scores.hypoxia_high[i],
            scores.inflammatory_high[i],
            scores.immune_suppression_high[i],
            scores.metabolic_suppression_high[i],
            scores.stromal_high[i],
            scores.microenv_stress_mode[i],
        ));
    }
    write_tsv_atomic(path, &lines)
}

fn parse_stage0_barcodes_path(path: &Path) -> Result<String> {
    let bytes = std::fs::read(path).map_err(|e| {
        KiraError::new(
            ErrorKind::Path,
            format!("failed to read {}: {e}", path.display()),
        )
    })?;
    let v: serde_json::Value = serde_json::from_slice(&bytes).map_err(|e| {
        KiraError::new(
            ErrorKind::Path,
            format!("failed to parse {}: {e}", path.display()),
        )
    })?;
    let arr = v
        .get("inputs")
        .and_then(|x| x.as_array())
        .ok_or_else(|| KiraError::new(ErrorKind::Path, "stage0 resolved.json missing inputs"))?;
    for item in arr {
        let kind = item
            .get("kind")
            .and_then(|x| x.as_str())
            .unwrap_or_default();
        if kind == "barcodes" {
            let p = item
                .get("path_abs")
                .and_then(|x| x.as_str())
                .ok_or_else(|| {
                    KiraError::new(ErrorKind::Path, "barcodes input missing path_abs")
                })?;
            return Ok(p.to_string());
        }
    }
    Err(KiraError::new(
        ErrorKind::Path,
        "stage0 resolved.json missing barcodes input",
    ))
}

fn parse_barcodes(path: &Path) -> Result<Vec<String>> {
    let input = read_tsv_input(path)?;
    let mut rdr = TsvReader::new(input.as_bytes());
    let mut row = Vec::new();
    let mut out = Vec::new();
    let mut first = true;
    let mut cell_col = 0usize;

    while let Some(rec) = rdr.next_record(&mut row) {
        rec?;
        if row.is_empty() {
            continue;
        }
        if first {
            first = false;
            if row.len() == 1 {
                if row[0] == "cell_id" {
                    continue;
                }
                out.push(row[0].clone());
                continue;
            }
            if let Some(idx) = row.iter().position(|x| x == "cell_id") {
                cell_col = idx;
                continue;
            }
            return Err(KiraError::new(
                ErrorKind::TsvHeader,
                format!("barcodes header missing cell_id in {}", path.display()),
            ));
        }
        if row.len() <= cell_col {
            return Err(KiraError::new(
                ErrorKind::TsvParse,
                format!("invalid barcodes row in {}", path.display()),
            ));
        }
        out.push(row[cell_col].clone());
    }
    Ok(out)
}

fn parse_groups_normalized(path: &Path) -> Result<BTreeMap<String, String>> {
    let input = read_tsv_input(path)?;
    let mut rdr = TsvReader::new(input.as_bytes());
    let mut row = Vec::new();
    let mut out = BTreeMap::new();
    let mut first = true;
    let mut cell_col = 0usize;
    let mut group_col = 1usize;

    while let Some(rec) = rdr.next_record(&mut row) {
        rec?;
        if row.is_empty() {
            continue;
        }
        if first {
            first = false;
            if row == ["cell_id", "group"] {
                continue;
            }
            cell_col = row.iter().position(|x| x == "cell_id").ok_or_else(|| {
                KiraError::new(
                    ErrorKind::TsvHeader,
                    format!("groups header missing cell_id in {}", path.display()),
                )
            })?;
            group_col = row.iter().position(|x| x == "group").ok_or_else(|| {
                KiraError::new(
                    ErrorKind::TsvHeader,
                    format!("groups header missing group in {}", path.display()),
                )
            })?;
            continue;
        }

        if row.len() <= cell_col || row.len() <= group_col {
            return Err(KiraError::new(
                ErrorKind::TsvParse,
                format!("invalid groups row in {}", path.display()),
            ));
        }
        out.entry(row[cell_col].clone())
            .or_insert(row[group_col].clone());
    }
    Ok(out)
}

pub mod anti;
pub mod assign;
pub mod coarse;
pub mod fine;
pub mod hierarchy;
pub mod markers;
pub mod output;
pub mod summary;

use crate::auto_groups::anti::load_anti_markers;
use crate::auto_groups::assign::{AssignConfig, assign_cells};
use crate::auto_groups::coarse::load_coarse_groups;
use crate::auto_groups::fine::load_fine_groups;
use crate::auto_groups::markers::load_marker_groups;
use crate::auto_groups::output::{write_groups, write_scores, write_summary};
use crate::auto_groups::summary::AutoGroupsSummary;
use crate::error::{ErrorKind, KiraError, Result};
use crate::expr::reader::ExprReader;
use crate::io::input::read_tsv_input;
use crate::io::tsv::TsvReader;
use crate::resources::aliases::{AliasResolver, parse_aliases};
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AutoGroupsMode {
    Flat,
    Hierarchical,
}

impl AutoGroupsMode {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Flat => "flat",
            Self::Hierarchical => "hierarchical",
        }
    }
}

#[derive(Debug, Clone)]
pub struct AutoGroupsConfig {
    pub expr: PathBuf,
    pub barcodes: PathBuf,
    pub mode: AutoGroupsMode,
    pub markers_flat: Option<PathBuf>,
    pub markers_coarse: Option<PathBuf>,
    pub markers_fine: Option<PathBuf>,
    pub anti_markers: Option<PathBuf>,
    pub resources_dir: PathBuf,
    pub out_dir: PathBuf,
    pub eps: f32,
    pub min_delta: f32,
    pub unknown_label: String,
    pub emit_scores: bool,
}

#[derive(Debug, Clone)]
pub struct AutoGroupsResult {
    pub groups_path: PathBuf,
    pub summary: AutoGroupsSummary,
    pub mode: AutoGroupsMode,
    pub markers_file: Option<String>,
    pub coarse_markers: Option<String>,
    pub fine_markers: Option<String>,
    pub anti_markers: Option<String>,
    pub eps: f32,
    pub min_delta: f32,
    pub unknown_label: String,
}

pub fn run_auto_groups(cfg: &AutoGroupsConfig) -> Result<AutoGroupsResult> {
    crate::logging::info(format!("Auto-groups: mode={}", cfg.mode.as_str()));
    let expr = ExprReader::open(&cfg.expr)?;
    let barcodes = parse_barcodes(&cfg.barcodes)?;
    crate::logging::info(format!(
        "Auto-groups: loaded expr cells={} genes={} barcodes={}",
        expr.n_cells(),
        expr.n_genes(),
        barcodes.len()
    ));

    let alias_path = cfg.resources_dir.join("gene_alias.tsv");
    let aliases = if alias_path.exists() {
        AliasResolver::from_entries(parse_aliases(&alias_path)?)
    } else {
        AliasResolver::default()
    };

    let anti_markers = if let Some(path) = &cfg.anti_markers {
        if !path.exists() {
            return Err(KiraError::new(
                ErrorKind::AutoGroupsMarkersMissing,
                format!("anti-markers missing: {}", path.display()),
            ));
        }
        load_anti_markers(path, &expr, &aliases)?
    } else {
        BTreeMap::new()
    };
    if !anti_markers.is_empty() {
        crate::logging::info(format!(
            "Auto-groups: anti-marker groups={}",
            anti_markers.len()
        ));
    }

    let (fine_groups, coarse_groups) = match cfg.mode {
        AutoGroupsMode::Flat => {
            let markers = cfg.markers_flat.as_ref().ok_or_else(|| {
                KiraError::new(
                    ErrorKind::AutoGroupsConfigInvalid,
                    "missing --auto-groups for flat mode",
                )
            })?;
            if !markers.exists() {
                return Err(KiraError::new(
                    ErrorKind::AutoGroupsMarkersMissing,
                    format!("marker panels missing: {}", markers.display()),
                ));
            }
            (load_marker_groups(markers, &expr, &aliases)?, None)
        }
        AutoGroupsMode::Hierarchical => {
            let coarse_path = cfg.markers_coarse.as_ref().ok_or_else(|| {
                KiraError::new(
                    ErrorKind::AutoGroupsConfigInvalid,
                    "missing --auto-groups-coarse for hierarchical mode",
                )
            })?;
            let fine_path = cfg.markers_fine.as_ref().ok_or_else(|| {
                KiraError::new(
                    ErrorKind::AutoGroupsConfigInvalid,
                    "missing --auto-groups-fine for hierarchical mode",
                )
            })?;
            if !coarse_path.exists() {
                return Err(KiraError::new(
                    ErrorKind::AutoGroupsMarkersMissing,
                    format!("coarse marker panels missing: {}", coarse_path.display()),
                ));
            }
            if !fine_path.exists() {
                return Err(KiraError::new(
                    ErrorKind::AutoGroupsMarkersMissing,
                    format!("fine marker panels missing: {}", fine_path.display()),
                ));
            }
            (
                load_fine_groups(fine_path, &expr, &aliases)?,
                Some(load_coarse_groups(coarse_path, &expr, &aliases)?),
            )
        }
    };

    if fine_groups.is_empty() {
        return Err(KiraError::new(
            ErrorKind::AutoGroupsParseError,
            "no groups parsed from marker panels",
        ));
    }
    if let Some(coarse) = &coarse_groups
        && coarse.is_empty()
    {
        return Err(KiraError::new(
            ErrorKind::AutoGroupsParseError,
            "no groups parsed from coarse marker panels",
        ));
    }

    let assign = assign_cells(
        &expr,
        &barcodes,
        &fine_groups,
        coarse_groups.as_deref(),
        &anti_markers,
        &AssignConfig {
            mode: cfg.mode,
            eps: cfg.eps,
            min_delta: cfg.min_delta,
            unknown_label: cfg.unknown_label.clone(),
        },
    )?;
    crate::logging::info(format!(
        "Auto-groups: assigned cells={}, unknown={}, coarse_only={}, fine={}",
        assign.rows.len(),
        assign.counts.n_unknown,
        assign.counts.n_coarse_only,
        assign.counts.n_fine
    ));

    if assign.rows.is_empty() {
        return Err(KiraError::new(
            ErrorKind::AutoGroupsEmptyAssignment,
            "auto-groups produced no rows",
        ));
    }

    let summary = AutoGroupsSummary {
        mode: cfg.mode.as_str().to_string(),
        counts: crate::auto_groups::summary::Counts {
            n_cells: barcodes.len(),
            n_unknown: assign.counts.n_unknown,
            n_coarse_only: assign.counts.n_coarse_only,
            n_fine: assign.counts.n_fine,
        },
        per_group: assign.per_group.clone(),
        anti_markers_used: !anti_markers.is_empty(),
    };

    let out_dir = cfg
        .out_dir
        .join("kira-microenvironment")
        .join("auto_groups");
    let groups_path = out_dir.join("groups.tsv");
    write_groups(&groups_path, &assign.rows)?;
    write_summary(&out_dir.join("summary.json"), &summary)?;
    if cfg.emit_scores {
        write_scores(&out_dir.join("cell_scores.tsv"), &assign.cell_scores)?;
    }
    crate::logging::info("Auto-groups: outputs written");

    Ok(AutoGroupsResult {
        groups_path,
        summary,
        mode: cfg.mode,
        markers_file: cfg
            .markers_flat
            .as_ref()
            .map(|p| p.to_string_lossy().into_owned()),
        coarse_markers: cfg
            .markers_coarse
            .as_ref()
            .map(|p| p.to_string_lossy().into_owned()),
        fine_markers: cfg
            .markers_fine
            .as_ref()
            .map(|p| p.to_string_lossy().into_owned()),
        anti_markers: cfg
            .anti_markers
            .as_ref()
            .map(|p| p.to_string_lossy().into_owned()),
        eps: cfg.eps,
        min_delta: cfg.min_delta,
        unknown_label: cfg.unknown_label.clone(),
    })
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

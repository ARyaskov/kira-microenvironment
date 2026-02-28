use crate::determinism::{cmp_opt_str, stable_sort_by};
use crate::error::{ErrorKind, KiraError, Result};
use crate::io::atomic::{write_json_atomic, write_tsv_atomic};
use crate::io::input::read_tsv_input;
use crate::io::tsv::TsvReader;
use crate::paths::{absolutize, metadata_for_path, must_exist_file, must_exist_path};
use crate::resources::{ResourcePaths, load_resources, resolve_lr_path};
use crate::simd::dispatch::SimdLevel;
use crate::version::{EDITION, RUST_FLOOR, TOOL_NAME, TOOL_VERSION};
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

#[derive(Debug, Clone)]
pub struct Stage0Config {
    pub expr: PathBuf,
    pub barcodes: PathBuf,
    pub groups: PathBuf,
    pub resources_dir: PathBuf,
    pub secretion: Option<PathBuf>,
    pub out_dir: PathBuf,
    pub validate_only: bool,
    pub lr_profile: String,
    pub simd_level: SimdLevel,
}

#[derive(Debug, Clone, Serialize)]
pub struct Stage0Counts {
    pub n_cells_barcodes: usize,
    pub n_rows_groups: usize,
    pub n_groups: usize,
    pub n_lr_pairs_raw: usize,
}

#[derive(Debug, Clone, Serialize)]
pub struct Stage0Resolved {
    pub tool: ToolInfo,
    pub simd: SimdInfo,
    pub inputs: Vec<crate::paths::PathMeta>,
    pub resources: Vec<crate::paths::PathMeta>,
    pub config: ConfigInfo,
    pub counts: Stage0Counts,
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize)]
pub struct ToolInfo {
    pub name: String,
    pub version: String,
    pub rust: String,
    pub edition: String,
}

#[derive(Debug, Clone, Serialize)]
pub struct SimdInfo {
    pub level: String,
}

#[derive(Debug, Clone, Serialize)]
pub struct ConfigInfo {
    pub lr_profile: String,
    pub validate_only: bool,
}

#[derive(Debug, Clone)]
pub struct Stage0Result {
    pub out_stage_dir: PathBuf,
    pub resolved: Stage0Resolved,
}

#[derive(Debug, Clone)]
struct GroupRow {
    cell_id: String,
    group: String,
    row_idx: usize,
}

#[derive(Debug, Clone)]
struct GeneRequest {
    role: String,
    query_symbol: String,
    complex_id: Option<String>,
    source: String,
}

pub fn run_stage0(cfg: Stage0Config) -> Result<Stage0Result> {
    crate::logging::info("Stage0: validating input/resource paths");
    let expr = must_exist_file(&cfg.expr, "expr")?;
    let barcodes = must_exist_file(&cfg.barcodes, "barcodes")?;
    let groups = must_exist_file(&cfg.groups, "groups")?;
    let resources_dir = must_exist_path(&cfg.resources_dir, "resources")?;
    let secretion = match &cfg.secretion {
        Some(p) => Some(must_exist_path(p, "secretion")?),
        None => None,
    };

    let lr_path = absolutize(&resolve_lr_path(&resources_dir, &cfg.lr_profile)?)?;
    let resource_paths = ResourcePaths {
        lr_pairs: lr_path.clone(),
        cofactors: resources_dir.join("cofactors.tsv"),
        aliases: resources_dir.join("gene_alias.tsv"),
        labels: resources_dir.join("labels.tsv"),
    };
    let resources = load_resources(&resource_paths)?;
    crate::logging::info("Stage0: loading barcodes and groups");

    let barcodes_vec = parse_barcodes(&barcodes)?;
    let barcodes_set: BTreeSet<String> = barcodes_vec.iter().cloned().collect();
    let (groups_rows, n_rows_groups, mut warnings) =
        parse_and_normalize_groups(&groups, &barcodes_set)?;

    let group_sizes = compute_group_sizes(&groups_rows);
    let n_groups = group_sizes.len();
    let missing_in_groups = barcodes_set
        .iter()
        .filter(|cell_id| !groups_rows.iter().any(|r| &r.cell_id == *cell_id))
        .count();
    if missing_in_groups > 0 {
        warnings.push(format!(
            "{missing_in_groups} barcodes missing group assignment"
        ));
    }

    let genes_requested = build_genes_requested(
        &resources.lr_pairs,
        &resources.cofactors,
        &resources.aliases,
    );

    let stage_dir = cfg
        .out_dir
        .join("kira-microenvironment")
        .join("stage0_resolve");
    crate::logging::info("Stage0: writing resolved outputs");

    let groups_lines = render_groups_normalized(&groups_rows);
    let group_sizes_lines = render_group_sizes(&group_sizes);
    let genes_lines = render_genes_requested(&genes_requested);

    write_tsv_atomic(&stage_dir.join("groups_normalized.tsv"), &groups_lines)?;
    write_tsv_atomic(&stage_dir.join("group_sizes.tsv"), &group_sizes_lines)?;
    write_tsv_atomic(&stage_dir.join("genes_requested.tsv"), &genes_lines)?;

    let mut inputs = vec![
        metadata_for_path("expr", &expr)?,
        metadata_for_path("barcodes", &barcodes)?,
        metadata_for_path("groups", &groups)?,
    ];
    if let Some(sec) = &secretion {
        inputs.push(metadata_for_path("secretion", sec)?);
    }

    let resolved = Stage0Resolved {
        tool: ToolInfo {
            name: TOOL_NAME.to_string(),
            version: TOOL_VERSION.to_string(),
            rust: RUST_FLOOR.to_string(),
            edition: EDITION.to_string(),
        },
        simd: SimdInfo {
            level: cfg.simd_level.as_str().to_string(),
        },
        inputs,
        resources: resources.metadata,
        config: ConfigInfo {
            lr_profile: lr_path.to_string_lossy().into_owned(),
            validate_only: cfg.validate_only,
        },
        counts: Stage0Counts {
            n_cells_barcodes: barcodes_vec.len(),
            n_rows_groups,
            n_groups,
            n_lr_pairs_raw: resources.lr_pairs.len(),
        },
        warnings,
    };

    write_json_atomic(&stage_dir.join("resolved.json"), &resolved)?;
    crate::logging::info(format!(
        "Stage0: counts cells={} groups={} lr_pairs={}",
        resolved.counts.n_cells_barcodes, resolved.counts.n_groups, resolved.counts.n_lr_pairs_raw
    ));

    Ok(Stage0Result {
        out_stage_dir: stage_dir,
        resolved,
    })
}

fn parse_barcodes(path: &Path) -> Result<Vec<String>> {
    let input = read_tsv_input(path)?;
    let mut rdr = TsvReader::new(input.as_bytes());
    let mut row = Vec::new();
    let mut out = Vec::new();
    let mut first_row_seen = false;
    let mut cell_col: usize = 0;

    while let Some(rec) = rdr.next_record(&mut row) {
        rec?;
        if row.is_empty() {
            continue;
        }
        if !first_row_seen {
            first_row_seen = true;
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
                format!("barcodes header must include cell_id in {}", path.display()),
            ));
        }

        if row.len() <= cell_col {
            return Err(KiraError::new(
                ErrorKind::TsvParse,
                format!("barcodes line {} missing cell_id column", rdr.line_no()),
            ));
        }
        out.push(row[cell_col].clone());
    }

    Ok(out)
}

fn parse_and_normalize_groups(
    path: &Path,
    barcodes: &BTreeSet<String>,
) -> Result<(Vec<GroupRow>, usize, Vec<String>)> {
    let input = read_tsv_input(path)?;
    let mut rdr = TsvReader::new(input.as_bytes());
    let mut row = Vec::new();
    let mut seen_header = false;
    let mut cell_col = None;
    let mut group_col = None;
    let mut rows = Vec::new();
    let mut warnings = Vec::new();
    let mut by_cell_first: BTreeMap<String, GroupRow> = BTreeMap::new();
    let mut n_rows = 0usize;

    while let Some(rec) = rdr.next_record(&mut row) {
        rec?;
        if row.is_empty() {
            continue;
        }
        if !seen_header {
            seen_header = true;
            cell_col = row.iter().position(|x| x == "cell_id");
            if cell_col.is_none() {
                return Err(KiraError::new(
                    ErrorKind::GroupsFormat,
                    format!("groups header missing cell_id in {}", path.display()),
                ));
            }
            group_col = row
                .iter()
                .position(|x| x == "group")
                .or_else(|| row.iter().position(|x| x == "cluster_id"))
                .or_else(|| row.iter().position(|x| x == "cell_type"));
            if group_col.is_none() {
                return Err(KiraError::new(
                    ErrorKind::GroupsFormat,
                    format!(
                        "groups header missing one of group/cluster_id/cell_type in {}",
                        path.display()
                    ),
                ));
            }
            continue;
        }
        n_rows += 1;

        let cidx = cell_col.unwrap_or(0);
        let gidx = group_col.unwrap_or(1);
        if row.len() <= cidx || row.len() <= gidx {
            return Err(KiraError::new(
                ErrorKind::TsvParse,
                format!("groups line {} has too few columns", rdr.line_no()),
            ));
        }

        let cell_id = row[cidx].clone();
        let group = row[gidx].clone();
        if by_cell_first.contains_key(&cell_id) {
            warnings.push(format!(
                "duplicate cell_id in groups: {cell_id} (kept first occurrence)"
            ));
            continue;
        }

        let rowobj = GroupRow {
            cell_id: cell_id.clone(),
            group,
            row_idx: n_rows,
        };
        by_cell_first.insert(cell_id.clone(), rowobj.clone());
        rows.push(rowobj);
    }

    let not_in_barcodes = rows
        .iter()
        .filter(|x| !barcodes.contains(&x.cell_id))
        .count();
    if not_in_barcodes > 0 {
        warnings.push(format!(
            "{not_in_barcodes} group rows reference unknown barcodes"
        ));
    }

    stable_sort_by(&mut rows, |a, b| {
        a.group
            .cmp(&b.group)
            .then(a.cell_id.cmp(&b.cell_id))
            .then(a.row_idx.cmp(&b.row_idx))
    });

    Ok((rows, n_rows, warnings))
}

fn compute_group_sizes(rows: &[GroupRow]) -> Vec<(String, usize)> {
    let mut map: BTreeMap<String, usize> = BTreeMap::new();
    for row in rows {
        *map.entry(row.group.clone()).or_insert(0) += 1;
    }
    map.into_iter().collect()
}

fn build_genes_requested(
    lr_pairs: &[crate::resources::lr_pairs::LrPair],
    cofactors: &[crate::resources::cofactors::Cofactor],
    aliases: &crate::resources::aliases::AliasResolver,
) -> Vec<GeneRequest> {
    let mut dedup: BTreeMap<(String, String, Option<String>), GeneRequest> = BTreeMap::new();

    for lr in lr_pairs {
        let ligand = aliases.resolve(&lr.ligand_symbol).to_string();
        let receptor = aliases.resolve(&lr.receptor_symbol).to_string();

        let lig = GeneRequest {
            role: "ligand".to_string(),
            query_symbol: ligand,
            complex_id: None,
            source: "lr_pairs".to_string(),
        };
        let rec = GeneRequest {
            role: "receptor".to_string(),
            query_symbol: receptor,
            complex_id: None,
            source: "lr_pairs".to_string(),
        };

        for item in [lig, rec] {
            let key = (
                item.role.clone(),
                item.query_symbol.clone(),
                item.complex_id.clone(),
            );
            dedup
                .entry(key)
                .and_modify(|existing| {
                    if item.source < existing.source {
                        *existing = item.clone();
                    }
                })
                .or_insert(item);
        }
    }

    for cf in cofactors {
        let sym = aliases.resolve(&cf.subunit_symbol).to_string();
        let item = GeneRequest {
            role: "subunit".to_string(),
            query_symbol: sym,
            complex_id: Some(cf.complex_id.clone()),
            source: "cofactors".to_string(),
        };
        let key = (
            item.role.clone(),
            item.query_symbol.clone(),
            item.complex_id.clone(),
        );
        dedup
            .entry(key)
            .and_modify(|existing| {
                if item.source < existing.source {
                    *existing = item.clone();
                }
            })
            .or_insert(item);
    }

    let mut out: Vec<GeneRequest> = dedup.into_values().collect();
    stable_sort_by(&mut out, |a, b| {
        a.role
            .cmp(&b.role)
            .then(a.query_symbol.cmp(&b.query_symbol))
            .then(cmp_opt_str(&a.complex_id, &b.complex_id))
            .then(a.source.cmp(&b.source))
    });
    out
}

fn render_groups_normalized(rows: &[GroupRow]) -> Vec<String> {
    let mut lines = Vec::with_capacity(rows.len() + 1);
    lines.push("cell_id\tgroup".to_string());
    for row in rows {
        lines.push(format!("{}\t{}", row.cell_id, row.group));
    }
    lines
}

fn render_group_sizes(group_sizes: &[(String, usize)]) -> Vec<String> {
    let mut lines = Vec::with_capacity(group_sizes.len() + 1);
    lines.push("group\tn_cells".to_string());
    for (group, n) in group_sizes {
        lines.push(format!("{}\t{}", group, n));
    }
    lines
}

fn render_genes_requested(genes: &[GeneRequest]) -> Vec<String> {
    let mut lines = Vec::with_capacity(genes.len() + 1);
    lines.push("role\tquery_symbol\tcomplex_id\tsource".to_string());
    for g in genes {
        let complex = g.complex_id.clone().unwrap_or_default();
        lines.push(format!(
            "{}\t{}\t{}\t{}",
            g.role, g.query_symbol, complex, g.source
        ));
    }
    lines
}

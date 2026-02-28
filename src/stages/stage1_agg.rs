use crate::cli::{AggModeArg, CapModeArg};
use crate::error::{ErrorKind, KiraError, Result};
use crate::expr::reader::ExprReader;
use crate::io::atomic::{write_bytes_atomic, write_json_atomic, write_tsv_atomic};
use crate::io::format::fmt_f32;
use crate::io::input::read_tsv_input;
use crate::io::mmap::MmapFile;
use crate::io::tsv::TsvReader;
use crate::select::{median_in_place, quantile_in_place, trimmed_mean_in_place};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

const CACHE_MAGIC: [u8; 8] = *b"KIRAEAGG";
const CACHE_VERSION: u32 = 1;

#[derive(Debug, Clone)]
pub struct Stage1Config {
    pub expr: PathBuf,
    pub out_dir: PathBuf,
    pub agg_mode: AggModeArg,
    pub trim: f32,
    pub eps: f32,
    pub cap_mode: CapModeArg,
    pub cap_p: f32,
    pub cap_fixed: Option<f32>,
}

#[derive(Debug, Clone, Serialize)]
pub struct Stage1Summary {
    pub counts: SummaryCounts,
    pub agg: SummaryAgg,
    pub eps: f32,
    pub cap: SummaryCap,
    pub missing_genes: usize,
}

#[derive(Debug, Clone, Serialize)]
pub struct SummaryCounts {
    pub n_groups: usize,
    pub n_genes_requested: usize,
    pub n_genes_found: usize,
}

#[derive(Debug, Clone, Serialize)]
pub struct SummaryAgg {
    pub mode: String,
    pub trim: f32,
}

#[derive(Debug, Clone, Serialize)]
pub struct SummaryCap {
    pub mode: String,
    pub value: f32,
}

#[derive(Debug, Clone)]
pub struct Stage1Result {
    pub out_stage_dir: PathBuf,
    pub summary: Stage1Summary,
}

#[derive(Debug, Clone)]
struct GroupRow {
    cell_id: String,
    group: String,
}

#[derive(Debug, Clone)]
struct GeneSpec {
    symbol: String,
    found_expr_idx: Option<u32>,
}

#[derive(Debug, Clone, Deserialize)]
struct Stage0ResolvedInput {
    inputs: Vec<Stage0InputPath>,
}

#[derive(Debug, Clone, Deserialize)]
struct Stage0InputPath {
    kind: String,
    path_abs: String,
}

pub fn run_stage1(cfg: Stage1Config) -> Result<Stage1Result> {
    crate::logging::info("Stage1: reading stage0 outputs and expression cache");
    if matches!(cfg.cap_mode, CapModeArg::Fixed) && cfg.cap_fixed.is_none() {
        return Err(KiraError::new(
            ErrorKind::TsvParse,
            "cap-fixed is required when --cap-mode=fixed",
        ));
    }

    let stage0_dir = cfg
        .out_dir
        .join("kira-microenvironment")
        .join("stage0_resolve");
    let stage1_dir = cfg.out_dir.join("kira-microenvironment").join("stage1_agg");

    let genes_requested_path = stage0_dir.join("genes_requested.tsv");
    let groups_normalized_path = stage0_dir.join("groups_normalized.tsv");
    let group_sizes_path = stage0_dir.join("group_sizes.tsv");
    let resolved_json_path = stage0_dir.join("resolved.json");

    for p in [
        &genes_requested_path,
        &groups_normalized_path,
        &group_sizes_path,
        &resolved_json_path,
    ] {
        if !p.exists() {
            return Err(KiraError::new(
                ErrorKind::InputMissing,
                format!("missing stage0 output: {}", p.display()),
            ));
        }
    }

    let reader = ExprReader::open(&cfg.expr)?;
    let barcodes = parse_stage0_barcodes_path(&resolved_json_path)?;
    let cell_to_index = load_barcodes_index(Path::new(&barcodes), reader.n_cells())?;

    let groups_rows = parse_groups_normalized(&groups_normalized_path)?;
    let group_sizes = parse_group_sizes(&group_sizes_path)?;
    let mut groups_sorted: Vec<String> = group_sizes.keys().cloned().collect();
    groups_sorted.sort();

    let genes_requested = parse_genes_requested(&genes_requested_path)?;
    let mut gene_symbols: Vec<String> = genes_requested.into_iter().collect();
    gene_symbols.sort();

    let mut genes = Vec::with_capacity(gene_symbols.len());
    let mut n_genes_found = 0usize;
    for symbol in &gene_symbols {
        let found = reader.gene_index(symbol);
        if found.is_some() {
            n_genes_found += 1;
        }
        genes.push(GeneSpec {
            symbol: symbol.clone(),
            found_expr_idx: found,
        });
    }

    let missing_genes = genes.len().saturating_sub(n_genes_found);

    let needed_expr_indices: Vec<u32> = genes.iter().filter_map(|g| g.found_expr_idx).collect();
    let max_gene_idx = needed_expr_indices.iter().copied().max().unwrap_or(0) as usize;
    let mut needed_lookup = vec![i32::MIN; max_gene_idx.saturating_add(1)];
    for (local_idx, g) in genes.iter().enumerate() {
        if let Some(idx) = g.found_expr_idx {
            needed_lookup[idx as usize] = local_idx as i32;
        }
    }

    let mut group_to_cells: BTreeMap<String, Vec<u32>> = BTreeMap::new();
    for row in &groups_rows {
        if let Some(idx) = cell_to_index.get(&row.cell_id) {
            group_to_cells
                .entry(row.group.clone())
                .or_default()
                .push(*idx as u32);
        }
    }
    for cells in group_to_cells.values_mut() {
        cells.sort();
    }

    let n_groups = groups_sorted.len();
    let n_genes = genes.len();
    let matrix_len = n_groups.saturating_mul(n_genes);
    let mut expr_matrix = vec![0.0f32; matrix_len];
    let mut cov_matrix = vec![0.0f32; matrix_len];

    let mut buffers: Vec<Vec<f32>> = (0..n_genes).map(|_| Vec::new()).collect();
    let mut positive_counts = vec![0usize; n_genes];
    crate::logging::info(format!(
        "Stage1: aggregating {} groups x {} genes",
        n_groups, n_genes
    ));

    for (group_idx, group_name) in groups_sorted.iter().enumerate() {
        let empty: Vec<u32> = Vec::new();
        let cells = group_to_cells.get(group_name).unwrap_or(&empty);
        let n_cells = *group_sizes.get(group_name).unwrap_or(&cells.len());

        for values in &mut buffers {
            values.clear();
        }
        positive_counts.fill(0);

        for &cell_idx in cells {
            for (gene_idx, value) in reader.iter_cell(cell_idx) {
                let gi = gene_idx as usize;
                if gi >= needed_lookup.len() {
                    continue;
                }
                let local = needed_lookup[gi];
                if local < 0 {
                    continue;
                }
                let local_idx = local as usize;
                buffers[local_idx].push(value);
                if value > cfg.eps {
                    positive_counts[local_idx] += 1;
                }
            }
        }

        for local_idx in 0..n_genes {
            let idx = group_idx * n_genes + local_idx;
            let gene = &genes[local_idx];
            if gene.found_expr_idx.is_none() || n_cells == 0 {
                expr_matrix[idx] = 0.0;
                cov_matrix[idx] = 0.0;
                continue;
            }

            if buffers[local_idx].len() < n_cells {
                buffers[local_idx].resize(n_cells, 0.0);
            } else if buffers[local_idx].len() > n_cells {
                buffers[local_idx].truncate(n_cells);
            }

            let expr = match cfg.agg_mode {
                AggModeArg::Median => median_in_place(&mut buffers[local_idx]),
                AggModeArg::TrimmedMean => trimmed_mean_in_place(&mut buffers[local_idx], cfg.trim),
            };
            let cov = (positive_counts[local_idx] as f32) / (n_cells as f32);
            expr_matrix[idx] = expr;
            cov_matrix[idx] = cov;
        }
    }

    let cap_value = match cfg.cap_mode {
        CapModeArg::Fixed => cfg.cap_fixed.unwrap_or(0.0),
        CapModeArg::P99 => {
            let mut values = expr_matrix.clone();
            quantile_in_place(&mut values, cfg.cap_p)
        }
    };

    write_group_gene_agg(
        &stage1_dir.join("group_gene_agg.tsv"),
        &groups_sorted,
        &genes,
        &group_sizes,
        &expr_matrix,
        &cov_matrix,
    )?;

    write_gene_stats(
        &stage1_dir.join("gene_stats.tsv"),
        &genes,
        n_groups,
        &expr_matrix,
        cap_value,
    )?;

    write_cache_bin(
        &stage1_dir.join("cache.bin"),
        &groups_sorted,
        &genes,
        &expr_matrix,
        &cov_matrix,
    )?;

    let summary = Stage1Summary {
        counts: SummaryCounts {
            n_groups,
            n_genes_requested: n_genes,
            n_genes_found,
        },
        agg: SummaryAgg {
            mode: match cfg.agg_mode {
                AggModeArg::Median => "median",
                AggModeArg::TrimmedMean => "trimmed_mean",
            }
            .to_string(),
            trim: cfg.trim,
        },
        eps: cfg.eps,
        cap: SummaryCap {
            mode: match cfg.cap_mode {
                CapModeArg::Fixed => "fixed",
                CapModeArg::P99 => "p99",
            }
            .to_string(),
            value: cap_value,
        },
        missing_genes,
    };

    write_json_atomic(&stage1_dir.join("stage1_summary.json"), &summary)?;
    crate::logging::info(format!(
        "Stage1: done, found_genes={} missing_genes={} cap={}",
        n_genes_found,
        missing_genes,
        fmt_f32(cap_value)
    ));

    Ok(Stage1Result {
        out_stage_dir: stage1_dir,
        summary,
    })
}

fn parse_stage0_barcodes_path(path: &Path) -> Result<String> {
    let bytes = std::fs::read(path).map_err(|e| {
        KiraError::new(
            ErrorKind::Path,
            format!("failed to read {}: {e}", path.display()),
        )
    })?;
    let parsed: Stage0ResolvedInput = serde_json::from_slice(&bytes).map_err(|e| {
        KiraError::new(
            ErrorKind::TsvParse,
            format!("failed to parse {}: {e}", path.display()),
        )
    })?;
    parsed
        .inputs
        .into_iter()
        .find(|x| x.kind == "barcodes")
        .map(|x| x.path_abs)
        .ok_or_else(|| {
            KiraError::new(
                ErrorKind::TsvParse,
                "stage0 resolved.json missing barcodes input",
            )
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

fn load_barcodes_index(path: &Path, n_cells: usize) -> Result<BTreeMap<String, usize>> {
    let barcodes = parse_barcodes(path)?;
    let mut map = BTreeMap::new();
    for (idx, cell_id) in barcodes.into_iter().enumerate() {
        if idx >= n_cells {
            break;
        }
        map.entry(cell_id).or_insert(idx);
    }
    Ok(map)
}

fn parse_groups_normalized(path: &Path) -> Result<Vec<GroupRow>> {
    let mmap = MmapFile::open(path)?;
    let mut rdr = TsvReader::new(mmap.as_bytes());
    let mut row = Vec::new();
    let mut out = Vec::new();
    let mut header_done = false;

    while let Some(rec) = rdr.next_record(&mut row) {
        rec?;
        if row.is_empty() {
            continue;
        }
        if !header_done {
            header_done = true;
            if row != ["cell_id", "group"] {
                return Err(KiraError::new(
                    ErrorKind::TsvHeader,
                    format!("groups_normalized header mismatch in {}", path.display()),
                ));
            }
            continue;
        }
        if row.len() != 2 {
            return Err(KiraError::new(
                ErrorKind::TsvParse,
                format!("invalid groups_normalized row in {}", path.display()),
            ));
        }
        out.push(GroupRow {
            cell_id: row[0].clone(),
            group: row[1].clone(),
        });
    }

    if !header_done {
        return Err(KiraError::new(
            ErrorKind::TsvHeader,
            format!("missing header in {}", path.display()),
        ));
    }

    Ok(out)
}

fn parse_group_sizes(path: &Path) -> Result<BTreeMap<String, usize>> {
    let mmap = MmapFile::open(path)?;
    let mut rdr = TsvReader::new(mmap.as_bytes());
    let mut row = Vec::new();
    let mut out = BTreeMap::new();
    let mut header_done = false;

    while let Some(rec) = rdr.next_record(&mut row) {
        rec?;
        if row.is_empty() {
            continue;
        }
        if !header_done {
            header_done = true;
            if row != ["group", "n_cells"] {
                return Err(KiraError::new(
                    ErrorKind::TsvHeader,
                    format!("group_sizes header mismatch in {}", path.display()),
                ));
            }
            continue;
        }
        if row.len() != 2 {
            return Err(KiraError::new(
                ErrorKind::TsvParse,
                format!("invalid group_sizes row in {}", path.display()),
            ));
        }
        let n_cells: usize = row[1].parse().map_err(|e| {
            KiraError::new(
                ErrorKind::TsvParse,
                format!("invalid n_cells in {}: {e}", path.display()),
            )
        })?;
        out.insert(row[0].clone(), n_cells);
    }

    if !header_done {
        return Err(KiraError::new(
            ErrorKind::TsvHeader,
            format!("missing header in {}", path.display()),
        ));
    }

    Ok(out)
}

fn parse_genes_requested(path: &Path) -> Result<BTreeSet<String>> {
    let mmap = MmapFile::open(path)?;
    let mut rdr = TsvReader::new(mmap.as_bytes());
    let mut row = Vec::new();
    let mut out = BTreeSet::new();
    let mut header_done = false;
    let mut query_col = None;

    while let Some(rec) = rdr.next_record(&mut row) {
        rec?;
        if row.is_empty() {
            continue;
        }
        if !header_done {
            header_done = true;
            query_col = row.iter().position(|x| x == "query_symbol");
            if query_col.is_none() {
                return Err(KiraError::new(
                    ErrorKind::TsvHeader,
                    format!("genes_requested missing query_symbol in {}", path.display()),
                ));
            }
            continue;
        }
        let idx = query_col.unwrap_or(1);
        if row.len() <= idx {
            return Err(KiraError::new(
                ErrorKind::TsvParse,
                format!("invalid genes_requested row in {}", path.display()),
            ));
        }
        out.insert(row[idx].clone());
    }

    if !header_done {
        return Err(KiraError::new(
            ErrorKind::TsvHeader,
            format!("missing header in {}", path.display()),
        ));
    }

    Ok(out)
}

fn write_group_gene_agg(
    path: &Path,
    groups: &[String],
    genes: &[GeneSpec],
    group_sizes: &BTreeMap<String, usize>,
    expr: &[f32],
    cov: &[f32],
) -> Result<()> {
    let n_genes = genes.len();
    let mut lines = Vec::with_capacity(groups.len().saturating_mul(n_genes) + 1);
    lines.push("group\tgene\texpr\tcov\tn_cells".to_string());

    for (gidx, group) in groups.iter().enumerate() {
        let n_cells = *group_sizes.get(group).unwrap_or(&0);
        for (lidx, gene) in genes.iter().enumerate() {
            let idx = gidx * n_genes + lidx;
            lines.push(format!(
                "{}\t{}\t{}\t{}\t{}",
                group,
                gene.symbol,
                fmt_f32(expr[idx]),
                fmt_f32(cov[idx]),
                n_cells
            ));
        }
    }

    write_tsv_atomic(path, &lines)
}

fn write_gene_stats(
    path: &Path,
    genes: &[GeneSpec],
    n_groups: usize,
    expr: &[f32],
    cap_used: f32,
) -> Result<()> {
    let mut lines = Vec::with_capacity(genes.len() + 1);
    lines.push("gene\tmean_over_groups\tmissing_fraction\tcap_used".to_string());
    let n_genes = genes.len();

    for (lidx, gene) in genes.iter().enumerate() {
        let mut sum = 0.0f32;
        for gidx in 0..n_groups {
            sum += expr[gidx * n_genes + lidx];
        }
        let mean = if n_groups == 0 {
            0.0
        } else {
            sum / (n_groups as f32)
        };
        let missing_fraction = if gene.found_expr_idx.is_some() {
            0.0
        } else {
            1.0
        };
        lines.push(format!(
            "{}\t{}\t{}\t{}",
            gene.symbol,
            fmt_f32(mean),
            fmt_f32(missing_fraction),
            fmt_f32(cap_used)
        ));
    }

    write_tsv_atomic(path, &lines)
}

fn write_cache_bin(
    path: &Path,
    groups: &[String],
    genes: &[GeneSpec],
    expr: &[f32],
    cov: &[f32],
) -> Result<()> {
    let mut out = Vec::new();
    out.extend_from_slice(&CACHE_MAGIC);
    out.extend_from_slice(&CACHE_VERSION.to_le_bytes());
    out.extend_from_slice(&(groups.len() as u32).to_le_bytes());
    out.extend_from_slice(&(genes.len() as u32).to_le_bytes());

    for group in groups {
        let b = group.as_bytes();
        out.extend_from_slice(&(b.len() as u32).to_le_bytes());
        out.extend_from_slice(b);
    }
    for gene in genes {
        let b = gene.symbol.as_bytes();
        out.extend_from_slice(&(b.len() as u32).to_le_bytes());
        out.extend_from_slice(b);
    }

    for &v in expr {
        out.extend_from_slice(&v.to_le_bytes());
    }
    for &v in cov {
        out.extend_from_slice(&v.to_le_bytes());
    }

    write_bytes_atomic(path, &out)
}

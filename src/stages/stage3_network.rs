use crate::error::{ErrorKind, KiraError, Result};
use crate::io::atomic::{write_json_atomic, write_tsv_atomic};
use crate::io::format::fmt_f64;
use crate::io::mmap::MmapFile;
use crate::io::tsv::TsvReader;
use serde::Serialize;
use serde_json::{Map, Value, json};
use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

#[derive(Debug, Clone)]
pub struct Stage3Config {
    pub out_dir: PathBuf,
}

#[derive(Debug, Clone)]
pub struct Stage3Result {
    pub out_root_dir: PathBuf,
}

#[derive(Debug, Clone)]
struct Edge {
    source_group: String,
    target_group: String,
    ligand: String,
    receptor: String,
    score: f64,
    l_expr: f64,
    r_expr: f64,
    cov_l: f64,
    cov_r: f64,
    spec_l: f64,
    spec_r: f64,
    flags: String,
}

#[derive(Debug, Clone, Serialize)]
struct TopEdge {
    source_group: String,
    target_group: String,
    ligand: String,
    receptor: String,
    score: f64,
}

pub fn run_stage3(cfg: Stage3Config) -> Result<Stage3Result> {
    crate::logging::info("Stage3: loading Stage2 edges");
    let root = cfg.out_dir.join("kira-microenvironment");
    let stage2_dir = root.join("stage2_score");

    let edges_raw_path = stage2_dir.join("edges_raw.tsv");
    let summary2_path = stage2_dir.join("stage2_summary.json");
    for p in [&edges_raw_path, &summary2_path] {
        if !p.exists() {
            return Err(KiraError::new(
                ErrorKind::InputMissing,
                format!("missing stage2 output: {}", p.display()),
            ));
        }
    }

    let edges_in_file_order = parse_edges_raw(&edges_raw_path)?;
    crate::logging::info(format!(
        "Stage3: building network outputs for {} edges",
        edges_in_file_order.len()
    ));

    let mut edges_sorted = edges_in_file_order.clone();
    edges_sorted.sort_by(|a, b| {
        b.score
            .total_cmp(&a.score)
            .then(a.ligand.cmp(&b.ligand))
            .then(a.receptor.cmp(&b.receptor))
            .then(a.source_group.cmp(&b.source_group))
            .then(a.target_group.cmp(&b.target_group))
    });

    write_edges_tsv(&root.join("edges.tsv"), &edges_sorted)?;

    let group_strength_rows = compute_group_strength(&edges_in_file_order);
    write_group_strength(&root.join("group_strength.tsv"), &group_strength_rows)?;

    let top_pairs = compute_top_pairs(&edges_in_file_order);
    write_top_pairs(&root.join("top_pairs.tsv"), &top_pairs)?;

    let stage2_summary = read_json_value(&summary2_path)?;
    let summary = build_summary(stage2_summary, &edges_sorted, &group_strength_rows)?;
    write_json_atomic(&root.join("summary.json"), &summary)?;
    crate::logging::info("Stage3: done");

    Ok(Stage3Result { out_root_dir: root })
}

fn parse_edges_raw(path: &Path) -> Result<Vec<Edge>> {
    let mmap = MmapFile::open(path)?;
    let mut rdr = TsvReader::new(mmap.as_bytes());
    let mut row = Vec::new();
    let mut out = Vec::new();
    let mut header_seen = false;

    while let Some(rec) = rdr.next_record(&mut row) {
        rec?;
        if row.is_empty() {
            continue;
        }
        if !header_seen {
            header_seen = true;
            if row
                != [
                    "source_group",
                    "target_group",
                    "ligand",
                    "receptor",
                    "score",
                    "L_expr",
                    "R_expr",
                    "cov_L",
                    "cov_R",
                    "spec_L",
                    "spec_R",
                    "flags",
                ]
            {
                return Err(KiraError::new(
                    ErrorKind::TsvHeader,
                    format!("edges_raw header mismatch in {}", path.display()),
                ));
            }
            continue;
        }
        if row.len() != 12 {
            return Err(KiraError::new(
                ErrorKind::TsvParse,
                format!("invalid edges_raw row in {}", path.display()),
            ));
        }

        out.push(Edge {
            source_group: row[0].clone(),
            target_group: row[1].clone(),
            ligand: row[2].clone(),
            receptor: row[3].clone(),
            score: parse_f64(&row[4], "score", path)?,
            l_expr: parse_f64(&row[5], "L_expr", path)?,
            r_expr: parse_f64(&row[6], "R_expr", path)?,
            cov_l: parse_f64(&row[7], "cov_L", path)?,
            cov_r: parse_f64(&row[8], "cov_R", path)?,
            spec_l: parse_f64(&row[9], "spec_L", path)?,
            spec_r: parse_f64(&row[10], "spec_R", path)?,
            flags: row[11].clone(),
        });
    }

    if !header_seen {
        return Err(KiraError::new(
            ErrorKind::TsvHeader,
            format!("missing header in {}", path.display()),
        ));
    }

    Ok(out)
}

fn parse_f64(s: &str, name: &str, path: &Path) -> Result<f64> {
    s.parse::<f64>().map_err(|e| {
        KiraError::new(
            ErrorKind::TsvParse,
            format!("invalid {name} in {}: {e}", path.display()),
        )
    })
}

fn write_edges_tsv(path: &Path, edges: &[Edge]) -> Result<()> {
    let mut lines = Vec::with_capacity(edges.len() + 1);
    lines.push(
        "source_group\ttarget_group\tligand\treceptor\tscore\tL_expr\tR_expr\tcov_L\tcov_R\tspec_L\tspec_R\tflags"
            .to_string(),
    );
    for e in edges {
        lines.push(format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            e.source_group,
            e.target_group,
            e.ligand,
            e.receptor,
            fmt_f64(e.score),
            fmt_f64(e.l_expr),
            fmt_f64(e.r_expr),
            fmt_f64(e.cov_l),
            fmt_f64(e.cov_r),
            fmt_f64(e.spec_l),
            fmt_f64(e.spec_r),
            e.flags
        ));
    }
    write_tsv_atomic(path, &lines)
}

#[derive(Debug, Clone)]
struct GroupStrengthRow {
    group: String,
    out_strength: f64,
    in_strength: f64,
    top_targets: String,
    top_sources: String,
}

fn compute_group_strength(edges: &[Edge]) -> Vec<GroupStrengthRow> {
    let mut groups: BTreeSet<String> = BTreeSet::new();
    let mut out_strength: BTreeMap<String, f64> = BTreeMap::new();
    let mut in_strength: BTreeMap<String, f64> = BTreeMap::new();
    let mut out_to_targets: BTreeMap<String, BTreeMap<String, f64>> = BTreeMap::new();
    let mut in_from_sources: BTreeMap<String, BTreeMap<String, f64>> = BTreeMap::new();

    for e in edges {
        groups.insert(e.source_group.clone());
        groups.insert(e.target_group.clone());

        *out_strength.entry(e.source_group.clone()).or_insert(0.0) += e.score;
        *in_strength.entry(e.target_group.clone()).or_insert(0.0) += e.score;

        *out_to_targets
            .entry(e.source_group.clone())
            .or_default()
            .entry(e.target_group.clone())
            .or_insert(0.0) += e.score;

        *in_from_sources
            .entry(e.target_group.clone())
            .or_default()
            .entry(e.source_group.clone())
            .or_insert(0.0) += e.score;
    }

    let mut rows = Vec::new();
    for group in groups {
        let out_map = out_to_targets.get(&group).cloned().unwrap_or_default();
        let in_map = in_from_sources.get(&group).cloned().unwrap_or_default();

        rows.push(GroupStrengthRow {
            group: group.clone(),
            out_strength: *out_strength.get(&group).unwrap_or(&0.0),
            in_strength: *in_strength.get(&group).unwrap_or(&0.0),
            top_targets: format_top_list(&out_map),
            top_sources: format_top_list(&in_map),
        });
    }
    rows
}

fn format_top_list(scores: &BTreeMap<String, f64>) -> String {
    let mut v: Vec<(String, f64)> = scores.iter().map(|(k, v)| (k.clone(), *v)).collect();
    v.sort_by(|a, b| b.1.total_cmp(&a.1).then(a.0.cmp(&b.0)));
    if v.len() > 5 {
        v.truncate(5);
    }
    v.into_iter()
        .map(|(name, score)| format!("{}:{}", name, fmt_f64(score)))
        .collect::<Vec<_>>()
        .join(";")
}

fn write_group_strength(path: &Path, rows: &[GroupStrengthRow]) -> Result<()> {
    let mut lines = Vec::with_capacity(rows.len() + 1);
    lines.push("group\tout_strength\tin_strength\ttop_targets\ttop_sources".to_string());
    for r in rows {
        lines.push(format!(
            "{}\t{}\t{}\t{}\t{}",
            r.group,
            fmt_f64(r.out_strength),
            fmt_f64(r.in_strength),
            r.top_targets,
            r.top_sources
        ));
    }
    write_tsv_atomic(path, &lines)
}

#[derive(Debug, Clone)]
struct TopPairRow {
    ligand: String,
    receptor: String,
    top_source_group: String,
    top_target_group: String,
    score: f64,
}

fn compute_top_pairs(edges: &[Edge]) -> Vec<TopPairRow> {
    let mut pair_map: BTreeMap<(String, String), TopPairRow> = BTreeMap::new();

    for e in edges {
        pair_map
            .entry((e.ligand.clone(), e.receptor.clone()))
            .and_modify(|cur| {
                if e.score > cur.score
                    || (e.score == cur.score
                        && (e.source_group.as_str(), e.target_group.as_str())
                            < (cur.top_source_group.as_str(), cur.top_target_group.as_str()))
                {
                    cur.top_source_group = e.source_group.clone();
                    cur.top_target_group = e.target_group.clone();
                    cur.score = e.score;
                }
            })
            .or_insert(TopPairRow {
                ligand: e.ligand.clone(),
                receptor: e.receptor.clone(),
                top_source_group: e.source_group.clone(),
                top_target_group: e.target_group.clone(),
                score: e.score,
            });
    }

    let mut rows: Vec<TopPairRow> = pair_map.into_values().collect();
    rows.sort_by(|a, b| {
        b.score
            .total_cmp(&a.score)
            .then(a.ligand.cmp(&b.ligand))
            .then(a.receptor.cmp(&b.receptor))
    });
    rows
}

fn write_top_pairs(path: &Path, rows: &[TopPairRow]) -> Result<()> {
    let mut lines = Vec::with_capacity(rows.len() + 1);
    lines.push(
        "global_rank\tligand\treceptor\ttop_source_group\ttop_target_group\tscore".to_string(),
    );
    for (i, r) in rows.iter().enumerate() {
        lines.push(format!(
            "{}\t{}\t{}\t{}\t{}\t{}",
            i + 1,
            r.ligand,
            r.receptor,
            r.top_source_group,
            r.top_target_group,
            fmt_f64(r.score)
        ));
    }
    write_tsv_atomic(path, &lines)
}

fn read_json_value(path: &Path) -> Result<Value> {
    let bytes = std::fs::read(path).map_err(|e| {
        KiraError::new(
            ErrorKind::Path,
            format!("failed to read {}: {e}", path.display()),
        )
    })?;
    serde_json::from_slice::<Value>(&bytes).map_err(|e| {
        KiraError::new(
            ErrorKind::TsvParse,
            format!("failed to parse {}: {e}", path.display()),
        )
    })
}

fn build_summary(
    stage2_summary: Value,
    edges_sorted: &[Edge],
    group_rows: &[GroupStrengthRow],
) -> Result<Value> {
    let mut obj: Map<String, Value> = match stage2_summary {
        Value::Object(m) => m,
        _ => {
            return Err(KiraError::new(
                ErrorKind::TsvParse,
                "stage2_summary.json is not a JSON object",
            ));
        }
    };

    let n_edges = edges_sorted.len();
    let n_groups = group_rows.len();

    let mut out_hubs: Vec<(&str, f64)> = group_rows
        .iter()
        .map(|r| (r.group.as_str(), r.out_strength))
        .collect();
    out_hubs.sort_by(|a, b| b.1.total_cmp(&a.1).then(a.0.cmp(b.0)));
    if out_hubs.len() > 5 {
        out_hubs.truncate(5);
    }

    let mut in_hubs: Vec<(&str, f64)> = group_rows
        .iter()
        .map(|r| (r.group.as_str(), r.in_strength))
        .collect();
    in_hubs.sort_by(|a, b| b.1.total_cmp(&a.1).then(a.0.cmp(b.0)));
    if in_hubs.len() > 5 {
        in_hubs.truncate(5);
    }

    let top_edges_by_score: Vec<TopEdge> = edges_sorted
        .iter()
        .take(10)
        .map(|e| TopEdge {
            source_group: e.source_group.clone(),
            target_group: e.target_group.clone(),
            ligand: e.ligand.clone(),
            receptor: e.receptor.clone(),
            score: e.score,
        })
        .collect();

    obj.insert(
        "network".to_string(),
        json!({
            "n_edges": n_edges,
            "n_groups": n_groups,
            "out_strength_hubs": out_hubs
                .iter()
                .map(|(g, s)| json!({"group": g, "strength": s}))
                .collect::<Vec<_>>(),
            "in_strength_hubs": in_hubs
                .iter()
                .map(|(g, s)| json!({"group": g, "strength": s}))
                .collect::<Vec<_>>(),
        }),
    );
    obj.insert(
        "top_edges_by_score".to_string(),
        serde_json::to_value(top_edges_by_score).unwrap_or_else(|_| Value::Array(Vec::new())),
    );

    Ok(Value::Object(obj))
}

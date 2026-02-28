use crate::error::{ErrorKind, KiraError, Result};
use crate::io::atomic::{write_json_atomic, write_tsv_atomic};
use crate::io::format::fmt_f64;
use crate::io::mmap::MmapFile;
use crate::io::tsv::TsvReader;
use crate::resources::regime_map::{RegimeClass, parse_regime_map};
use serde::Serialize;
use serde_json::{Map, Value, json};
use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

#[derive(Debug, Clone)]
pub struct Stage4Config {
    pub out_dir: PathBuf,
    pub resources_dir: PathBuf,
    pub secretion_dir: PathBuf,
    pub regime_map_path: Option<PathBuf>,
    pub loud_thresh: f64,
    pub silent_thresh: f64,
    pub min_support: usize,
    pub eps: f64,
}

#[derive(Debug, Clone)]
pub struct Stage4Result {
    pub out_stage_dir: PathBuf,
}

#[derive(Debug, Clone)]
struct Edge {
    source_group: String,
    target_group: String,
    ligand: String,
    receptor: String,
    score: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum GroupClass {
    Loud,
    Silent,
    Mixed,
}

#[derive(Debug, Clone)]
struct PairAgg {
    sum_loud: f64,
    n_loud: usize,
    sum_silent: f64,
    n_silent: usize,
    ex_loud: Option<(String, String, f64)>,
    ex_silent: Option<(String, String, f64)>,
}

#[derive(Debug, Clone, Serialize)]
struct DriverRow {
    rank: usize,
    ligand: String,
    receptor: String,
    fold: f64,
    mean_primary: f64,
    mean_other: f64,
    n_primary: usize,
    n_other: usize,
    exemplar_source_group: String,
    exemplar_target_group: String,
    exemplar_score: f64,
}

#[derive(Debug, Clone, Serialize)]
struct LinkCounts {
    n_groups_total: usize,
    n_groups_loud: usize,
    n_groups_silent: usize,
    n_groups_mixed: usize,
    n_pairs_scored: usize,
}

#[derive(Debug, Clone, Serialize)]
struct LinkConfigOut {
    loud_thresh: f64,
    silent_thresh: f64,
    min_support: usize,
    eps: f64,
}

#[derive(Debug, Clone, Serialize)]
struct LinkSummary {
    config: LinkConfigOut,
    counts: LinkCounts,
    unknown_regimes: Vec<String>,
    top_loud_pairs: Vec<DriverRow>,
    top_silent_pairs: Vec<DriverRow>,
    top_loud_edges: Vec<TopEdge>,
    top_silent_edges: Vec<TopEdge>,
}

#[derive(Debug, Clone, Serialize)]
struct TopEdge {
    source_group: String,
    target_group: String,
    ligand: String,
    receptor: String,
    score: f64,
}

pub fn run_stage4(cfg: Stage4Config) -> Result<Stage4Result> {
    crate::logging::info("Stage4: loading edges, groups and secretion inputs");
    let root = cfg.out_dir.join("kira-microenvironment");
    let stage0_dir = root.join("stage0_resolve");
    let stage4_dir = root.join("stage4_link");

    let edges_path = root.join("edges.tsv");
    let root_summary_path = root.join("summary.json");
    let groups_norm_path = stage0_dir.join("groups_normalized.tsv");

    for p in [&edges_path, &root_summary_path, &groups_norm_path] {
        if !p.exists() {
            return Err(KiraError::new(
                ErrorKind::InputMissing,
                format!("missing required input: {}", p.display()),
            ));
        }
    }

    let edges = parse_edges(&edges_path)?;
    let groups_from_edges = collect_groups_from_edges(&edges);

    let regime_map_path = cfg
        .regime_map_path
        .clone()
        .unwrap_or_else(|| cfg.resources_dir.join("regime_map.tsv"));
    let regime_map = parse_regime_map(&regime_map_path)?;

    let group_fractions = load_group_regime_fractions(&cfg.secretion_dir, &groups_norm_path)?;
    crate::logging::info(format!(
        "Stage4: groups_with_regime_fractions={}",
        group_fractions.len()
    ));

    let (group_classes, unknown_regimes) = classify_groups(
        &groups_from_edges,
        &group_fractions,
        &regime_map,
        cfg.loud_thresh,
        cfg.silent_thresh,
    );

    let mut pair_map: BTreeMap<(String, String), PairAgg> = BTreeMap::new();
    let mut loud_edges = Vec::new();
    let mut silent_edges = Vec::new();

    for e in &edges {
        match group_classes
            .get(&e.source_group)
            .copied()
            .unwrap_or(GroupClass::Mixed)
        {
            GroupClass::Loud => {
                loud_edges.push(e.clone());
                let ent = pair_map
                    .entry((e.ligand.clone(), e.receptor.clone()))
                    .or_insert(PairAgg {
                        sum_loud: 0.0,
                        n_loud: 0,
                        sum_silent: 0.0,
                        n_silent: 0,
                        ex_loud: None,
                        ex_silent: None,
                    });
                ent.sum_loud += e.score;
                ent.n_loud += 1;
                update_exemplar(&mut ent.ex_loud, e);
            }
            GroupClass::Silent => {
                silent_edges.push(e.clone());
                let ent = pair_map
                    .entry((e.ligand.clone(), e.receptor.clone()))
                    .or_insert(PairAgg {
                        sum_loud: 0.0,
                        n_loud: 0,
                        sum_silent: 0.0,
                        n_silent: 0,
                        ex_loud: None,
                        ex_silent: None,
                    });
                ent.sum_silent += e.score;
                ent.n_silent += 1;
                update_exemplar(&mut ent.ex_silent, e);
            }
            GroupClass::Mixed => {}
        }
    }

    let mut loud_rows = Vec::new();
    let mut silent_rows = Vec::new();

    for ((lig, rec), agg) in &pair_map {
        let support = agg.n_loud + agg.n_silent;
        if support < cfg.min_support {
            continue;
        }

        let mean_loud = if agg.n_loud == 0 {
            0.0
        } else {
            agg.sum_loud / (agg.n_loud as f64)
        };
        let mean_silent = if agg.n_silent == 0 {
            0.0
        } else {
            agg.sum_silent / (agg.n_silent as f64)
        };

        let fold_loud = (mean_loud + cfg.eps) / (mean_silent + cfg.eps);
        let fold_silent = (mean_silent + cfg.eps) / (mean_loud + cfg.eps);

        let (ls, lt, lscore) = agg
            .ex_loud
            .clone()
            .unwrap_or_else(|| ("".to_string(), "".to_string(), 0.0));
        loud_rows.push(DriverRow {
            rank: 0,
            ligand: lig.clone(),
            receptor: rec.clone(),
            fold: fold_loud,
            mean_primary: mean_loud,
            mean_other: mean_silent,
            n_primary: agg.n_loud,
            n_other: agg.n_silent,
            exemplar_source_group: ls,
            exemplar_target_group: lt,
            exemplar_score: lscore,
        });

        let (ss, st, sscore) = agg
            .ex_silent
            .clone()
            .unwrap_or_else(|| ("".to_string(), "".to_string(), 0.0));
        silent_rows.push(DriverRow {
            rank: 0,
            ligand: lig.clone(),
            receptor: rec.clone(),
            fold: fold_silent,
            mean_primary: mean_silent,
            mean_other: mean_loud,
            n_primary: agg.n_silent,
            n_other: agg.n_loud,
            exemplar_source_group: ss,
            exemplar_target_group: st,
            exemplar_score: sscore,
        });
    }

    loud_rows.sort_by(|a, b| {
        b.fold
            .total_cmp(&a.fold)
            .then(b.mean_primary.total_cmp(&a.mean_primary))
            .then(a.ligand.cmp(&b.ligand))
            .then(a.receptor.cmp(&b.receptor))
    });
    for (i, row) in loud_rows.iter_mut().enumerate() {
        row.rank = i + 1;
    }

    silent_rows.sort_by(|a, b| {
        b.fold
            .total_cmp(&a.fold)
            .then(b.mean_primary.total_cmp(&a.mean_primary))
            .then(a.ligand.cmp(&b.ligand))
            .then(a.receptor.cmp(&b.receptor))
    });
    for (i, row) in silent_rows.iter_mut().enumerate() {
        row.rank = i + 1;
    }

    write_loudness_tsv(&stage4_dir.join("loudness_drivers.tsv"), &loud_rows)?;
    write_silence_tsv(&stage4_dir.join("silence_drivers.tsv"), &silent_rows)?;

    loud_edges.sort_by(edge_sort_key_desc);
    silent_edges.sort_by(edge_sort_key_desc);
    if loud_edges.len() > 50 {
        loud_edges.truncate(50);
    }
    if silent_edges.len() > 50 {
        silent_edges.truncate(50);
    }

    let counts = summarize_counts(&group_classes, loud_rows.len());
    crate::logging::info(format!(
        "Stage4: class_counts loud={} silent={} mixed={}",
        counts.n_groups_loud, counts.n_groups_silent, counts.n_groups_mixed
    ));

    let link_summary = LinkSummary {
        config: LinkConfigOut {
            loud_thresh: cfg.loud_thresh,
            silent_thresh: cfg.silent_thresh,
            min_support: cfg.min_support,
            eps: cfg.eps,
        },
        counts: counts.clone(),
        unknown_regimes,
        top_loud_pairs: loud_rows.iter().take(20).cloned().collect(),
        top_silent_pairs: silent_rows.iter().take(20).cloned().collect(),
        top_loud_edges: loud_edges
            .iter()
            .map(|e| TopEdge {
                source_group: e.source_group.clone(),
                target_group: e.target_group.clone(),
                ligand: e.ligand.clone(),
                receptor: e.receptor.clone(),
                score: e.score,
            })
            .collect(),
        top_silent_edges: silent_edges
            .iter()
            .map(|e| TopEdge {
                source_group: e.source_group.clone(),
                target_group: e.target_group.clone(),
                ligand: e.ligand.clone(),
                receptor: e.receptor.clone(),
                score: e.score,
            })
            .collect(),
    };

    write_json_atomic(&stage4_dir.join("link_summary.json"), &link_summary)?;

    let mut root_summary = read_json_object(&root_summary_path)?;
    root_summary.insert(
        "linking".to_string(),
        json!({
            "counts": counts,
            "top_loud_pairs": loud_rows.iter().take(10).cloned().collect::<Vec<_>>(),
            "top_silent_pairs": silent_rows.iter().take(10).cloned().collect::<Vec<_>>(),
            "note": "communication potential; not causal"
        }),
    );
    write_json_atomic(&root_summary_path, &Value::Object(root_summary))?;
    crate::logging::info("Stage4: wrote linking outputs");

    Ok(Stage4Result {
        out_stage_dir: stage4_dir,
    })
}

fn edge_sort_key_desc(a: &Edge, b: &Edge) -> std::cmp::Ordering {
    b.score
        .total_cmp(&a.score)
        .then(a.ligand.cmp(&b.ligand))
        .then(a.receptor.cmp(&b.receptor))
        .then(a.source_group.cmp(&b.source_group))
        .then(a.target_group.cmp(&b.target_group))
}

fn summarize_counts(
    group_classes: &BTreeMap<String, GroupClass>,
    n_pairs_scored: usize,
) -> LinkCounts {
    let mut n_loud = 0usize;
    let mut n_silent = 0usize;
    let mut n_mixed = 0usize;
    for class in group_classes.values() {
        match class {
            GroupClass::Loud => n_loud += 1,
            GroupClass::Silent => n_silent += 1,
            GroupClass::Mixed => n_mixed += 1,
        }
    }
    LinkCounts {
        n_groups_total: group_classes.len(),
        n_groups_loud: n_loud,
        n_groups_silent: n_silent,
        n_groups_mixed: n_mixed,
        n_pairs_scored,
    }
}

fn update_exemplar(slot: &mut Option<(String, String, f64)>, edge: &Edge) {
    match slot {
        Some((src, tgt, score)) => {
            if edge.score > *score
                || (edge.score == *score
                    && (edge.source_group.as_str(), edge.target_group.as_str())
                        < (src.as_str(), tgt.as_str()))
            {
                *src = edge.source_group.clone();
                *tgt = edge.target_group.clone();
                *score = edge.score;
            }
        }
        None => {
            *slot = Some((
                edge.source_group.clone(),
                edge.target_group.clone(),
                edge.score,
            ));
        }
    }
}

fn write_loudness_tsv(path: &Path, rows: &[DriverRow]) -> Result<()> {
    let mut lines = Vec::with_capacity(rows.len() + 1);
    lines.push("rank\tligand\treceptor\tfold\tmean_loud\tmean_silent\tn_edges_loud\tn_edges_silent\texemplar_source_group\texemplar_target_group\texemplar_score".to_string());
    for r in rows {
        lines.push(format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            r.rank,
            r.ligand,
            r.receptor,
            fmt_f64(r.fold),
            fmt_f64(r.mean_primary),
            fmt_f64(r.mean_other),
            r.n_primary,
            r.n_other,
            r.exemplar_source_group,
            r.exemplar_target_group,
            fmt_f64(r.exemplar_score)
        ));
    }
    write_tsv_atomic(path, &lines)
}

fn write_silence_tsv(path: &Path, rows: &[DriverRow]) -> Result<()> {
    let mut lines = Vec::with_capacity(rows.len() + 1);
    lines.push("rank\tligand\treceptor\tfold\tmean_silent\tmean_loud\tn_edges_silent\tn_edges_loud\texemplar_source_group\texemplar_target_group\texemplar_score".to_string());
    for r in rows {
        lines.push(format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            r.rank,
            r.ligand,
            r.receptor,
            fmt_f64(r.fold),
            fmt_f64(r.mean_primary),
            fmt_f64(r.mean_other),
            r.n_primary,
            r.n_other,
            r.exemplar_source_group,
            r.exemplar_target_group,
            fmt_f64(r.exemplar_score)
        ));
    }
    write_tsv_atomic(path, &lines)
}

fn read_json_object(path: &Path) -> Result<Map<String, Value>> {
    let bytes = std::fs::read(path).map_err(|e| {
        KiraError::new(
            ErrorKind::Path,
            format!("failed to read {}: {e}", path.display()),
        )
    })?;
    let v: Value = serde_json::from_slice(&bytes).map_err(|e| {
        KiraError::new(
            ErrorKind::TsvParse,
            format!("failed to parse {}: {e}", path.display()),
        )
    })?;
    match v {
        Value::Object(m) => Ok(m),
        _ => Err(KiraError::new(
            ErrorKind::TsvParse,
            format!("{} is not a JSON object", path.display()),
        )),
    }
}

fn parse_edges(path: &Path) -> Result<Vec<Edge>> {
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
                    format!("edges header mismatch in {}", path.display()),
                ));
            }
            continue;
        }
        if row.len() != 12 {
            return Err(KiraError::new(
                ErrorKind::TsvParse,
                format!("invalid edges row in {}", path.display()),
            ));
        }
        let score: f64 = row[4].parse().map_err(|e| {
            KiraError::new(
                ErrorKind::TsvParse,
                format!("invalid edge score in {}: {e}", path.display()),
            )
        })?;
        out.push(Edge {
            source_group: row[0].clone(),
            target_group: row[1].clone(),
            ligand: row[2].clone(),
            receptor: row[3].clone(),
            score,
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

fn collect_groups_from_edges(edges: &[Edge]) -> BTreeSet<String> {
    let mut groups = BTreeSet::new();
    for e in edges {
        groups.insert(e.source_group.clone());
        groups.insert(e.target_group.clone());
    }
    groups
}

fn classify_groups(
    groups_from_edges: &BTreeSet<String>,
    group_fractions: &BTreeMap<String, BTreeMap<String, f64>>,
    regime_map: &BTreeMap<String, RegimeClass>,
    loud_thresh: f64,
    silent_thresh: f64,
) -> (BTreeMap<String, GroupClass>, Vec<String>) {
    let mut unknown = BTreeSet::new();
    let mut classes = BTreeMap::new();

    let all_groups: BTreeSet<String> = groups_from_edges
        .iter()
        .cloned()
        .chain(group_fractions.keys().cloned())
        .collect();

    for group in all_groups {
        let regimes = group_fractions.get(&group);
        let mut loud = 0.0f64;
        let mut silent = 0.0f64;

        if let Some(map) = regimes {
            for (regime, frac) in map {
                match regime_map.get(regime).copied().unwrap_or_else(|| {
                    unknown.insert(regime.clone());
                    RegimeClass::Ignore
                }) {
                    RegimeClass::Loud => loud += *frac,
                    RegimeClass::Silent => silent += *frac,
                    RegimeClass::Mixed | RegimeClass::Ignore => {}
                }
            }
        }

        let class = if loud >= loud_thresh {
            GroupClass::Loud
        } else if silent >= silent_thresh {
            GroupClass::Silent
        } else {
            GroupClass::Mixed
        };

        classes.insert(group, class);
    }

    (classes, unknown.into_iter().collect())
}

fn load_group_regime_fractions(
    secretion_dir: &Path,
    groups_normalized_path: &Path,
) -> Result<BTreeMap<String, BTreeMap<String, f64>>> {
    let group_candidates = ["secretion_groups.tsv", "groups.tsv", "regime_fractions.tsv"];
    let cell_candidates = [
        "secretion_cells.tsv",
        "cells.tsv",
        "stage6_classify.tsv",
        "classify.tsv",
        "secretion.tsv",
    ];

    let group_file = find_candidate_file(secretion_dir, &group_candidates)
        .or_else(|| find_tsv_by_header(secretion_dir, &["group", "regime", "fraction"]));
    let cell_file = find_candidate_file(secretion_dir, &cell_candidates)
        .or_else(|| find_tsv_with_required_columns(secretion_dir, &["cell_id", "regime"]));

    if group_file.is_none() && cell_file.is_none() {
        return Err(KiraError::new(
            ErrorKind::SecretionInputMissing,
            format!(
                "no supported secretion input files found in {}",
                secretion_dir.display()
            ),
        ));
    }

    if let Some(p) = group_file {
        return parse_group_fractions_file(&p);
    }

    let cell_file = cell_file.unwrap();
    let cell_regime = parse_cell_regimes_file(&cell_file)?;
    let group_by_cell = parse_groups_normalized(groups_normalized_path)?;

    let mut counts: BTreeMap<String, BTreeMap<String, usize>> = BTreeMap::new();
    let mut totals: BTreeMap<String, usize> = BTreeMap::new();

    for (cell_id, group) in group_by_cell {
        let Some(regime) = cell_regime.get(&cell_id) else {
            continue;
        };
        *counts
            .entry(group.clone())
            .or_default()
            .entry(regime.clone())
            .or_insert(0) += 1;
        *totals.entry(group).or_insert(0) += 1;
    }

    let mut out: BTreeMap<String, BTreeMap<String, f64>> = BTreeMap::new();
    for (group, regimes) in counts {
        let total = *totals.get(&group).unwrap_or(&0);
        let mut fracs = BTreeMap::new();
        for (regime, n) in regimes {
            let frac = if total == 0 {
                0.0
            } else {
                (n as f64) / (total as f64)
            };
            fracs.insert(regime, frac);
        }
        out.insert(group, fracs);
    }

    Ok(out)
}

fn find_candidate_file(base: &Path, names: &[&str]) -> Option<PathBuf> {
    let search_roots = [base.to_path_buf(), base.join("out")];

    // Priority 1: direct files in base, then base/out
    for root in &search_roots {
        for name in names {
            let p = root.join(name);
            if p.exists() && p.is_file() {
                return Some(p);
            }
        }
    }

    // Priority 2: one level deep in base, then base/out (lexicographic)
    for root in &search_roots {
        if !root.exists() || !root.is_dir() {
            continue;
        }
        let mut children: Vec<PathBuf> = match std::fs::read_dir(root) {
            Ok(rd) => rd
                .filter_map(|x| x.ok().map(|e| e.path()))
                .filter(|p| p.is_dir())
                .collect(),
            Err(_) => Vec::new(),
        };
        children.sort();
        for child in children {
            for name in names {
                let p = child.join(name);
                if p.exists() && p.is_file() {
                    return Some(p);
                }
            }
        }
    }

    None
}

fn find_tsv_by_header(base: &Path, expected: &[&str]) -> Option<PathBuf> {
    for p in list_candidate_tsv_paths(base) {
        let Ok(mmap) = MmapFile::open(&p) else {
            continue;
        };
        let mut rdr = TsvReader::new(mmap.as_bytes());
        let mut row = Vec::new();
        while let Some(rec) = rdr.next_record(&mut row) {
            if rec.is_err() {
                break;
            }
            if row.is_empty() {
                continue;
            }
            if row.len() == expected.len() && row.iter().zip(expected.iter()).all(|(a, b)| a == b) {
                return Some(p);
            }
            break;
        }
    }
    None
}

fn find_tsv_with_required_columns(base: &Path, required: &[&str]) -> Option<PathBuf> {
    for p in list_candidate_tsv_paths(base) {
        let Ok(mmap) = MmapFile::open(&p) else {
            continue;
        };
        let mut rdr = TsvReader::new(mmap.as_bytes());
        let mut row = Vec::new();
        while let Some(rec) = rdr.next_record(&mut row) {
            if rec.is_err() {
                break;
            }
            if row.is_empty() {
                continue;
            }
            if required.iter().all(|col| row.iter().any(|h| h == col)) {
                return Some(p);
            }
            break;
        }
    }
    None
}

fn list_candidate_tsv_paths(base: &Path) -> Vec<PathBuf> {
    let mut paths = Vec::new();
    let search_roots = [base.to_path_buf(), base.join("out")];
    for root in &search_roots {
        if !root.exists() || !root.is_dir() {
            continue;
        }

        // root/*.tsv
        if let Ok(rd) = std::fs::read_dir(root) {
            let mut entries: Vec<PathBuf> = rd.filter_map(|e| e.ok().map(|x| x.path())).collect();
            entries.sort();
            for p in entries {
                if p.is_file() && p.extension().map(|x| x == "tsv").unwrap_or(false) {
                    paths.push(p);
                }
            }
        }

        // root/*/*.tsv
        if let Ok(rd) = std::fs::read_dir(root) {
            let mut dirs: Vec<PathBuf> = rd
                .filter_map(|e| e.ok().map(|x| x.path()))
                .filter(|p| p.is_dir())
                .collect();
            dirs.sort();
            for d in dirs {
                if let Ok(rd2) = std::fs::read_dir(&d) {
                    let mut files: Vec<PathBuf> =
                        rd2.filter_map(|e| e.ok().map(|x| x.path())).collect();
                    files.sort();
                    for p in files {
                        if p.is_file() && p.extension().map(|x| x == "tsv").unwrap_or(false) {
                            paths.push(p);
                        }
                    }
                }
            }
        }
    }
    paths
}

fn parse_group_fractions_file(path: &Path) -> Result<BTreeMap<String, BTreeMap<String, f64>>> {
    let mmap = MmapFile::open(path)?;
    let mut rdr = TsvReader::new(mmap.as_bytes());
    let mut row = Vec::new();
    let mut header_seen = false;
    let mut out: BTreeMap<String, BTreeMap<String, f64>> = BTreeMap::new();

    while let Some(rec) = rdr.next_record(&mut row) {
        rec?;
        if row.is_empty() {
            continue;
        }
        if !header_seen {
            header_seen = true;
            if row != ["group", "regime", "fraction"] {
                return Err(KiraError::new(
                    ErrorKind::TsvHeader,
                    format!("group fractions header mismatch in {}", path.display()),
                ));
            }
            continue;
        }
        if row.len() != 3 {
            return Err(KiraError::new(
                ErrorKind::TsvParse,
                format!("invalid group fractions row in {}", path.display()),
            ));
        }
        let frac: f64 = row[2].parse().map_err(|e| {
            KiraError::new(
                ErrorKind::TsvParse,
                format!("invalid fraction in {}: {e}", path.display()),
            )
        })?;
        *out.entry(row[0].clone())
            .or_default()
            .entry(row[1].clone())
            .or_insert(0.0) += frac;
    }

    if !header_seen {
        return Err(KiraError::new(
            ErrorKind::TsvHeader,
            format!("missing header in {}", path.display()),
        ));
    }

    Ok(out)
}

fn parse_cell_regimes_file(path: &Path) -> Result<BTreeMap<String, String>> {
    let mmap = MmapFile::open(path)?;
    let mut rdr = TsvReader::new(mmap.as_bytes());
    let mut row = Vec::new();
    let mut header_seen = false;
    let mut cell_col = None;
    let mut regime_col = None;
    let mut out = BTreeMap::new();

    while let Some(rec) = rdr.next_record(&mut row) {
        rec?;
        if row.is_empty() {
            continue;
        }
        if !header_seen {
            header_seen = true;
            cell_col = row.iter().position(|x| x == "cell_id");
            regime_col = row.iter().position(|x| x == "regime");
            if cell_col.is_none() || regime_col.is_none() {
                return Err(KiraError::new(
                    ErrorKind::TsvHeader,
                    format!("cell secretion header mismatch in {}", path.display()),
                ));
            }
            continue;
        }

        let cidx = cell_col.unwrap_or(0);
        let ridx = regime_col.unwrap_or(1);
        if row.len() <= cidx || row.len() <= ridx {
            return Err(KiraError::new(
                ErrorKind::TsvParse,
                format!("invalid cell secretion row in {}", path.display()),
            ));
        }
        out.entry(row[cidx].clone()).or_insert(row[ridx].clone());
    }

    if !header_seen {
        return Err(KiraError::new(
            ErrorKind::TsvHeader,
            format!("missing header in {}", path.display()),
        ));
    }

    Ok(out)
}

fn parse_groups_normalized(path: &Path) -> Result<Vec<(String, String)>> {
    let mmap = MmapFile::open(path)?;
    let mut rdr = TsvReader::new(mmap.as_bytes());
    let mut row = Vec::new();
    let mut header_seen = false;
    let mut out = Vec::new();

    while let Some(rec) = rdr.next_record(&mut row) {
        rec?;
        if row.is_empty() {
            continue;
        }
        if !header_seen {
            header_seen = true;
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
        out.push((row[0].clone(), row[1].clone()));
    }

    if !header_seen {
        return Err(KiraError::new(
            ErrorKind::TsvHeader,
            format!("missing header in {}", path.display()),
        ));
    }

    Ok(out)
}

use crate::agg::reader::AggReader;
use crate::error::{ErrorKind, KiraError, Result};
use crate::io::atomic::{write_json_atomic, write_tsv_atomic};
use crate::io::format::fmt_f32;
use crate::io::mmap::MmapFile;
use crate::io::tsv::TsvReader;
use crate::resources::resolve_lr_path;
use serde::Serialize;
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

const COV_WARN: f32 = 0.20;

#[derive(Debug, Clone)]
pub struct Stage2Config {
    pub out_dir: PathBuf,
    pub resources_dir: PathBuf,
    pub lr_profile: String,
    pub eps: f32,
    pub cov_min: f32,
    pub expr_min: f32,
    pub spec_on: bool,
    pub spec_cap: f32,
    pub top_n_per_pair: usize,
    pub top_n_per_source: usize,
}

#[derive(Debug, Clone)]
pub struct Stage2Result {
    pub out_stage_dir: PathBuf,
    pub summary: Stage2Summary,
}

#[derive(Debug, Clone, Serialize)]
pub struct Stage2Summary {
    pub counts: Stage2Counts,
    pub thresholds: Stage2Thresholds,
    pub cap_used: f32,
    pub skipped: Stage2Skipped,
}

#[derive(Debug, Clone, Serialize)]
pub struct Stage2Counts {
    pub n_pairs_raw: usize,
    pub n_pairs_expanded: usize,
    pub n_edges_before_filter: usize,
    pub n_edges_after_filter: usize,
}

#[derive(Debug, Clone, Serialize)]
pub struct Stage2Thresholds {
    pub cov_min: f32,
    pub expr_min: f32,
    pub spec_on: bool,
}

#[derive(Debug, Clone, Serialize)]
pub struct Stage2Skipped {
    pub missing_components: usize,
    pub missing_genes: usize,
}

#[derive(Debug, Clone)]
struct LrPairRaw {
    ligand: String,
    receptor: String,
    weight: f32,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum Logic {
    And,
    Or,
}

#[derive(Debug, Clone)]
struct CofactorEntry {
    role: String,
    complex_id: String,
    subunit_symbol: String,
    required: bool,
    logic: Logic,
}

#[derive(Debug, Clone)]
struct LabelRule {
    ligand_pat: String,
    receptor_pat: String,
    label: String,
    weight: f32,
}

#[derive(Debug, Clone)]
struct GeneStat {
    mean_over_groups: f32,
}

#[derive(Debug, Clone)]
struct ComponentRef {
    gene_idx: usize,
    mean: f32,
    required: bool,
    logic: Logic,
}

#[derive(Debug, Clone)]
enum Entity {
    Single { gene_idx: usize, mean: f32 },
    Complex { components: Vec<ComponentRef> },
}

#[derive(Debug, Clone)]
struct ResolvedPair {
    ligand: String,
    receptor: String,
    weight_pair: f32,
    weight_labels: f32,
    label_flags: Vec<String>,
    ligand_entity: Entity,
    receptor_entity: Entity,
}

#[derive(Debug, Clone)]
struct Edge {
    source_group: String,
    target_group: String,
    ligand: String,
    receptor: String,
    score: f32,
    l_expr: f32,
    r_expr: f32,
    cov_l: f32,
    cov_r: f32,
    spec_l: f32,
    spec_r: f32,
    flags: String,
}

pub fn run_stage2(cfg: Stage2Config) -> Result<Stage2Result> {
    crate::logging::info("Stage2: loading Stage1 cache and resource tables");
    let stage1_dir = cfg.out_dir.join("kira-microenvironment").join("stage1_agg");
    let stage2_dir = cfg
        .out_dir
        .join("kira-microenvironment")
        .join("stage2_score");

    let cache_path = stage1_dir.join("cache.bin");
    let gene_stats_path = stage1_dir.join("gene_stats.tsv");
    let summary1_path = stage1_dir.join("stage1_summary.json");
    for p in [&cache_path, &gene_stats_path, &summary1_path] {
        if !p.exists() {
            return Err(KiraError::new(
                ErrorKind::InputMissing,
                format!("missing stage1 output: {}", p.display()),
            ));
        }
    }

    // Ensure stage1 summary is present and parseable.
    let _: serde_json::Value =
        serde_json::from_slice(&std::fs::read(&summary1_path).map_err(|e| {
            KiraError::new(
                ErrorKind::TsvParse,
                format!("failed to read {}: {e}", summary1_path.display()),
            )
        })?)
        .map_err(|e| {
            KiraError::new(
                ErrorKind::TsvParse,
                format!("failed to parse {}: {e}", summary1_path.display()),
            )
        })?;

    let agg = AggReader::open(&cache_path)?;
    let (gene_stats, cap_used) = parse_gene_stats(&gene_stats_path)?;
    let lr_path = resolve_lr_path(&cfg.resources_dir, &cfg.lr_profile)?;
    let lr_pairs = parse_lr_pairs(&lr_path)?;
    let n_pairs_raw = lr_pairs.len();
    let cofactors = parse_cofactors(&cfg.resources_dir.join("cofactors.tsv"))?;
    let labels = parse_labels(&cfg.resources_dir.join("labels.tsv"))?;
    crate::logging::info(format!(
        "Stage2: lr_pairs={} cofactors={} labels={}",
        n_pairs_raw,
        cofactors.len(),
        labels.len()
    ));

    let mut group_names = agg.groups().to_vec();
    group_names.sort();
    let group_idx_by_name: BTreeMap<String, usize> = agg
        .groups()
        .iter()
        .enumerate()
        .map(|(i, g)| (g.clone(), i))
        .collect();

    let gene_idx_by_name: BTreeMap<String, usize> = agg
        .genes()
        .iter()
        .enumerate()
        .map(|(i, g)| (g.clone(), i))
        .collect();

    let complex_map = build_complex_map(&cofactors);

    let mut pairs_sorted = lr_pairs;
    pairs_sorted.sort_by(|a, b| a.ligand.cmp(&b.ligand).then(a.receptor.cmp(&b.receptor)));

    let mut missing_components = 0usize;
    let mut missing_genes = 0usize;
    let mut resolved_pairs = Vec::new();

    for pair in pairs_sorted {
        let (weight_labels, label_flags) = label_effects(&pair.ligand, &pair.receptor, &labels);

        let ligand_entity = match resolve_entity(
            &pair.ligand,
            "ligand",
            &complex_map,
            &gene_idx_by_name,
            &gene_stats,
            &mut missing_components,
            &mut missing_genes,
        ) {
            Some(x) => x,
            None => continue,
        };

        let receptor_entity = match resolve_entity(
            &pair.receptor,
            "receptor",
            &complex_map,
            &gene_idx_by_name,
            &gene_stats,
            &mut missing_components,
            &mut missing_genes,
        ) {
            Some(x) => x,
            None => continue,
        };

        resolved_pairs.push(ResolvedPair {
            ligand: pair.ligand,
            receptor: pair.receptor,
            weight_pair: pair.weight,
            weight_labels,
            label_flags,
            ligand_entity,
            receptor_entity,
        });
    }
    crate::logging::info(format!(
        "Stage2: resolved_pairs={} skipped_missing_components={} skipped_missing_genes={}",
        resolved_pairs.len(),
        missing_components,
        missing_genes
    ));

    let n_pairs_expanded = resolved_pairs.len();

    let mut edges_before_filter = 0usize;
    let mut filtered_edges = Vec::new();

    for rp in &resolved_pairs {
        for src in &group_names {
            let src_idx = *group_idx_by_name.get(src).unwrap_or(&0);
            let (l_expr, cov_l, mean_l) = eval_entity(&rp.ligand_entity, src_idx, &agg);

            for tgt in &group_names {
                let tgt_idx = *group_idx_by_name.get(tgt).unwrap_or(&0);
                let (r_expr, cov_r, mean_r) = eval_entity(&rp.receptor_entity, tgt_idx, &agg);

                edges_before_filter += 1;

                if cov_l < cfg.cov_min
                    || cov_r < cfg.cov_min
                    || l_expr < cfg.expr_min
                    || r_expr < cfg.expr_min
                {
                    continue;
                }

                let fl = clamp_expr(l_expr, cap_used);
                let fr = clamp_expr(r_expr, cap_used);
                let score0 = fl * fr * rp.weight_pair * rp.weight_labels;

                let (spec_l, spec_r, spec_boost) = if cfg.spec_on {
                    let sl = clamp_ratio(l_expr, mean_l, cfg.eps, cfg.spec_cap);
                    let sr = clamp_ratio(r_expr, mean_r, cfg.eps, cfg.spec_cap);
                    (sl, sr, (sl * sr).sqrt())
                } else {
                    (1.0, 1.0, 1.0)
                };

                let final_score = score0 * spec_boost;

                let mut flags = Vec::new();
                if cov_l.min(cov_r) < COV_WARN {
                    flags.push("LOW_COVERAGE_EDGE".to_string());
                }
                if cfg.spec_on && spec_boost >= 3.0 {
                    flags.push("HIGH_SPECIFICITY".to_string());
                }
                flags.extend(rp.label_flags.iter().cloned());
                flags.sort();
                flags.dedup();

                filtered_edges.push(Edge {
                    source_group: src.clone(),
                    target_group: tgt.clone(),
                    ligand: rp.ligand.clone(),
                    receptor: rp.receptor.clone(),
                    score: final_score,
                    l_expr,
                    r_expr,
                    cov_l,
                    cov_r,
                    spec_l,
                    spec_r,
                    flags: flags.join(","),
                });
            }
        }
    }

    // Step 1 limit: per (ligand, receptor, source_group)
    let mut by_pair_source: BTreeMap<(String, String, String), Vec<Edge>> = BTreeMap::new();
    for edge in filtered_edges {
        by_pair_source
            .entry((
                edge.ligand.clone(),
                edge.receptor.clone(),
                edge.source_group.clone(),
            ))
            .or_default()
            .push(edge);
    }

    let mut after_pair_limit = Vec::new();
    for (_, mut edges) in by_pair_source {
        edges.sort_by(|a, b| {
            b.score
                .total_cmp(&a.score)
                .then(a.target_group.cmp(&b.target_group))
        });
        if edges.len() > cfg.top_n_per_pair {
            edges.truncate(cfg.top_n_per_pair);
        }
        after_pair_limit.extend(edges);
    }

    // Step 2 limit: per source_group
    let mut by_source: BTreeMap<String, Vec<Edge>> = BTreeMap::new();
    for edge in after_pair_limit {
        by_source
            .entry(edge.source_group.clone())
            .or_default()
            .push(edge);
    }

    let mut final_edges = Vec::new();
    for (_, mut edges) in by_source {
        edges.sort_by(|a, b| {
            b.score
                .total_cmp(&a.score)
                .then(a.ligand.cmp(&b.ligand))
                .then(a.receptor.cmp(&b.receptor))
                .then(a.target_group.cmp(&b.target_group))
        });
        if edges.len() > cfg.top_n_per_source {
            edges.truncate(cfg.top_n_per_source);
        }
        final_edges.extend(edges);
    }

    final_edges.sort_by(|a, b| {
        b.score
            .total_cmp(&a.score)
            .then(a.ligand.cmp(&b.ligand))
            .then(a.receptor.cmp(&b.receptor))
            .then(a.source_group.cmp(&b.source_group))
            .then(a.target_group.cmp(&b.target_group))
    });

    write_edges_raw(&stage2_dir.join("edges_raw.tsv"), &final_edges)?;
    write_pairs_stats(&stage2_dir.join("pairs_stats.tsv"), &final_edges)?;

    let summary = Stage2Summary {
        counts: Stage2Counts {
            n_pairs_raw,
            n_pairs_expanded,
            n_edges_before_filter: edges_before_filter,
            n_edges_after_filter: final_edges.len(),
        },
        thresholds: Stage2Thresholds {
            cov_min: cfg.cov_min,
            expr_min: cfg.expr_min,
            spec_on: cfg.spec_on,
        },
        cap_used,
        skipped: Stage2Skipped {
            missing_components,
            missing_genes,
        },
    };

    write_json_atomic(&stage2_dir.join("stage2_summary.json"), &summary)?;
    crate::logging::info(format!(
        "Stage2: done, edges_before_filter={} edges_after_filter={}",
        summary.counts.n_edges_before_filter, summary.counts.n_edges_after_filter
    ));

    Ok(Stage2Result {
        out_stage_dir: stage2_dir,
        summary,
    })
}

fn clamp_expr(x: f32, cap: f32) -> f32 {
    x.max(0.0).min(cap)
}

fn clamp_ratio(x: f32, mean: f32, eps: f32, cap: f32) -> f32 {
    (x / (mean + eps)).clamp(0.0, cap)
}

fn eval_entity(entity: &Entity, group_idx: usize, agg: &AggReader) -> (f32, f32, f32) {
    match entity {
        Entity::Single { gene_idx, mean } => (
            agg.expr(group_idx, *gene_idx),
            agg.cov(group_idx, *gene_idx),
            *mean,
        ),
        Entity::Complex { components } => {
            let mut and_expr: Option<f32> = None;
            let mut and_cov: Option<f32> = None;
            let mut and_mean: Option<f32> = None;
            let mut or_expr: Option<f32> = None;
            let mut or_cov: Option<f32> = None;
            let mut or_mean: Option<f32> = None;

            for c in components {
                let e = agg.expr(group_idx, c.gene_idx);
                let v = agg.cov(group_idx, c.gene_idx);
                let m = c.mean;
                match c.logic {
                    Logic::And => {
                        and_expr = Some(match and_expr {
                            Some(cur) => cur.min(e),
                            None => e,
                        });
                        and_cov = Some(match and_cov {
                            Some(cur) => cur.min(v),
                            None => v,
                        });
                        and_mean = Some(match and_mean {
                            Some(cur) => cur.min(m),
                            None => m,
                        });
                    }
                    Logic::Or => {
                        or_expr = Some(match or_expr {
                            Some(cur) => cur.max(e),
                            None => e,
                        });
                        or_cov = Some(match or_cov {
                            Some(cur) => cur.max(v),
                            None => v,
                        });
                        or_mean = Some(match or_mean {
                            Some(cur) => cur.max(m),
                            None => m,
                        });
                    }
                }
            }

            // Required AND components dominate; OR contributes only when AND set is absent.
            match (and_expr, and_cov, and_mean, or_expr, or_cov, or_mean) {
                (Some(e), Some(c), Some(m), _, _, _) => (e, c, m),
                (None, None, None, Some(e), Some(c), Some(m)) => (e, c, m),
                _ => (0.0, 0.0, 0.0),
            }
        }
    }
}

fn resolve_entity(
    symbol: &str,
    role: &str,
    complex_map: &BTreeMap<(String, String), Vec<CofactorEntry>>,
    gene_idx_by_name: &BTreeMap<String, usize>,
    gene_stats: &BTreeMap<String, GeneStat>,
    missing_components: &mut usize,
    missing_genes: &mut usize,
) -> Option<Entity> {
    if let Some(entries) = complex_map.get(&(role.to_string(), symbol.to_string())) {
        let mut components = Vec::new();
        for c in entries {
            let Some(&gene_idx) = gene_idx_by_name.get(&c.subunit_symbol) else {
                if c.required {
                    *missing_components += 1;
                    return None;
                }
                *missing_genes += 1;
                continue;
            };
            let mean = gene_stats
                .get(&c.subunit_symbol)
                .map(|x| x.mean_over_groups)
                .unwrap_or(0.0);
            components.push(ComponentRef {
                gene_idx,
                mean,
                required: c.required,
                logic: c.logic,
            });
        }
        components.sort_by(|a, b| {
            a.gene_idx
                .cmp(&b.gene_idx)
                .then(a.required.cmp(&b.required))
                .then(a.logic.cmp(&b.logic))
        });

        if components.is_empty() {
            *missing_components += 1;
            return None;
        }
        Some(Entity::Complex { components })
    } else {
        let Some(&gene_idx) = gene_idx_by_name.get(symbol) else {
            *missing_genes += 1;
            return None;
        };
        let mean = gene_stats
            .get(symbol)
            .map(|x| x.mean_over_groups)
            .unwrap_or(0.0);
        Some(Entity::Single { gene_idx, mean })
    }
}

fn build_complex_map(
    cofactors: &[CofactorEntry],
) -> BTreeMap<(String, String), Vec<CofactorEntry>> {
    let mut map: BTreeMap<(String, String), Vec<CofactorEntry>> = BTreeMap::new();
    for c in cofactors {
        map.entry((c.role.clone(), c.complex_id.clone()))
            .or_default()
            .push(c.clone());
    }
    for entries in map.values_mut() {
        entries.sort_by(|a, b| {
            a.subunit_symbol
                .cmp(&b.subunit_symbol)
                .then(a.required.cmp(&b.required))
                .then(a.logic.cmp(&b.logic))
        });
    }
    map
}

fn label_effects(ligand: &str, receptor: &str, labels: &[LabelRule]) -> (f32, Vec<String>) {
    let mut matched: Vec<&LabelRule> = labels
        .iter()
        .filter(|r| {
            pattern_match(&r.ligand_pat, ligand) && pattern_match(&r.receptor_pat, receptor)
        })
        .collect();
    matched.sort_by(|a, b| {
        a.label
            .cmp(&b.label)
            .then(a.ligand_pat.cmp(&b.ligand_pat))
            .then(a.receptor_pat.cmp(&b.receptor_pat))
    });

    let mut weight = 1.0f32;
    let mut flags = Vec::new();
    for m in matched {
        weight *= m.weight;
        flags.push(m.label.clone());
    }
    (weight, flags)
}

fn pattern_match(pat: &str, value: &str) -> bool {
    pat == "*" || pat == value
}

fn write_edges_raw(path: &Path, edges: &[Edge]) -> Result<()> {
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
            fmt_f32(e.score),
            fmt_f32(e.l_expr),
            fmt_f32(e.r_expr),
            fmt_f32(e.cov_l),
            fmt_f32(e.cov_r),
            fmt_f32(e.spec_l),
            fmt_f32(e.spec_r),
            e.flags
        ));
    }
    write_tsv_atomic(path, &lines)
}

fn write_pairs_stats(path: &Path, edges: &[Edge]) -> Result<()> {
    let mut by_pair: BTreeMap<(String, String), (usize, f32)> = BTreeMap::new();
    for e in edges {
        by_pair
            .entry((e.ligand.clone(), e.receptor.clone()))
            .and_modify(|v| {
                v.0 += 1;
                if e.score > v.1 {
                    v.1 = e.score;
                }
            })
            .or_insert((1, e.score));
    }

    let mut lines = Vec::with_capacity(by_pair.len() + 1);
    lines.push("ligand\treceptor\tn_edges_kept\ttop_score".to_string());
    for ((lig, rec), (n, top)) in by_pair {
        lines.push(format!("{}\t{}\t{}\t{}", lig, rec, n, fmt_f32(top)));
    }
    write_tsv_atomic(path, &lines)
}

fn parse_gene_stats(path: &Path) -> Result<(BTreeMap<String, GeneStat>, f32)> {
    let mmap = MmapFile::open(path)?;
    let mut rdr = TsvReader::new(mmap.as_bytes());
    let mut row = Vec::new();
    let mut header_seen = false;
    let mut out = BTreeMap::new();
    let mut cap_used: Option<f32> = None;

    while let Some(rec) = rdr.next_record(&mut row) {
        rec?;
        if row.is_empty() || row.iter().all(|c| c.trim().is_empty()) {
            continue;
        }
        if row[0].trim_start().starts_with('#') {
            continue;
        }
        if !header_seen {
            header_seen = true;
            if row != ["gene", "mean_over_groups", "missing_fraction", "cap_used"] {
                return Err(KiraError::new(
                    ErrorKind::TsvHeader,
                    format!("gene_stats header mismatch in {}", path.display()),
                ));
            }
            continue;
        }
        if row.len() != 4 {
            return Err(KiraError::new(
                ErrorKind::TsvParse,
                format!("invalid row in {}", path.display()),
            ));
        }
        let mean: f32 = row[1].parse().map_err(|e| {
            KiraError::new(
                ErrorKind::TsvParse,
                format!("invalid mean_over_groups in {}: {e}", path.display()),
            )
        })?;
        let cap: f32 = row[3].parse().map_err(|e| {
            KiraError::new(
                ErrorKind::TsvParse,
                format!("invalid cap_used in {}: {e}", path.display()),
            )
        })?;
        if cap_used.is_none() {
            cap_used = Some(cap);
        }
        out.insert(
            row[0].clone(),
            GeneStat {
                mean_over_groups: mean,
            },
        );
    }

    if !header_seen {
        return Err(KiraError::new(
            ErrorKind::TsvHeader,
            format!("missing header in {}", path.display()),
        ));
    }

    Ok((out, cap_used.unwrap_or(0.0)))
}

fn parse_lr_pairs(path: &Path) -> Result<Vec<LrPairRaw>> {
    let mmap = MmapFile::open(path)?;
    let mut rdr = TsvReader::new(mmap.as_bytes());
    let mut row = Vec::new();
    let mut out = Vec::new();
    let mut header_seen = false;

    while let Some(rec) = rdr.next_record(&mut row) {
        rec?;
        if row.is_empty() || row.iter().all(|c| c.trim().is_empty()) {
            continue;
        }
        if row[0].trim_start().starts_with('#') {
            continue;
        }
        if !header_seen {
            header_seen = true;
            let ok =
                row == [
                    "ligand_symbol",
                    "receptor_symbol",
                    "family",
                    "directionality",
                ] || row
                    == [
                        "ligand_symbol",
                        "receptor_symbol",
                        "family",
                        "directionality",
                        "weight",
                    ]
                    || row
                        == [
                            "ligand_symbol",
                            "receptor_symbol",
                            "family",
                            "directionality",
                            "weight",
                            "notes",
                        ]
                    || row
                        == [
                            "ligand_symbol",
                            "receptor_symbol",
                            "family",
                            "directionality",
                            "notes",
                        ];
            if !ok {
                return Err(KiraError::new(
                    ErrorKind::TsvHeader,
                    format!("lr_pairs header mismatch in {}", path.display()),
                ));
            }
            continue;
        }

        if row.len() < 4 {
            return Err(KiraError::new(
                ErrorKind::TsvParse,
                format!("invalid lr_pairs row in {}", path.display()),
            ));
        }
        let weight = parse_lr_weight(&row).ok_or_else(|| {
            KiraError::new(
                ErrorKind::TsvParse,
                format!("invalid lr_pairs weight in {}", path.display()),
            )
        })?;
        out.push(LrPairRaw {
            ligand: row[0].clone(),
            receptor: row[1].clone(),
            weight,
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

fn parse_lr_weight(row: &[String]) -> Option<f32> {
    if row.len() < 5 {
        return Some(1.0);
    }
    for col in row.iter().rev() {
        let raw = col
            .split('#')
            .next()
            .unwrap_or("")
            .split_whitespace()
            .next()
            .unwrap_or("");
        if raw.is_empty() {
            continue;
        }
        if let Ok(v) = raw.parse::<f32>() {
            return Some(v);
        }
    }
    Some(1.0)
}

fn parse_cofactors(path: &Path) -> Result<Vec<CofactorEntry>> {
    let mmap = MmapFile::open(path)?;
    let mut rdr = TsvReader::new(mmap.as_bytes());
    let mut row = Vec::new();
    let mut out = Vec::new();
    let mut header_seen = false;

    while let Some(rec) = rdr.next_record(&mut row) {
        rec?;
        if row.is_empty() || row.iter().all(|c| c.trim().is_empty()) {
            continue;
        }
        if row[0].trim_start().starts_with('#') {
            continue;
        }
        if !header_seen {
            header_seen = true;
            if row != ["complex_id", "role", "subunit_symbol", "required", "logic"] {
                return Err(KiraError::new(
                    ErrorKind::TsvHeader,
                    format!("cofactors header mismatch in {}", path.display()),
                ));
            }
            continue;
        }
        if row.len() != 5 {
            return Err(KiraError::new(
                ErrorKind::TsvParse,
                format!("invalid cofactors row in {}", path.display()),
            ));
        }
        let required = match row[3].as_str() {
            "1" => true,
            "0" => false,
            _ => {
                return Err(KiraError::new(
                    ErrorKind::TsvParse,
                    format!("invalid required value in {}", path.display()),
                ));
            }
        };
        let logic = match row[4].as_str() {
            "AND" => Logic::And,
            "OR" => Logic::Or,
            _ => {
                return Err(KiraError::new(
                    ErrorKind::TsvParse,
                    format!("invalid logic value in {}", path.display()),
                ));
            }
        };
        out.push(CofactorEntry {
            role: row[1].clone(),
            complex_id: row[0].clone(),
            subunit_symbol: row[2].clone(),
            required,
            logic,
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

fn parse_labels(path: &Path) -> Result<Vec<LabelRule>> {
    let mmap = MmapFile::open(path)?;
    let mut rdr = TsvReader::new(mmap.as_bytes());
    let mut row = Vec::new();
    let mut out = Vec::new();
    let mut header_seen = false;

    while let Some(rec) = rdr.next_record(&mut row) {
        rec?;
        if row.is_empty() || row.iter().all(|c| c.trim().is_empty()) {
            continue;
        }
        if row[0].trim_start().starts_with('#') {
            continue;
        }
        if !header_seen {
            header_seen = true;
            let ok =
                row == [
                    "ligand_symbol",
                    "receptor_symbol",
                    "label",
                    "weight",
                    "notes",
                ] || row == ["ligand_symbol", "receptor_symbol", "label", "weight"]
                    || row == ["ligand_symbol", "receptor_symbol", "label", "notes"]
                    || row == ["ligand_symbol", "receptor_symbol", "label"];
            if !ok {
                return Err(KiraError::new(
                    ErrorKind::TsvHeader,
                    format!("labels header mismatch in {}", path.display()),
                ));
            }
            continue;
        }

        if row.len() < 3 {
            return Err(KiraError::new(
                ErrorKind::TsvParse,
                format!(
                    "invalid labels row in {} at line {}",
                    path.display(),
                    rdr.line_no()
                ),
            ));
        }
        let weight = if row.len() >= 4 {
            let raw = row[3]
                .split('#')
                .next()
                .unwrap_or("")
                .split_whitespace()
                .next()
                .unwrap_or("1.0");
            raw.parse().map_err(|e| {
                KiraError::new(
                    ErrorKind::TsvParse,
                    format!("invalid label weight in {}: {e}", path.display()),
                )
            })?
        } else {
            1.0
        };
        out.push(LabelRule {
            ligand_pat: row[0].clone(),
            receptor_pat: row[1].clone(),
            label: row[2].clone(),
            weight,
        });
    }

    // stable deterministic order for rule evaluation
    out.sort_by(|a, b| {
        a.ligand_pat
            .cmp(&b.ligand_pat)
            .then(a.receptor_pat.cmp(&b.receptor_pat))
            .then(a.label.cmp(&b.label))
    });

    if !header_seen {
        return Err(KiraError::new(
            ErrorKind::TsvHeader,
            format!("missing header in {}", path.display()),
        ));
    }

    Ok(out)
}

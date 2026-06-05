use crate::agg::reader::AggReader;
use crate::error::{ErrorKind, KiraError, Result};
use crate::io::atomic::{write_json_atomic, write_tsv_atomic};
use crate::io::format::fmt_f32;
use crate::io::mmap::MmapFile;
use crate::io::tsv::TsvReader;
use crate::resources::resolve_lr_path;
use ahash::AHashMap;
use rayon::prelude::*;
use serde::Serialize;
use std::cmp::Reverse;
use std::collections::{BTreeMap, BinaryHeap};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering as AtomicOrdering};

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

/// Compact edge representation. String names live in shared `Vec<String>`
/// pools; the edge carries only indices and floats (≤ 40 bytes vs the 200+
/// bytes of the old String-heavy struct). The `flags` field is short and
/// retained as a String because labels are dynamic.
#[derive(Debug, Clone)]
struct Edge {
    source_idx: u32,
    target_idx: u32,
    pair_idx: u32,
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

    // Build sorted group ordering (by name) but keep agg-internal indices
    // for fast lookups. `groups_sorted[i]` is the i-th group in name order;
    // `group_agg_idx[i]` is its index into the AggReader cache.
    let agg_groups = agg.groups();
    let mut order: Vec<usize> = (0..agg_groups.len()).collect();
    order.sort_by(|&a, &b| agg_groups[a].cmp(&agg_groups[b]));
    let groups_sorted: Vec<String> = order.iter().map(|&i| agg_groups[i].clone()).collect();
    let group_agg_idx: Vec<usize> = order;

    let gene_idx_by_name: AHashMap<&str, usize> = agg
        .genes()
        .iter()
        .enumerate()
        .map(|(i, g)| (g.as_str(), i))
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
    let n_groups = groups_sorted.len();

    // Parallel scan over LR pairs. Each pair is independent: its bucket keys
    // are `(pair_idx, src_idx)` which are unique to that pair, so we can merge
    // the per-pair AHashMaps without conflict. `edges_before_filter` is an
    // exact counter and accumulates atomically.
    let edges_before_filter_atomic = AtomicUsize::new(0);

    let per_pair_buckets: Vec<Vec<(u64, TopK)>> = resolved_pairs
        .par_iter()
        .enumerate()
        .map(|(pair_idx, rp)| {
            // Hoist per-source ligand and per-target receptor evaluations.
            let l_data: Vec<(f32, f32, f32)> = group_agg_idx
                .iter()
                .map(|&aix| eval_entity(&rp.ligand_entity, aix, &agg))
                .collect();
            let r_data: Vec<(f32, f32, f32)> = group_agg_idx
                .iter()
                .map(|&aix| eval_entity(&rp.receptor_entity, aix, &agg))
                .collect();

            let mut local_buckets: AHashMap<u64, TopK> = AHashMap::new();
            let mut local_before: usize = 0;

            for src in 0..n_groups {
                let (l_expr, cov_l, mean_l) = l_data[src];
                if cov_l < cfg.cov_min || l_expr < cfg.expr_min {
                    local_before += n_groups;
                    continue;
                }

                for tgt in 0..n_groups {
                    let (r_expr, cov_r, mean_r) = r_data[tgt];
                    local_before += 1;

                    if cov_r < cfg.cov_min || r_expr < cfg.expr_min {
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

                    let flags_str = build_flags_string(
                        cov_l,
                        cov_r,
                        cfg.spec_on,
                        spec_boost,
                        &rp.label_flags,
                    );

                    let edge = Edge {
                        source_idx: src as u32,
                        target_idx: tgt as u32,
                        pair_idx: pair_idx as u32,
                        score: final_score,
                        l_expr,
                        r_expr,
                        cov_l,
                        cov_r,
                        spec_l,
                        spec_r,
                        flags: flags_str,
                    };

                    let key = pack_u32_pair(pair_idx as u32, src as u32);
                    let topk = local_buckets
                        .entry(key)
                        .or_insert_with(|| TopK::new(cfg.top_n_per_pair));
                    topk.push(edge);
                }
            }

            edges_before_filter_atomic.fetch_add(local_before, AtomicOrdering::Relaxed);
            local_buckets.into_iter().collect()
        })
        .collect();

    let mut by_pair_source: AHashMap<u64, TopK> =
        AHashMap::with_capacity(per_pair_buckets.iter().map(|v| v.len()).sum());
    for bucket_list in per_pair_buckets {
        for (k, v) in bucket_list {
            // Per-pair keys are disjoint across pairs (pair_idx is in the key),
            // so insert is unconditional and never collides.
            by_pair_source.insert(k, v);
        }
    }
    let edges_before_filter = edges_before_filter_atomic.load(AtomicOrdering::Relaxed);

    // Drain (pair, source) buckets, sort each by score-desc with tie-break by
    // target name; then group by source_idx and apply top_n_per_source.
    let mut by_source: AHashMap<u32, TopK> = AHashMap::new();
    let target_cmp = |a: &Edge, b: &Edge| {
        b.score
            .total_cmp(&a.score)
            .then_with(|| groups_sorted[a.target_idx as usize].cmp(
                &groups_sorted[b.target_idx as usize],
            ))
    };
    for (_, topk) in by_pair_source {
        for edge in topk.into_sorted_vec(target_cmp) {
            let src = edge.source_idx;
            let topk2 = by_source
                .entry(src)
                .or_insert_with(|| TopK::new(cfg.top_n_per_source));
            topk2.push(edge);
        }
    }

    // Flatten with final global sort.
    let mut final_edges: Vec<Edge> = Vec::new();
    let per_source_cmp = |a: &Edge, b: &Edge| {
        b.score
            .total_cmp(&a.score)
            .then_with(|| resolved_pairs[a.pair_idx as usize].ligand.cmp(
                &resolved_pairs[b.pair_idx as usize].ligand,
            ))
            .then_with(|| resolved_pairs[a.pair_idx as usize].receptor.cmp(
                &resolved_pairs[b.pair_idx as usize].receptor,
            ))
            .then_with(|| groups_sorted[a.target_idx as usize].cmp(
                &groups_sorted[b.target_idx as usize],
            ))
    };
    for (_, topk) in by_source {
        final_edges.extend(topk.into_sorted_vec(per_source_cmp));
    }
    final_edges.sort_by(|a, b| {
        b.score
            .total_cmp(&a.score)
            .then_with(|| resolved_pairs[a.pair_idx as usize].ligand.cmp(
                &resolved_pairs[b.pair_idx as usize].ligand,
            ))
            .then_with(|| resolved_pairs[a.pair_idx as usize].receptor.cmp(
                &resolved_pairs[b.pair_idx as usize].receptor,
            ))
            .then_with(|| groups_sorted[a.source_idx as usize].cmp(
                &groups_sorted[b.source_idx as usize],
            ))
            .then_with(|| groups_sorted[a.target_idx as usize].cmp(
                &groups_sorted[b.target_idx as usize],
            ))
    });

    write_edges_raw(
        &stage2_dir.join("edges_raw.tsv"),
        &final_edges,
        &resolved_pairs,
        &groups_sorted,
    )?;
    write_pairs_stats(
        &stage2_dir.join("pairs_stats.tsv"),
        &final_edges,
        &resolved_pairs,
    )?;

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

#[inline]
fn pack_u32_pair(hi: u32, lo: u32) -> u64 {
    ((hi as u64) << 32) | (lo as u64)
}

/// Deterministic, allocation-light flag-string builder.
/// Built-in flags (`LOW_COVERAGE_EDGE`, `HIGH_SPECIFICITY`) are inserted at
/// their sorted position. Label flags are already deduplicated in stage2 setup
/// and sorted lexicographically here.
fn build_flags_string(
    cov_l: f32,
    cov_r: f32,
    spec_on: bool,
    spec_boost: f32,
    label_flags: &[String],
) -> String {
    let low_cov = cov_l.min(cov_r) < COV_WARN;
    let high_spec = spec_on && spec_boost >= 3.0;

    // Pre-size: each flag is short; estimate.
    let mut flags: Vec<&str> = Vec::with_capacity(2 + label_flags.len());
    if high_spec {
        flags.push("HIGH_SPECIFICITY");
    }
    if low_cov {
        flags.push("LOW_COVERAGE_EDGE");
    }
    for lf in label_flags {
        flags.push(lf.as_str());
    }
    flags.sort();
    flags.dedup();
    flags.join(",")
}

/// Bounded top-K accumulator for `Edge`, scored by `Edge::score` (total order
/// via `f32::total_cmp`). Implemented as a min-heap on the negated score: the
/// element currently retained with the SMALLEST score sits at the heap root,
/// so when capacity is reached we compare a new candidate against root.score
/// and only insert if better — O(log K) per insert, O(N log K) overall.
///
/// Drain via `into_sorted_vec(cmp)` applies the final ordering (with tie-breaks)
/// once at the end.
struct TopK {
    cap: usize,
    heap: BinaryHeap<HeapNode>,
}

#[derive(Clone)]
struct HeapNode {
    // Reverse-ordered so BinaryHeap (max-heap) yields the smallest score at top.
    score: Reverse<F32Total>,
    edge: Edge,
}

#[derive(Clone, Copy)]
struct F32Total(f32);

impl PartialEq for F32Total {
    fn eq(&self, other: &Self) -> bool {
        self.0.total_cmp(&other.0).is_eq()
    }
}
impl Eq for F32Total {}
impl PartialOrd for F32Total {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for F32Total {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.0.total_cmp(&other.0)
    }
}

impl PartialEq for HeapNode {
    fn eq(&self, other: &Self) -> bool {
        self.score == other.score
    }
}
impl Eq for HeapNode {}
impl PartialOrd for HeapNode {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for HeapNode {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.score.cmp(&other.score)
    }
}

impl TopK {
    fn new(cap: usize) -> Self {
        Self {
            cap,
            heap: BinaryHeap::with_capacity(cap.saturating_add(1)),
        }
    }

    fn push(&mut self, edge: Edge) {
        if self.cap == 0 {
            return;
        }
        if self.heap.len() < self.cap {
            self.heap.push(HeapNode {
                score: Reverse(F32Total(edge.score)),
                edge,
            });
            return;
        }
        // Replace root if the new candidate is strictly better.
        if let Some(top) = self.heap.peek()
            && edge.score.total_cmp(&top.score.0.0).is_gt()
        {
            self.heap.pop();
            self.heap.push(HeapNode {
                score: Reverse(F32Total(edge.score)),
                edge,
            });
        }
    }

    fn into_sorted_vec<F>(self, cmp: F) -> Vec<Edge>
    where
        F: FnMut(&Edge, &Edge) -> std::cmp::Ordering,
    {
        let mut v: Vec<Edge> = self.heap.into_iter().map(|n| n.edge).collect();
        v.sort_by(cmp);
        v
    }
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
    gene_idx_by_name: &AHashMap<&str, usize>,
    gene_stats: &BTreeMap<String, GeneStat>,
    missing_components: &mut usize,
    missing_genes: &mut usize,
) -> Option<Entity> {
    if let Some(entries) = complex_map.get(&(role.to_string(), symbol.to_string())) {
        let mut components = Vec::new();
        for c in entries {
            let Some(&gene_idx) = gene_idx_by_name.get(c.subunit_symbol.as_str()) else {
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

fn write_edges_raw(
    path: &Path,
    edges: &[Edge],
    resolved_pairs: &[ResolvedPair],
    groups_sorted: &[String],
) -> Result<()> {
    let mut lines = Vec::with_capacity(edges.len() + 1);
    lines.push(
        "source_group\ttarget_group\tligand\treceptor\tscore\tL_expr\tR_expr\tcov_L\tcov_R\tspec_L\tspec_R\tflags"
            .to_string(),
    );
    for e in edges {
        let rp = &resolved_pairs[e.pair_idx as usize];
        lines.push(format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            groups_sorted[e.source_idx as usize],
            groups_sorted[e.target_idx as usize],
            rp.ligand,
            rp.receptor,
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

fn write_pairs_stats(
    path: &Path,
    edges: &[Edge],
    resolved_pairs: &[ResolvedPair],
) -> Result<()> {
    let mut by_pair: BTreeMap<(String, String), (usize, f32)> = BTreeMap::new();
    for e in edges {
        let rp = &resolved_pairs[e.pair_idx as usize];
        by_pair
            .entry((rp.ligand.clone(), rp.receptor.clone()))
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

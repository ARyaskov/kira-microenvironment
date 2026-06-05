use crate::auto_groups::AutoGroupsMode;
use crate::auto_groups::anti::AntiMarkerGene;
use crate::auto_groups::hierarchy::allowed_fine_groups;
use crate::auto_groups::markers::MarkerGroup;
use crate::error::Result;
use crate::expr::reader::ExprReader;
use rayon::prelude::*;
use std::cell::RefCell;
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone)]
pub struct AssignConfig {
    pub mode: AutoGroupsMode,
    pub eps: f32,
    pub min_delta: f32,
    pub unknown_label: String,
}

#[derive(Debug, Clone)]
pub struct AssignRow {
    pub cell_id: String,
    pub group: String,
}

#[derive(Debug, Clone)]
pub struct CellScoreRow {
    pub cell_id: String,
    pub group_name: String,
    pub score: f32,
}

#[derive(Debug, Clone, Default)]
pub struct AssignCounts {
    pub n_unknown: usize,
    pub n_coarse_only: usize,
    pub n_fine: usize,
}

#[derive(Debug, Clone)]
pub struct AssignResult {
    pub rows: Vec<AssignRow>,
    pub cell_scores: Vec<CellScoreRow>,
    pub counts: AssignCounts,
    pub per_group: BTreeMap<String, usize>,
}

/// Per-cell output produced inside the parallel map; aggregated afterwards
/// so the final `AssignResult` ordering matches the input barcode order.
enum CellOutcome {
    /// Hierarchical with neither coarse nor fine → unknown.
    Unknown,
    /// Flat or hierarchical with a coarse-only fallback / fine assignment.
    /// `is_fine` indicates whether the fine bucket (hierarchical) was hit.
    Assigned { group: String, is_fine: bool },
}

pub fn assign_cells(
    expr: &ExprReader,
    barcodes: &[String],
    marker_groups: &[MarkerGroup],
    coarse_groups: Option<&[MarkerGroup]>,
    anti_markers: &BTreeMap<String, Vec<AntiMarkerGene>>,
    cfg: &AssignConfig,
) -> Result<AssignResult> {
    let n_genes = expr.n_genes();

    // Per-thread scratch: (dense_values, touched_indices, scored_vec, per_cell_score_rows_local).
    thread_local! {
        static SCRATCH: RefCell<(Vec<f32>, Vec<u32>, Vec<(String, f32)>)> =
            const { RefCell::new((Vec::new(), Vec::new(), Vec::new())) };
    }

    // Parallel per-cell map. Each item returns (CellOutcome, score rows for this cell).
    let per_cell: Vec<(CellOutcome, Vec<CellScoreRow>)> = barcodes
        .par_iter()
        .enumerate()
        .map(|(cell_i, cell_id)| {
            SCRATCH.with(|cell| {
                let mut scratch = cell.borrow_mut();
                let (values, touched, scored) = &mut *scratch;
                if values.len() != n_genes {
                    values.clear();
                    values.resize(n_genes, 0.0);
                    touched.clear();
                } else {
                    for &g in touched.iter() {
                        values[g as usize] = 0.0;
                    }
                    touched.clear();
                }

                if cell_i < expr.n_cells() {
                    for (gidx, v) in expr.iter_cell(cell_i as u32) {
                        let gi = gidx as usize;
                        if gi < n_genes {
                            values[gi] = v;
                            touched.push(gidx);
                        }
                    }
                }

                let mut local_scores: Vec<CellScoreRow> = Vec::new();

                let outcome = match cfg.mode {
                    AutoGroupsMode::Flat => {
                        compute_scores(marker_groups, None, values, anti_markers, scored);
                        record_scores(cell_id, scored, &mut local_scores);
                        match pick_group(scored, cfg) {
                            Some(group) => CellOutcome::Assigned { group, is_fine: true },
                            None => CellOutcome::Assigned {
                                group: cfg.unknown_label.clone(),
                                is_fine: false,
                            },
                        }
                    }
                    AutoGroupsMode::Hierarchical => {
                        let coarse = coarse_groups.unwrap_or(&[]);
                        compute_scores(coarse, None, values, anti_markers, scored);
                        record_scores(cell_id, scored, &mut local_scores);
                        let Some(coarse_group) = pick_group(scored, cfg) else {
                            return (CellOutcome::Unknown, local_scores);
                        };

                        let allowed: BTreeSet<&str> =
                            allowed_fine_groups(&coarse_group).iter().copied().collect();
                        compute_scores(
                            marker_groups,
                            Some(&allowed),
                            values,
                            anti_markers,
                            scored,
                        );
                        record_scores(cell_id, scored, &mut local_scores);

                        if let Some(fine_group) = pick_group(scored, cfg) {
                            CellOutcome::Assigned {
                                group: fine_group,
                                is_fine: true,
                            }
                        } else {
                            CellOutcome::Assigned {
                                group: coarse_group,
                                is_fine: false,
                            }
                        }
                    }
                };
                (outcome, local_scores)
            })
        })
        .collect();

    // Sequential aggregation (cheap) preserves deterministic ordering.
    let mut rows = Vec::with_capacity(barcodes.len());
    let mut cell_scores = Vec::with_capacity(barcodes.len() * marker_groups.len());
    let mut counts = AssignCounts::default();
    let mut per_group: BTreeMap<String, usize> = BTreeMap::new();

    for ((outcome, local_scores), cell_id) in per_cell.into_iter().zip(barcodes.iter()) {
        cell_scores.extend(local_scores);
        let group = match outcome {
            CellOutcome::Unknown => {
                counts.n_unknown += 1;
                cfg.unknown_label.clone()
            }
            CellOutcome::Assigned { group, is_fine } => {
                if group == cfg.unknown_label {
                    counts.n_unknown += 1;
                } else if is_fine {
                    // Flat-mode "assigned" is always is_fine=true (set in map).
                    // Hierarchical fine match also sets is_fine=true.
                    counts.n_fine += 1;
                } else {
                    // Hierarchical fallback to coarse.
                    counts.n_coarse_only += 1;
                }
                group
            }
        };
        *per_group.entry(group.clone()).or_insert(0) += 1;
        rows.push(AssignRow {
            cell_id: cell_id.clone(),
            group,
        });
    }

    Ok(AssignResult {
        rows,
        cell_scores,
        counts,
        per_group,
    })
}

fn compute_scores(
    groups: &[MarkerGroup],
    allowed: Option<&BTreeSet<&str>>,
    values: &[f32],
    anti_markers: &BTreeMap<String, Vec<AntiMarkerGene>>,
    out: &mut Vec<(String, f32)>,
) {
    out.clear();
    let n_genes = values.len();
    for group in groups {
        if let Some(allowed_groups) = allowed
            && !allowed_groups.contains(group.name.as_str())
        {
            continue;
        }

        let mut sum = 0.0f32;
        for m in &group.genes {
            if let Some(idx) = m.gene_idx {
                let gi = idx as usize;
                if gi < n_genes {
                    sum += values[gi] * m.weight;
                }
            }
        }
        let mut score = sum / (group.genes.len().max(1) as f32);
        if let Some(anti) = anti_markers.get(&group.name) {
            let mut penalty = 0.0f32;
            for a in anti {
                if let Some(idx) = a.gene_idx {
                    let gi = idx as usize;
                    if gi < n_genes {
                        penalty += values[gi] * a.penalty;
                    }
                }
            }
            penalty /= anti.len().max(1) as f32;
            score = (score - penalty).max(0.0);
        }
        out.push((group.name.clone(), score));
    }
}

fn record_scores(cell_id: &str, scored: &[(String, f32)], out: &mut Vec<CellScoreRow>) {
    for (group_name, score) in scored {
        out.push(CellScoreRow {
            cell_id: cell_id.to_string(),
            group_name: group_name.clone(),
            score: *score,
        });
    }
}

fn pick_group(scored: &mut [(String, f32)], cfg: &AssignConfig) -> Option<String> {
    if scored.is_empty() {
        return None;
    }
    scored.sort_by(|a, b| b.1.total_cmp(&a.1).then(a.0.cmp(&b.0)));

    let best = &scored[0];
    let second_score = scored.get(1).map(|x| x.1).unwrap_or(0.0);
    if best.1 < cfg.eps || (best.1 - second_score) < cfg.min_delta {
        None
    } else {
        Some(best.0.clone())
    }
}

use crate::auto_groups::AutoGroupsMode;
use crate::auto_groups::anti::AntiMarkerGene;
use crate::auto_groups::hierarchy::allowed_fine_groups;
use crate::auto_groups::markers::MarkerGroup;
use crate::error::Result;
use crate::expr::reader::ExprReader;
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

pub fn assign_cells(
    expr: &ExprReader,
    barcodes: &[String],
    marker_groups: &[MarkerGroup],
    coarse_groups: Option<&[MarkerGroup]>,
    anti_markers: &BTreeMap<String, Vec<AntiMarkerGene>>,
    cfg: &AssignConfig,
) -> Result<AssignResult> {
    let mut rows = Vec::with_capacity(barcodes.len());
    let mut cell_scores = Vec::new();
    let mut counts = AssignCounts::default();
    let mut per_group: BTreeMap<String, usize> = BTreeMap::new();

    for (cell_i, cell_id) in barcodes.iter().enumerate() {
        let mut values: BTreeMap<u32, f32> = BTreeMap::new();
        if cell_i < expr.n_cells() {
            for (gidx, v) in expr.iter_cell(cell_i as u32) {
                values.insert(gidx, v);
            }
        }

        let assigned = match cfg.mode {
            AutoGroupsMode::Flat => {
                let scored = compute_scores(marker_groups, None, &values, anti_markers);
                record_scores(cell_id, &scored, &mut cell_scores);
                pick_group(&scored, cfg).unwrap_or_else(|| cfg.unknown_label.clone())
            }
            AutoGroupsMode::Hierarchical => {
                let coarse = coarse_groups.unwrap_or(&[]);
                let coarse_scores = compute_scores(coarse, None, &values, anti_markers);
                record_scores(cell_id, &coarse_scores, &mut cell_scores);
                let Some(coarse_group) = pick_group(&coarse_scores, cfg) else {
                    counts.n_unknown += 1;
                    rows.push(AssignRow {
                        cell_id: cell_id.clone(),
                        group: cfg.unknown_label.clone(),
                    });
                    *per_group.entry(cfg.unknown_label.clone()).or_insert(0) += 1;
                    continue;
                };

                let allowed: BTreeSet<String> = allowed_fine_groups(&coarse_group)
                    .iter()
                    .map(|x| x.to_string())
                    .collect();
                let fine_scores =
                    compute_scores(marker_groups, Some(&allowed), &values, anti_markers);
                record_scores(cell_id, &fine_scores, &mut cell_scores);

                if let Some(fine_group) = pick_group(&fine_scores, cfg) {
                    counts.n_fine += 1;
                    fine_group
                } else {
                    counts.n_coarse_only += 1;
                    coarse_group
                }
            }
        };

        if assigned == cfg.unknown_label {
            counts.n_unknown += 1;
        } else if matches!(cfg.mode, AutoGroupsMode::Flat) {
            counts.n_fine += 1;
        }
        *per_group.entry(assigned.clone()).or_insert(0) += 1;
        rows.push(AssignRow {
            cell_id: cell_id.clone(),
            group: assigned,
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
    allowed: Option<&BTreeSet<String>>,
    values: &BTreeMap<u32, f32>,
    anti_markers: &BTreeMap<String, Vec<AntiMarkerGene>>,
) -> Vec<(String, f32)> {
    let mut scored = Vec::new();
    for group in groups {
        if let Some(allowed_groups) = allowed
            && !allowed_groups.contains(&group.name)
        {
            continue;
        }

        let mut sum = 0.0f32;
        for m in &group.genes {
            if let Some(idx) = m.gene_idx {
                sum += values.get(&idx).copied().unwrap_or(0.0) * m.weight;
            }
        }
        let mut score = sum / (group.genes.len().max(1) as f32);
        if let Some(anti) = anti_markers.get(&group.name) {
            let mut penalty = 0.0f32;
            for a in anti {
                if let Some(idx) = a.gene_idx {
                    penalty += values.get(&idx).copied().unwrap_or(0.0) * a.penalty;
                }
            }
            penalty /= anti.len().max(1) as f32;
            score = (score - penalty).max(0.0);
        }
        scored.push((group.name.clone(), score));
    }
    scored
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

fn pick_group(scored: &[(String, f32)], cfg: &AssignConfig) -> Option<String> {
    if scored.is_empty() {
        return None;
    }
    let mut ordered = scored.to_vec();
    ordered.sort_by(|a, b| b.1.total_cmp(&a.1).then(a.0.cmp(&b.0)));

    let best = &ordered[0];
    let second_score = ordered.get(1).map(|x| x.1).unwrap_or(0.0);
    if best.1 < cfg.eps || (best.1 - second_score) < cfg.min_delta {
        None
    } else {
        Some(best.0.clone())
    }
}

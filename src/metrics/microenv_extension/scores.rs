use crate::expr::reader::ExprReader;
use crate::metrics::microenv_extension::panels::{
    MIN_PANEL_GENES, PANEL_TRIM_FRACTION, PanelKind, panel_specs,
};
use crate::select::{median_in_place, trimmed_mean_in_place};
use serde::Serialize;
use std::collections::BTreeSet;

const ROBUST_Z_EPS: f32 = 1e-6;

pub const HSI_HYPOXIA_HIGH: f32 = 2.0;
pub const IAS_INFLAMMATORY_HIGH: f32 = 2.0;
pub const ISS_IMMUNE_SUPPRESSION_HIGH: f32 = 1.5;
pub const MIO_METABOLIC_SUPPRESSION_HIGH: f32 = 1.5;
pub const SII_STROMAL_HIGH: f32 = 2.0;
pub const MSM_STRESS_MODE_HIGH: f32 = 2.0;

#[derive(Debug, Clone, Serialize)]
pub struct MicroenvThresholds {
    pub hypoxia_high: f32,
    pub inflammatory_high: f32,
    pub immune_suppression_high: f32,
    pub metabolic_suppression_high: f32,
    pub stromal_high: f32,
    pub microenv_stress_mode: f32,
}

impl Default for MicroenvThresholds {
    fn default() -> Self {
        Self {
            hypoxia_high: HSI_HYPOXIA_HIGH,
            inflammatory_high: IAS_INFLAMMATORY_HIGH,
            immune_suppression_high: ISS_IMMUNE_SUPPRESSION_HIGH,
            metabolic_suppression_high: MIO_METABOLIC_SUPPRESSION_HIGH,
            stromal_high: SII_STROMAL_HIGH,
            microenv_stress_mode: MSM_STRESS_MODE_HIGH,
        }
    }
}

#[derive(Debug, Clone, Serialize)]
pub struct PanelCoverage {
    pub panel: String,
    pub genes_total: usize,
    pub genes_found: usize,
    pub found_genes: Vec<String>,
    pub missing_genes: Vec<String>,
    pub available: bool,
}

#[derive(Debug, Clone)]
pub struct ResolvedPanel {
    pub kind: PanelKind,
    pub name: &'static str,
    pub slots: Vec<usize>,
    pub found_gene_count: usize,
}

#[derive(Debug, Clone)]
pub struct ResolvedPanels {
    pub panels: Vec<ResolvedPanel>,
    pub gene_to_slot: Vec<i32>,
    pub n_union_slots: usize,
    pub coverage: Vec<PanelCoverage>,
}

#[derive(Debug, Clone)]
pub struct CellMicroenvScores {
    pub hyp_core: Vec<f32>,
    pub nfkb_core: Vec<f32>,
    pub ifn_core: Vec<f32>,
    pub checkpoint_core: Vec<f32>,
    pub adenosine_core: Vec<f32>,
    pub stromal_core: Vec<f32>,
    pub hsi: Vec<f32>,
    pub ias: Vec<f32>,
    pub iss: Vec<f32>,
    pub mio: Vec<f32>,
    pub sii: Vec<f32>,
    pub msm: Vec<f32>,
    pub hypoxia_high: Vec<bool>,
    pub inflammatory_high: Vec<bool>,
    pub immune_suppression_high: Vec<bool>,
    pub metabolic_suppression_high: Vec<bool>,
    pub stromal_high: Vec<bool>,
    pub microenv_stress_mode: Vec<bool>,
}

impl CellMicroenvScores {
    pub fn n_cells(&self) -> usize {
        self.hsi.len()
    }
}

pub fn resolve_panels(expr: &ExprReader) -> ResolvedPanels {
    let specs = panel_specs();
    let mut gene_to_slot = vec![-1i32; expr.n_genes()];
    let mut next_slot = 0usize;
    let mut panels = Vec::with_capacity(specs.len());
    let mut coverage = Vec::with_capacity(specs.len());

    for spec in specs {
        let mut slots = Vec::with_capacity(spec.genes.len());
        let mut found_queries = Vec::new();
        let mut missing_queries = Vec::new();
        let mut seen_gene_idx = BTreeSet::new();
        for &gene in spec.genes {
            if let Some(gidx) = expr.gene_index(gene) {
                if seen_gene_idx.insert(gidx) {
                    let gi = gidx as usize;
                    if gene_to_slot[gi] < 0 {
                        gene_to_slot[gi] = next_slot as i32;
                        next_slot += 1;
                    }
                    slots.push(gene_to_slot[gi] as usize);
                    found_queries.push(gene.to_string());
                }
            } else {
                missing_queries.push(gene.to_string());
            }
        }

        let found_gene_count = slots.len();
        coverage.push(PanelCoverage {
            panel: spec.name.to_string(),
            genes_total: spec.genes.len(),
            genes_found: found_gene_count,
            found_genes: found_queries,
            missing_genes: missing_queries,
            available: found_gene_count >= MIN_PANEL_GENES,
        });
        panels.push(ResolvedPanel {
            kind: spec.kind,
            name: spec.name,
            slots,
            found_gene_count,
        });
    }

    ResolvedPanels {
        panels,
        gene_to_slot,
        n_union_slots: next_slot,
        coverage,
    }
}

pub fn compute_microenv_scores(expr: &ExprReader, resolved: &ResolvedPanels) -> CellMicroenvScores {
    let n_cells = expr.n_cells();
    let mut hyp_core = vec![f32::NAN; n_cells];
    let mut nfkb_core = vec![f32::NAN; n_cells];
    let mut ifn_core = vec![f32::NAN; n_cells];
    let mut checkpoint_core = vec![f32::NAN; n_cells];
    let mut adenosine_core = vec![f32::NAN; n_cells];
    let mut stromal_core = vec![f32::NAN; n_cells];

    let mut union_values = vec![0.0f32; resolved.n_union_slots];
    let mut panel_buf = Vec::<f32>::new();

    for cell_idx in 0..n_cells {
        union_values.fill(0.0);
        for (gene_idx, value) in expr.iter_cell(cell_idx as u32) {
            let slot = resolved.gene_to_slot[gene_idx as usize];
            if slot >= 0 {
                union_values[slot as usize] = value;
            }
        }

        for panel in &resolved.panels {
            let core = if panel.found_gene_count < MIN_PANEL_GENES {
                f32::NAN
            } else {
                panel_buf.clear();
                for &slot in &panel.slots {
                    panel_buf.push(union_values[slot]);
                }
                trimmed_mean_in_place(&mut panel_buf, PANEL_TRIM_FRACTION)
            };

            match panel.kind {
                PanelKind::Hypoxia => hyp_core[cell_idx] = core,
                PanelKind::Inflammation => nfkb_core[cell_idx] = core,
                PanelKind::Interferon => ifn_core[cell_idx] = core,
                PanelKind::Checkpoint => checkpoint_core[cell_idx] = core,
                PanelKind::Adenosine => adenosine_core[cell_idx] = core,
                PanelKind::Stromal => stromal_core[cell_idx] = core,
            }
        }
    }

    let hsi = robust_zscore(&hyp_core);
    let znf = robust_zscore(&nfkb_core);
    let zif = robust_zscore(&ifn_core);
    let iss = robust_zscore(&checkpoint_core);
    let zad = robust_zscore(&adenosine_core);
    let sii = robust_zscore(&stromal_core);

    let mut ias = vec![f32::NAN; n_cells];
    let mut mio = vec![f32::NAN; n_cells];
    let mut msm = vec![f32::NAN; n_cells];

    let thresholds = MicroenvThresholds::default();
    let mut hypoxia_high = vec![false; n_cells];
    let mut inflammatory_high = vec![false; n_cells];
    let mut immune_suppression_high = vec![false; n_cells];
    let mut metabolic_suppression_high = vec![false; n_cells];
    let mut stromal_high = vec![false; n_cells];
    let mut microenv_stress_mode = vec![false; n_cells];

    for i in 0..n_cells {
        if znf[i].is_finite() && zif[i].is_finite() {
            ias[i] = 0.6 * znf[i] + 0.4 * zif[i];
        }
        if zad[i].is_finite() && hsi[i].is_finite() {
            mio[i] = 0.5 * zad[i] + 0.5 * hsi[i].max(0.0);
        }
        if hsi[i].is_finite() && iss[i].is_finite() && mio[i].is_finite() {
            msm[i] =
                (0.4 * hsi[i].max(0.0) + 0.3 * iss[i].max(0.0) + 0.3 * mio[i].max(0.0)).max(0.0);
        }

        if hsi[i].is_finite() {
            hypoxia_high[i] = hsi[i] >= thresholds.hypoxia_high;
        }
        if ias[i].is_finite() {
            inflammatory_high[i] = ias[i] >= thresholds.inflammatory_high;
        }
        if iss[i].is_finite() {
            immune_suppression_high[i] = iss[i] >= thresholds.immune_suppression_high;
        }
        if mio[i].is_finite() {
            metabolic_suppression_high[i] = mio[i] >= thresholds.metabolic_suppression_high;
        }
        if sii[i].is_finite() {
            stromal_high[i] = sii[i] >= thresholds.stromal_high;
        }
        if msm[i].is_finite() {
            microenv_stress_mode[i] = msm[i] >= thresholds.microenv_stress_mode;
        }
    }

    CellMicroenvScores {
        hyp_core,
        nfkb_core,
        ifn_core,
        checkpoint_core,
        adenosine_core,
        stromal_core,
        hsi,
        ias,
        iss,
        mio,
        sii,
        msm,
        hypoxia_high,
        inflammatory_high,
        immune_suppression_high,
        metabolic_suppression_high,
        stromal_high,
        microenv_stress_mode,
    }
}

fn robust_zscore(raw: &[f32]) -> Vec<f32> {
    let mut finite = raw
        .iter()
        .copied()
        .filter(|x| x.is_finite())
        .collect::<Vec<_>>();
    if finite.is_empty() {
        return vec![f32::NAN; raw.len()];
    }
    let median = median_in_place(&mut finite);
    let mut abs_dev = raw
        .iter()
        .copied()
        .filter(|x| x.is_finite())
        .map(|x| (x - median).abs())
        .collect::<Vec<_>>();
    let mad = median_in_place(&mut abs_dev);

    let mut out = vec![f32::NAN; raw.len()];
    if mad == 0.0 {
        for (i, x) in raw.iter().enumerate() {
            if x.is_finite() {
                out[i] = 0.0;
            }
        }
        return out;
    }

    let denom = 1.4826 * mad + ROBUST_Z_EPS;
    for (i, x) in raw.iter().enumerate() {
        if x.is_finite() {
            out[i] = (*x - median) / denom;
        }
    }
    out
}

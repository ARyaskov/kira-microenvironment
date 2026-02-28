use serde::Serialize;
use std::collections::BTreeMap;

#[derive(Debug, Clone, Serialize)]
pub struct AutoGroupsSummary {
    pub mode: String,
    pub counts: Counts,
    pub per_group: BTreeMap<String, usize>,
    pub anti_markers_used: bool,
}

#[derive(Debug, Clone, Serialize)]
pub struct Counts {
    pub n_cells: usize,
    pub n_unknown: usize,
    pub n_coarse_only: usize,
    pub n_fine: usize,
}

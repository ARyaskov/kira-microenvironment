pub fn allowed_fine_groups(coarse_group: &str) -> &'static [&'static str] {
    match coarse_group {
        "tumor" => &["tumor_epithelial", "Neuroendocrine"],
        "immune" => &[
            "T_cell",
            "T_cell_CD4",
            "T_cell_CD8",
            "T_cell_CD8_cytotoxic",
            "T_cell_CD8_exhausted",
            "Treg",
            "NK",
            "B_cell",
            "Plasma_cell",
            "Myeloid",
            "Monocyte",
            "Macrophage",
            "Macrophage_M2",
            "DC",
            "pDC",
            "Neutrophil",
        ],
        "stromal" => &[
            "Fibroblast",
            "CAF",
            "Endothelial",
            "Lymphatic_endothelial",
            "Pericyte",
            "Smooth_muscle",
        ],
        _ => &[],
    }
}

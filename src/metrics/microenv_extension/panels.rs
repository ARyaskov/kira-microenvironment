pub const MICROENV_EXTENSION_PANEL_V1: &str = "MICROENV_EXTENSION_PANEL_V1";
pub const PANEL_TRIM_FRACTION: f32 = 0.10;
pub const MIN_PANEL_GENES: usize = 2;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum PanelKind {
    Hypoxia,
    Inflammation,
    Interferon,
    Checkpoint,
    Adenosine,
    Stromal,
}

#[derive(Debug, Clone, Copy)]
pub struct PanelSpec {
    pub kind: PanelKind,
    pub name: &'static str,
    pub genes: &'static [&'static str],
}

const HYPOXIA_PANEL: &[&str] = &["HIF1A", "VEGFA", "SLC2A1", "LDHA", "CA9", "BNIP3"];
const INFLAMMATION_PANEL: &[&str] = &["NFKB1", "RELA", "IL1B", "TNF", "CXCL8", "IL8", "CCL2"];
const INTERFERON_PANEL: &[&str] = &["STAT1", "IRF1", "ISG15", "IFIT1", "MX1"];
const CHECKPOINT_PANEL: &[&str] = &["CD274", "PDCD1LG2", "TGFB1", "IDO1", "LAG3"];
const ADENOSINE_PANEL: &[&str] = &["NT5E", "ENTPD1", "ADORA2A"];
const STROMAL_PANEL: &[&str] = &["COL1A1", "COL3A1", "FN1", "SPARC", "ACTA2", "MMP2"];

const PANEL_SPECS: &[PanelSpec] = &[
    PanelSpec {
        kind: PanelKind::Hypoxia,
        name: "hypoxia_hif",
        genes: HYPOXIA_PANEL,
    },
    PanelSpec {
        kind: PanelKind::Inflammation,
        name: "inflammation_nfkb",
        genes: INFLAMMATION_PANEL,
    },
    PanelSpec {
        kind: PanelKind::Interferon,
        name: "interferon_response",
        genes: INTERFERON_PANEL,
    },
    PanelSpec {
        kind: PanelKind::Checkpoint,
        name: "immune_checkpoint_suppression",
        genes: CHECKPOINT_PANEL,
    },
    PanelSpec {
        kind: PanelKind::Adenosine,
        name: "adenosine_pathway",
        genes: ADENOSINE_PANEL,
    },
    PanelSpec {
        kind: PanelKind::Stromal,
        name: "ecm_stromal_interaction",
        genes: STROMAL_PANEL,
    },
];

pub fn panel_specs() -> &'static [PanelSpec] {
    PANEL_SPECS
}

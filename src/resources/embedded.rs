use crate::error::{ErrorKind, KiraError, Result};
use crate::io::atomic::write_bytes_atomic;
use std::path::{Path, PathBuf};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EmbeddedProfile {
    Onco,
    Immune,
}

impl EmbeddedProfile {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Onco => "onco",
            Self::Immune => "immune",
        }
    }
}

const EMBEDDED_LR_PAIRS_ONCO: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/resources/lr_pairs_mvp_onco.tsv"
));
const EMBEDDED_LR_PAIRS_IMMUNE: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/resources/lr_pairs_mvp_immune.tsv"
));
const EMBEDDED_COFACTORS_ONCO: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/resources_onco/cofactors.tsv"
));
const EMBEDDED_COFACTORS_IMMUNE: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/resources_immune/cofactors.tsv"
));
const EMBEDDED_GENE_ALIAS_ONCO: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/resources_onco/gene_alias.tsv"
));
const EMBEDDED_GENE_ALIAS_IMMUNE: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/resources_immune/gene_alias.tsv"
));
const EMBEDDED_LABELS_ONCO: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/resources_onco/labels.tsv"
));
const EMBEDDED_LABELS_IMMUNE: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/resources_immune/labels.tsv"
));
const EMBEDDED_REGIME_MAP_ONCO: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/resources_onco/regime_map.tsv"
));

pub fn materialize_embedded_resources(
    out_root: &Path,
    profile: EmbeddedProfile,
) -> Result<PathBuf> {
    let target = out_root.join(".kira-microenvironment-embedded-resources");
    std::fs::create_dir_all(&target).map_err(|e| {
        KiraError::new(
            ErrorKind::Path,
            format!(
                "failed creating embedded resources dir {}: {e}",
                target.display()
            ),
        )
    })?;

    let (lr_pairs, cofactors, aliases, labels) = match profile {
        EmbeddedProfile::Onco => (
            EMBEDDED_LR_PAIRS_ONCO,
            EMBEDDED_COFACTORS_ONCO,
            EMBEDDED_GENE_ALIAS_ONCO,
            EMBEDDED_LABELS_ONCO,
        ),
        EmbeddedProfile::Immune => (
            EMBEDDED_LR_PAIRS_IMMUNE,
            EMBEDDED_COFACTORS_IMMUNE,
            EMBEDDED_GENE_ALIAS_IMMUNE,
            EMBEDDED_LABELS_IMMUNE,
        ),
    };

    write_resource(&target.join("lr_pairs.tsv"), lr_pairs)?;
    write_resource(&target.join("lr_pairs_mvp.tsv"), lr_pairs)?;
    write_resource(&target.join("cofactors.tsv"), cofactors)?;
    write_resource(&target.join("gene_alias.tsv"), aliases)?;
    write_resource(&target.join("labels.tsv"), labels)?;
    write_resource(&target.join("regime_map.tsv"), EMBEDDED_REGIME_MAP_ONCO)?;
    let profile_marker = format!("{}\n", profile.as_str());
    write_resource(&target.join("embedded_profile.txt"), &profile_marker)?;

    Ok(target)
}

fn write_resource(path: &Path, content: &str) -> Result<()> {
    write_bytes_atomic(path, content.as_bytes()).map_err(|e| {
        KiraError::new(
            ErrorKind::Path,
            format!("failed writing embedded resource {}: {e}", path.display()),
        )
    })
}

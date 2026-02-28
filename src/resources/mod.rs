pub mod aliases;
pub mod cofactors;
pub mod embedded;
pub mod labels;
pub mod lr_pairs;
pub mod regime_map;

use crate::error::{ErrorKind, KiraError, Result};
use crate::paths::{PathMeta, metadata_for_path};
use aliases::{AliasResolver, parse_aliases};
use cofactors::{Cofactor, parse_cofactors};
use labels::validate_labels;
use lr_pairs::{LrPair, parse_lr_pairs};
use std::path::{Path, PathBuf};

#[derive(Debug, Clone)]
pub struct LoadedResources {
    pub lr_pairs: Vec<LrPair>,
    pub cofactors: Vec<Cofactor>,
    pub aliases: AliasResolver,
    pub metadata: Vec<PathMeta>,
}

#[derive(Debug, Clone)]
pub struct ResourcePaths {
    pub lr_pairs: PathBuf,
    pub cofactors: PathBuf,
    pub aliases: PathBuf,
    pub labels: PathBuf,
}

pub fn resolve_lr_path(resources_dir: &Path, lr_profile: &str) -> Result<PathBuf> {
    match lr_profile {
        "mvp" => {
            let p = resources_dir.join("lr_pairs_mvp.tsv");
            if !p.exists() {
                return Err(KiraError::new(
                    ErrorKind::ResourcesIncomplete,
                    format!("lr profile mvp missing {}", p.display()),
                ));
            }
            Ok(p)
        }
        "full" => Ok(resources_dir.join("lr_pairs.tsv")),
        other => {
            let raw = PathBuf::from(other);
            if raw.is_absolute() {
                Ok(raw)
            } else {
                Ok(std::env::current_dir()
                    .map_err(|e| {
                        KiraError::new(ErrorKind::Path, format!("current_dir failed: {e}"))
                    })?
                    .join(raw))
            }
        }
    }
}

pub fn load_resources(paths: &ResourcePaths) -> Result<LoadedResources> {
    for (name, p) in [
        ("lr_pairs", &paths.lr_pairs),
        ("cofactors", &paths.cofactors),
        ("gene_alias", &paths.aliases),
        ("labels", &paths.labels),
    ] {
        if !p.exists() {
            return Err(KiraError::new(
                ErrorKind::ResourcesIncomplete,
                format!("missing resource {name}: {}", p.display()),
            ));
        }
    }

    let lr_pairs = parse_lr_pairs(&paths.lr_pairs)?;
    let cofactors = parse_cofactors(&paths.cofactors)?;
    let aliases = AliasResolver::from_entries(parse_aliases(&paths.aliases)?);
    validate_labels(&paths.labels)?;

    let metadata = vec![
        metadata_for_path("lr_pairs", &paths.lr_pairs)?,
        metadata_for_path("cofactors", &paths.cofactors)?,
        metadata_for_path("gene_alias", &paths.aliases)?,
        metadata_for_path("labels", &paths.labels)?,
    ];

    Ok(LoadedResources {
        lr_pairs,
        cofactors,
        aliases,
        metadata,
    })
}

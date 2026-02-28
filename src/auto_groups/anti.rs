use crate::error::{ErrorKind, KiraError, Result};
use crate::expr::reader::ExprReader;
use crate::resources::aliases::AliasResolver;
use std::collections::BTreeMap;
use std::path::Path;

#[derive(Debug, Clone)]
pub struct AntiMarkerGene {
    pub penalty: f32,
    pub gene_idx: Option<u32>,
}

pub fn load_anti_markers(
    path: &Path,
    expr: &ExprReader,
    aliases: &AliasResolver,
) -> Result<BTreeMap<String, Vec<AntiMarkerGene>>> {
    let content = std::fs::read_to_string(path).map_err(|e| {
        KiraError::new(
            ErrorKind::AutoGroupsAntiParse,
            format!("failed reading {}: {e}", path.display()),
        )
    })?;

    let mut out: BTreeMap<String, Vec<AntiMarkerGene>> = BTreeMap::new();
    let mut header_checked = false;
    for (i, line) in content.lines().enumerate() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = trimmed.split('\t').collect();
        if !header_checked {
            header_checked = true;
            if cols == ["group_name", "gene_symbol", "penalty"] {
                continue;
            }
        }
        if cols.len() != 3 {
            return Err(KiraError::new(
                ErrorKind::AutoGroupsAntiParse,
                format!("invalid anti-marker line {} in {}", i + 1, path.display()),
            ));
        }
        let penalty = cols[2].parse::<f32>().map_err(|e| {
            KiraError::new(
                ErrorKind::AutoGroupsAntiParse,
                format!("invalid penalty line {}: {e}", i + 1),
            )
        })?;
        let symbol = aliases.resolve(cols[1]).to_string();
        out.entry(cols[0].to_string())
            .or_default()
            .push(AntiMarkerGene {
                penalty,
                gene_idx: expr.gene_index(&symbol),
            });
    }

    for genes in out.values_mut() {
        genes.sort_by(|a, b| a.penalty.total_cmp(&b.penalty));
    }
    Ok(out)
}

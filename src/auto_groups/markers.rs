use crate::error::{ErrorKind, KiraError, Result};
use crate::expr::reader::ExprReader;
use crate::resources::aliases::AliasResolver;
use std::collections::BTreeMap;
use std::path::Path;

#[derive(Debug, Clone)]
pub struct MarkerGene {
    pub symbol: String,
    pub weight: f32,
    pub gene_idx: Option<u32>,
}

#[derive(Debug, Clone)]
pub struct MarkerGroup {
    pub name: String,
    pub genes: Vec<MarkerGene>,
}

pub fn load_marker_groups(
    path: &Path,
    expr: &ExprReader,
    aliases: &AliasResolver,
) -> Result<Vec<MarkerGroup>> {
    let content = std::fs::read_to_string(path).map_err(|e| {
        KiraError::new(
            ErrorKind::AutoGroupsParseError,
            format!("failed reading {}: {e}", path.display()),
        )
    })?;

    let mut by_group: BTreeMap<String, Vec<(String, f32)>> = BTreeMap::new();
    for (i, line) in content.lines().enumerate() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let line_wo_inline_comment = if let Some((left, _)) = trimmed.split_once("\t#") {
            left
        } else {
            trimmed
        };
        let cols: Vec<&str> = line_wo_inline_comment.split('\t').collect();
        if i == 0 && cols.len() >= 2 && cols[0] == "group_name" && cols[1] == "gene_symbol" {
            continue;
        }
        if cols.len() < 2 || cols.len() > 3 {
            return Err(KiraError::new(
                ErrorKind::AutoGroupsParseError,
                format!("invalid marker line {} in {}", i + 1, path.display()),
            ));
        }
        let weight = if cols.len() == 3 {
            let raw = cols[2]
                .split('#')
                .next()
                .unwrap_or("")
                .split_whitespace()
                .next()
                .unwrap_or("");
            raw.parse::<f32>().map_err(|e| {
                KiraError::new(
                    ErrorKind::AutoGroupsParseError,
                    format!("invalid marker weight line {}: {e}", i + 1),
                )
            })?
        } else {
            1.0
        };

        by_group
            .entry(cols[0].to_string())
            .or_default()
            .push((cols[1].to_string(), weight));
    }

    let mut groups = Vec::new();
    for (group, mut markers) in by_group {
        markers.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.total_cmp(&b.1)));
        let genes = markers
            .into_iter()
            .map(|(sym, w)| {
                let resolved = aliases.resolve(&sym).to_string();
                MarkerGene {
                    symbol: resolved.clone(),
                    weight: w,
                    gene_idx: expr.gene_index(&resolved),
                }
            })
            .collect();
        groups.push(MarkerGroup { name: group, genes });
    }
    groups.sort_by(|a, b| a.name.cmp(&b.name));

    Ok(groups)
}

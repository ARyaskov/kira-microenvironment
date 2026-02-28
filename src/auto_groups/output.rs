use crate::auto_groups::assign::{AssignRow, CellScoreRow};
use crate::auto_groups::summary::AutoGroupsSummary;
use crate::error::Result;
use crate::io::atomic::{write_json_atomic, write_tsv_atomic};
use crate::io::format::fmt_f32;
use std::path::Path;

pub fn write_groups(path: &Path, rows: &[AssignRow]) -> Result<()> {
    let mut lines = Vec::with_capacity(rows.len() + 1);
    lines.push("cell_id\tgroup".to_string());
    for r in rows {
        lines.push(format!("{}\t{}", r.cell_id, r.group));
    }
    write_tsv_atomic(path, &lines)
}

pub fn write_scores(path: &Path, rows: &[CellScoreRow]) -> Result<()> {
    let mut lines = Vec::with_capacity(rows.len() + 1);
    lines.push("cell_id\tgroup_name\tscore".to_string());
    for r in rows {
        lines.push(format!(
            "{}\t{}\t{}",
            r.cell_id,
            r.group_name,
            fmt_f32(r.score)
        ));
    }
    write_tsv_atomic(path, &lines)
}

pub fn write_summary(path: &Path, summary: &AutoGroupsSummary) -> Result<()> {
    write_json_atomic(path, summary)
}

use crate::error::{ErrorKind, KiraError, Result};
use crate::io::mmap::MmapFile;
use crate::io::tsv::{TsvReader, expect_header_exact};
use std::path::Path;

#[derive(Debug, Clone)]
pub struct LrPair {
    pub ligand_symbol: String,
    pub receptor_symbol: String,
}

pub fn parse_lr_pairs(path: &Path) -> Result<Vec<LrPair>> {
    let mmap = MmapFile::open(path)?;
    let mut rdr = TsvReader::new(mmap.as_bytes());
    let mut row = Vec::new();

    let mut out = Vec::new();
    let mut header_done = false;
    while let Some(rec) = rdr.next_record(&mut row) {
        rec?;
        if row.is_empty() || row.iter().all(|c| c.trim().is_empty()) {
            continue;
        }
        if row[0].trim_start().starts_with('#') {
            continue;
        }
        if !header_done {
            if row
                == [
                    "ligand_symbol",
                    "receptor_symbol",
                    "family",
                    "directionality",
                ]
                || row
                    == [
                        "ligand_symbol",
                        "receptor_symbol",
                        "family",
                        "directionality",
                        "weight",
                    ]
                || row
                    == [
                        "ligand_symbol",
                        "receptor_symbol",
                        "family",
                        "directionality",
                        "weight",
                        "notes",
                    ]
                || row
                    == [
                        "ligand_symbol",
                        "receptor_symbol",
                        "family",
                        "directionality",
                        "notes",
                    ]
            {
                header_done = true;
                continue;
            }
            expect_header_exact(
                &row,
                &[
                    "ligand_symbol",
                    "receptor_symbol",
                    "family",
                    "directionality",
                ],
                &path.display().to_string(),
            )?;
        }
        if row.len() < 4 {
            return Err(KiraError::new(
                ErrorKind::TsvParse,
                format!("{} line {} has <4 columns", path.display(), rdr.line_no()),
            ));
        }
        out.push(LrPair {
            ligand_symbol: row[0].clone(),
            receptor_symbol: row[1].clone(),
        });
    }

    if !header_done {
        return Err(KiraError::new(
            ErrorKind::TsvHeader,
            format!("missing header in {}", path.display()),
        ));
    }

    Ok(out)
}

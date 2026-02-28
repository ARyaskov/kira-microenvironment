use crate::error::{ErrorKind, KiraError, Result};
use crate::io::mmap::MmapFile;
use crate::io::tsv::TsvReader;
use std::path::Path;

pub fn validate_labels(path: &Path) -> Result<()> {
    let mmap = MmapFile::open(path)?;
    let mut rdr = TsvReader::new(mmap.as_bytes());
    let mut row = Vec::new();
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
                != [
                    "ligand_symbol",
                    "receptor_symbol",
                    "label",
                    "weight",
                    "notes",
                ]
                && row != ["ligand_symbol", "receptor_symbol", "label", "weight"]
                && row != ["ligand_symbol", "receptor_symbol", "label", "notes"]
                && row != ["ligand_symbol", "receptor_symbol", "label"]
            {
                return Err(KiraError::new(
                    ErrorKind::TsvHeader,
                    format!("header mismatch in {}", path.display()),
                ));
            }
            header_done = true;
            continue;
        }
        if row.len() < 3 {
            return Err(KiraError::new(
                ErrorKind::TsvParse,
                format!("{} line {} has <3 columns", path.display(), rdr.line_no()),
            ));
        }
    }

    if !header_done {
        return Err(KiraError::new(
            ErrorKind::TsvHeader,
            format!("missing header in {}", path.display()),
        ));
    }

    Ok(())
}

use crate::error::{ErrorKind, KiraError, Result};
use crate::io::mmap::MmapFile;
use crate::io::tsv::TsvReader;
use std::path::Path;

#[derive(Debug, Clone)]
pub struct Cofactor {
    pub complex_id: String,
    pub role: String,
    pub subunit_symbol: String,
}

pub fn parse_cofactors(path: &Path) -> Result<Vec<Cofactor>> {
    let mmap = MmapFile::open(path)?;
    let mut rdr = TsvReader::new(mmap.as_bytes());
    let mut row = Vec::new();
    let mut header_done = false;
    let mut out = Vec::new();

    while let Some(rec) = rdr.next_record(&mut row) {
        rec?;
        if row.is_empty() {
            continue;
        }
        if !header_done {
            if row != ["complex_id", "role", "subunit_symbol", "required", "logic"] {
                return Err(KiraError::new(
                    ErrorKind::TsvHeader,
                    format!("header mismatch in {}", path.display()),
                ));
            }
            header_done = true;
            continue;
        }
        if row.len() != 5 {
            return Err(KiraError::new(
                ErrorKind::TsvParse,
                format!(
                    "{} line {} has {} columns",
                    path.display(),
                    rdr.line_no(),
                    row.len()
                ),
            ));
        }
        if row[1] != "ligand" && row[1] != "receptor" {
            return Err(KiraError::new(
                ErrorKind::TsvParse,
                format!(
                    "{} line {} invalid role {}",
                    path.display(),
                    rdr.line_no(),
                    row[1]
                ),
            ));
        }
        if row[3] != "0" && row[3] != "1" {
            return Err(KiraError::new(
                ErrorKind::TsvParse,
                format!(
                    "{} line {} invalid required {}",
                    path.display(),
                    rdr.line_no(),
                    row[3]
                ),
            ));
        }
        if row[4] != "AND" && row[4] != "OR" {
            return Err(KiraError::new(
                ErrorKind::TsvParse,
                format!(
                    "{} line {} invalid logic {}",
                    path.display(),
                    rdr.line_no(),
                    row[4]
                ),
            ));
        }

        out.push(Cofactor {
            complex_id: row[0].clone(),
            role: row[1].clone(),
            subunit_symbol: row[2].clone(),
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

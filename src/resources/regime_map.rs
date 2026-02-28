use crate::error::{ErrorKind, KiraError, Result};
use crate::io::mmap::MmapFile;
use crate::io::tsv::TsvReader;
use std::collections::BTreeMap;
use std::path::Path;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RegimeClass {
    Loud,
    Silent,
    Mixed,
    Ignore,
}

pub fn parse_regime_map(path: &Path) -> Result<BTreeMap<String, RegimeClass>> {
    if !path.exists() {
        return Err(KiraError::new(
            ErrorKind::RegimeMapMissing,
            format!("missing regime map: {}", path.display()),
        ));
    }

    let mmap = MmapFile::open(path)?;
    let mut rdr = TsvReader::new(mmap.as_bytes());
    let mut row = Vec::new();
    let mut header_seen = false;
    let mut out = BTreeMap::new();

    while let Some(rec) = rdr.next_record(&mut row) {
        rec?;
        if row.is_empty() || row.iter().all(|c| c.trim().is_empty()) {
            continue;
        }
        if row[0].trim_start().starts_with('#') {
            continue;
        }
        if !header_seen {
            header_seen = true;
            if row != ["regime", "class", "notes"] && row != ["regime", "class"] {
                return Err(KiraError::new(
                    ErrorKind::RegimeMapParse,
                    format!("regime_map header mismatch in {}", path.display()),
                ));
            }
            continue;
        }
        if row.len() < 2 {
            return Err(KiraError::new(
                ErrorKind::RegimeMapParse,
                format!("invalid regime_map row in {}", path.display()),
            ));
        }

        let class = match row[1].as_str() {
            "LOUD" => RegimeClass::Loud,
            "SILENT" => RegimeClass::Silent,
            "MIXED" => RegimeClass::Mixed,
            "IGNORE" => RegimeClass::Ignore,
            other => {
                return Err(KiraError::new(
                    ErrorKind::RegimeMapParse,
                    format!("unknown regime class '{other}' in {}", path.display()),
                ));
            }
        };
        out.insert(row[0].clone(), class);
    }

    if !header_seen {
        return Err(KiraError::new(
            ErrorKind::RegimeMapParse,
            format!("missing header in {}", path.display()),
        ));
    }

    Ok(out)
}

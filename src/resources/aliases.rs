use crate::error::{ErrorKind, KiraError, Result};
use crate::io::mmap::MmapFile;
use crate::io::tsv::TsvReader;
use std::collections::BTreeMap;
use std::path::Path;

#[derive(Debug, Clone)]
pub struct AliasEntry {
    pub alias: String,
    pub symbol: String,
    pub priority: i32,
}

#[derive(Debug, Clone, Default)]
pub struct AliasResolver {
    best: BTreeMap<String, AliasEntry>,
}

fn normalize_alias_key(raw: &str) -> String {
    raw.trim().to_ascii_uppercase()
}

impl AliasResolver {
    pub fn from_entries(entries: Vec<AliasEntry>) -> Self {
        let mut best: BTreeMap<String, AliasEntry> = BTreeMap::new();
        for entry in entries {
            let key = normalize_alias_key(&entry.alias);
            best.entry(key)
                .and_modify(|cur| {
                    if entry.priority < cur.priority
                        || (entry.priority == cur.priority && entry.symbol < cur.symbol)
                    {
                        *cur = entry.clone();
                    }
                })
                .or_insert(entry);
        }
        Self { best }
    }

    pub fn resolve<'a>(&'a self, symbol_or_alias: &'a str) -> &'a str {
        let key = normalize_alias_key(symbol_or_alias);
        self.best
            .get(&key)
            .map(|x| x.symbol.as_str())
            .unwrap_or(symbol_or_alias)
    }
}

pub fn parse_aliases(path: &Path) -> Result<Vec<AliasEntry>> {
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
        if row[0].trim_start().starts_with('#') {
            continue;
        }
        if !header_done {
            if row != ["alias", "symbol", "priority", "notes"]
                && row != ["alias", "symbol", "priority"]
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
        let priority: i32 = row[2].parse().map_err(|e| {
            KiraError::new(
                ErrorKind::TsvParse,
                format!(
                    "{} line {} invalid priority: {e}",
                    path.display(),
                    rdr.line_no()
                ),
            )
        })?;
        out.push(AliasEntry {
            alias: row[0].clone(),
            symbol: row[1].clone(),
            priority,
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

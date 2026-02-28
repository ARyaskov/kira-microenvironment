use crate::error::{ErrorKind, KiraError, Result};
use crate::io::mmap::MmapFile;
use std::path::Path;

const CACHE_MAGIC: [u8; 8] = *b"KIRAEAGG";
const CACHE_VERSION: u32 = 1;

pub struct AggReader {
    _mmap: MmapFile,
    groups: Vec<String>,
    genes: Vec<String>,
    expr: Vec<f32>,
    cov: Vec<f32>,
}

impl AggReader {
    pub fn open(path: &Path) -> Result<Self> {
        let mmap = MmapFile::open(path)?;
        let data = mmap.as_bytes();
        if data.len() < 20 {
            return Err(KiraError::new(
                ErrorKind::TsvParse,
                format!("cache too small: {}", path.display()),
            ));
        }

        if data[0..8] != CACHE_MAGIC {
            return Err(KiraError::new(
                ErrorKind::TsvParse,
                format!("invalid cache magic: {}", path.display()),
            ));
        }
        let version = read_u32(data, 8)?;
        if version != CACHE_VERSION {
            return Err(KiraError::new(
                ErrorKind::TsvParse,
                format!("unsupported cache version {version}: {}", path.display()),
            ));
        }

        let n_groups = read_u32(data, 12)? as usize;
        let n_genes = read_u32(data, 16)? as usize;

        let mut pos = 20usize;
        let mut groups = Vec::with_capacity(n_groups);
        for _ in 0..n_groups {
            let len = read_u32(data, pos)? as usize;
            pos += 4;
            let end = pos.checked_add(len).ok_or_else(|| {
                KiraError::new(ErrorKind::TsvParse, "cache groups length overflow")
            })?;
            if end > data.len() {
                return Err(KiraError::new(
                    ErrorKind::TsvParse,
                    "cache groups out of bounds",
                ));
            }
            let s = std::str::from_utf8(&data[pos..end])
                .map_err(|e| KiraError::new(ErrorKind::TsvParse, format!("utf8 error: {e}")))?;
            groups.push(s.to_string());
            pos = end;
        }

        let mut genes = Vec::with_capacity(n_genes);
        for _ in 0..n_genes {
            let len = read_u32(data, pos)? as usize;
            pos += 4;
            let end = pos.checked_add(len).ok_or_else(|| {
                KiraError::new(ErrorKind::TsvParse, "cache genes length overflow")
            })?;
            if end > data.len() {
                return Err(KiraError::new(
                    ErrorKind::TsvParse,
                    "cache genes out of bounds",
                ));
            }
            let s = std::str::from_utf8(&data[pos..end])
                .map_err(|e| KiraError::new(ErrorKind::TsvParse, format!("utf8 error: {e}")))?;
            genes.push(s.to_string());
            pos = end;
        }

        let n = n_groups
            .checked_mul(n_genes)
            .ok_or_else(|| KiraError::new(ErrorKind::TsvParse, "cache matrix size overflow"))?;

        let expr_bytes = n
            .checked_mul(4)
            .ok_or_else(|| KiraError::new(ErrorKind::TsvParse, "cache expr bytes overflow"))?;
        let cov_bytes = expr_bytes;
        if pos + expr_bytes + cov_bytes > data.len() {
            return Err(KiraError::new(
                ErrorKind::TsvParse,
                "cache matrix arrays out of bounds",
            ));
        }

        let mut expr = Vec::with_capacity(n);
        for _ in 0..n {
            expr.push(read_f32(data, pos)?);
            pos += 4;
        }

        let mut cov = Vec::with_capacity(n);
        for _ in 0..n {
            cov.push(read_f32(data, pos)?);
            pos += 4;
        }

        Ok(Self {
            _mmap: mmap,
            groups,
            genes,
            expr,
            cov,
        })
    }

    pub fn n_groups(&self) -> usize {
        self.groups.len()
    }

    pub fn n_genes(&self) -> usize {
        self.genes.len()
    }

    pub fn groups(&self) -> &[String] {
        &self.groups
    }

    pub fn genes(&self) -> &[String] {
        &self.genes
    }

    pub fn expr(&self, group_idx: usize, gene_idx: usize) -> f32 {
        let idx = group_idx
            .checked_mul(self.n_genes())
            .and_then(|x| x.checked_add(gene_idx));
        match idx {
            Some(i) if i < self.expr.len() => self.expr[i],
            _ => 0.0,
        }
    }

    pub fn cov(&self, group_idx: usize, gene_idx: usize) -> f32 {
        let idx = group_idx
            .checked_mul(self.n_genes())
            .and_then(|x| x.checked_add(gene_idx));
        match idx {
            Some(i) if i < self.cov.len() => self.cov[i],
            _ => 0.0,
        }
    }
}

fn read_u32(data: &[u8], offset: usize) -> Result<u32> {
    let end = offset + 4;
    if end > data.len() {
        return Err(KiraError::new(
            ErrorKind::TsvParse,
            "read_u32 out of bounds",
        ));
    }
    Ok(u32::from_le_bytes(data[offset..end].try_into().unwrap()))
}

fn read_f32(data: &[u8], offset: usize) -> Result<f32> {
    let end = offset + 4;
    if end > data.len() {
        return Err(KiraError::new(
            ErrorKind::TsvParse,
            "read_f32 out of bounds",
        ));
    }
    Ok(f32::from_le_bytes(data[offset..end].try_into().unwrap()))
}

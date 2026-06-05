use crate::error::{ErrorKind, KiraError, Result};
use crate::io::mmap::MmapFile;
use std::path::Path;

const CACHE_MAGIC: [u8; 8] = *b"KIRAEAGG";
const CACHE_VERSION: u32 = 1;

/// Reader for KIRAEAGG v1 (`stage1_agg/cache.bin`).
///
/// Holds the mmap alive and exposes the `expr`/`cov` matrices either as
/// borrowed slices of f32 (when alignment permits) or as owned Vec<f32>
/// (fallback). Lookups are O(1) by direct slice indexing.
pub struct AggReader {
    // Field drop order matters: storages may borrow from `_mmap`, so the mmap
    // MUST be dropped last. In Rust, fields drop in declaration order, so keep
    // `_mmap` at the bottom.
    groups: Vec<String>,
    genes: Vec<String>,
    expr_storage: MatrixStorage,
    cov_storage: MatrixStorage,
    n_genes: usize,
    _mmap: MmapFile,
}

enum MatrixStorage {
    /// Slice into the mmap (zero-copy fast path).
    Borrowed(&'static [f32]),
    /// Owned copy (fallback when alignment doesn't match f32).
    Owned(Vec<f32>),
}

impl MatrixStorage {
    #[inline]
    fn as_slice(&self) -> &[f32] {
        match self {
            Self::Borrowed(s) => s,
            Self::Owned(v) => v.as_slice(),
        }
    }
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

        let bytes = n
            .checked_mul(4)
            .ok_or_else(|| KiraError::new(ErrorKind::TsvParse, "cache expr bytes overflow"))?;
        if pos + bytes * 2 > data.len() {
            return Err(KiraError::new(
                ErrorKind::TsvParse,
                "cache matrix arrays out of bounds",
            ));
        }

        let expr_bytes = &data[pos..pos + bytes];
        let cov_bytes = &data[pos + bytes..pos + bytes * 2];

        // SAFETY: We extend the slice lifetime to 'static; this is sound only
        // because `_mmap` is stored in the same struct and is never dropped
        // before the slices are dropped (struct field drop order: declaration order).
        // Note: declaration order matters — `_mmap` is below the storages, so
        // they are dropped first. To be safe we keep `_mmap` LAST in field order.
        let expr_storage = make_storage(expr_bytes, n)?;
        let cov_storage = make_storage(cov_bytes, n)?;

        Ok(Self {
            groups,
            genes,
            expr_storage,
            cov_storage,
            n_genes,
            _mmap: mmap,
        })
    }

    pub fn n_groups(&self) -> usize {
        self.groups.len()
    }

    pub fn n_genes(&self) -> usize {
        self.n_genes
    }

    pub fn groups(&self) -> &[String] {
        &self.groups
    }

    pub fn genes(&self) -> &[String] {
        &self.genes
    }

    #[inline]
    pub fn expr(&self, group_idx: usize, gene_idx: usize) -> f32 {
        let s = self.expr_storage.as_slice();
        let n = self.n_genes;
        if gene_idx >= n {
            return 0.0;
        }
        let idx = group_idx * n + gene_idx;
        if idx < s.len() { s[idx] } else { 0.0 }
    }

    #[inline]
    pub fn cov(&self, group_idx: usize, gene_idx: usize) -> f32 {
        let s = self.cov_storage.as_slice();
        let n = self.n_genes;
        if gene_idx >= n {
            return 0.0;
        }
        let idx = group_idx * n + gene_idx;
        if idx < s.len() { s[idx] } else { 0.0 }
    }
}

fn make_storage(bytes: &[u8], n: usize) -> Result<MatrixStorage> {
    debug_assert_eq!(bytes.len(), n * 4);
    // Fast path: if mmap is f32-aligned, reinterpret directly.
    if bytes.as_ptr() as usize % std::mem::align_of::<f32>() == 0 {
        let slice: &[f32] = bytemuck::cast_slice(bytes);
        // SAFETY: see AggReader::open — `_mmap` outlives this slice via field order.
        let extended: &'static [f32] = unsafe { std::mem::transmute(slice) };
        Ok(MatrixStorage::Borrowed(extended))
    } else {
        // Fallback path: allocate and bulk-copy via from_le_bytes per chunk.
        let mut out = Vec::with_capacity(n);
        for chunk in bytes.chunks_exact(4) {
            out.push(f32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]));
        }
        Ok(MatrixStorage::Owned(out))
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

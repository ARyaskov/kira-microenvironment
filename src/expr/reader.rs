use crate::error::{ErrorKind, KiraError, Result};
use crate::expr::format::{
    EXPR_HEADER_LEN, EXPR_MAGIC, EXPR_VERSION, MTX_HEADER_LEN, MTX_MAGIC, MTX_VERSION,
};
use crate::io::input::read_tsv_input;
use crate::io::mmap::MmapFile;
use crate::io::tsv::TsvReader;
use std::collections::BTreeMap;
use std::path::Path;

pub struct ExprReader {
    _mmap: MmapFile,
    n_cells: u32,
    n_genes: u32,
    gene_symbols: Vec<String>,
    gene_lookup: BTreeMap<String, u32>,
    cell_offsets: Vec<u64>,
    gene_idx: Vec<u32>,
    values: Vec<f32>,
}

fn normalize_gene_key(raw: &str) -> String {
    let trimmed = raw.trim();
    let no_version = if trimmed.starts_with("ENS") {
        trimmed.split('.').next().unwrap_or(trimmed)
    } else {
        trimmed
    };
    no_version.to_ascii_uppercase()
}

fn add_gene_lookup_key(lookup: &mut BTreeMap<String, u32>, key: &str, idx: u32) {
    lookup.entry(key.to_string()).or_insert(idx);
    let normalized = normalize_gene_key(key);
    lookup.entry(normalized).or_insert(idx);
}

impl ExprReader {
    pub fn open(path: &Path) -> Result<Self> {
        let mmap = MmapFile::open(path)?;
        let data = mmap.as_bytes();
        if data.len() < 8 {
            return Err(KiraError::new(
                ErrorKind::TsvParse,
                format!("expr.bin too small: {}", path.display()),
            ));
        }

        let magic: [u8; 8] = data[0..8].try_into().map_err(|_| {
            KiraError::new(
                ErrorKind::TsvParse,
                format!("invalid expr magic in {}", path.display()),
            )
        })?;
        if magic == EXPR_MAGIC {
            return Self::open_kiraexpr(path, mmap);
        }
        if magic == MTX_MAGIC {
            return Self::open_kiramtx(path, mmap);
        }
        Err(KiraError::new(
            ErrorKind::TsvParse,
            format!("unexpected expr magic in {}", path.display()),
        ))
    }

    fn open_kiraexpr(path: &Path, mmap: MmapFile) -> Result<Self> {
        let data = mmap.as_bytes();
        if data.len() < EXPR_HEADER_LEN {
            return Err(KiraError::new(
                ErrorKind::TsvParse,
                format!("expr.bin too small: {}", path.display()),
            ));
        }

        let version = read_u32(data, 8)?;
        if version != EXPR_VERSION {
            return Err(KiraError::new(
                ErrorKind::TsvParse,
                format!("unsupported expr version {version} in {}", path.display()),
            ));
        }

        let n_cells = read_u32(data, 12)?;
        let n_genes = read_u32(data, 16)?;
        let gene_symbols_offset = read_u64(data, 20)? as usize;
        let cell_index_offset = read_u64(data, 28)? as usize;
        let nnz_offset = read_u64(data, 36)? as usize;

        if !(gene_symbols_offset <= cell_index_offset
            && cell_index_offset <= nnz_offset
            && nnz_offset <= data.len())
        {
            return Self::open_kiraexpr_csc_by_gene(path, mmap);
        }

        let mut pos = gene_symbols_offset;
        let mut gene_symbols = Vec::with_capacity(n_genes as usize);
        for _ in 0..n_genes {
            let len = read_u32(data, pos)? as usize;
            pos += 4;
            let end = pos + len;
            if end > cell_index_offset {
                return Err(KiraError::new(
                    ErrorKind::TsvParse,
                    format!("gene symbol block overflow in {}", path.display()),
                ));
            }
            let s = std::str::from_utf8(&data[pos..end]).map_err(|e| {
                KiraError::new(
                    ErrorKind::TsvParse,
                    format!("gene symbol utf8 error in {}: {e}", path.display()),
                )
            })?;
            gene_symbols.push(s.to_string());
            pos = end;
        }

        if cell_index_offset + ((n_cells as usize + 1) * 8) > nnz_offset {
            return Err(KiraError::new(
                ErrorKind::TsvParse,
                format!("cell index block overflow in {}", path.display()),
            ));
        }

        let mut cell_offsets = Vec::with_capacity(n_cells as usize + 1);
        let mut idx_pos = cell_index_offset;
        for _ in 0..=n_cells {
            let off = read_u64(data, idx_pos)?;
            cell_offsets.push(off);
            idx_pos += 8;
        }

        let nnz = *cell_offsets.last().unwrap_or(&0) as usize;
        if !cell_offsets.windows(2).all(|w| w[0] <= w[1]) {
            return Err(KiraError::new(
                ErrorKind::TsvParse,
                format!("non-monotonic cell offsets in {}", path.display()),
            ));
        }

        let nnz_gene_bytes = nnz.checked_mul(4).ok_or_else(|| {
            KiraError::new(
                ErrorKind::TsvParse,
                format!("nnz overflow in {}", path.display()),
            )
        })?;
        let nnz_val_bytes = nnz_gene_bytes;
        if nnz_offset + nnz_gene_bytes + nnz_val_bytes > data.len() {
            return Err(KiraError::new(
                ErrorKind::TsvParse,
                format!("nnz arrays overflow in {}", path.display()),
            ));
        }

        let mut gene_idx = Vec::with_capacity(nnz);
        let mut p = nnz_offset;
        for _ in 0..nnz {
            let g = read_u32(data, p)?;
            if g >= n_genes {
                return Err(KiraError::new(
                    ErrorKind::TsvParse,
                    format!("gene_idx out of range in {}", path.display()),
                ));
            }
            gene_idx.push(g);
            p += 4;
        }

        let mut values = Vec::with_capacity(nnz);
        for _ in 0..nnz {
            let v = read_f32(data, p)?;
            values.push(v);
            p += 4;
        }

        let mut gene_lookup = BTreeMap::new();
        for (idx, symbol) in gene_symbols.iter().enumerate() {
            add_gene_lookup_key(&mut gene_lookup, symbol, idx as u32);
        }

        Ok(Self {
            _mmap: mmap,
            n_cells,
            n_genes,
            gene_symbols,
            gene_lookup,
            cell_offsets,
            gene_idx,
            values,
        })
    }

    fn open_kiraexpr_csc_by_gene(path: &Path, mmap: MmapFile) -> Result<Self> {
        let data = mmap.as_bytes();
        const HEADER_LEN: usize = 32;
        if data.len() < HEADER_LEN {
            return Err(KiraError::new(
                ErrorKind::TsvParse,
                format!("expr header too small in {}", path.display()),
            ));
        }

        let n_genes = read_u32(data, 12)? as usize;
        let n_cells = read_u32(data, 16)? as usize;
        let nnz = read_u64(data, 20)? as usize;
        let layout = read_u32(data, 28)?;
        if layout != 1 {
            return Err(KiraError::new(
                ErrorKind::TsvParse,
                format!("unsupported expr layout {layout} in {}", path.display()),
            ));
        }

        let gene_ptr_len = n_genes.checked_add(1).ok_or_else(|| {
            KiraError::new(
                ErrorKind::TsvParse,
                format!("gene_ptr length overflow in {}", path.display()),
            )
        })?;
        let gene_ptr_bytes = gene_ptr_len.checked_mul(8).ok_or_else(|| {
            KiraError::new(
                ErrorKind::TsvParse,
                format!("gene_ptr bytes overflow in {}", path.display()),
            )
        })?;
        let cell_idx_bytes = nnz.checked_mul(4).ok_or_else(|| {
            KiraError::new(
                ErrorKind::TsvParse,
                format!("cell_idx bytes overflow in {}", path.display()),
            )
        })?;
        let values_bytes = nnz.checked_mul(4).ok_or_else(|| {
            KiraError::new(
                ErrorKind::TsvParse,
                format!("values bytes overflow in {}", path.display()),
            )
        })?;
        let expected_len = HEADER_LEN
            .checked_add(gene_ptr_bytes)
            .and_then(|x| x.checked_add(cell_idx_bytes))
            .and_then(|x| x.checked_add(values_bytes))
            .ok_or_else(|| {
                KiraError::new(
                    ErrorKind::TsvParse,
                    format!("expr size overflow in {}", path.display()),
                )
            })?;
        if expected_len != data.len() {
            return Err(KiraError::new(
                ErrorKind::TsvParse,
                format!(
                    "expr size mismatch in {}: expected {}, got {}",
                    path.display(),
                    expected_len,
                    data.len()
                ),
            ));
        }

        let mut gene_ptr = Vec::with_capacity(gene_ptr_len);
        let mut pos = HEADER_LEN;
        for _ in 0..gene_ptr_len {
            gene_ptr.push(read_u64(data, pos)? as usize);
            pos += 8;
        }
        if !gene_ptr.windows(2).all(|w| w[0] <= w[1]) {
            return Err(KiraError::new(
                ErrorKind::TsvParse,
                format!("non-monotonic gene_ptr in {}", path.display()),
            ));
        }
        if gene_ptr.last().copied().unwrap_or(0) != nnz {
            return Err(KiraError::new(
                ErrorKind::TsvParse,
                format!("gene_ptr end != nnz in {}", path.display()),
            ));
        }

        let mut cell_idx = Vec::with_capacity(nnz);
        for _ in 0..nnz {
            let c = read_u32(data, pos)?;
            if c as usize >= n_cells {
                return Err(KiraError::new(
                    ErrorKind::TsvParse,
                    format!("cell_idx out of range in {}", path.display()),
                ));
            }
            cell_idx.push(c as usize);
            pos += 4;
        }

        let mut src_values = Vec::with_capacity(nnz);
        for _ in 0..nnz {
            src_values.push(read_f32(data, pos)?);
            pos += 4;
        }

        let mut per_cell_counts = vec![0usize; n_cells];
        for &c in &cell_idx {
            per_cell_counts[c] = per_cell_counts[c].saturating_add(1);
        }

        let mut cell_offsets: Vec<u64> = Vec::with_capacity(n_cells + 1);
        cell_offsets.push(0);
        for &count in &per_cell_counts {
            let next = cell_offsets
                .last()
                .copied()
                .unwrap_or(0)
                .saturating_add(count as u64);
            cell_offsets.push(next);
        }

        let mut fill = cell_offsets[..n_cells]
            .iter()
            .map(|x| *x as usize)
            .collect::<Vec<_>>();
        let mut gene_idx = vec![0u32; nnz];
        let mut values = vec![0f32; nnz];
        for gene in 0..n_genes {
            let start = gene_ptr[gene];
            let end = gene_ptr[gene + 1];
            for idx in start..end {
                let c = cell_idx[idx];
                let out = fill[c];
                if out >= nnz {
                    return Err(KiraError::new(
                        ErrorKind::TsvParse,
                        format!("cell fill overflow in {}", path.display()),
                    ));
                }
                gene_idx[out] = gene as u32;
                values[out] = src_values[idx];
                fill[c] = out + 1;
            }
        }

        let gene_symbols = load_gene_symbols_for_kiramtx(path, n_genes)?;
        let mut gene_lookup = BTreeMap::new();
        for (idx, symbol) in gene_symbols.iter().enumerate() {
            add_gene_lookup_key(&mut gene_lookup, symbol, idx as u32);
        }

        Ok(Self {
            _mmap: mmap,
            n_cells: n_cells as u32,
            n_genes: n_genes as u32,
            gene_symbols,
            gene_lookup,
            cell_offsets,
            gene_idx,
            values,
        })
    }

    fn open_kiramtx(path: &Path, mmap: MmapFile) -> Result<Self> {
        let data = mmap.as_bytes();
        if data.len() < MTX_HEADER_LEN {
            return Err(KiraError::new(
                ErrorKind::TsvParse,
                format!("KIRAMTX header too small in {}", path.display()),
            ));
        }

        let version = read_u32(data, 8)?;
        if version != MTX_VERSION {
            return Err(KiraError::new(
                ErrorKind::TsvParse,
                format!(
                    "unsupported KIRAMTX version {version} in {}",
                    path.display()
                ),
            ));
        }

        let n_genes = read_u32(data, 12)?;
        let n_cells = read_u32(data, 16)?;
        let _flags = read_u32(data, 20)?;

        let dense_len = (n_genes as usize)
            .checked_mul(n_cells as usize)
            .ok_or_else(|| {
                KiraError::new(
                    ErrorKind::TsvParse,
                    format!("dense size overflow in {}", path.display()),
                )
            })?;
        let expected_bytes = MTX_HEADER_LEN
            .checked_add(dense_len.checked_mul(4).ok_or_else(|| {
                KiraError::new(
                    ErrorKind::TsvParse,
                    format!("dense byte size overflow in {}", path.display()),
                )
            })?)
            .ok_or_else(|| {
                KiraError::new(
                    ErrorKind::TsvParse,
                    format!("dense total size overflow in {}", path.display()),
                )
            })?;
        if data.len() < expected_bytes {
            return Err(KiraError::new(
                ErrorKind::TsvParse,
                format!(
                    "KIRAMTX truncated in {}: expected at least {expected_bytes} bytes, got {}",
                    path.display(),
                    data.len()
                ),
            ));
        }

        let symbols = load_gene_symbols_for_kiramtx(path, n_genes as usize)?;
        let mut gene_lookup = BTreeMap::new();
        for (idx, symbol) in symbols.iter().enumerate() {
            add_gene_lookup_key(&mut gene_lookup, symbol, idx as u32);
        }

        let mut cell_offsets = Vec::with_capacity(n_cells as usize + 1);
        let mut gene_idx = Vec::new();
        let mut values = Vec::new();
        cell_offsets.push(0);

        for cell in 0..(n_cells as usize) {
            for gene in 0..(n_genes as usize) {
                let dense_idx = gene
                    .checked_mul(n_cells as usize)
                    .and_then(|x| x.checked_add(cell))
                    .ok_or_else(|| {
                        KiraError::new(
                            ErrorKind::TsvParse,
                            format!("dense index overflow in {}", path.display()),
                        )
                    })?;
                let offset = MTX_HEADER_LEN + dense_idx * 4;
                let v = read_f32(data, offset)?;
                if v != 0.0 {
                    gene_idx.push(gene as u32);
                    values.push(v);
                }
            }
            cell_offsets.push(gene_idx.len() as u64);
        }

        Ok(Self {
            _mmap: mmap,
            n_cells,
            n_genes,
            gene_symbols: symbols,
            gene_lookup,
            cell_offsets,
            gene_idx,
            values,
        })
    }

    pub fn n_cells(&self) -> usize {
        self.n_cells as usize
    }

    pub fn n_genes(&self) -> usize {
        self.n_genes as usize
    }

    pub fn gene_index(&self, symbol: &str) -> Option<u32> {
        if let Some(v) = self.gene_lookup.get(symbol) {
            return Some(*v);
        }
        let normalized = normalize_gene_key(symbol);
        self.gene_lookup.get(&normalized).copied()
    }

    pub fn gene_symbol(&self, gene_idx: u32) -> &str {
        self.gene_symbols
            .get(gene_idx as usize)
            .map(|s| s.as_str())
            .unwrap_or("")
    }

    pub fn iter_cell(&self, cell_idx: u32) -> impl Iterator<Item = (u32, f32)> + '_ {
        let idx = cell_idx as usize;
        let start = self.cell_offsets.get(idx).copied().unwrap_or(0) as usize;
        let end = self
            .cell_offsets
            .get(idx + 1)
            .copied()
            .unwrap_or(start as u64) as usize;
        let end = end.min(self.gene_idx.len()).min(self.values.len());
        let start = start.min(end);
        self.gene_idx[start..end]
            .iter()
            .copied()
            .zip(self.values[start..end].iter().copied())
    }
}

fn load_gene_symbols_for_kiramtx(expr_path: &Path, n_genes: usize) -> Result<Vec<String>> {
    let dir = expr_path.parent().ok_or_else(|| {
        KiraError::new(
            ErrorKind::Path,
            format!(
                "failed to resolve parent directory for {}",
                expr_path.display()
            ),
        )
    })?;
    let candidates = [
        dir.join("features.tsv"),
        dir.join("features.tsv.gz"),
        dir.join("genes.tsv"),
        dir.join("genes.tsv.gz"),
    ];

    for path in candidates {
        if !path.exists() {
            continue;
        }
        if let Ok(symbols) = parse_feature_symbols(&path, n_genes)
            && symbols.len() == n_genes
        {
            return Ok(symbols);
        }
    }

    // Fallback for caches without a colocated features table.
    Ok((0..n_genes).map(|i| format!("GENE_{i}")).collect())
}

fn parse_feature_symbols(path: &Path, n_genes: usize) -> Result<Vec<String>> {
    let input = read_tsv_input(path)?;
    let mut reader = TsvReader::new(input.as_bytes());
    let mut row = Vec::new();
    let mut out = Vec::with_capacity(n_genes);
    let mut first = true;
    let mut symbol_col = 0usize;
    let mut header = false;

    while let Some(rec) = reader.next_record(&mut row) {
        rec?;
        if row.is_empty() {
            continue;
        }
        if first {
            first = false;
            if looks_like_feature_header(&row) {
                header = true;
                symbol_col = row
                    .iter()
                    .position(|v| {
                        let low = v.to_ascii_lowercase();
                        low == "gene_symbol" || low == "symbol" || low == "gene"
                    })
                    .unwrap_or(if row.len() >= 2 { 1 } else { 0 });
                continue;
            } else {
                symbol_col = if row.len() >= 2 { 1 } else { 0 };
            }
        }

        if row.len() <= symbol_col {
            return Err(KiraError::new(
                ErrorKind::TsvParse,
                format!("feature row missing symbol column in {}", path.display()),
            ));
        }
        out.push(row[symbol_col].clone());
        if out.len() == n_genes {
            break;
        }
    }

    if !header && first {
        return Err(KiraError::new(
            ErrorKind::TsvParse,
            format!("empty features table {}", path.display()),
        ));
    }
    if out.len() != n_genes {
        return Err(KiraError::new(
            ErrorKind::TsvParse,
            format!(
                "features count mismatch in {}: expected {}, got {}",
                path.display(),
                n_genes,
                out.len()
            ),
        ));
    }

    Ok(out)
}

fn looks_like_feature_header(row: &[String]) -> bool {
    if row.is_empty() {
        return false;
    }
    let c0 = row[0].to_ascii_lowercase();
    let c1 = row
        .get(1)
        .map(|x| x.to_ascii_lowercase())
        .unwrap_or_default();
    matches!(
        c0.as_str(),
        "gene_id" | "feature_id" | "id" | "gene" | "symbol"
    ) || matches!(c1.as_str(), "gene_symbol" | "symbol" | "gene" | "gene_name")
}

fn read_u32(data: &[u8], offset: usize) -> Result<u32> {
    let end = offset + 4;
    if end > data.len() {
        return Err(KiraError::new(
            ErrorKind::TsvParse,
            "expr read_u32 out of bounds",
        ));
    }
    Ok(u32::from_le_bytes(data[offset..end].try_into().unwrap()))
}

fn read_u64(data: &[u8], offset: usize) -> Result<u64> {
    let end = offset + 8;
    if end > data.len() {
        return Err(KiraError::new(
            ErrorKind::TsvParse,
            "expr read_u64 out of bounds",
        ));
    }
    Ok(u64::from_le_bytes(data[offset..end].try_into().unwrap()))
}

fn read_f32(data: &[u8], offset: usize) -> Result<f32> {
    let end = offset + 4;
    if end > data.len() {
        return Err(KiraError::new(
            ErrorKind::TsvParse,
            "expr read_f32 out of bounds",
        ));
    }
    Ok(f32::from_le_bytes(data[offset..end].try_into().unwrap()))
}

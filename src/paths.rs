use crate::error::{ErrorKind, KiraError, Result};
use sha2::{Digest, Sha256};
use std::fs;
use std::path::{Path, PathBuf};
use std::time::UNIX_EPOCH;

#[derive(Debug, Clone, serde::Serialize)]
pub struct PathMeta {
    pub kind: String,
    pub path_abs: String,
    pub size_bytes: u64,
    pub mtime_unix: i64,
    pub sha256_8: String,
}

pub fn absolutize(path: &Path) -> Result<PathBuf> {
    let abs = if path.is_absolute() {
        path.to_path_buf()
    } else {
        std::env::current_dir()
            .map_err(|e| KiraError::new(ErrorKind::Path, format!("failed current_dir: {e}")))?
            .join(path)
    };
    if abs.exists() {
        fs::canonicalize(&abs).map_err(|e| {
            KiraError::new(
                ErrorKind::Path,
                format!("failed to canonicalize {}: {e}", abs.display()),
            )
        })
    } else {
        Ok(abs)
    }
}

pub fn must_exist_file(path: &Path, kind: &str) -> Result<PathBuf> {
    let abs = absolutize(path)?;
    if !abs.exists() {
        return Err(KiraError::new(
            ErrorKind::InputMissing,
            format!("missing {kind}: {}", abs.display()),
        ));
    }
    if !abs.is_file() {
        return Err(KiraError::new(
            ErrorKind::InputMissing,
            format!("expected file for {kind}: {}", abs.display()),
        ));
    }
    Ok(abs)
}

pub fn must_exist_path(path: &Path, kind: &str) -> Result<PathBuf> {
    let abs = absolutize(path)?;
    if !abs.exists() {
        return Err(KiraError::new(
            ErrorKind::InputMissing,
            format!("missing {kind}: {}", abs.display()),
        ));
    }
    Ok(abs)
}

pub fn metadata_for_path(kind: impl Into<String>, path: &Path) -> Result<PathMeta> {
    let meta = fs::metadata(path).map_err(|e| {
        KiraError::new(
            ErrorKind::Path,
            format!("metadata failed {}: {e}", path.display()),
        )
    })?;
    let size = meta.len();
    let mtime_unix = meta
        .modified()
        .ok()
        .and_then(|t| t.duration_since(UNIX_EPOCH).ok())
        .map(|d| d.as_secs() as i64)
        .unwrap_or(0);
    let sha256_8 = sha256_8(path)?;
    Ok(PathMeta {
        kind: kind.into(),
        path_abs: path.to_string_lossy().into_owned(),
        size_bytes: size,
        mtime_unix,
        sha256_8,
    })
}

fn sha256_8(path: &Path) -> Result<String> {
    let data = if path.is_file() {
        fs::read(path).map_err(|e| {
            KiraError::new(
                ErrorKind::Path,
                format!("read failed {}: {e}", path.display()),
            )
        })?
    } else {
        Vec::new()
    };
    let mut hasher = Sha256::new();
    hasher.update(&data);
    let out = hasher.finalize();
    let hex = format!("{:x}", out);
    Ok(hex.chars().take(8).collect())
}

use crate::error::{ErrorKind, KiraError, Result};
use crate::io::json::to_pretty_json_bytes;
use serde::Serialize;
#[cfg(feature = "strict-fsync")]
use std::fs::File;
use std::fs::{self, OpenOptions};
use std::io::Write;
use std::path::Path;
use std::time::{SystemTime, UNIX_EPOCH};

pub fn write_bytes_atomic(path: &Path, bytes: &[u8]) -> Result<()> {
    let parent = path.parent().ok_or_else(|| {
        KiraError::new(ErrorKind::Path, format!("no parent for {}", path.display()))
    })?;
    fs::create_dir_all(parent).map_err(|e| {
        KiraError::new(
            ErrorKind::Path,
            format!("mkdir failed {}: {e}", parent.display()),
        )
    })?;

    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_nanos())
        .unwrap_or(0);
    let tmp = parent.join(format!(
        ".{}.tmp.{}.{}",
        path.file_name().and_then(|x| x.to_str()).unwrap_or("out"),
        std::process::id(),
        nanos
    ));

    let mut file = OpenOptions::new()
        .create_new(true)
        .write(true)
        .open(&tmp)
        .map_err(|e| {
            KiraError::new(
                ErrorKind::Path,
                format!("temp open failed {}: {e}", tmp.display()),
            )
        })?;
    file.write_all(bytes).map_err(|e| {
        KiraError::new(
            ErrorKind::Path,
            format!("write failed {}: {e}", tmp.display()),
        )
    })?;

    #[cfg(feature = "strict-fsync")]
    {
        file.sync_all().map_err(|e| {
            KiraError::new(
                ErrorKind::Path,
                format!("fsync temp failed {}: {e}", tmp.display()),
            )
        })?;
    }

    drop(file);
    fs::rename(&tmp, path).map_err(|e| {
        KiraError::new(
            ErrorKind::Path,
            format!("rename failed {} -> {}: {e}", tmp.display(), path.display()),
        )
    })?;

    #[cfg(feature = "strict-fsync")]
    {
        if let Ok(dirf) = File::open(parent) {
            let _ = dirf.sync_all();
        }
    }

    Ok(())
}

pub fn write_json_atomic<T: Serialize>(path: &Path, value: &T) -> Result<()> {
    let bytes = to_pretty_json_bytes(value)?;
    write_bytes_atomic(path, &bytes)
}

pub fn write_tsv_atomic(path: &Path, lines: &[String]) -> Result<()> {
    let mut out = Vec::new();
    for line in lines {
        out.extend_from_slice(line.as_bytes());
        out.push(b'\n');
    }
    write_bytes_atomic(path, &out)
}

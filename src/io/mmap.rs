use crate::error::{ErrorKind, KiraError, Result};
use memmap2::Mmap;
use std::fs::File;
use std::path::Path;

pub struct MmapFile {
    _file: File,
    map: Mmap,
}

impl MmapFile {
    pub fn open(path: &Path) -> Result<Self> {
        let file = File::open(path).map_err(|e| {
            KiraError::new(
                ErrorKind::InputMissing,
                format!("open failed {}: {e}", path.display()),
            )
        })?;
        // SAFETY: File is kept alive in struct for lifetime of mapping and used read-only.
        let map = unsafe { Mmap::map(&file) }.map_err(|e| {
            KiraError::new(
                ErrorKind::Path,
                format!("mmap failed {}: {e}", path.display()),
            )
        })?;
        Ok(Self { _file: file, map })
    }

    pub fn as_bytes(&self) -> &[u8] {
        &self.map
    }
}

use crate::error::{ErrorKind, KiraError, Result};
use crate::io::mmap::MmapFile;
use flate2::read::GzDecoder;
use std::io::Read;
use std::path::Path;

pub enum InputBytes {
    Mmap(MmapFile),
    Owned(Vec<u8>),
}

impl InputBytes {
    pub fn as_bytes(&self) -> &[u8] {
        match self {
            Self::Mmap(m) => m.as_bytes(),
            Self::Owned(v) => v.as_slice(),
        }
    }
}

pub fn read_tsv_input(path: &Path) -> Result<InputBytes> {
    if path.extension().map(|x| x == "gz").unwrap_or(false) {
        let file = std::fs::File::open(path).map_err(|e| {
            KiraError::new(
                ErrorKind::InputMissing,
                format!("open failed {}: {e}", path.display()),
            )
        })?;
        let mut decoder = GzDecoder::new(file);
        let mut out = Vec::new();
        decoder.read_to_end(&mut out).map_err(|e| {
            KiraError::new(
                ErrorKind::TsvParse,
                format!("gzip decode failed {}: {e}", path.display()),
            )
        })?;
        Ok(InputBytes::Owned(out))
    } else {
        Ok(InputBytes::Mmap(MmapFile::open(path)?))
    }
}

use crate::error::{ErrorKind, KiraError, Result};
use serde::Serialize;

pub fn to_pretty_json_bytes<T: Serialize>(value: &T) -> Result<Vec<u8>> {
    let mut out = serde_json::to_vec_pretty(value)
        .map_err(|e| KiraError::new(ErrorKind::Path, format!("json serialization failed: {e}")))?;
    out.push(b'\n');
    Ok(out)
}

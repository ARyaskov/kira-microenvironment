use crate::error::{ErrorKind, KiraError, Result};

pub struct TsvReader<'a> {
    data: &'a [u8],
    pos: usize,
    line_no: usize,
}

impl<'a> TsvReader<'a> {
    pub fn new(data: &'a [u8]) -> Self {
        Self {
            data,
            pos: 0,
            line_no: 0,
        }
    }

    pub fn line_no(&self) -> usize {
        self.line_no
    }

    pub fn next_record(&mut self, fields: &mut Vec<String>) -> Option<Result<()>> {
        if self.pos >= self.data.len() {
            return None;
        }
        fields.clear();

        let start = self.pos;
        while self.pos < self.data.len() && self.data[self.pos] != b'\n' {
            self.pos += 1;
        }
        let mut end = self.pos;
        if self.pos < self.data.len() && self.data[self.pos] == b'\n' {
            self.pos += 1;
        }
        if end > start && self.data[end - 1] == b'\r' {
            end -= 1;
        }
        self.line_no += 1;

        let line = &self.data[start..end];
        if line.is_empty() {
            return Some(Ok(()));
        }

        let mut field_start = 0usize;
        for i in 0..=line.len() {
            if i == line.len() || line[i] == b'\t' {
                let raw = &line[field_start..i];
                match std::str::from_utf8(raw) {
                    Ok(s) => fields.push(s.to_string()),
                    Err(e) => {
                        return Some(Err(KiraError::new(
                            ErrorKind::TsvParse,
                            format!("invalid utf-8 at line {}: {e}", self.line_no),
                        )));
                    }
                }
                field_start = i + 1;
            }
        }
        Some(Ok(()))
    }
}

pub fn expect_header_exact(found: &[String], expected: &[&str], file: &str) -> Result<()> {
    let exact =
        found.len() == expected.len() && found.iter().zip(expected.iter()).all(|(l, r)| l == r);
    if exact {
        return Ok(());
    }
    Err(KiraError::new(
        ErrorKind::TsvHeader,
        format!(
            "header mismatch in {file}: expected [{}], found [{}]",
            expected.join(","),
            found.join(",")
        ),
    ))
}

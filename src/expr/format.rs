pub const EXPR_MAGIC: [u8; 8] = *b"KIRAEXPR";
pub const EXPR_VERSION: u32 = 1;

pub const EXPR_HEADER_LEN: usize = 8 + 4 + 4 + 4 + 8 + 8 + 8;

pub const MTX_MAGIC: [u8; 8] = *b"KIRAMTX\0";
pub const MTX_VERSION: u32 = 1;
pub const MTX_HEADER_LEN: usize = 8 + 4 + 4 + 4 + 4;

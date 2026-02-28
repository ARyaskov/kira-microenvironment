use crate::auto_groups::markers::{MarkerGroup, load_marker_groups};
use crate::error::Result;
use crate::expr::reader::ExprReader;
use crate::resources::aliases::AliasResolver;
use std::path::Path;

pub fn load_fine_groups(
    path: &Path,
    expr: &ExprReader,
    aliases: &AliasResolver,
) -> Result<Vec<MarkerGroup>> {
    load_marker_groups(path, expr, aliases)
}

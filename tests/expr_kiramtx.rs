use kira_environment::expr::reader::ExprReader;
use std::fs;
use std::path::Path;
use tempfile::tempdir;

#[test]
fn expr_reader_supports_kiramtx_with_features_sidecar() {
    let td = tempdir().expect("tempdir");
    let root = td.path();
    let expr = root.join("expr.bin");
    let features = root.join("features.tsv");

    fs::write(&features, "ENSG1\tEPCAM\tGene\nENSG2\tCD3D\tGene\n").expect("features");
    write_kiramtx_expr(&expr).expect("expr");

    let reader = ExprReader::open(&expr).expect("open expr");
    assert_eq!(reader.n_genes(), 2);
    assert_eq!(reader.n_cells(), 3);
    assert_eq!(reader.gene_index("EPCAM"), Some(0));
    assert_eq!(reader.gene_index("CD3D"), Some(1));

    let c0: Vec<(u32, f32)> = reader.iter_cell(0).collect();
    assert_eq!(c0, vec![(0, 1.0), (1, 2.0)]);
    let c1: Vec<(u32, f32)> = reader.iter_cell(1).collect();
    assert_eq!(c1, vec![(1, 3.0)]);
    let c2: Vec<(u32, f32)> = reader.iter_cell(2).collect();
    assert_eq!(c2, vec![(0, 4.0)]);
}

fn write_kiramtx_expr(path: &Path) -> Result<(), std::io::Error> {
    let mut out = Vec::new();
    out.extend_from_slice(b"KIRAMTX\0");
    out.extend_from_slice(&1u32.to_le_bytes()); // version
    out.extend_from_slice(&2u32.to_le_bytes()); // genes
    out.extend_from_slice(&3u32.to_le_bytes()); // samples(cells)
    out.extend_from_slice(&3u32.to_le_bytes()); // flags

    // Row-major [gene][sample]
    let values: [f32; 6] = [
        1.0, 0.0, 4.0, // gene0 across samples
        2.0, 3.0, 0.0, // gene1 across samples
    ];
    for v in values {
        out.extend_from_slice(&v.to_le_bytes());
    }
    fs::write(path, out)
}

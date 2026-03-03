use kira_environment::stages::stage5_microenv_extension::{Stage5Config, run_stage5};
use std::fs;
use std::path::Path;
use tempfile::tempdir;

#[test]
fn stage5_microenv_outputs_are_deterministic() {
    let td = tempdir().expect("tempdir");
    let root = td.path();
    let out = root.join("out");
    let km = out.join("kira-microenvironment");
    let s0 = km.join("stage0_resolve");
    fs::create_dir_all(&s0).expect("mkdir stage0");

    let expr = root.join("expr.bin");
    write_expr_bin(&expr).expect("expr");
    let barcodes = root.join("barcodes.tsv");
    fs::write(&barcodes, "cell_id\nC1\nC2\nC3\nC4\n").expect("barcodes");
    fs::write(
        s0.join("groups_normalized.tsv"),
        "cell_id\tgroup\nC1\tA\nC2\tA\nC3\tB\nC4\tB\n",
    )
    .expect("groups");
    fs::write(
        s0.join("resolved.json"),
        format!(
            "{{\"inputs\":[{{\"kind\":\"barcodes\",\"path_abs\":\"{}\"}}]}}\n",
            barcodes.to_string_lossy()
        ),
    )
    .expect("resolved");

    let cfg = Stage5Config {
        expr: expr.clone(),
        out_dir: out.clone(),
    };

    let r1 = run_stage5(cfg.clone()).expect("stage5 #1");
    let m1 = fs::read(km.join("metrics.tsv")).expect("metrics1");
    let s1 = fs::read(r1.out_stage_dir.join("stage5_summary.json")).expect("summary1");

    let r2 = run_stage5(cfg).expect("stage5 #2");
    let m2 = fs::read(km.join("metrics.tsv")).expect("metrics2");
    let s2 = fs::read(r2.out_stage_dir.join("stage5_summary.json")).expect("summary2");

    assert_eq!(m1, m2);
    assert_eq!(s1, s2);

    let metrics = String::from_utf8(m2).expect("metrics utf8");
    assert!(
        metrics
            .lines()
            .next()
            .unwrap_or("")
            .contains("microenv_stress_mode")
    );
    assert_eq!(metrics.lines().count(), 5);

    let summary = String::from_utf8(s2).expect("summary utf8");
    assert!(summary.contains("\"panel_version\": \"MICROENV_EXTENSION_PANEL_V1\""));
    assert!(summary.contains("\"thresholds\""));
    assert!(summary.contains("\"cluster_stats\""));
    assert!(summary.contains("\"missingness\""));
}

fn write_expr_bin(path: &Path) -> Result<(), std::io::Error> {
    let magic = *b"KIRAEXPR";
    let version = 1u32;
    let n_cells = 4u32;
    let n_genes = 4u32;
    let genes = ["HIF1A", "VEGFA", "NT5E", "ENTPD1"];

    let mut gene_block = Vec::new();
    for g in genes {
        gene_block.extend_from_slice(&(g.len() as u32).to_le_bytes());
        gene_block.extend_from_slice(g.as_bytes());
    }

    let cell_offsets: [u64; 5] = [0, 2, 4, 6, 8];
    let gene_idx: [u32; 8] = [0, 2, 0, 2, 1, 3, 1, 3];
    let values: [f32; 8] = [2.0, 1.0, 1.5, 1.2, 1.8, 1.4, 1.0, 0.8];

    let header_len = 8 + 4 + 4 + 4 + 8 + 8 + 8;
    let gene_symbols_offset = header_len as u64;
    let cell_index_offset = gene_symbols_offset + gene_block.len() as u64;
    let nnz_offset = cell_index_offset + ((n_cells as usize + 1) * 8) as u64;

    let mut out = Vec::new();
    out.extend_from_slice(&magic);
    out.extend_from_slice(&version.to_le_bytes());
    out.extend_from_slice(&n_cells.to_le_bytes());
    out.extend_from_slice(&n_genes.to_le_bytes());
    out.extend_from_slice(&gene_symbols_offset.to_le_bytes());
    out.extend_from_slice(&cell_index_offset.to_le_bytes());
    out.extend_from_slice(&nnz_offset.to_le_bytes());
    out.extend_from_slice(&gene_block);
    for off in cell_offsets {
        out.extend_from_slice(&off.to_le_bytes());
    }
    for g in gene_idx {
        out.extend_from_slice(&g.to_le_bytes());
    }
    for v in values {
        out.extend_from_slice(&v.to_le_bytes());
    }
    fs::write(path, out)
}

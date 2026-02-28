use kira_environment::cli::{AggModeArg, CapModeArg};
use kira_environment::stages::stage1_agg::{Stage1Config, run_stage1};
use std::fs;
use std::path::Path;
use tempfile::tempdir;

#[test]
fn stage1_agg_outputs_and_cache_are_deterministic() {
    let td = tempdir().expect("tempdir");
    let root = td.path();

    let out = root.join("out");
    let s0 = out.join("kira-microenvironment").join("stage0_resolve");
    fs::create_dir_all(&s0).expect("mkdir stage0");

    let barcodes_path = root.join("barcodes.tsv");
    fs::write(&barcodes_path, "cell_id\nC1\nC2\nC3\nC4\n").expect("barcodes");

    fs::write(
        s0.join("groups_normalized.tsv"),
        "cell_id\tgroup\nC1\tG1\nC2\tG1\nC3\tG2\nC4\tG2\n",
    )
    .expect("groups_normalized");

    fs::write(s0.join("group_sizes.tsv"), "group\tn_cells\nG1\t2\nG2\t2\n").expect("group_sizes");

    fs::write(
        s0.join("genes_requested.tsv"),
        "role\tquery_symbol\tcomplex_id\tsource\nligand\tG1\t\tlr_pairs\nreceptor\tG2\t\tlr_pairs\nsubunit\tG3\tCX\tcofactors\n",
    )
    .expect("genes_requested");

    fs::write(
        s0.join("resolved.json"),
        format!(
            "{{\n  \"inputs\": [{{\"kind\": \"barcodes\", \"path_abs\": \"{}\"}}]\n}}\n",
            barcodes_path.to_string_lossy()
        ),
    )
    .expect("resolved");

    let expr_path = root.join("expr.bin");
    write_expr_bin(&expr_path).expect("expr");

    let cfg = Stage1Config {
        expr: expr_path.clone(),
        out_dir: out.clone(),
        agg_mode: AggModeArg::Median,
        trim: 0.05,
        eps: 1e-6,
        cap_mode: CapModeArg::P99,
        cap_p: 0.99,
        cap_fixed: None,
    };

    let r1 = run_stage1(cfg.clone()).expect("stage1 first");
    let gga1 = fs::read(r1.out_stage_dir.join("group_gene_agg.tsv")).expect("gga1");
    let gst1 = fs::read(r1.out_stage_dir.join("gene_stats.tsv")).expect("gst1");
    let c1 = fs::read(r1.out_stage_dir.join("cache.bin")).expect("cache1");

    let r2 = run_stage1(cfg).expect("stage1 second");
    let gga2 = fs::read(r2.out_stage_dir.join("group_gene_agg.tsv")).expect("gga2");
    let gst2 = fs::read(r2.out_stage_dir.join("gene_stats.tsv")).expect("gst2");
    let c2 = fs::read(r2.out_stage_dir.join("cache.bin")).expect("cache2");

    assert_eq!(gga1, gga2);
    assert_eq!(gst1, gst2);
    assert_eq!(c1, c2);

    let gga = String::from_utf8(gga1).expect("utf8 gga");
    let expected = concat!(
        "group\tgene\texpr\tcov\tn_cells\n",
        "G1\tG1\t1.000000\t0.500000\t2\n",
        "G1\tG2\t1.500000\t1.000000\t2\n",
        "G1\tG3\t0.500000\t0.500000\t2\n",
        "G2\tG1\t1.500000\t0.500000\t2\n",
        "G2\tG2\t0.000000\t0.000000\t2\n",
        "G2\tG3\t1.000000\t0.500000\t2\n",
    );
    assert_eq!(gga, expected);

    let (groups, genes, expr, cov) = read_cache_bin(&c1).expect("read cache");
    assert_eq!(groups, vec!["G1", "G2"]);
    assert_eq!(genes, vec!["G1", "G2", "G3"]);
    assert_eq!(expr, vec![1.0, 1.5, 0.5, 1.5, 0.0, 1.0]);
    assert_eq!(cov, vec![0.5, 1.0, 0.5, 0.5, 0.0, 0.5]);
}

fn write_expr_bin(path: &Path) -> Result<(), std::io::Error> {
    let magic = *b"KIRAEXPR";
    let version = 1u32;
    let n_cells = 4u32;
    let n_genes = 3u32;

    let genes = ["G1", "G2", "G3"];
    let mut gene_block = Vec::new();
    for g in genes {
        gene_block.extend_from_slice(&(g.len() as u32).to_le_bytes());
        gene_block.extend_from_slice(g.as_bytes());
    }

    // cell offsets for 4 cells with 8 nnz total
    let cell_offsets: [u64; 5] = [0, 2, 4, 6, 8];

    // C1: G1=2, G2=1
    // C2: G2=2, G3=1
    // C3: G1=3, G3=2
    // C4: G1=0, G2=0 (explicit zeros)
    let gene_idx: [u32; 8] = [0, 1, 1, 2, 0, 2, 0, 1];
    let values: [f32; 8] = [2.0, 1.0, 2.0, 1.0, 3.0, 2.0, 0.0, 0.0];

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

fn read_cache_bin(bytes: &[u8]) -> Result<(Vec<String>, Vec<String>, Vec<f32>, Vec<f32>), String> {
    if bytes.len() < 20 {
        return Err("cache too small".to_string());
    }
    if &bytes[0..8] != b"KIRAEAGG" {
        return Err("bad magic".to_string());
    }
    let version = u32::from_le_bytes(bytes[8..12].try_into().unwrap());
    if version != 1 {
        return Err("bad version".to_string());
    }
    let n_groups = u32::from_le_bytes(bytes[12..16].try_into().unwrap()) as usize;
    let n_genes = u32::from_le_bytes(bytes[16..20].try_into().unwrap()) as usize;

    let mut pos = 20usize;
    let mut groups = Vec::with_capacity(n_groups);
    for _ in 0..n_groups {
        let len = u32::from_le_bytes(bytes[pos..pos + 4].try_into().unwrap()) as usize;
        pos += 4;
        let s = std::str::from_utf8(&bytes[pos..pos + len]).map_err(|e| e.to_string())?;
        groups.push(s.to_string());
        pos += len;
    }

    let mut genes = Vec::with_capacity(n_genes);
    for _ in 0..n_genes {
        let len = u32::from_le_bytes(bytes[pos..pos + 4].try_into().unwrap()) as usize;
        pos += 4;
        let s = std::str::from_utf8(&bytes[pos..pos + len]).map_err(|e| e.to_string())?;
        genes.push(s.to_string());
        pos += len;
    }

    let n = n_groups * n_genes;
    let mut expr = Vec::with_capacity(n);
    for _ in 0..n {
        expr.push(f32::from_le_bytes(bytes[pos..pos + 4].try_into().unwrap()));
        pos += 4;
    }

    let mut cov = Vec::with_capacity(n);
    for _ in 0..n {
        cov.push(f32::from_le_bytes(bytes[pos..pos + 4].try_into().unwrap()));
        pos += 4;
    }

    Ok((groups, genes, expr, cov))
}

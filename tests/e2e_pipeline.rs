use kira_environment::cli::{AggModeArg, AutoGroupsModeArg, CapModeArg, SpecModeArg};
use kira_environment::simd::dispatch::SimdLevel;
use kira_environment::{RunConfig, run};
use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};
use tempfile::tempdir;

#[test]
fn e2e_pipeline_is_deterministic_and_schema_stable() {
    let td = tempdir().expect("tempdir");
    let root = td.path();

    let expr = root.join("expr.bin");
    write_expr_bin(&expr).expect("expr");

    let barcodes = root.join("barcodes.tsv");
    fs::write(&barcodes, "cell_id\nC1\nC2\nC3\nC4\n").expect("barcodes");

    let groups = root.join("groups.tsv");
    fs::write(&groups, "cell_id\tgroup\nC1\tA\nC2\tA\nC3\tB\nC4\tB\n").expect("groups");

    let resources = root.join("resources");
    fs::create_dir_all(&resources).expect("resources");
    fs::write(
        resources.join("lr_pairs.tsv"),
        "ligand_symbol\treceptor_symbol\tfamily\tdirectionality\tweight\nG1\tG2\tF\tbidirectional\t1.0\n",
    )
    .expect("lr_pairs");
    fs::write(
        resources.join("cofactors.tsv"),
        "complex_id\trole\tsubunit_symbol\trequired\tlogic\n",
    )
    .expect("cofactors");
    fs::write(
        resources.join("gene_alias.tsv"),
        "alias\tsymbol\tpriority\tnotes\nG1\tG1\t1\t-\nG2\tG2\t1\t-\n",
    )
    .expect("aliases");
    fs::write(
        resources.join("labels.tsv"),
        "ligand_symbol\treceptor_symbol\tlabel\tweight\tnotes\n*\t*\tBASE\t1.0\t-\n",
    )
    .expect("labels");
    fs::write(
        resources.join("regime_map.tsv"),
        "regime\tclass\tnotes\nR_LOUD\tLOUD\t-\nR_SILENT\tSILENT\t-\n",
    )
    .expect("regime map");

    let secretion = root.join("secretion");
    fs::create_dir_all(&secretion).expect("secretion");
    fs::write(
        secretion.join("secretion_groups.tsv"),
        "group\tregime\tfraction\nA\tR_LOUD\t1.0\nB\tR_SILENT\t1.0\n",
    )
    .expect("secretion groups");

    // Full run with Stage4
    let out_full = root.join("out_full");
    let cfg_full = base_cfg(
        expr.clone(),
        barcodes.clone(),
        groups.clone(),
        resources.clone(),
        out_full.clone(),
        Some(secretion.clone()),
    );

    run(cfg_full.clone()).expect("run full #1");
    let snap1 = snapshot_dir(&out_full.join("kira-microenvironment")).expect("snapshot1");
    run(cfg_full).expect("run full #2");
    let snap2 = snapshot_dir(&out_full.join("kira-microenvironment")).expect("snapshot2");
    assert_eq!(snap1, snap2);

    assert_file_header(
        &out_full.join("kira-microenvironment/stage2_score/edges_raw.tsv"),
        "source_group\ttarget_group\tligand\treceptor\tscore\tL_expr\tR_expr\tcov_L\tcov_R\tspec_L\tspec_R\tflags",
    );
    assert_file_header(
        &out_full.join("kira-microenvironment/edges.tsv"),
        "source_group\ttarget_group\tligand\treceptor\tscore\tL_expr\tR_expr\tcov_L\tcov_R\tspec_L\tspec_R\tflags",
    );
    assert_file_header(
        &out_full.join("kira-microenvironment/stage4_link/loudness_drivers.tsv"),
        "rank\tligand\treceptor\tfold\tmean_loud\tmean_silent\tn_edges_loud\tn_edges_silent\texemplar_source_group\texemplar_target_group\texemplar_score",
    );

    let summary_full = fs::read_to_string(out_full.join("kira-microenvironment/summary.json"))
        .expect("summary full");
    assert!(summary_full.contains("\"tool\""));
    assert!(summary_full.contains("\"inputs\""));
    assert!(summary_full.contains("\"resources\""));
    assert!(summary_full.contains("\"stage0\""));
    assert!(summary_full.contains("\"stage1\""));
    assert!(summary_full.contains("\"stage2\""));
    assert!(summary_full.contains("\"stage3\""));
    assert!(summary_full.contains("\"stage4\""));
    assert!(summary_full.contains("\"network\""));
    assert!(summary_full.contains("\"linking\""));

    // Run without secretion (Stage4 skipped)
    let out_no_link = root.join("out_no_link");
    let cfg_no_link = base_cfg(expr, barcodes, groups, resources, out_no_link.clone(), None);
    run(cfg_no_link).expect("run without secretion");

    assert!(out_no_link.join("kira-microenvironment/edges.tsv").exists());
    assert!(
        !out_no_link
            .join("kira-microenvironment/stage4_link/loudness_drivers.tsv")
            .exists()
    );
}

fn base_cfg(
    expr: PathBuf,
    barcodes: PathBuf,
    groups: PathBuf,
    resources: PathBuf,
    out: PathBuf,
    secretion: Option<PathBuf>,
) -> RunConfig {
    RunConfig {
        expr,
        barcodes,
        groups: Some(groups),
        auto_groups: None,
        auto_groups_coarse: None,
        auto_groups_fine: None,
        auto_groups_anti: None,
        auto_groups_mode: AutoGroupsModeArg::Flat,
        auto_groups_eps: 1e-6,
        auto_groups_min_delta: 0.1,
        auto_groups_unknown: "unknown".to_string(),
        auto_groups_emit_scores: false,
        resources,
        secretion,
        out,
        validate_only: false,
        lr_profile: "full".to_string(),
        simd_level: SimdLevel::Scalar,
        agg: AggModeArg::Median,
        trim: 0.05,
        eps: 1e-6,
        cap_mode: CapModeArg::P99,
        cap_p: 0.99,
        cap_fixed: None,
        cov_min: 0.10,
        expr_min: 0.05,
        spec: SpecModeArg::On,
        spec_cap: 10.0,
        top_n_per_pair: 200,
        top_n_per_source: 200,
        loud_thresh: 0.5,
        silent_thresh: 0.5,
        min_support: 1,
        regime_map: None,
    }
}

fn assert_file_header(path: &Path, expected: &str) {
    let content = fs::read_to_string(path).expect("read header file");
    let first = content.lines().next().unwrap_or("");
    assert_eq!(first, expected);
}

fn snapshot_dir(root: &Path) -> Result<BTreeMap<String, Vec<u8>>, std::io::Error> {
    let mut out = BTreeMap::new();
    walk(root, root, &mut out)?;
    Ok(out)
}

fn walk(
    root: &Path,
    dir: &Path,
    out: &mut BTreeMap<String, Vec<u8>>,
) -> Result<(), std::io::Error> {
    let mut entries = fs::read_dir(dir)?.collect::<Result<Vec<_>, _>>()?;
    entries.sort_by_key(|e| e.path());
    for e in entries {
        let p = e.path();
        if p.is_dir() {
            walk(root, &p, out)?;
        } else {
            let rel = p
                .strip_prefix(root)
                .unwrap_or(&p)
                .to_string_lossy()
                .to_string();
            out.insert(rel, fs::read(&p)?);
        }
    }
    Ok(())
}

fn write_expr_bin(path: &Path) -> Result<(), std::io::Error> {
    let magic = *b"KIRAEXPR";
    let version = 1u32;
    let n_cells = 4u32;
    let n_genes = 2u32;

    let genes = ["G1", "G2"];
    let mut gene_block = Vec::new();
    for g in genes {
        gene_block.extend_from_slice(&(g.len() as u32).to_le_bytes());
        gene_block.extend_from_slice(g.as_bytes());
    }

    // C1: G1=2, C2:G1=1, C3:G2=2, C4:G2=1
    let cell_offsets: [u64; 5] = [0, 1, 2, 3, 4];
    let gene_idx: [u32; 4] = [0, 0, 1, 1];
    let values: [f32; 4] = [2.0, 1.0, 2.0, 1.0];

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

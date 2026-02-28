use kira_environment::cli::{AggModeArg, AutoGroupsModeArg, CapModeArg, SpecModeArg};
use kira_environment::simd::dispatch::SimdLevel;
use kira_environment::{RunConfig, run};
use std::collections::BTreeMap;
use std::fs;
use std::path::Path;
use tempfile::tempdir;

#[test]
fn pipeline_runs_with_auto_groups_and_is_deterministic() {
    let td = tempdir().expect("tempdir");
    let root = td.path();

    let expr = root.join("expr.bin");
    write_expr_bin(&expr).expect("expr");

    let barcodes = root.join("barcodes.tsv");
    fs::write(&barcodes, "cell_id\nC1\nC2\nC3\nC4\n").expect("barcodes");

    let markers = root.join("marker_panels.tsv");
    fs::write(
        &markers,
        "group_name\tgene_symbol\tweight\nA\tG1\t1.0\nB\tG2\t1.0\n",
    )
    .expect("markers");

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

    let out = root.join("out");
    let cfg = RunConfig {
        expr,
        barcodes,
        groups: None,
        auto_groups: Some(markers),
        auto_groups_coarse: None,
        auto_groups_fine: None,
        auto_groups_anti: None,
        auto_groups_mode: AutoGroupsModeArg::Flat,
        auto_groups_eps: 1e-6,
        auto_groups_min_delta: 0.1,
        auto_groups_unknown: "unknown".to_string(),
        auto_groups_emit_scores: true,
        resources,
        embedded_profile: None,
        secretion: None,
        out: out.clone(),
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
        min_support: 3,
        regime_map: None,
    };

    run(cfg.clone()).expect("run #1");
    let snap1 = snapshot_dir(&out.join("kira-microenvironment")).expect("snapshot1");
    run(cfg).expect("run #2");
    let snap2 = snapshot_dir(&out.join("kira-microenvironment")).expect("snapshot2");
    assert_eq!(snap1, snap2);

    assert!(
        out.join("kira-microenvironment/auto_groups/groups.tsv")
            .exists()
    );
    assert!(
        out.join("kira-microenvironment/auto_groups/summary.json")
            .exists()
    );
    assert!(
        out.join("kira-microenvironment/auto_groups/cell_scores.tsv")
            .exists()
    );
    assert!(
        out.join("kira-microenvironment/stage0_resolve/groups_normalized.tsv")
            .exists()
    );
    assert!(
        out.join("kira-microenvironment/stage1_agg/group_gene_agg.tsv")
            .exists()
    );
    assert!(
        out.join("kira-microenvironment/stage2_score/edges_raw.tsv")
            .exists()
    );
    assert!(out.join("kira-microenvironment/edges.tsv").exists());

    let summary =
        fs::read_to_string(out.join("kira-microenvironment/summary.json")).expect("summary");
    assert!(summary.contains("\"auto_groups\""));
    assert!(summary.contains("\"enabled\": true"));
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

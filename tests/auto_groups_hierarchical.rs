use kira_environment::cli::{AggModeArg, AutoGroupsModeArg, CapModeArg, SpecModeArg};
use kira_environment::error::ErrorKind;
use kira_environment::simd::dispatch::SimdLevel;
use kira_environment::{RunConfig, run};
use std::collections::BTreeMap;
use std::fs;
use std::path::Path;
use tempfile::tempdir;

#[test]
fn hierarchical_assignment_fallback_and_determinism() {
    let td = tempdir().expect("tempdir");
    let root = td.path();
    write_common_inputs(root);

    let coarse = root.join("marker_panels_coarse.tsv");
    fs::write(
        &coarse,
        "group_name\tgene_symbol\tweight\ntumor\tEPCAM\t1.0\nimmune\tCD3D\t1.0\nstromal\tCOL1A1\t1.0\n",
    )
    .expect("coarse");
    let fine = root.join("marker_panels_fine.tsv");
    fs::write(
        &fine,
        "group_name\tgene_symbol\tweight\ntumor_epithelial\tEPCAM\t1.0\nT_cell\tCD3E\t1.0\nFibroblast\tCOL1A1\t1.0\n",
    )
    .expect("fine");

    let out = root.join("out_hier");
    let cfg = base_cfg(root, out.clone(), None);
    let cfg = RunConfig {
        auto_groups: None,
        auto_groups_coarse: Some(coarse),
        auto_groups_fine: Some(fine),
        auto_groups_anti: None,
        auto_groups_mode: AutoGroupsModeArg::Hierarchical,
        ..cfg
    };

    run(cfg.clone()).expect("run #1");
    let groups_path = out.join("kira-microenvironment/auto_groups/groups.tsv");
    let groups = fs::read_to_string(&groups_path).expect("groups");
    assert_eq!(
        groups,
        "cell_id\tgroup\nC1\ttumor_epithelial\nC2\timmune\nC3\tFibroblast\nC4\tunknown\n"
    );

    let auto_summary =
        fs::read_to_string(out.join("kira-microenvironment/auto_groups/summary.json"))
            .expect("summary");
    assert!(auto_summary.contains("\"mode\": \"hierarchical\""));
    assert!(auto_summary.contains("\"n_unknown\": 1"));
    assert!(auto_summary.contains("\"n_coarse_only\": 1"));
    assert!(auto_summary.contains("\"n_fine\": 2"));

    let snap1 = snapshot_dir(&out.join("kira-microenvironment/auto_groups")).expect("snapshot1");
    run(cfg).expect("run #2");
    let snap2 = snapshot_dir(&out.join("kira-microenvironment/auto_groups")).expect("snapshot2");
    assert_eq!(snap1, snap2);
}

#[test]
fn anti_markers_can_suppress_group_score() {
    let td = tempdir().expect("tempdir");
    let root = td.path();
    write_common_inputs(root);

    let markers = root.join("marker_panels.tsv");
    fs::write(
        &markers,
        "group_name\tgene_symbol\tweight\ntumor\tEPCAM\t1.0\nimmune\tCD3D\t1.0\n",
    )
    .expect("markers");
    let anti = root.join("anti_markers.tsv");
    fs::write(
        &anti,
        "group_name\tgene_symbol\tpenalty\nimmune\tEPCAM\t1.0\n",
    )
    .expect("anti");

    let out = root.join("out_anti");
    let cfg = RunConfig {
        auto_groups: Some(markers),
        auto_groups_coarse: None,
        auto_groups_fine: None,
        auto_groups_anti: Some(anti),
        auto_groups_mode: AutoGroupsModeArg::Flat,
        ..base_cfg(root, out.clone(), None)
    };
    run(cfg).expect("run");

    let groups = fs::read_to_string(out.join("kira-microenvironment/auto_groups/groups.tsv"))
        .expect("groups");
    assert!(groups.contains("C4\ttumor\n"));
}

#[test]
fn cli_validation_errors_are_stable() {
    let td = tempdir().expect("tempdir");
    let root = td.path();
    write_common_inputs(root);

    let markers = root.join("marker_panels.tsv");
    fs::write(
        &markers,
        "group_name\tgene_symbol\tweight\ntumor\tEPCAM\t1.0\nimmune\tCD3D\t1.0\n",
    )
    .expect("markers");

    let err = run(RunConfig {
        groups: Some(root.join("groups.tsv")),
        auto_groups: Some(markers.clone()),
        ..base_cfg(root, root.join("out_err1"), None)
    })
    .expect_err("expected error");
    assert_eq!(err.kind, ErrorKind::AutoGroupsConfigInvalid);

    let err = run(RunConfig {
        groups: None,
        auto_groups: None,
        auto_groups_coarse: Some(markers),
        auto_groups_fine: None,
        auto_groups_mode: AutoGroupsModeArg::Hierarchical,
        ..base_cfg(root, root.join("out_err2"), None)
    })
    .expect_err("expected error");
    assert_eq!(err.kind, ErrorKind::AutoGroupsConfigInvalid);
}

fn base_cfg(
    root: &Path,
    out: std::path::PathBuf,
    secretion: Option<std::path::PathBuf>,
) -> RunConfig {
    RunConfig {
        expr: root.join("expr.bin"),
        barcodes: root.join("barcodes.tsv"),
        groups: None,
        auto_groups: Some(root.join("marker_panels.tsv")),
        auto_groups_coarse: None,
        auto_groups_fine: None,
        auto_groups_anti: None,
        auto_groups_mode: AutoGroupsModeArg::Flat,
        auto_groups_eps: 0.05,
        auto_groups_min_delta: 0.2,
        auto_groups_unknown: "unknown".to_string(),
        auto_groups_emit_scores: false,
        resources: root.join("resources"),
        embedded_profile: None,
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
        min_support: 3,
        regime_map: None,
    }
}

fn write_common_inputs(root: &Path) {
    write_expr_bin(&root.join("expr.bin")).expect("expr");
    fs::write(root.join("barcodes.tsv"), "cell_id\nC1\nC2\nC3\nC4\n").expect("barcodes");
    fs::write(
        root.join("groups.tsv"),
        "cell_id\tgroup\nC1\tA\nC2\tA\nC3\tB\nC4\tB\n",
    )
    .expect("groups");

    let resources = root.join("resources");
    fs::create_dir_all(&resources).expect("resources");
    fs::write(
        resources.join("lr_pairs.tsv"),
        "ligand_symbol\treceptor_symbol\tfamily\tdirectionality\tweight\nEPCAM\tCOL1A1\tF\tbidirectional\t1.0\n",
    )
    .expect("lr_pairs");
    fs::write(
        resources.join("cofactors.tsv"),
        "complex_id\trole\tsubunit_symbol\trequired\tlogic\n",
    )
    .expect("cofactors");
    fs::write(
        resources.join("gene_alias.tsv"),
        "alias\tsymbol\tpriority\tnotes\nEPCAM\tEPCAM\t1\t-\nCD3D\tCD3D\t1\t-\nCOL1A1\tCOL1A1\t1\t-\n",
    )
    .expect("aliases");
    fs::write(
        resources.join("labels.tsv"),
        "ligand_symbol\treceptor_symbol\tlabel\tweight\tnotes\n*\t*\tBASE\t1.0\t-\n",
    )
    .expect("labels");
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
    let n_genes = 3u32;

    let genes = ["EPCAM", "CD3D", "COL1A1"];
    let mut gene_block = Vec::new();
    for g in genes {
        gene_block.extend_from_slice(&(g.len() as u32).to_le_bytes());
        gene_block.extend_from_slice(g.as_bytes());
    }

    // C1: EPCAM=3; C2: CD3D=3; C3: COL1A1=3; C4: EPCAM=2,CD3D=2
    let cell_offsets: [u64; 5] = [0, 1, 2, 3, 5];
    let gene_idx: [u32; 5] = [0, 1, 2, 0, 1];
    let values: [f32; 5] = [3.0, 3.0, 3.0, 2.0, 2.0];

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

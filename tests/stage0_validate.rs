use kira_environment::cli::{AggModeArg, AutoGroupsModeArg, CapModeArg, SpecModeArg};
use kira_environment::simd::dispatch::SimdLevel;
use kira_environment::{RunConfig, run};
use std::fs;
use std::io::Write;
use tempfile::tempdir;

#[test]
fn stage0_validate_outputs_are_stable() {
    let td = tempdir().expect("tempdir");
    let root = td.path();

    let resources = root.join("resources");
    fs::create_dir_all(&resources).expect("resources");

    fs::write(
        resources.join("lr_pairs.tsv"),
        "ligand_symbol\treceptor_symbol\tfamily\tdirectionality\nLIG1\tREC1\tF1\tbidirectional\nLIG2\tREC2\tF1\tbidirectional\n",
    )
    .expect("lr_pairs");

    fs::write(
        resources.join("cofactors.tsv"),
        "complex_id\trole\tsubunit_symbol\trequired\tlogic\nCX1\tligand\tSUBA\t1\tAND\nCX1\tligand\tSUBB\t1\tAND\n",
    )
    .expect("cofactors");

    fs::write(
        resources.join("gene_alias.tsv"),
        "alias\tsymbol\tpriority\tnotes\nSUBA\tSUBA\t1\t-\n",
    )
    .expect("aliases");

    fs::write(
        resources.join("labels.tsv"),
        "ligand_symbol\treceptor_symbol\tlabel\tweight\tnotes\n*\t*\tdefault\t1\t-\n",
    )
    .expect("labels");

    fs::write(root.join("barcodes.tsv"), "cell_id\nC1\nC2\nC3\nC4\n").expect("barcodes");
    fs::write(
        root.join("groups.tsv"),
        "cell_id\tgroup\nC1\tG2\nC2\tG1\nC3\tG1\nC4\tG2\n",
    )
    .expect("groups");
    fs::write(root.join("expr.bin"), b"").expect("expr");

    let out = root.join("out");

    let result = run(RunConfig {
        expr: root.join("expr.bin"),
        barcodes: root.join("barcodes.tsv"),
        groups: Some(root.join("groups.tsv")),
        auto_groups: None,
        auto_groups_coarse: None,
        auto_groups_fine: None,
        auto_groups_anti: None,
        auto_groups_mode: AutoGroupsModeArg::Flat,
        auto_groups_eps: 1e-6,
        auto_groups_min_delta: 0.1,
        auto_groups_unknown: "unknown".to_string(),
        auto_groups_emit_scores: false,
        resources: resources.clone(),
        embedded_profile: None,
        secretion: None,
        out: out.clone(),
        validate_only: true,
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
    })
    .expect("run stage0");

    let stage_dir = result.stage0.out_stage_dir;
    assert!(stage_dir.join("resolved.json").exists());
    assert!(stage_dir.join("genes_requested.tsv").exists());
    assert!(stage_dir.join("groups_normalized.tsv").exists());
    assert!(stage_dir.join("group_sizes.tsv").exists());

    let groups_norm =
        fs::read_to_string(stage_dir.join("groups_normalized.tsv")).expect("read groups");
    let expected_groups = "cell_id\tgroup\nC2\tG1\nC3\tG1\nC1\tG2\nC4\tG2\n";
    assert_eq!(groups_norm, expected_groups);

    let group_sizes = fs::read_to_string(stage_dir.join("group_sizes.tsv")).expect("read sizes");
    let expected_sizes = "group\tn_cells\nG1\t2\nG2\t2\n";
    assert_eq!(group_sizes, expected_sizes);

    let genes = fs::read_to_string(stage_dir.join("genes_requested.tsv")).expect("read genes");
    let expected_genes = concat!(
        "role\tquery_symbol\tcomplex_id\tsource\n",
        "ligand\tLIG1\t\tlr_pairs\n",
        "ligand\tLIG2\t\tlr_pairs\n",
        "receptor\tREC1\t\tlr_pairs\n",
        "receptor\tREC2\t\tlr_pairs\n",
        "subunit\tSUBA\tCX1\tcofactors\n",
        "subunit\tSUBB\tCX1\tcofactors\n",
    );
    assert_eq!(genes, expected_genes);

    let resolved = fs::read_to_string(stage_dir.join("resolved.json")).expect("read resolved");
    assert!(resolved.contains("\"level\": \"scalar\""));
    assert!(resolved.contains("\"n_cells_barcodes\": 4"));
    assert!(resolved.contains("\"n_rows_groups\": 4"));
    assert!(resolved.contains("\"n_groups\": 2"));
    assert!(resolved.contains("\"n_lr_pairs_raw\": 2"));
}

#[test]
fn stage0_accepts_gz_barcodes_and_groups() {
    let td = tempdir().expect("tempdir");
    let root = td.path();

    let resources = root.join("resources");
    fs::create_dir_all(&resources).expect("resources");
    fs::write(
        resources.join("lr_pairs.tsv"),
        "ligand_symbol\treceptor_symbol\tfamily\tdirectionality\nLIG1\tREC1\tF1\tbidirectional\n",
    )
    .expect("lr_pairs");
    fs::write(
        resources.join("cofactors.tsv"),
        "complex_id\trole\tsubunit_symbol\trequired\tlogic\n",
    )
    .expect("cofactors");
    fs::write(
        resources.join("gene_alias.tsv"),
        "alias\tsymbol\tpriority\tnotes\nLIG1\tLIG1\t1\t-\nREC1\tREC1\t1\t-\n",
    )
    .expect("aliases");
    fs::write(
        resources.join("labels.tsv"),
        "ligand_symbol\treceptor_symbol\tlabel\tweight\tnotes\n*\t*\tdefault\t1\t-\n",
    )
    .expect("labels");

    write_gz(&root.join("barcodes.tsv.gz"), "cell_id\nC1\nC2\nC3\nC4\n");
    write_gz(
        &root.join("groups.tsv.gz"),
        "cell_id\tgroup\nC1\tG1\nC2\tG1\nC3\tG2\nC4\tG2\n",
    );
    fs::write(root.join("expr.bin"), b"").expect("expr");

    let out = root.join("out");
    let result = run(RunConfig {
        expr: root.join("expr.bin"),
        barcodes: root.join("barcodes.tsv.gz"),
        groups: Some(root.join("groups.tsv.gz")),
        auto_groups: None,
        auto_groups_coarse: None,
        auto_groups_fine: None,
        auto_groups_anti: None,
        auto_groups_mode: AutoGroupsModeArg::Flat,
        auto_groups_eps: 1e-6,
        auto_groups_min_delta: 0.1,
        auto_groups_unknown: "unknown".to_string(),
        auto_groups_emit_scores: false,
        resources: resources.clone(),
        embedded_profile: None,
        secretion: None,
        out: out.clone(),
        validate_only: true,
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
    })
    .expect("run stage0 gz");

    assert!(
        result
            .stage0
            .out_stage_dir
            .join("groups_normalized.tsv")
            .exists()
    );
}

fn write_gz(path: &std::path::Path, content: &str) {
    let file = fs::File::create(path).expect("create gz");
    let mut gz = flate2::write::GzEncoder::new(file, flate2::Compression::default());
    gz.write_all(content.as_bytes()).expect("write gz");
    gz.finish().expect("finish gz");
}

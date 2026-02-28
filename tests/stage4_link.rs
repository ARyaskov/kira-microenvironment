use kira_environment::stages::stage4_link::{Stage4Config, run_stage4};
use std::fs;
use tempfile::tempdir;

#[test]
fn stage4_link_outputs_are_deterministic() {
    let td = tempdir().expect("tempdir");
    let root = td.path();
    let out = root.join("out");

    let root_out = out.join("kira-microenvironment");
    let s0 = root_out.join("stage0_resolve");
    fs::create_dir_all(&s0).expect("mkdir s0");

    fs::write(
        root_out.join("edges.tsv"),
        concat!(
            "source_group\ttarget_group\tligand\treceptor\tscore\tL_expr\tR_expr\tcov_L\tcov_R\tspec_L\tspec_R\tflags\n",
            "A\tB\tL1\tR1\t10.000000\t1\t1\t1\t1\t1\t1\t\n",
            "B\tA\tL1\tR1\t2.000000\t1\t1\t1\t1\t1\t1\t\n",
            "A\tC\tL2\tR2\t1.000000\t1\t1\t1\t1\t1\t1\t\n",
            "B\tA\tL2\tR2\t8.000000\t1\t1\t1\t1\t1\t1\t\n",
            "C\tA\tL1\tR1\t9.000000\t1\t1\t1\t1\t1\t1\t\n",
        ),
    )
    .expect("edges");

    fs::write(
        root_out.join("summary.json"),
        "{\"counts\":{\"n_pairs_raw\":2}}\n",
    )
    .expect("summary");

    fs::write(
        s0.join("groups_normalized.tsv"),
        "cell_id\tgroup\nca\tA\ncb\tB\ncc\tC\n",
    )
    .expect("groups normalized");

    let secretion = root.join("secretion");
    fs::create_dir_all(&secretion).expect("mkdir secretion");
    fs::write(
        secretion.join("secretion_groups.tsv"),
        concat!(
            "group\tregime\tfraction\n",
            "A\tR_LOUD\t1.0\n",
            "B\tR_SILENT\t1.0\n",
            "C\tR_LOUD\t0.4\n",
            "C\tR_SILENT\t0.4\n",
            "C\tMYSTERY\t0.2\n",
        ),
    )
    .expect("secretion_groups");

    let resources = root.join("resources");
    fs::create_dir_all(&resources).expect("mkdir resources");
    fs::write(
        resources.join("regime_map.tsv"),
        concat!(
            "regime\tclass\tnotes\n",
            "R_LOUD\tLOUD\t-\n",
            "R_SILENT\tSILENT\t-\n",
        ),
    )
    .expect("regime map");

    let cfg = Stage4Config {
        out_dir: out.clone(),
        resources_dir: resources,
        secretion_dir: secretion,
        regime_map_path: None,
        loud_thresh: 0.5,
        silent_thresh: 0.5,
        min_support: 2,
        eps: 1e-6,
    };

    let r1 = run_stage4(cfg.clone()).expect("run stage4 #1");
    let loud1 = fs::read(r1.out_stage_dir.join("loudness_drivers.tsv")).expect("loud1");
    let silent1 = fs::read(r1.out_stage_dir.join("silence_drivers.tsv")).expect("silent1");
    let sum1 = fs::read(r1.out_stage_dir.join("link_summary.json")).expect("sum1");
    let root_sum1 = fs::read(root_out.join("summary.json")).expect("root sum1");

    let r2 = run_stage4(cfg).expect("run stage4 #2");
    let loud2 = fs::read(r2.out_stage_dir.join("loudness_drivers.tsv")).expect("loud2");
    let silent2 = fs::read(r2.out_stage_dir.join("silence_drivers.tsv")).expect("silent2");
    let sum2 = fs::read(r2.out_stage_dir.join("link_summary.json")).expect("sum2");
    let root_sum2 = fs::read(root_out.join("summary.json")).expect("root sum2");

    assert_eq!(loud1, loud2);
    assert_eq!(silent1, silent2);
    assert_eq!(sum1, sum2);
    assert_eq!(root_sum1, root_sum2);

    let loud = String::from_utf8(loud1).expect("loud utf8");
    let expected_loud = concat!(
        "rank\tligand\treceptor\tfold\tmean_loud\tmean_silent\tn_edges_loud\tn_edges_silent\texemplar_source_group\texemplar_target_group\texemplar_score\n",
        "1\tL1\tR1\t4.999998\t10.000000\t2.000000\t1\t1\tA\tB\t10.000000\n",
        "2\tL2\tR2\t0.125000\t1.000000\t8.000000\t1\t1\tA\tC\t1.000000\n",
    );
    assert_eq!(loud, expected_loud);

    let silent = String::from_utf8(silent1).expect("silent utf8");
    let expected_silent = concat!(
        "rank\tligand\treceptor\tfold\tmean_silent\tmean_loud\tn_edges_silent\tn_edges_loud\texemplar_source_group\texemplar_target_group\texemplar_score\n",
        "1\tL2\tR2\t7.999993\t8.000000\t1.000000\t1\t1\tB\tA\t8.000000\n",
        "2\tL1\tR1\t0.200000\t2.000000\t10.000000\t1\t1\tB\tA\t2.000000\n",
    );
    assert_eq!(silent, expected_silent);

    let root_summary = String::from_utf8(root_sum2).expect("root summary utf8");
    assert!(root_summary.contains("\"linking\""));
    assert!(root_summary.contains("\"n_groups_loud\": 1"));
    assert!(root_summary.contains("\"n_groups_silent\": 1"));
    assert!(root_summary.contains("\"note\": \"communication potential; not causal\""));
}

#[test]
fn stage4_accepts_kira_secretion_classify_tsv_in_out_dir() {
    let td = tempdir().expect("tempdir");
    let root = td.path();
    let out = root.join("out");

    let root_out = out.join("kira-microenvironment");
    let s0 = root_out.join("stage0_resolve");
    fs::create_dir_all(&s0).expect("mkdir s0");

    fs::write(
        root_out.join("edges.tsv"),
        concat!(
            "source_group\ttarget_group\tligand\treceptor\tscore\tL_expr\tR_expr\tcov_L\tcov_R\tspec_L\tspec_R\tflags\n",
            "A\tB\tL1\tR1\t10.000000\t1\t1\t1\t1\t1\t1\t\n",
            "B\tA\tL1\tR1\t2.000000\t1\t1\t1\t1\t1\t1\t\n",
        ),
    )
    .expect("edges");
    fs::write(root_out.join("summary.json"), "{\"counts\":{}}\n").expect("summary");
    fs::write(
        s0.join("groups_normalized.tsv"),
        "cell_id\tgroup\nc1\tA\nc2\tA\nc3\tB\nc4\tB\n",
    )
    .expect("groups normalized");

    let secretion_repo = root.join("kira-secretion");
    let secretion_out = secretion_repo.join("out");
    fs::create_dir_all(&secretion_out).expect("mkdir secretion out");
    fs::write(
        secretion_out.join("classify.tsv"),
        "cell_id\tregime\trule_id\tflags\nc1\tR_LOUD\tr\t.\nc2\tR_LOUD\tr\t.\nc3\tR_SILENT\tr\t.\nc4\tR_SILENT\tr\t.\n",
    )
    .expect("classify");

    let resources = root.join("resources");
    fs::create_dir_all(&resources).expect("mkdir resources");
    fs::write(
        resources.join("regime_map.tsv"),
        "regime\tclass\tnotes\nR_LOUD\tLOUD\t-\nR_SILENT\tSILENT\t-\n",
    )
    .expect("regime map");

    run_stage4(Stage4Config {
        out_dir: out,
        resources_dir: resources,
        secretion_dir: secretion_repo,
        regime_map_path: None,
        loud_thresh: 0.5,
        silent_thresh: 0.5,
        min_support: 1,
        eps: 1e-6,
    })
    .expect("run stage4");

    assert!(root_out.join("stage4_link/loudness_drivers.tsv").exists());
    assert!(root_out.join("stage4_link/silence_drivers.tsv").exists());
}

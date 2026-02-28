use kira_environment::stages::stage3_network::{Stage3Config, run_stage3};
use std::fs;
use tempfile::tempdir;

#[test]
fn stage3_network_outputs_are_deterministic() {
    let td = tempdir().expect("tempdir");
    let root = td.path();
    let out = root.join("out");
    let stage2 = out.join("kira-microenvironment").join("stage2_score");
    fs::create_dir_all(&stage2).expect("mkdir stage2");

    fs::write(
        stage2.join("edges_raw.tsv"),
        concat!(
            "source_group\ttarget_group\tligand\treceptor\tscore\tL_expr\tR_expr\tcov_L\tcov_R\tspec_L\tspec_R\tflags\n",
            "B\tC\tL2\tR2\t4.000000\t1.000000\t1.000000\t0.500000\t0.500000\t1.000000\t1.000000\t\n",
            "A\tB\tL1\tR1\t5.000000\t1.000000\t1.000000\t0.500000\t0.500000\t1.000000\t1.000000\t\n",
            "C\tA\tL2\tR2\t2.000000\t1.000000\t1.000000\t0.500000\t0.500000\t1.000000\t1.000000\t\n",
            "A\tC\tL1\tR1\t1.000000\t1.000000\t1.000000\t0.500000\t0.500000\t1.000000\t1.000000\t\n",
            "B\tA\tL1\tR1\t3.000000\t1.000000\t1.000000\t0.500000\t0.500000\t1.000000\t1.000000\t\n",
        ),
    )
    .expect("edges_raw");

    fs::write(
        stage2.join("stage2_summary.json"),
        "{\n  \"counts\": {\"n_pairs_raw\": 2}\n}\n",
    )
    .expect("summary");

    let cfg = Stage3Config {
        out_dir: out.clone(),
    };

    let r1 = run_stage3(cfg.clone()).expect("stage3 first");
    let edges1 = fs::read(r1.out_root_dir.join("edges.tsv")).expect("edges1");
    let gs1 = fs::read(r1.out_root_dir.join("group_strength.tsv")).expect("gs1");
    let tp1 = fs::read(r1.out_root_dir.join("top_pairs.tsv")).expect("tp1");
    let s1 = fs::read(r1.out_root_dir.join("summary.json")).expect("s1");

    let r2 = run_stage3(cfg).expect("stage3 second");
    let edges2 = fs::read(r2.out_root_dir.join("edges.tsv")).expect("edges2");
    let gs2 = fs::read(r2.out_root_dir.join("group_strength.tsv")).expect("gs2");
    let tp2 = fs::read(r2.out_root_dir.join("top_pairs.tsv")).expect("tp2");
    let s2 = fs::read(r2.out_root_dir.join("summary.json")).expect("s2");

    assert_eq!(edges1, edges2);
    assert_eq!(gs1, gs2);
    assert_eq!(tp1, tp2);
    assert_eq!(s1, s2);

    let edges = String::from_utf8(edges1).expect("edges utf8");
    let expected_edges = concat!(
        "source_group\ttarget_group\tligand\treceptor\tscore\tL_expr\tR_expr\tcov_L\tcov_R\tspec_L\tspec_R\tflags\n",
        "A\tB\tL1\tR1\t5.000000\t1.000000\t1.000000\t0.500000\t0.500000\t1.000000\t1.000000\t\n",
        "B\tC\tL2\tR2\t4.000000\t1.000000\t1.000000\t0.500000\t0.500000\t1.000000\t1.000000\t\n",
        "B\tA\tL1\tR1\t3.000000\t1.000000\t1.000000\t0.500000\t0.500000\t1.000000\t1.000000\t\n",
        "C\tA\tL2\tR2\t2.000000\t1.000000\t1.000000\t0.500000\t0.500000\t1.000000\t1.000000\t\n",
        "A\tC\tL1\tR1\t1.000000\t1.000000\t1.000000\t0.500000\t0.500000\t1.000000\t1.000000\t\n",
    );
    assert_eq!(edges, expected_edges);

    let group_strength = String::from_utf8(gs1).expect("gs utf8");
    let expected_group_strength = concat!(
        "group\tout_strength\tin_strength\ttop_targets\ttop_sources\n",
        "A\t6.000000\t5.000000\tB:5.000000;C:1.000000\tB:3.000000;C:2.000000\n",
        "B\t7.000000\t5.000000\tC:4.000000;A:3.000000\tA:5.000000\n",
        "C\t2.000000\t5.000000\tA:2.000000\tB:4.000000;A:1.000000\n",
    );
    assert_eq!(group_strength, expected_group_strength);

    let top_pairs = String::from_utf8(tp1).expect("tp utf8");
    let expected_top_pairs = concat!(
        "global_rank\tligand\treceptor\ttop_source_group\ttop_target_group\tscore\n",
        "1\tL1\tR1\tA\tB\t5.000000\n",
        "2\tL2\tR2\tB\tC\t4.000000\n",
    );
    assert_eq!(top_pairs, expected_top_pairs);

    let summary = String::from_utf8(s2).expect("summary utf8");
    assert!(summary.contains("\"n_pairs_raw\": 2"));
    assert!(summary.contains("\"network\""));
    assert!(summary.contains("\"n_edges\": 5"));
    assert!(summary.contains("\"n_groups\": 3"));
    assert!(summary.contains("\"out_strength_hubs\""));
    assert!(summary.contains("\"in_strength_hubs\""));
    assert!(summary.contains("\"top_edges_by_score\""));
    assert!(summary.contains("\"source_group\": \"A\""));
    assert!(summary.contains("\"ligand\": \"L1\""));
}

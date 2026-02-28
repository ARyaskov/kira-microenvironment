use kira_environment::stages::stage2_score::{Stage2Config, run_stage2};
use std::fs;
use tempfile::tempdir;

#[test]
fn stage2_scores_are_correct_and_deterministic() {
    let td = tempdir().expect("tempdir");
    let root = td.path();

    let out = root.join("out");
    let stage1 = out.join("kira-microenvironment").join("stage1_agg");
    fs::create_dir_all(&stage1).expect("mkdir stage1");

    write_stage1_cache(&stage1.join("cache.bin")).expect("write cache");
    fs::write(
        stage1.join("gene_stats.tsv"),
        concat!(
            "gene\tmean_over_groups\tmissing_fraction\tcap_used\n",
            "G1\t1.250000\t0.000000\t1.500000\n",
            "G2\t0.750000\t0.000000\t1.500000\n",
            "G3\t0.750000\t0.000000\t1.500000\n",
        ),
    )
    .expect("gene_stats");
    fs::write(stage1.join("stage1_summary.json"), "{}\n").expect("summary");

    let resources = root.join("resources");
    fs::create_dir_all(&resources).expect("mkdir resources");

    fs::write(
        resources.join("lr_pairs.tsv"),
        concat!(
            "ligand_symbol\treceptor_symbol\tfamily\tdirectionality\tweight\n",
            "G1\tG2\tF\tbidirectional\t1.0\n",
            "G1\tCX1\tF\tbidirectional\t1.0\n",
        ),
    )
    .expect("lr_pairs");

    fs::write(
        resources.join("cofactors.tsv"),
        concat!(
            "complex_id\trole\tsubunit_symbol\trequired\tlogic\n",
            "CX1\treceptor\tG2\t1\tAND\n",
            "CX1\treceptor\tG3\t1\tAND\n",
        ),
    )
    .expect("cofactors");

    fs::write(
        resources.join("labels.tsv"),
        concat!(
            "ligand_symbol\treceptor_symbol\tlabel\tweight\tnotes\n",
            "*\tG2\tALERT\t2.0\t-\n",
            "G1\tG2\tIMMUNE\t1.5\t-\n",
        ),
    )
    .expect("labels");

    let cfg = Stage2Config {
        out_dir: out.clone(),
        resources_dir: resources,
        lr_profile: "full".to_string(),
        eps: 1e-6,
        cov_min: 0.10,
        expr_min: 0.05,
        spec_on: true,
        spec_cap: 10.0,
        top_n_per_pair: 200,
        top_n_per_source: 200,
    };

    let r1 = run_stage2(cfg.clone()).expect("stage2 first");
    let edges1 = fs::read(r1.out_stage_dir.join("edges_raw.tsv")).expect("edges1");
    let pairs1 = fs::read(r1.out_stage_dir.join("pairs_stats.tsv")).expect("pairs1");

    let r2 = run_stage2(cfg).expect("stage2 second");
    let edges2 = fs::read(r2.out_stage_dir.join("edges_raw.tsv")).expect("edges2");
    let pairs2 = fs::read(r2.out_stage_dir.join("pairs_stats.tsv")).expect("pairs2");

    assert_eq!(edges1, edges2);
    assert_eq!(pairs1, pairs2);

    let edges = String::from_utf8(edges1).expect("edges utf8");
    let expected_edges = concat!(
        "source_group\ttarget_group\tligand\treceptor\tscore\tL_expr\tR_expr\tcov_L\tcov_R\tspec_L\tspec_R\tflags\n",
        "G2\tG1\tG1\tG2\t10.457044\t1.500000\t1.500000\t0.500000\t1.000000\t1.199999\t1.999997\tALERT,IMMUNE\n",
        "G1\tG1\tG1\tG2\t5.692094\t1.000000\t1.500000\t0.500000\t1.000000\t0.799999\t1.999997\tALERT,IMMUNE\n",
        "G2\tG1\tG1\tCX1\t0.670820\t1.500000\t0.500000\t0.500000\t0.500000\t1.199999\t0.666666\t\n",
        "G1\tG1\tG1\tCX1\t0.365148\t1.000000\t0.500000\t0.500000\t0.500000\t0.799999\t0.666666\t\n",
    );
    assert_eq!(edges, expected_edges);

    let pairs = String::from_utf8(pairs1).expect("pairs utf8");
    let expected_pairs = concat!(
        "ligand\treceptor\tn_edges_kept\ttop_score\n",
        "G1\tCX1\t2\t0.670820\n",
        "G1\tG2\t2\t10.457044\n",
    );
    assert_eq!(pairs, expected_pairs);

    let summary =
        fs::read_to_string(r2.out_stage_dir.join("stage2_summary.json")).expect("summary");
    assert!(summary.contains("\"n_pairs_raw\": 2"));
    assert!(summary.contains("\"n_pairs_expanded\": 2"));
    assert!(summary.contains("\"n_edges_before_filter\": 8"));
    assert!(summary.contains("\"n_edges_after_filter\": 4"));
    assert!(summary.contains("\"cap_used\": 1.5"));
    assert!(summary.contains("\"missing_components\": 0"));
    assert!(summary.contains("\"missing_genes\": 0"));
}

fn write_stage1_cache(path: &std::path::Path) -> Result<(), std::io::Error> {
    let mut out = Vec::new();
    out.extend_from_slice(b"KIRAEAGG");
    out.extend_from_slice(&1u32.to_le_bytes());
    out.extend_from_slice(&2u32.to_le_bytes());
    out.extend_from_slice(&3u32.to_le_bytes());

    for g in ["G1", "G2"] {
        out.extend_from_slice(&(g.len() as u32).to_le_bytes());
        out.extend_from_slice(g.as_bytes());
    }
    for gene in ["G1", "G2", "G3"] {
        out.extend_from_slice(&(gene.len() as u32).to_le_bytes());
        out.extend_from_slice(gene.as_bytes());
    }

    let expr: [f32; 6] = [1.0, 1.5, 0.5, 1.5, 0.0, 1.0];
    let cov: [f32; 6] = [0.5, 1.0, 0.5, 0.5, 0.0, 0.5];

    for v in expr {
        out.extend_from_slice(&v.to_le_bytes());
    }
    for v in cov {
        out.extend_from_slice(&v.to_le_bytes());
    }

    fs::write(path, out)
}

use crate::metrics::microenv_extension::scores::{CellMicroenvScores, PanelCoverage};
use crate::select::{median_in_place, quantile_in_place};
use serde::Serialize;
use std::collections::BTreeMap;

#[derive(Debug, Clone, Serialize)]
pub struct MetricSummary {
    pub hsi: Option<f32>,
    pub ias: Option<f32>,
    pub iss: Option<f32>,
    pub mio: Option<f32>,
    pub sii: Option<f32>,
    pub msm: Option<f32>,
}

#[derive(Debug, Clone, Serialize)]
pub struct MetricQuantiles {
    pub median: Option<f32>,
    pub p10: Option<f32>,
    pub p90: Option<f32>,
}

#[derive(Debug, Clone, Serialize)]
pub struct ClusterMetricStats {
    pub hsi: MetricQuantiles,
    pub ias: MetricQuantiles,
    pub iss: MetricQuantiles,
    pub mio: MetricQuantiles,
    pub sii: MetricQuantiles,
    pub msm: MetricQuantiles,
}

#[derive(Debug, Clone, Serialize)]
pub struct FlagFractions {
    pub hypoxia_high: f32,
    pub inflammatory_high: f32,
    pub immune_suppression_high: f32,
    pub metabolic_suppression_high: f32,
    pub stromal_high: f32,
    pub microenv_stress_mode: f32,
}

#[derive(Debug, Clone, Serialize)]
pub struct ClusterStatsRow {
    pub cluster: String,
    pub n_cells: usize,
    pub metrics: ClusterMetricStats,
    pub flags: FlagFractions,
}

#[derive(Debug, Clone, Serialize)]
pub struct TopClusterByStress {
    pub cluster: String,
    pub n_cells: usize,
    pub median_msm: Option<f32>,
    pub fraction_microenv_stress_mode: f32,
}

#[derive(Debug, Clone, Serialize)]
pub struct GlobalStats {
    pub medians: MetricSummary,
    pub mads: MetricSummary,
    pub top_clusters_by_microenv_stress_mode: Vec<TopClusterByStress>,
}

#[derive(Debug, Clone, Serialize)]
pub struct MetricMissingness {
    pub metric: String,
    pub n_nan: usize,
    pub n_total: usize,
    pub fraction_nan: f32,
}

#[derive(Debug, Clone, Serialize)]
pub struct Missingness {
    pub panels: Vec<PanelCoverage>,
    pub metrics: Vec<MetricMissingness>,
}

pub fn build_global_stats(
    scores: &CellMicroenvScores,
    cluster_stats: &[ClusterStatsRow],
) -> GlobalStats {
    let medians = MetricSummary {
        hsi: median_finite(&scores.hsi),
        ias: median_finite(&scores.ias),
        iss: median_finite(&scores.iss),
        mio: median_finite(&scores.mio),
        sii: median_finite(&scores.sii),
        msm: median_finite(&scores.msm),
    };
    let mads = MetricSummary {
        hsi: mad_finite(&scores.hsi),
        ias: mad_finite(&scores.ias),
        iss: mad_finite(&scores.iss),
        mio: mad_finite(&scores.mio),
        sii: mad_finite(&scores.sii),
        msm: mad_finite(&scores.msm),
    };

    let mut top = cluster_stats
        .iter()
        .map(|c| TopClusterByStress {
            cluster: c.cluster.clone(),
            n_cells: c.n_cells,
            median_msm: c.metrics.msm.median,
            fraction_microenv_stress_mode: c.flags.microenv_stress_mode,
        })
        .collect::<Vec<_>>();
    top.sort_by(|a, b| {
        b.median_msm
            .unwrap_or(f32::NEG_INFINITY)
            .total_cmp(&a.median_msm.unwrap_or(f32::NEG_INFINITY))
            .then(a.cluster.cmp(&b.cluster))
    });
    if top.len() > 10 {
        top.truncate(10);
    }

    GlobalStats {
        medians,
        mads,
        top_clusters_by_microenv_stress_mode: top,
    }
}

pub fn build_cluster_stats(
    scores: &CellMicroenvScores,
    clusters: &[String],
) -> Vec<ClusterStatsRow> {
    let mut by_cluster: BTreeMap<String, Vec<usize>> = BTreeMap::new();
    for (idx, cluster) in clusters.iter().enumerate() {
        by_cluster.entry(cluster.clone()).or_default().push(idx);
    }

    let mut out = Vec::with_capacity(by_cluster.len());
    for (cluster, idxs) in by_cluster {
        let n = idxs.len().max(1);
        let metrics = ClusterMetricStats {
            hsi: quantiles_for(&scores.hsi, &idxs),
            ias: quantiles_for(&scores.ias, &idxs),
            iss: quantiles_for(&scores.iss, &idxs),
            mio: quantiles_for(&scores.mio, &idxs),
            sii: quantiles_for(&scores.sii, &idxs),
            msm: quantiles_for(&scores.msm, &idxs),
        };
        let flags = FlagFractions {
            hypoxia_high: flag_fraction(&scores.hypoxia_high, &idxs, n),
            inflammatory_high: flag_fraction(&scores.inflammatory_high, &idxs, n),
            immune_suppression_high: flag_fraction(&scores.immune_suppression_high, &idxs, n),
            metabolic_suppression_high: flag_fraction(&scores.metabolic_suppression_high, &idxs, n),
            stromal_high: flag_fraction(&scores.stromal_high, &idxs, n),
            microenv_stress_mode: flag_fraction(&scores.microenv_stress_mode, &idxs, n),
        };
        out.push(ClusterStatsRow {
            cluster,
            n_cells: idxs.len(),
            metrics,
            flags,
        });
    }
    out
}

pub fn build_missingness(
    scores: &CellMicroenvScores,
    panel_coverage: Vec<PanelCoverage>,
) -> Missingness {
    let n_total = scores.n_cells();
    let metrics = vec![
        missing_metric("hyp_core", &scores.hyp_core, n_total),
        missing_metric("nfkb_core", &scores.nfkb_core, n_total),
        missing_metric("ifn_core", &scores.ifn_core, n_total),
        missing_metric("checkpoint_core", &scores.checkpoint_core, n_total),
        missing_metric("adenosine_core", &scores.adenosine_core, n_total),
        missing_metric("stromal_core", &scores.stromal_core, n_total),
        missing_metric("HSI", &scores.hsi, n_total),
        missing_metric("IAS", &scores.ias, n_total),
        missing_metric("ISS", &scores.iss, n_total),
        missing_metric("MIO", &scores.mio, n_total),
        missing_metric("SII", &scores.sii, n_total),
        missing_metric("MSM", &scores.msm, n_total),
    ];

    Missingness {
        panels: panel_coverage,
        metrics,
    }
}

fn quantiles_for(values: &[f32], idxs: &[usize]) -> MetricQuantiles {
    let v = idxs
        .iter()
        .map(|&i| values[i])
        .filter(|x| x.is_finite())
        .collect::<Vec<_>>();
    if v.is_empty() {
        return MetricQuantiles {
            median: None,
            p10: None,
            p90: None,
        };
    }

    let mut v_median = v.clone();
    let median = median_in_place(&mut v_median);
    let mut v_p10 = v.clone();
    let p10 = quantile_in_place(&mut v_p10, 0.10);
    let mut v_p90 = v;
    let p90 = quantile_in_place(&mut v_p90, 0.90);
    MetricQuantiles {
        median: Some(median),
        p10: Some(p10),
        p90: Some(p90),
    }
}

fn median_finite(values: &[f32]) -> Option<f32> {
    let mut v = values
        .iter()
        .copied()
        .filter(|x| x.is_finite())
        .collect::<Vec<_>>();
    if v.is_empty() {
        None
    } else {
        Some(median_in_place(&mut v))
    }
}

fn mad_finite(values: &[f32]) -> Option<f32> {
    let mut finite = values
        .iter()
        .copied()
        .filter(|x| x.is_finite())
        .collect::<Vec<_>>();
    if finite.is_empty() {
        return None;
    }
    let med = median_in_place(&mut finite);
    let mut dev = values
        .iter()
        .copied()
        .filter(|x| x.is_finite())
        .map(|x| (x - med).abs())
        .collect::<Vec<_>>();
    if dev.is_empty() {
        None
    } else {
        Some(median_in_place(&mut dev))
    }
}

fn flag_fraction(flags: &[bool], idxs: &[usize], denom: usize) -> f32 {
    if denom == 0 {
        return 0.0;
    }
    let count = idxs.iter().filter(|&&i| flags[i]).count();
    (count as f32) / (denom as f32)
}

fn missing_metric(name: &str, values: &[f32], n_total: usize) -> MetricMissingness {
    let n_nan = values.iter().filter(|x| !x.is_finite()).count();
    MetricMissingness {
        metric: name.to_string(),
        n_nan,
        n_total,
        fraction_nan: if n_total == 0 {
            0.0
        } else {
            (n_nan as f32) / (n_total as f32)
        },
    }
}

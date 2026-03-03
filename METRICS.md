# kira-microenvironment Metrics Specification

This document defines canonical metrics, formulas, and constants used by `kira-microenvironment`.

Scope:
- optional Auto-groups assignment metrics
- Stage1 aggregation metrics
- Stage5 microenvironment extension metrics (transcriptional proxies)
- Stage2 LR scoring metrics
- Stage3 network metrics
- Stage4 secretion-linking metrics

## Canonical Conventions

1. Determinism
- No stochastic steps.
- Stable tie-breaking and sort keys are explicit.
- TSV numeric output formatting uses fixed precision: `6` decimals.

2. Axis semantics
- Stage1/Stage2 scores are computed per:
  - `group` (for expression/covariate metrics), and
  - ordered pair `(source_group, target_group)` (for edge metrics).
- Stage4 classifies source groups by regime fractions.

3. Value domains
- Expression values are non-negative in scoring (`clamp_expr`).
- Coverage metrics (`cov`) are in `[0, 1]`.
- Specificity metrics (`spec_L`, `spec_R`) are in `[0, spec_cap]` when enabled.
- Fold metrics in Stage4 are strictly positive due to epsilon smoothing.

## Default Constants (CLI)

From `src/cli.rs`:
- `agg=median`
- `trim=0.05` (used only for `trimmed_mean`)
- `eps=1e-6`
- `cap_mode=p99`
- `cap_p=0.99`
- `cov_min=0.10`
- `expr_min=0.05`
- `spec=on`
- `spec_cap=10.0`
- `top_n_per_pair=200`
- `top_n_per_source=200`
- `loud_thresh=0.5`
- `silent_thresh=0.5`
- `min_support=3`

Internal constants:
- Stage2 warning threshold: `COV_WARN=0.20`

## Notation

Let:
- `E[g, x]` be aggregated expression of gene/entity `x` in group `g` (Stage1 output).
- `Cov[g, x]` be coverage fraction in group `g`: cells with `expr > eps`.
- `Mean[x]` be `mean_over_groups` from `gene_stats.tsv`.
- `cap` be `cap_used` from Stage1.
- `src`, `tgt` be source/target groups.

For LR pair `p=(L,R)`:
- `w_pair` = LR pair weight from `lr_pairs.tsv` (`1.0` if absent/unparseable).
- `w_label` = product of matched label weights from `labels.tsv`.

Helper functions:
- `clamp_expr(x, cap) = min(max(x, 0), cap)`
- `clamp_ratio(x, m, eps, cap) = clamp(x / (m + eps), 0, cap)`

## Auto-groups Metrics (Optional)

Used only with `--auto-groups*`.

Per-cell group score:
- marker contribution:
  - `score_raw(group) = sum_i(expr_i * marker_weight_i) / n_markers`
- anti-marker penalty (if configured):
  - `penalty(group) = sum_j(expr_j * anti_penalty_j) / n_anti_markers`
- final:
  - `score(group) = max(score_raw - penalty, 0)`

Assignment rule:
- sort group scores descending (tie -> group name ascending).
- assign best group only if:
  - `best_score >= auto_groups_eps`, and
  - `(best_score - second_score) >= auto_groups_min_delta`.
- otherwise assign `unknown_label`.

Defaults:
- `auto_groups_eps=1e-6`
- `auto_groups_min_delta=0.1`
- `auto_groups_unknown="unknown"`

Summary counts:
- `n_unknown`, `n_coarse_only`, `n_fine`, `n_cells`.

## Stage1 Aggregation Metrics

For each `(group, gene)`:
- collect one value per cell (missing values padded by `0` to `n_cells`).
- aggregation mode:
  - `median`: median of group cell values.
  - `trimmed_mean`: after sorting, trim `k=floor(n*trim)` values from each side where `trim` is clamped to `[0, 0.499999]`, then average remaining.
- coverage:
  - `cov = (#cells with value > eps) / n_cells`

Cap metric:
- if `cap_mode=fixed`: `cap = cap_fixed`
- if `cap_mode=p99`: `cap = quantile(expr_matrix, cap_p)` with rank:
  - `rank = ceil(cap_p * n) - 1`, clamped to `[0, n-1]`

Gene statistics:
- `mean_over_groups(gene) = average_g E[g, gene]`
- `missing_fraction = 1.0` if gene absent in expression cache, else `0.0`

## Stage5 Microenvironment Extension Metrics

All Stage5 metrics are transcriptional proxies computed from one normalized scRNA-seq matrix.
They do not measure oxygen, cytokines, checkpoint protein levels, or metabolites directly.

Panel version:
- `MICROENV_EXTENSION_PANEL_V1`

Panel cores:
- For panel `P` and cell `c`, core is `TM(P,c)`:
  - collect panel-gene expression values for `c` (zeros for absent entries in sparse storage)
  - compute trimmed mean with trim fraction `0.10`
  - if resolved genes in panel `< MIN_PANEL_GENES (2)`, core is `NaN` for all cells

Robust normalization:
- For core series `S(c)`, compute:
  - `Z(c) = (S(c) - median(S)) / (1.4826 * MAD(S) + eps)`
  - `eps = 1e-6`
  - if `MAD == 0`, all finite `Z(c)` are set to `0`

Core metrics:
- `hyp_core`: HIF/hypoxia panel core
- `nfkb_core`: inflammatory (NF-kB) panel core
- `ifn_core`: interferon response panel core
- `checkpoint_core`: immune suppression/checkpoint panel core
- `adenosine_core`: adenosine metabolic suppression panel core
- `stromal_core`: ECM/stromal interaction panel core

Derived scores:
- `HSI = Z(hyp_core)`
- `IAS = 0.6 * Z(nfkb_core) + 0.4 * Z(ifn_core)`
- `ISS = Z(checkpoint_core)`
- `MIO = 0.5 * Z(adenosine_core) + 0.5 * max(0, HSI)`
- `SII = Z(stromal_core)`
- `MSM = max(0, 0.4 * max(0, HSI) + 0.3 * max(0, ISS) + 0.3 * max(0, MIO))`

Flags:
- `hypoxia_high`: `HSI >= 2.0`
- `inflammatory_high`: `IAS >= 2.0`
- `immune_suppression_high`: `ISS >= 1.5`
- `metabolic_suppression_high`: `MIO >= 1.5`
- `stromal_high`: `SII >= 2.0`
- `microenv_stress_mode`: `MSM >= 2.0`

NaN handling:
- if any required core/score input is `NaN`, downstream score is `NaN`
- flags are `false` for `NaN` scores
- missingness is recorded in Stage5 summary (`panels`, per-metric `fraction_nan`)

## Stage2 LR Edge Scoring Metrics

### Entity evaluation

Single-gene entity:
- `expr = E[g, gene]`
- `cov = Cov[g, gene]`
- `mean = Mean[gene]`

Complex entity (from `cofactors.tsv`):
- `AND` components aggregate by `min`
- `OR` components aggregate by `max`
- if any `AND` components exist: only `AND` aggregate is used
- otherwise `OR` aggregate is used
- if unresolved/incomplete in a required way -> pair skipped (`missing_components` / `missing_genes`)

### Pair label effects

Label rule matching uses exact symbol or wildcard `*`.

For matched rules:
- `w_label = Π weight_i`
- all matched labels are appended to `flags`

### Pre-filter metrics

For each resolved LR pair and each ordered `(src, tgt)`:
- ligand side at source: `(L_expr, cov_L, mean_L)`
- receptor side at target: `(R_expr, cov_R, mean_R)`

Edge is kept only if:
- `cov_L >= cov_min`
- `cov_R >= cov_min`
- `L_expr >= expr_min`
- `R_expr >= expr_min`

### Score formula

Base score:
- `L_cap = clamp_expr(L_expr, cap)`
- `R_cap = clamp_expr(R_expr, cap)`
- `score0 = L_cap * R_cap * w_pair * w_label`

Specificity (if `spec=on`):
- `spec_L = clamp_ratio(L_expr, mean_L, eps, spec_cap)`
- `spec_R = clamp_ratio(R_expr, mean_R, eps, spec_cap)`
- `spec_boost = sqrt(spec_L * spec_R)`

If `spec=off`: `spec_L=1`, `spec_R=1`, `spec_boost=1`.

Final edge score:
- `score = score0 * spec_boost`

### Edge flags

- add `LOW_COVERAGE_EDGE` when `min(cov_L, cov_R) < 0.20`
- add `HIGH_SPECIFICITY` when `spec=on` and `spec_boost >= 3.0`
- add matched label names
- deduplicate + sort flag tokens, then join with comma

### Ranking and pruning

1. Per `(ligand, receptor, source_group)` keep top `top_n_per_pair` by:
- score desc, then `target_group` asc.

2. Per `source_group` keep top `top_n_per_source` by:
- score desc, then `ligand` asc, `receptor` asc, `target_group` asc.

3. Global sort in output `edges_raw.tsv`:
- score desc, then `ligand`, `receptor`, `source_group`, `target_group` asc.

### Pair summary metrics (`pairs_stats.tsv`)

For each `(ligand, receptor)`:
- `n_edges_kept` = number of final kept edges.
- `top_score` = maximum kept edge score.

## Stage3 Network Metrics

Input: Stage2 `edges_raw.tsv`; output root-level network tables.

### Group strength

For each group `g`:
- `out_strength(g) = Σ score(e)` for edges with `source_group=g`
- `in_strength(g) = Σ score(e)` for edges with `target_group=g`

Top partner lists:
- `top_targets`: top 5 targets by summed score from `g`
- `top_sources`: top 5 sources by summed score into `g`
- sorting: summed score desc, then group name asc
- string format: `group:score;group:score;...`

### Top pairs

For each `(ligand, receptor)` choose one exemplar edge:
- highest score; tie -> lexicographically smallest `(source_group, target_group)`.

`top_pairs.tsv` sorting:
- score desc, then ligand asc, receptor asc.
- `global_rank` is 1-based index in that order.

### Root summary additions

In `summary.json`:
- `network.n_edges`
- `network.n_groups`
- `network.out_strength_hubs`: top 5 groups by `out_strength`
- `network.in_strength_hubs`: top 5 groups by `in_strength`
- `top_edges_by_score`: top 10 edges by Stage3 global edge sort

## Stage4 Secretion Linking Metrics (Optional)

Executed only with `--secretion`.

### Group classification from regime fractions

Given per-group regime fractions and `regime_map.tsv` (`LOUD|SILENT|MIXED|IGNORE`):
- `loud_frac(group) = Σ fractions of regimes mapped to LOUD`
- `silent_frac(group) = Σ fractions of regimes mapped to SILENT`

Class rule:
- `Loud` if `loud_frac >= loud_thresh`
- else `Silent` if `silent_frac >= silent_thresh`
- else `Mixed`

Unknown regime names (not present in map) are collected into `unknown_regimes` and treated as `IGNORE`.

### Pair driver metrics

Aggregation by LR pair over source-group class:
- `mean_loud = sum_loud / n_loud` (or `0` if `n_loud=0`)
- `mean_silent = sum_silent / n_silent` (or `0` if `n_silent=0`)

Support filter:
- only keep pairs with `n_loud + n_silent >= min_support`

Fold metrics with smoothing epsilon:
- `fold_loud = (mean_loud + eps) / (mean_silent + eps)`
- `fold_silent = (mean_silent + eps) / (mean_loud + eps)`

Ranking in both driver tables:
- fold desc, then primary mean desc, then ligand asc, receptor asc.
- ranks are 1-based.

Exemplar edge per class:
- highest edge score for that LR pair/class; tie -> smallest `(source_group, target_group)`.

### Output truncation in summary

`link_summary.json`:
- `top_loud_pairs`: top 20
- `top_silent_pairs`: top 20
- `top_loud_edges`: top 50 source-loud edges
- `top_silent_edges`: top 50 source-silent edges

Root `summary.json` linking section:
- same counts
- top 10 loud pairs
- top 10 silent pairs
- note: `"communication potential; not causal"`

# Pipeline Specification

## Stages

1. Stage0 `resolve`
2. Stage1 `agg`
3. Stage5 `microenv_extension` (per-cell transcriptional proxy metrics)
4. Stage2 `score`
5. Stage3 `network`
6. Stage4 `link` (optional; only when `--secretion` is provided)

All outputs are under `out/kira-microenvironment/`.

## Stage0 Resolve

Inputs:
- CLI inputs: `expr.bin`, `barcodes`, `groups`, `resources`, optional `secretion` path.

Outputs: `out/kira-microenvironment/stage0_resolve/`
- `resolved.json`
  - `tool`: `{name, version, rust, edition}`
  - `simd`: `{level}`
  - `inputs`: `[{kind, path_abs, size_bytes, mtime_unix, sha256_8}]`
  - `resources`: `[{kind/name, path_abs, size_bytes, mtime_unix, sha256_8}]`
  - `config`: `{lr_profile, validate_only}`
  - `counts`: `{n_cells_barcodes, n_rows_groups, n_groups, n_lr_pairs_raw}`
  - `warnings`: `[string]`
- `genes_requested.tsv`
  - columns: `role, query_symbol, complex_id, source`
- `groups_normalized.tsv`
  - columns: `cell_id, group`
- `group_sizes.tsv`
  - columns: `group, n_cells`

## Stage1 Agg

Inputs:
- Stage0: `genes_requested.tsv`, `groups_normalized.tsv`, `group_sizes.tsv`, `resolved.json`
- `expr.bin` (KIRAEXPR v1)

Outputs: `out/kira-microenvironment/stage1_agg/`
- `group_gene_agg.tsv`
  - columns: `group, gene, expr, cov, n_cells`
- `gene_stats.tsv`
  - columns: `gene, mean_over_groups, missing_fraction, cap_used`
- `cache.bin` (KIRAEAGG v1)
- `stage1_summary.json`
  - `counts: {n_groups, n_genes_requested, n_genes_found}`
  - `agg: {mode, trim}`
  - `eps`
  - `cap: {mode, value}`
  - `missing_genes`

## Stage2 Score

Inputs:
- Stage1: `cache.bin`, `gene_stats.tsv`, `stage1_summary.json`
- Resources: `lr_pairs.tsv` (or selected profile), `cofactors.tsv`, `labels.tsv`

Outputs: `out/kira-microenvironment/stage2_score/`
- `edges_raw.tsv`
  - columns: `source_group, target_group, ligand, receptor, score, L_expr, R_expr, cov_L, cov_R, spec_L, spec_R, flags`
- `pairs_stats.tsv`
  - columns: `ligand, receptor, n_edges_kept, top_score`
- `stage2_summary.json`
  - `counts: {n_pairs_raw, n_pairs_expanded, n_edges_before_filter, n_edges_after_filter}`
  - `thresholds: {cov_min, expr_min, spec_on}`
  - `cap_used`
  - `skipped: {missing_components, missing_genes}`

## Stage5 Microenvironment Extension

Inputs:
- Stage0: `groups_normalized.tsv`, `resolved.json`
- `expr.bin` (same matrix as Stage1)

Outputs:
- Root: `out/kira-microenvironment/metrics.tsv`
  - columns:
    - `cell_id, group`
    - `hyp_core, nfkb_core, ifn_core, checkpoint_core, adenosine_core, stromal_core`
    - `HSI, IAS, ISS, MIO, SII, MSM`
    - `hypoxia_high, inflammatory_high, immune_suppression_high, metabolic_suppression_high, stromal_high, microenv_stress_mode`
- Stage5: `out/kira-microenvironment/stage5_microenv_extension/stage5_summary.json`
  - `panel_version`
  - `thresholds`
  - `global_stats`
  - `cluster_stats`
  - `missingness`

Root `summary.json` additive block:
- `microenvironment_extension` (same payload shape as `stage5_summary.json`)

## Stage3 Network

Inputs:
- Stage2: `edges_raw.tsv`, `stage2_summary.json` (`pairs_stats.tsv` optional informational)

Outputs (root): `out/kira-microenvironment/`
- `edges.tsv`
  - same columns as `edges_raw.tsv`
- `group_strength.tsv`
  - columns: `group, out_strength, in_strength, top_targets, top_sources`
  - `top_targets` / `top_sources` format: `group:score;group:score;...`
- `top_pairs.tsv`
  - columns: `global_rank, ligand, receptor, top_source_group, top_target_group, score`
- `summary.json`
  - carries Stage2 summary content and Stage3 additions (`network`, `top_edges_by_score`)

## Stage4 Link (Optional)

Inputs:
- Root: `edges.tsv`, `summary.json`
- Stage0: `stage0_resolve/groups_normalized.tsv`
- Secretion dir:
  - per-group candidate files (priority):
    - `secretion_groups.tsv`
    - `groups.tsv`
    - `regime_fractions.tsv`
  - per-cell candidate files (priority):
    - `secretion_cells.tsv`
    - `cells.tsv`
    - `stage6_classify.tsv`
- `resources/regime_map.tsv` (or `--regime-map` override)

Outputs: `out/kira-microenvironment/stage4_link/`
- `loudness_drivers.tsv`
  - columns: `rank, ligand, receptor, fold, mean_loud, mean_silent, n_edges_loud, n_edges_silent, exemplar_source_group, exemplar_target_group, exemplar_score`
- `silence_drivers.tsv`
  - columns: `rank, ligand, receptor, fold, mean_silent, mean_loud, n_edges_silent, n_edges_loud, exemplar_source_group, exemplar_target_group, exemplar_score`
- `link_summary.json`
  - `config: {loud_thresh, silent_thresh, min_support, eps}`
  - `counts: {n_groups_total, n_groups_loud, n_groups_silent, n_groups_mixed, n_pairs_scored}`
  - `unknown_regimes: [string]`
  - `top_loud_pairs` (first 20)
  - `top_silent_pairs` (first 20)
  - `top_loud_edges` (top 50)
  - `top_silent_edges` (top 50)

Also updates root `summary.json` key `linking`:
- `counts`
- `top_loud_pairs` (top 10)
- `top_silent_pairs` (top 10)
- `note: "communication potential; not causal"`

## Binary Formats

### KIRAEXPR v1 (`expr.bin`)
Header:
- `magic[8] = "KIRAEXPR"`
- `u32 version`
- `u32 n_cells`
- `u32 n_genes`
- `u64 gene_symbols_offset`
- `u64 cell_index_offset`
- `u64 nnz_offset`

Blocks:
- Gene symbols: `n_genes * (u32 len + utf8 bytes)`
- Cell index: `(n_cells + 1) * u64` offsets into NNZ arrays
- NNZ arrays:
  - `gene_idx[u32; nnz]`
  - `value[f32; nnz]`

### KIRAEAGG v1 (`stage1_agg/cache.bin`)
Header:
- `magic[8] = "KIRAEAGG"`
- `u32 version`
- `u32 n_groups`
- `u32 n_genes`

String tables:
- groups: `n_groups * (u32 len + utf8)`
- genes: `n_genes * (u32 len + utf8)`

Arrays (row-major, contiguous):
- `expr[f32; n_groups*n_genes]`
- `cov[f32; n_groups*n_genes]`

Index: `idx = group_idx * n_genes + gene_idx`.

## Formatting and Validation

- All emitted TSV floats use fixed precision formatting (`6` decimals) via shared helpers (`src/io/format.rs`).
- TSV parsers validate expected headers and row column counts.
- Output file writes are atomic across stages (`src/io/atomic.rs`).

## Error Codes / Exit Status

- `E_INPUT_MISSING` -> `2`
- `E_TSV_HEADER` -> `3`
- `E_TSV_PARSE` -> `4`
- `E_PATH` -> `5`
- `E_RESOURCES_INCOMPLETE` -> `6`
- `E_GROUPS_FORMAT` -> `7`
- `E_SECRETION_INPUT_MISSING` -> `8`
- `E_REGIME_MAP_MISSING` -> `9`
- `E_REGIME_MAP_PARSE` -> `10`

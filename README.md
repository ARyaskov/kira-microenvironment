# kira-microenvironment

`kira-microenvironment` builds a deterministic directed interaction-potential graph over groups: `source_group -> ligand -> receptor -> target_group`. It resolves LR resources, aggregates expression by group, scores pairwise group interactions, summarizes network structure, and can optionally link with `kira-secretion` regimes. Reported scores are communication potential, not causality.

## Build requirements

- Rust >= 1.95

## Install

Install from crates.io:

```bash
cargo install kira-microenvironment
```

## Inputs

- `expr.bin` (KIRAEXPR v1 binary): expression cache used by Stage1 mmap reader. In normal workflows this is produced by upstream expression cache generation.
- `barcodes.tsv` / `cell_ids.tsv`: cell IDs.
- `groups.tsv`: `cell_id` plus one of `group | cluster_id | cell_type`.
- `resources/` directory:
  - `lr_pairs.tsv`
  - `cofactors.tsv`
  - `gene_alias.tsv`
  - `labels.tsv`
  - `regime_map.tsv` (for Stage4 linking)
- If `--resources` path is missing, the binary falls back to an embedded default
  resource bundle and materializes it under `--out/.kira-microenvironment-embedded-resources`.
- Embedded fallback profile can be selected with `--embedded-profile onco|immune`
  or `KIRA_MICROENV_EMBEDDED_PROFILE=onco|immune` (CLI has priority).
- Optional `--secretion <dir>` for Stage4 linking with kira-secretion outputs.

## Usage examples

Base usage (auto-groups + secretion linking):

```bash
cargo run -- run --expr data/br-t/expr.bin --barcodes data/br-t/barcodes.tsv.gz --auto-groups-coarse marker_panels_coarse.tsv --auto-groups-fine marker_panels.tsv --auto-groups-anti anti_markers.tsv --auto-groups-mode hierarchical --auto-groups-eps 0.05 --auto-groups-min-delta 0.2 --resources ./resources_onco --lr-profile ./resources/lr_pairs_mvp_onco.tsv --out ./out_env --secretion ./data/br-t/
```

Validate-only (Stage0 only):

```bash
kira-microenvironment validate \
  --expr ./expr.bin \
  --barcodes ./barcodes.tsv \
  --groups ./groups.tsv \
  --resources ./resources \
  --out ./out
```

Run with MVP LR profile:

```bash
kira-microenvironment run \
  --expr ./expr.bin \
  --barcodes ./barcodes.tsv \
  --groups ./groups.tsv \
  --resources ./resources \
  --out ./out \
  --lr-profile mvp
```

Run with fixed cap:

```bash
kira-microenvironment run \
  --expr ./expr.bin \
  --barcodes ./barcodes.tsv \
  --groups ./groups.tsv \
  --resources ./resources \
  --out ./out \
  --cap-mode fixed \
  --cap-fixed 2.5
```

Run with secretion linking (Stage4):

```bash
kira-microenvironment run \
  --expr ./expr.bin \
  --barcodes ./barcodes.tsv \
  --groups ./groups.tsv \
  --resources ./resources \
  --secretion ./out/kira-secretion \
  --out ./out
```

## Outputs Overview

Root outputs in `out/kira-microenvironment/`:
- `edges.tsv`
- `group_strength.tsv`
- `top_pairs.tsv`
- `summary.json`
- Stage folders:
  - `stage0_resolve/`
  - `stage1_agg/`
  - `stage2_score/`
  - `stage4_link/` (only when `--secretion` is provided)

## Interpreting Outputs

- `score` (edges): interaction potential magnitude after expression, coverage, and optional specificity scaling.
- `cov_L` / `cov_R`: fraction of cells in group with expression above `eps` for ligand/receptor-side entities.
- `spec_L` / `spec_R`: relative expression vs mean-over-groups (if `--spec on`).
- Stage4 driver tables compare LOUD-vs-SILENT source-group means per LR pair:
  - `loudness_drivers.tsv`: large `fold` means louder in LOUD groups.
  - `silence_drivers.tsv`: large `fold` means louder in SILENT groups.

## Determinism Guarantees

- Same inputs and config produce identical outputs.
- Stable sort keys are explicit at every ranking step.
- Output writes are atomic (`write temp + rename`).
- Logging is stage-based and always includes UTC timestamps.

# PhageLysClass — Stage 08 homology-aware splitting

## v3 BioSieve compatibility correction

The Stage-08 adapter uses the current BioSieve `HomologyAwareSplitter`
parameter names observed from the installed CLI:

```yaml
homology_aware:
  mode: precomputed
  clusters_format: csv
  member_col: id
  cluster_col: cluster_id
  test_size: 0.30
  val_size: 0.0
  seed: 42
```

`random_state` is not used. MMseqs2 connected-component clustering is also run
with `--single-step-clustering 1`.


## Design

Stage 08 creates only an **outer 70/30 train/test split**.

There is no frozen validation set because the planned benchmark uses fixed/default
model hyperparameters. If an internal holdout is needed during model fitting, it
must be sampled later from the **training partition only** and must never alter
or expose the external test set.

For every Stage-07 dataset:

```text
H90 → 70/30 homology-aware train/test
H70 → 70/30 homology-aware train/test
H50 → 70/30 homology-aware train/test
H30 → 70/30 homology-aware train/test
```

## Hybrid homology strategy

The outer split combines:

1. **MMseqs2 full-task homology components**
   - same identity regime as the Stage-07 dataset;
   - `-c 0.8`;
   - `--cov-mode 0`;
   - `--cluster-mode 1` (connected components).

2. **BioSieve `homology_aware` splitting**
   - `mode: precomputed`;
   - entire MMseqs2 component stays in one partition;
   - target `test_size = 0.30`.

This is intentionally different from Stage 07. Stage 07 reduced redundancy
inside biological strata. Stage 08 clusters the **whole task together**, so
cross-class / cross-state homologs are also forced into the same outer split.

## Why connected components?

Using MMseqs2 connected components makes each detected homology component an
indivisible split group. Therefore a component cannot be shared between train
and test.

The pipeline fails if any precomputed homology cluster is present in both
partitions.

## Output

Example:

```text
data/splits/global_mixed/H70/
├── train.tsv.gz
├── test.tsv.gz
├── split_assignments.tsv.gz
├── homology_clusters.tsv.gz
├── split_summary.tsv
├── cluster_summary.tsv
└── manifest.json
```

Each task root also receives:

```text
split_overview.tsv
```

## Scripts

```text
_split_utils.py
08a_split_global_easy.py
08b_split_global_hard.py
08c_split_global_mixed.py
08d_split_subtype_binary.py
08e_split_multiclass_clean.py
08f_split_pu_global.py
```

Challenge sets are not split here.

## Internal model-development split

Later model-training code may take the outer training partition and create a
random/stratified internal holdout, for example:

```text
outer train
    ├── model-fit subset
    └── internal monitoring subset
```

That internal partition is **not** the final benchmark. The outer test partition
remains untouched until final evaluation.

No balancing, embeddings, one-hot encoding, or model training is performed by
Stage 08.

### v3 note

The BioSieve parameter file is serialized as JSON (which is valid YAML) rather
than assembled manually. This prevents literal `\n` escaping from corrupting
the configuration passed through `--params`.

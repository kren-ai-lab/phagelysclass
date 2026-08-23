# PhageLysClass training implementation — Stage 11

## Design

The trainer is representation-driven. The biological split is fixed upstream, and the same training code is reused for:

- `roxy_48`
- `onehot_1024`
- `ankh2_ext1`
- `esmc_300m`
- `esm2_t6_8m`
- `esm2_t12_35m`
- `esm2_t30_150m`
- `prot_bert`
- `prot_t5_uniref50`

One invocation corresponds to:

`task × redundancy threshold × split threshold × representation`

and internally executes LR, RF, KNN, XGBoost, and MLP unless one or more algorithms are selected explicitly.

## Important prerequisite

Before the full factorial benchmark, Stage 08 must contain all 16 combinations:

`R90/R70/R50/R30 × S90/S70/S50/S30`

The current diagonal splits can already be used for smoke tests.

## Example smoke test

```bash
python src/training_models/11a_train_supervised.py \
  --task global_easy \
  --task-type binary \
  --redundancy R90 \
  --split-threshold S90 \
  --split-file PATH_TO_R90_S90_SPLIT.tsv.gz \
  --representation roxy_48 \
  --algorithm LR
```

Then:

```bash
python src/training_models/11a_train_supervised.py \
  --task global_easy \
  --task-type binary \
  --redundancy R90 \
  --split-threshold S90 \
  --split-file PATH_TO_R90_S90_SPLIT.tsv.gz \
  --representation esm2_t6_8m \
  --algorithm all
```

And the same command later works for every PLM simply by changing `--representation`.

## Common eligibility

The trainer recomputes the frozen predictive eligibility directly from Stage 10A:

- canonical 20 AA;
- length <=1024 aa.

This guarantees that Roxy, one-hot, and PLMs are evaluated on the same biological entities.

## Frozen inner validation

For each task × R × S, the trainer creates once:

`data/training_splits/<task>/<R>_<S>/inner_seed42.tsv.gz`

It is reused by every representation and algorithm.

## Output

`results/training/<task>/<R>_<S>/<representation>/<algorithm>/`

containing validation metrics, selected configuration, final outer-test metrics,
predictions, confusion matrix, and manifest.

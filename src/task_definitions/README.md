# PhageLysClass task definitions

Copy all files in this directory to:

```text
src/task_definitions/
```

Each task script reads:

```text
data/processed/pnu_universes/pnu_sequence_master.tsv.gz
```

and creates its own folder under:

```text
data/tasks/
```

## Scripts

- `06a_global_easy.py` → `data/tasks/global_easy/`
- `06b_global_hard.py` → `data/tasks/global_hard/`
- `06c_global_mixed.py` → `data/tasks/global_mixed/`
- `06d_subtype_binary.py` → `data/tasks/subtype_binary/`
- `06e_multiclass_clean.py` → `data/tasks/multiclass_clean/`
- `06f_pu_global.py` → `data/tasks/pu_global/`
- `06g_challenge_sets.py` → `data/tasks/challenge_sets/`

Every task folder contains the materialized data, a `summary.tsv`, and a
`manifest.json`. Multi-regime scripts additionally create subfolders.

No script performs balancing, homology filtering, train/validation/test
splitting, representation calculation, or model training.

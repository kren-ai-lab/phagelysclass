# PhageLysClass — Stage 07 redundancy reduction

BioSieve/MMseqs2 redundancy reduction is applied **independently per task** and
**independently per biological stratum inside each task**.

Thresholds:

```text
H90
H70
H50
H30
```

Every threshold starts from the original Stage-06 task dataset. The reductions
are not chained.

## Why reduce within strata?

A positive sequence is not allowed to disappear merely because it clusters with
a negative sequence, and minority biological classes are not allowed to be
displaced by representatives from majority classes. Cross-stratum homology is
therefore intentionally retained at Stage 07 and is handled during the
subsequent leakage-aware splitting stage.

## Scripts

```text
_redundancy_utils.py
07a_reduce_global_easy.py
07b_reduce_global_hard.py
07c_reduce_global_mixed.py
07d_reduce_subtype_binary.py
07e_reduce_multiclass_clean.py
07f_reduce_pu_global.py
```

Challenge sets are not redundancy-reduced at this stage.

## Output layout

Example:

```text
data/redundancy/global_mixed/
├── H90/
│   ├── dataset.tsv.gz
│   ├── mapping.tsv.gz
│   ├── summary.tsv
│   └── manifest.json
├── H70/
├── H50/
├── H30/
└── redundancy_overview.tsv
```

Subtype and PU tasks have an additional subtype/regime directory.

Each `summary.tsv` reports before/after counts for each biological stratum and
for the whole task. Each `manifest.json` stores BioSieve/MMseqs2 versions,
threshold, input SHA-256, and the raw BioSieve reports returned for every
stratum.

## Preflight

The helper requires both:

```bash
biosieve
mmseqs
```

in `PATH`.

Because BioSieve strategy parameter names can differ between installed
versions, the helper probes the installed CLI once to identify the MMseqs2
identity-threshold key. The detected key is cached in each task output root.

If automatic detection is not possible, pass:

```bash
--threshold-param <name>
```

explicitly.

No balancing or train/validation/test split is performed in Stage 07.

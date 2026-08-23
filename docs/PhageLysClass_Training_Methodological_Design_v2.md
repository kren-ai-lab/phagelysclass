# PhageLysClass — Predictive Training and Factorial Analysis Methodological Contract

**Version:** 2.0  
**Document status:** implementation-ready methodological contract  
**Supersedes:** `PhageLysClass_Training_Methodological_Design_v1.md`  
**Scope:** supervised binary classification, multiclass classification, positive–unlabeled learning, uncertainty estimation, factorial performance analysis, challenge-set evaluation, and computational benchmarking  
**Project stage:** downstream of frozen evidence/task definitions, redundancy reduction, full factorial homology-aware outer splitting, and representation generation

---

## 0. Revision rationale and major changes from v1

Version 2.0 updates the predictive-training design after completion of the full redundancy × homology-split factorial stage and reframes the training stage as a reproducible experimental benchmark rather than a collection of independent model fits.

The most important changes are:

1. **Stage 08 is now fully materialized and frozen.** All 144 expected `task × R × S` cells are present and pass the split audit.
2. **Redundancy (`R`) and split stringency (`S`) remain independent experimental factors.** The complete 4 × 4 factorial surface is now the authoritative outer benchmark.
3. **Representation eligibility is separated from representation availability.** A fixed biological eligibility universe is defined independently of whether a particular representation has already been generated.
4. **Representation completeness becomes a hard gate.** Missing embeddings/features for otherwise eligible sequences are treated as representation failures, not as a reason to silently shrink the benchmark.
5. **The inner train/validation partition is now homology-group-aware.** The previous row-wise stratified 80/20 split is superseded because it could place members of the same S-level homology component in both inner train and inner validation.
6. **Experiment and representation registries become first-class artifacts.** The training engine consumes frozen registries and is not allowed to make scientific design decisions dynamically.
7. **Per-sequence predictions are mandatory outputs.** Final test probabilities/scores must be retained so that paired bootstrap, calibration analysis, common-test sensitivity analysis, error-consistency analysis, and later challenge/discovery work can be performed without retraining.
8. **Validation information is retained at configuration level.** Metrics and predictions for all predefined configurations are preserved, not only the selected winner.
9. **Paired and non-paired comparisons are explicitly separated.** Representation/algorithm comparisons within a fixed `task × R × S` cell are paired; comparisons across R, S, or alternative task definitions generally are not because the evaluated entities differ.
10. **A common-test sensitivity analysis is added.** R/S effects are evaluated both on their native frozen test sets and, when feasible, on intersections of common test entities.
11. **Runtime, memory, convergence, and failure state are part of the scientific output.**
12. **PU learning is kept as a separate analytical branch.** Its evaluation contract is strengthened to avoid treating U as a biological negative class.
13. **The paper-level framing is methodological.** The benchmark is designed to quantify how evidence definition, redundancy, homology-aware generalization, numerical representation, and learning strategy jointly affect functional recognition of phage lytic enzymes.

This document is the authoritative contract for implementation of Stages 10E onward.

---

# Part I — Scientific and data contract

## 1. Purpose

The objective of the predictive stage is not simply to identify the classifier with the highest scalar performance. The benchmark is designed to determine **why apparent performance changes**, which experimental decisions are responsible for those changes, and which modeling strategies remain reliable when sequence similarity and evidence structure become more challenging.

The central factorial experiment is:

`Task × R × S × Representation × Algorithm`

where:

- `Task` controls the biological/evidence formulation;
- `R` controls redundancy reduction before modeling;
- `S` controls homology-aware separation of outer train and test;
- `Representation` controls the numerical description of the sequence;
- `Algorithm` controls the downstream learning strategy.

The design therefore decomposes apparent predictive performance into independently inspectable components.

The broader methodological questions are:

1. How strongly does **negative-evidence definition** alter apparent functional-classification performance?
2. How much performance is lost when **sequence redundancy is progressively removed**?
3. How much performance is lost when train and test are forced to be **increasingly remote in sequence space**?
4. Do redundancy reduction and remote-generalization stringency **interact**?
5. Do protein language models outperform compact physicochemical descriptors and one-hot encoding under matched test conditions?
6. Are larger PLMs consistently better, or do gains depend on task, redundancy, and split stringency?
7. Are performance gains primarily **representation-driven or classifier-driven**?
8. Which modeling strategies remain useful under **Low-N subtype regimes**, especially VAL?
9. Can PU-specific methods rank known positives effectively without converting unlabeled proteins into false negatives?
10. What is the computational cost of each modeling strategy relative to its predictive gain?

---

## 2. Immutable upstream inputs

Training must never redefine:

- evidence state;
- P/N/U membership;
- task membership;
- subtype labels;
- redundancy membership;
- outer train/test membership;
- S-level homology components;
- challenge-set membership.

The authoritative upstream layers are:

- **Stage 05:** exact-sequence P/N/U universe;
- **Stage 06:** frozen task definitions;
- **Stage 07:** independently generated redundancy-reduced datasets;
- **Stage 08:** full factorial homology-aware outer train/test splitting;
- **Stage 10A:** exact unique-sequence representation universe;
- **Stage 10B:** Roxy descriptors;
- **Stage 10C:** Sylphy one-hot representation;
- **Stage 10D:** Sylphy PLM embeddings;
- **Stage 10E:** common predictive-eligibility projection, introduced in this document;
- **Stage 10F:** frozen group-aware inner train/validation assignment, introduced in this document.

No task label may ever be inferred from a representation table.

The exact sequence remains the identity unit. The preferred representation join key is the stable SHA-256-based `sequence_uid`. Upstream biological IDs such as `entity_id` must be retained for traceability but must not be silently assumed to be equivalent to `sequence_uid`.

---

## 3. Current Stage-08 freeze status

The full factorial outer-split stage is complete.

Current audit:

- expected factorial cells: **144**;
- observed factorial cells: **144**;
- PASS cells: **144**;
- missing/failed cells: **0**;
- final status: **PASS**;
- train/test exact-entity overlap: **0 in every cell**;
- shared S-level homology clusters across train/test: **0 in every cell**.

The observed Stage-07/08 dataset sizes are:

| Task | R90 | R70 | R50 | R30 |
|---|---:|---:|---:|---:|
| Global easy | 30,747 | 22,008 | 15,138 | 9,216 |
| Global hard | 6,907 | 5,083 | 3,661 | 2,565 |
| Global mixed | 32,641 | 23,437 | 16,244 | 10,076 |
| Endolysin | 31,760 | 22,748 | 15,691 | 9,637 |
| Depolymerase | 28,139 | 20,236 | 14,093 | 8,721 |
| VAL | 27,719 | 19,865 | 13,754 | 8,420 |
| Multiclass clean | 4,734 | 3,500 | 2,471 | 1,665 |
| PU conservative | 152,820 | 123,073 | 103,912 | 89,153 |
| PU broad | 222,235 | 175,249 | 143,736 | 117,481 |

Within a given task and R condition, `n_total` is fixed across S as expected. S changes the homology-component structure and train/test membership, not the Stage-07 entity universe.

The outer-test proportion remains approximately 0.30 across the complete factorial. Small deviations from exactly 30% are accepted because homology-component integrity has priority over exact row-level partition size.

---

## 4. Training tasks

### 4.1 Supervised global binary tasks

#### Global easy

- positive: `P_core`;
- negative: `N_easy`.

#### Global hard

- positive: `P_core`;
- negative: `N_hard`.

#### Global mixed

- positive: `P_core`;
- negative: `N_easy + N_hard`.

These three tasks form the principal negative-definition experiment.

They must not be collapsed into a single global benchmark because their biological interpretation differs.

### 4.2 Supervised subtype binary tasks

- **Endolysin**
- **Depolymerase**
- **Virion-associated lytic enzyme (VAL)**

For each subtype task:

- the positive class is the selected lytic subtype;
- the negative class is the frozen non-lytic background;
- other lytic classes are not converted into negatives.

VAL remains in the benchmark despite its Low-N positive support after stringent sequence control.

The Stage-08 reports show that VAL outer-test positive support can be as low as approximately **14–32 positives** depending on R/S. Consequently, all VAL analyses must expose absolute counts alongside rates.

### 4.3 Multiclass clean task

Classes:

- endolysin;
- depolymerase;
- virion-associated lytic enzyme.

The frozen `multiple_nonexclusive` and `other` lytic groups are excluded from clean multiclass training and retained as challenge sets.

VAL support in the multiclass test can also become extremely small and must always be reported explicitly.

### 4.4 Positive–unlabeled tasks

#### PU conservative

- P: `P_core`;
- U: `U_explicit`.

#### PU broad

- P: `P_core`;
- U: `U_explicit + U_other_unresolved`.

U means **unlabeled**. It must never be described as a confirmed biological negative class.

---

## 5. Full redundancy × split factorial design

### 5.1 Redundancy factor R

Four independently generated Stage-07 datasets are used:

- **R90**
- **R70**
- **R50**
- **R30**

Each R dataset was generated independently from the frozen task definition. R thresholds are not sequential filters.

### 5.2 Homology-aware split factor S

Each R dataset has been independently split using four S-level homology criteria:

- **S90**
- **S70**
- **S50**
- **S30**

The complete design is:

| Redundancy \ Split | S90 | S70 | S50 | S30 |
|---|---:|---:|---:|---:|
| **R90** | R90_S90 | R90_S70 | R90_S50 | R90_S30 |
| **R70** | R70_S90 | R70_S70 | R70_S50 | R70_S30 |
| **R50** | R50_S90 | R50_S70 | R50_S50 | R50_S30 |
| **R30** | R30_S90 | R30_S70 | R30_S50 | R30_S30 |

Every trainable task therefore has **16 R/S cells**.

### 5.3 Interpretation

#### Redundancy effect

At fixed S:

`R90 → R70 → R50 → R30`

This tests the effect of progressively removing sequence redundancy while holding the outer split stringency fixed.

#### Generalization effect

At fixed R:

`S90 → S70 → S50 → S30`

This tests increasingly remote train/test generalization while holding the Stage-07 entity universe fixed.

#### Interaction

The full response surface is primary.

A compact extreme-condition contrast may be reported as:

`I = [M(R30,S30) - M(R90,S30)] - [M(R30,S90) - M(R90,S90)]`

where `M` is a selected performance metric.

R and S are interpreted as **ordered categorical factors**. The identity thresholds must not be treated as if equal numeric differences implied equal biological distances.

---

# Part II — Predictive-universe and partition contract

## 6. Stage 10E — Common predictive-eligibility freeze

### 6.1 Rationale

Representation filtering must not be performed dynamically inside each training job.

If different representations silently remove different proteins, performance differences become confounded with changes in test-set composition.

The common predictive universe must therefore be materialized once before training.

### 6.2 Biological eligibility definition

The canonical predictive-eligibility rule is:

- exact sequence belongs to the frozen Stage-10A universe;
- canonical 20-amino-acid alphabet;
- sequence length `<= 1024 aa`.

This rule is **representation-independent**.

It defines which biological entities are eligible for the main representation benchmark.

### 6.3 Representation availability is not eligibility

A representation being missing for an otherwise eligible sequence does **not** make that sequence ineligible.

Instead:

- expected sequence exists + representation exists → `AVAILABLE`;
- expected sequence exists + representation missing → representation completeness failure;
- sequence is noncanonical or >1024 aa → biologically `INELIGIBLE_FOR_MAIN_BENCHMARK`.

This distinction is critical.

The benchmark must not use:

`intersection of whatever each representation happened to generate`

because that would allow an incomplete representation to redefine the test set.

### 6.4 Stage-10E output

For each of the 144 `task × R × S` cells materialize a compact eligibility table with at least:

- `entity_id`;
- `sequence_uid`;
- `cluster_id`;
- `task`;
- `redundancy_threshold`;
- `split_threshold`;
- `outer_split`;
- `label`;
- `class_name` or equivalent frozen role;
- `sequence_length`;
- `is_canonical20`;
- `is_main_predictive_eligible`;
- upstream source/split hashes.

Recommended canonical path:

`data/predictive_universe/<task>/<R>_<S>/eligible_entities.tsv.gz`

### 6.5 Stage-10E audit

Each cell must record:

- `n_before_eligibility`;
- `n_after_eligibility`;
- `n_dropped_noncanonical`;
- `n_dropped_length_gt_1024`;
- train count before/after;
- test count before/after;
- class/role counts before/after;
- train/test prevalence before/after;
- cluster counts before/after;
- train/test exact-entity overlap;
- shared train/test S-level cluster count.

Required closure:

- no entity can move from train to test or test to train;
- filtering can only remove entities;
- exact overlap remains 0;
- shared S-level clusters remain 0.

### 6.6 Representation-completeness matrix

Before a representation can enter training, generate:

`representation_completeness.tsv`

with at least:

- `representation_id`;
- `expected_eligible_uids`;
- `available_uids`;
- `missing_uids`;
- `unexpected_uids`;
- `duplicate_uids`;
- `dimension`;
- `storage_type`;
- `status`;
- file path;
- file SHA-256.

Main-benchmark status is `READY` only if:

- all expected eligible sequence_uids required by the intended tasks are available;
- there are no duplicate sequence_uids;
- representation dimensionality is internally consistent;
- no unexpected row-to-sequence mismatch exists.

Representation incompleteness must never be repaired by silently shrinking one model's test set.

---

## 7. Stage 10F — Frozen homology-group-aware inner validation split

### 7.1 Superseded v1 rule

The v1 design proposed a row-wise stratified 80/20 split of the outer-train partition.

That rule is superseded.

A row-wise split could place proteins from the same S-level homology component into inner train and inner validation, creating an avoidable validation leak even though the outer test remained protected.

### 7.2 New inner-split principle

For every `task × R × S` cell:

1. start from the **eligible outer-train entities only**;
2. retain the S-level `cluster_id` from Stage 08;
3. assign complete clusters to inner train or inner validation;
4. target approximately:
   - inner train: 80%;
   - inner validation: 20%;
5. use deterministic `seed = 42`;
6. optimize the cluster assignment to preserve class/role composition as closely as feasible.

The inner split therefore becomes:

`group-aware + label-aware`, not row-wise stratified.

### 7.3 Constraints

Hard constraints:

- a cluster must never appear in both inner train and inner validation;
- all eligible outer-train entities must be assigned exactly once;
- no outer-test entity may enter the inner split;
- binary tasks must contain both observed classes in inner train and validation whenever mathematically feasible;
- multiclass must retain every class in both partitions whenever feasible;
- PU inner partitions must retain both observed roles P and U whenever feasible.

Soft objectives:

- validation fraction close to 0.20;
- class/role prevalence close to the eligible outer-train prevalence;
- avoid pathological concentration of minority-positive clusters in one partition.

### 7.4 Reuse

For a fixed `task × R × S`, the same inner assignment is reused across:

- all representations;
- all algorithms;
- all predefined parameter configurations.

The inner split is therefore part of the **data contract**, not the model contract.

### 7.5 Output

Recommended path:

`data/training_splits/<task>/<R>_<S>/inner_seed42.tsv.gz`

Required columns:

- `entity_id`;
- `sequence_uid`;
- `cluster_id`;
- `outer_split`;
- `inner_split`;
- `label`;
- `class_name`/role;
- `task`;
- `R`;
- `S`.

A corresponding `inner_split_manifest.json` must store:

- input hashes;
- target and observed validation fraction;
- class/role counts;
- cluster counts;
- shared inner-train/validation clusters;
- seed;
- assignment algorithm/version;
- status.

Required closure:

`shared_inner_train_validation_clusters = 0`

---

# Part III — Representation and model registries

## 8. Representation registry

Representation loading must be registry-driven, never hard-coded into the training engine.

Recommended artifact:

`config/representation_registry.tsv`

Minimum fields:

- `representation_id`;
- `representation_family`;
- `source_tool`;
- `source_model`;
- `format`;
- `dimension`;
- `dtype`;
- `sparse`;
- `requires_scaling_linear`;
- `requires_scaling_knn`;
- `requires_scaling_mlp`;
- `file_path`;
- `index_path`;
- `sha256`;
- `expected_uid_count`;
- `available_uid_count`;
- `status`;
- `notes`.

### 8.1 Frozen main benchmark representations

Nine representations are included.

#### Roxy

- ID: `roxy_48`;
- family: physicochemical/descriptive;
- 48 descriptors;
- dense matrix.

#### One-hot

- ID: `onehot_1024`;
- family: classical sequence encoding;
- canonical 20-aa alphabet;
- max length: 1024 aa;
- sparse CSR;
- no PCA/SVD in the main benchmark.

#### PLMs

1. `ankh2_ext1`
2. `esmc_300m`
3. `esm2_t6_8m`
4. `esm2_t12_35m`
5. `esm2_t30_150m`
6. `protbert`
7. `prott5_xl_uniref50`

Frozen PLM extraction contract:

- canonical sequences only;
- max length 1024 aa;
- final hidden layer;
- mean sequence pooling;
- FP32;
- one row per `sequence_uid`;
- deterministic extraction settings;
- model/source metadata and hashes recorded.

### 8.2 Representation state

The registry controls execution:

- `READY`: may be trained;
- `PENDING`: expected representation not yet materialized;
- `INCOMPLETE`: materialized but missing expected eligible entities;
- `INVALID`: schema/dimension/index audit failed;
- `BLOCKED`: technically unavailable for a requested learner.

The training launcher executes only `READY` representations.

Adding a newly completed PLM must require only changing its registry status from `PENDING` to `READY`; no trainer code should change.

---

## 9. Learning algorithms

Five supervised algorithms remain frozen:

1. **Logistic Regression (LR)**
2. **Random Forest (RF)**
3. **K-Nearest Neighbors (KNN)**
4. **XGBoost**
5. **Multilayer Perceptron (MLP)**

The objective is a controlled, low-resolution parameter exploration, not exhaustive optimization.

Exactly three predefined supervised configurations are evaluated per algorithm whenever technically compatible.

---

## 10. Feature preprocessing

All preprocessing is fitted only on the current fitting partition.

During configuration selection:

`inner train → fit preprocessing → transform inner train/validation`

After selection:

`full eligible outer train → refit preprocessing → transform outer train/test`

No scaler, feature statistic, class weight, or data-derived preprocessing parameter may use outer-test observations.

### 10.1 Roxy

Use `StandardScaler` for:

- LR;
- KNN;
- MLP.

Use raw descriptors for:

- RF;
- XGBoost.

No main-benchmark feature selection.

### 10.2 PLMs

Use `StandardScaler` for:

- LR;
- KNN;
- MLP.

Use raw FP32 embeddings for:

- RF;
- XGBoost.

No PCA or supervised dimensionality reduction in the main benchmark.

A future L2-normalized sensitivity analysis is allowed only as a separately named representation variant. It must not replace the frozen raw/standardized contract silently.

### 10.3 One-hot

One-hot remains sparse.

Rules:

- no mean centering;
- no dense conversion unless explicitly demonstrated safe for the full target matrix;
- LR uses sparse input directly;
- RF uses sparse-compatible input;
- KNN uses brute-force sparse distance;
- XGBoost uses sparse input;
- MLP runs only if the installed implementation can safely consume the sparse matrix without uncontrolled densification.

If one-hot + MLP is unsafe:

`status = BLOCKED`

No PCA/SVD fallback is allowed because dimensionality reduction would define a different representation.

---

## 11. Fixed supervised parameter panels

### 11.1 Logistic Regression

| Config | Parameters |
|---|---|
| LR01 | `C=0.1`, L2 |
| LR02 | `C=1.0`, L2 |
| LR03 | `C=10.0`, L2 |

Common settings:

- `max_iter = 2000`;
- explicit deterministic solver choice frozen in implementation;
- balanced class weighting for binary supervised tasks;
- library-supported multinomial formulation for multiclass.

Convergence state and number of iterations must be recorded.

### 11.2 Random Forest

| Config | Parameters |
|---|---|
| RF01 | 300 trees, unlimited depth, leaf=1, `max_features=sqrt` |
| RF02 | 500 trees, unlimited depth, leaf=2, `max_features=sqrt` |
| RF03 | 500 trees, depth=30, leaf=1, `max_features=sqrt` |

Use balanced/balanced-subsample weighting where applicable.

Freeze:

- `random_state = 42`;
- `n_jobs` controlled by runtime configuration, not by scientific config.

### 11.3 KNN

| Config | Parameters |
|---|---|
| KNN01 | `k=3`, distance weighting, cosine distance |
| KNN02 | `k=5`, distance weighting, cosine distance |
| KNN03 | `k=11`, distance weighting, cosine distance |

Use brute-force search for sparse/high-dimensional representations.

Prediction time and peak memory are scientifically relevant outputs.

### 11.4 XGBoost

| Config | Parameters |
|---|---|
| XGB01 | 300 trees, depth=4, lr=0.10, subsample=0.8, colsample=0.8 |
| XGB02 | 500 trees, depth=6, lr=0.05, subsample=0.8, colsample=0.8 |
| XGB03 | 800 trees, depth=4, lr=0.03, subsample=0.8, colsample=0.8 |

Use histogram-based tree construction.

Class weighting is derived from the current fitting partition only.

Freeze explicit random seed and objective settings by task type.

### 11.5 MLP

| Config | Hidden architecture | alpha | learning rate |
|---|---|---:|---:|
| MLP01 | `(256,)` | 1e-4 | 1e-3 |
| MLP02 | `(512, 128)` | 1e-4 | 5e-4 |
| MLP03 | `(256, 128, 64)` | 1e-3 | 5e-4 |

Common settings to freeze in code before production execution:

- `random_state = 42`;
- fixed batch-size rule;
- fixed maximum iterations/epochs;
- fixed activation and optimizer;
- fixed stopping tolerance;
- deterministic early-stopping rule if early stopping is retained.

Early stopping, if used, may inspect only the frozen inner validation subset.

Outer-test-informed stopping is prohibited.

---

# Part IV — Experiment registry and supervised training

## 12. Experiment registry

Before model execution, materialize the complete experiment space.

Recommended file:

`config/supervised_experiment_registry.tsv`

One row per **experiment family**:

`Task × R × S × Representation × Algorithm`

Minimum fields:

- `experiment_id`;
- `task`;
- `task_type`;
- `R`;
- `S`;
- `representation_id`;
- `algorithm`;
- `outer_split_path`;
- `inner_split_path`;
- `representation_path`;
- `representation_status`;
- `expected_train_n`;
- `expected_validation_n`;
- `expected_test_n`;
- `expected_positive_train_n`;
- `expected_positive_validation_n`;
- `expected_positive_test_n`;
- `parameter_panel_id`;
- `seed`;
- `status`;
- `priority`;
- `notes`.

### 12.1 Human-readable experiment IDs

Use deterministic IDs such as:

`SUP__global_easy__R70__S30__esm2_t6_8m__RF`

Multiclass:

`SUP__multiclass_clean__R50__S30__roxy_48__XGB`

Avoid opaque UUID-only naming.

A short deterministic hash of the experiment contract may be stored in addition to the human-readable ID.

### 12.2 Registry status

Allowed values:

- `PENDING`;
- `READY`;
- `RUNNING`;
- `DONE`;
- `FAILED`;
- `BLOCKED`;
- `SKIP`.

No missing row may be interpreted as an implicit failed experiment.

---

## 13. Supervised training protocol

For each experiment family:

1. verify representation status is `READY`;
2. verify outer-split hash;
3. verify inner-split hash;
4. verify expected eligible sequence_uids are present in the representation;
5. load eligible outer train/test membership;
6. load frozen inner train/validation assignment;
7. load only required representation rows;
8. for each of the three predefined configurations:
   1. fit preprocessing on inner train only;
   2. fit model on inner train;
   3. predict probabilities/scores on inner train and validation;
   4. calculate validation metrics;
   5. record runtime, warnings, convergence, and resource use;
   6. persist validation predictions;
9. select exactly one configuration using the frozen task-specific selection rule;
10. refit preprocessing on the complete eligible outer train;
11. refit the selected model configuration on complete eligible outer train;
12. predict outer test exactly once;
13. calculate final test metrics;
14. persist per-sequence probabilities/scores;
15. persist class-level metrics/confusion information;
16. write a complete immutable manifest;
17. set registry status to `DONE`.

The outer-test result must never alter configuration choice, preprocessing, class weighting, stopping, threshold, or representation handling.

---

## 14. Configuration selection

### 14.1 Binary supervised

Primary:

1. MCC.

Secondary/tie support:

2. AP.

Recommended deterministic tie rule:

1. larger MCC;
2. if absolute MCC difference is below a frozen numerical tolerance, larger AP;
3. if still tied, choose the simpler configuration according to the predefined config order.

The tolerance and deterministic order must be stored in the implementation config.

### 14.2 Multiclass

Primary:

1. multiclass MCC.

Secondary:

2. macro-F1.

Use the same deterministic tie philosophy.

### 14.3 Outer test never breaks ties

Outer-test performance cannot be consulted under any circumstances during configuration selection.

---

# Part V — Mandatory result bundle

## 15. Output philosophy

The result bundle must preserve enough information to reconstruct every principal analysis **without retraining models**.

Scalar final metrics alone are insufficient.

The most important retained object is the per-sequence prediction table.

No serialized production model is required during the benchmark phase unless needed for a later discovery stage.

---

## 16. Required supervised output directory

Recommended path:

`results/training_supervised/<task>/<R>_<S>/<representation>/<algorithm>/`

Required files:

- `manifest.json`
- `validation_metrics.tsv`
- `validation_predictions.parquet`
- `selected_configuration.json`
- `outer_test_metrics.tsv`
- `outer_test_predictions.parquet`
- `per_class_metrics.tsv`
- `confusion_matrix.tsv`
- `runtime.json`
- `warnings.log`

Optional but useful:

- `training_metrics.tsv`;
- `feature_preprocessing_summary.json`;
- `calibration_bins.tsv`;
- `model_diagnostics.json`.

### 16.1 `validation_metrics.tsv`

One row per predefined configuration.

Minimum fields:

- `experiment_id`;
- `config_id`;
- all hyperparameters;
- `n_inner_train`;
- `n_validation`;
- class counts;
- prevalence;
- primary metric;
- secondary metric;
- all supporting metrics;
- fit time;
- validation prediction time;
- peak memory if available;
- convergence status;
- warning count;
- selected flag.

### 16.2 `validation_predictions.parquet`

Retain configuration-level validation predictions.

Binary minimum fields:

- `experiment_id`;
- `config_id`;
- `entity_id`;
- `sequence_uid`;
- `y_true`;
- `y_pred`;
- `score_positive`.

Multiclass minimum fields:

- `experiment_id`;
- `config_id`;
- `entity_id`;
- `sequence_uid`;
- `y_true`;
- `y_pred`;
- one probability column per class.

Use compressed Parquet, preferably Zstandard compression.

Retaining these predictions enables:

- post hoc selection audits;
- calibration analysis;
- threshold-sensitivity analysis;
- model-agreement analysis;
- examination of whether config choice is driven by a few difficult sequences.

### 16.3 `outer_test_predictions.parquet`

This file is mandatory.

Binary:

- `experiment_id`;
- `entity_id`;
- `sequence_uid`;
- `cluster_id`;
- `y_true`;
- `y_pred`;
- `score_positive`.

Multiclass:

- `experiment_id`;
- `entity_id`;
- `sequence_uid`;
- `cluster_id`;
- `y_true`;
- `y_pred`;
- `p_endolysin`;
- `p_depolymerase`;
- `p_val`.

The full biological metadata need not be duplicated in every prediction row because it can be joined from the frozen master universe, but stable IDs and cluster membership must be preserved.

### 16.4 `manifest.json`

Minimum content:

- experiment ID;
- complete factor values;
- input file paths and hashes;
- representation path/hash;
- outer split hash;
- inner split hash;
- software versions;
- algorithm/configuration;
- random state;
- preprocessing;
- number of entities by partition/class;
- dimensionality;
- sparse/dense status;
- start/end timestamps;
- hardware summary;
- final status;
- error/warning summary.

---

# Part VI — Supervised metrics

## 17. Binary metrics

### 17.1 Primary metric

**Matthews correlation coefficient (MCC)**

MCC is the primary configuration-selection and reporting metric.

### 17.2 Principal ranking metric

**Average Precision (AP)**

AP must always be stored with:

`pi = N_positive / N`

and:

`AP_gain = AP - pi`

because AP baseline depends on prevalence.

The Stage-08 factorial reports already show that train/test class composition can change modestly as S changes, so prevalence must be treated as an explicit analytical covariate rather than assumed constant.

### 17.3 Required supporting metrics

Store:

- ROC-AUC;
- F1;
- precision;
- recall/sensitivity;
- specificity;
- balanced accuracy;
- TP;
- TN;
- FP;
- FN;
- positive prevalence;
- predicted-positive fraction;
- Brier score when probability output is available;
- log loss when probability output is available.

Optional descriptive calibration metrics may include:

- expected calibration error;
- calibration slope/intercept;
- reliability-bin summaries.

These are diagnostics, not primary model-selection criteria.

### 17.4 Decision threshold

Canonical threshold-dependent metrics use:

`threshold = 0.5`

No outer-test threshold optimization.

Because raw scores are retained, threshold-sensitivity curves may be calculated later as secondary descriptive analyses without changing the canonical benchmark.

---

## 18. Multiclass metrics

### 18.1 Primary

**Multiclass MCC**

### 18.2 Secondary

**Macro-F1**

### 18.3 Supporting

Store:

- balanced accuracy;
- macro precision;
- macro recall;
- weighted F1;
- multiclass log loss;
- per-class precision;
- per-class recall;
- per-class F1;
- per-class support;
- confusion matrix;
- class probability distributions.

VAL support in train/validation/test must always be displayed as absolute counts.

---

# Part VII — Uncertainty and statistical comparison contract

## 19. Bootstrap uncertainty

No repeated-training stability-seed experiment is used.

Uncertainty in final outer-test metrics is estimated from frozen predictions.

Default:

- 1,000 bootstrap replicates;
- percentile 95% confidence intervals;
- class-stratified resampling for supervised binary/multiclass metrics where mathematically feasible;
- no model retraining during bootstrap.

Bootstrap objects are generated in an analysis stage after training.

For Low-N VAL conditions:

- always show absolute positive support;
- bootstrap distributions may be unstable or discrete;
- intervals are interpreted as finite-test-set uncertainty, not training-seed uncertainty.

---

## 20. Paired versus non-paired comparisons

This distinction is mandatory.

### 20.1 Valid paired comparisons

Within the same:

`task × R × S`

all representations and algorithms evaluate the same eligible outer-test entities.

Therefore the following can use paired bootstrap differences:

- PLM vs Roxy;
- PLM vs one-hot;
- PLM vs PLM;
- algorithm vs algorithm;
- ESM2 size ladder;
- selected model vs another selected model in the same cell.

Store paired deltas such as:

- `ΔMCC`;
- `ΔAP`;
- `Δbalanced_accuracy`;
- 95% paired-bootstrap CI.

### 20.2 Comparisons that are not automatically paired

The following generally evaluate different entity sets:

- R90 vs R70/R50/R30;
- S90 vs S70/S50/S30;
- Global easy vs Global hard vs Global mixed;
- subtype task vs global task.

Therefore native-test comparisons across these factors must **not** be called paired sequence-level comparisons.

They are matched-factor contrasts or factorial response-surface comparisons.

### 20.3 Negative-definition comparison

Global easy, hard, and mixed share a conceptual positive definition but use different negative universes and independently generated splits.

Primary analysis:

- compare performance surfaces while holding R/S/representation/algorithm constant;
- report metric differences as matched experimental contrasts;
- compare prevalence and FP/FN structure;
- do not use entity-level paired bootstrap unless an explicitly defined common entity subset exists.

---

# Part VIII — Factorial performance analysis

## 21. Native-test R × S analysis

For each:

`task × representation × algorithm`

construct the complete 4 × 4 surface:

`M(R,S)`

Primary binary surfaces:

- MCC;
- AP;
- AP gain.

Multiclass:

- MCC;
- macro-F1.

Supporting surfaces may include:

- recall;
- specificity;
- balanced accuracy;
- runtime.

### 21.1 Redundancy contrast

At fixed S:

`ΔR(S) = M(R30,S) - M(R90,S)`

### 21.2 Generalization contrast

At fixed R:

`ΔS(R) = M(R,S30) - M(R,S90)`

### 21.3 Extreme interaction contrast

`I = [M(R30,S30)-M(R90,S30)] - [M(R30,S90)-M(R90,S90)]`

The full surface remains primary; these contrasts are summaries.

### 21.4 Marginal summaries

For descriptive synthesis compute:

- marginal median/mean performance by R;
- marginal median/mean performance by S;
- rank distributions;
- worst-case performance;
- best-case performance;
- performance range across S;
- performance range across R.

Do not imply linearity in R or S.

---

## 22. Common-test sensitivity analysis

### 22.1 Motivation

Native R/S test sets differ in composition.

A performance drop from S90 to S30 can reflect both:

1. a genuinely harder train/test generalization condition;
2. a different set of proteins assigned to test.

Because per-sequence predictions are retained, this can be decomposed.

### 22.2 Fixed-R common-S analysis

For each:

`task × R`

identify:

`T_common = T_S90 ∩ T_S70 ∩ T_S50 ∩ T_S30`

Recalculate every model's metrics on this common entity subset.

Report:

- common-test size;
- class supports;
- fraction of each native test retained;
- common-test metrics.

This produces a sensitivity analysis in which the evaluated proteins are held constant across S.

### 22.3 Pairwise common-test analysis

If the four-way intersection is too small or lacks adequate class support, compute pairwise intersections:

- S90 ∩ S70;
- S90 ∩ S50;
- S90 ∩ S30;
- S70 ∩ S50;
- S70 ∩ S30;
- S50 ∩ S30.

All intersection sizes and class supports must be reported.

No common-test result is considered primary if its class support is clearly insufficient for meaningful interpretation.

### 22.4 R comparisons

R changes the biological entity universe because Stage-07 redundancy filtering removes sequences.

A common-entity sensitivity analysis across R may be performed on intersections, but it is secondary because increasingly stringent R conditions intentionally define smaller biological datasets.

The native R effect therefore remains the primary interpretation.

---

## 23. Representation analyses

Questions:

1. Do PLMs outperform Roxy?
2. Do PLMs outperform one-hot?
3. Does the advantage depend on R or S?
4. Which representations degrade least at S30?
5. Which representations remain useful at R30?
6. Is representation rank stable across algorithms?
7. What predictive gain is obtained per unit of computational cost?

Primary matched analysis unit:

`task × R × S × algorithm`

Within that cell, representation comparisons are paired.

Recommended summaries:

- paired ΔMCC/ΔAP;
- win/tie/loss counts;
- rank distributions;
- generalization robustness;
- redundancy robustness;
- worst-case MCC;
- compute-adjusted performance.

---

## 24. ESM2 size ladder

Compare:

- ESM2 8M;
- ESM2 35M;
- ESM2 150M.

The comparison is especially valuable because extraction/training protocol is otherwise matched.

Questions:

- does capacity improve median performance?
- does capacity improve S30 performance?
- does capacity improve Low-N subtype performance?
- are gains consistent enough to justify compute cost?

Use paired comparisons within fixed cells.

---

## 25. Algorithm analyses

Within each fixed representation and cell compare:

- LR;
- RF;
- KNN;
- XGBoost;
- MLP.

Store and analyze:

- predictive metrics;
- fit time;
- prediction time;
- peak memory;
- convergence/failure;
- sensitivity to R;
- sensitivity to S;
- rank stability;
- frequency of algorithm wins;
- interaction between representation family and classifier family.

A high-performing algorithm that repeatedly fails or becomes computationally prohibitive under large PU/one-hot conditions must not be presented as universally preferable.

---

## 26. Robustness summaries

### 26.1 Generalization robustness

At fixed R:

`G = M(R,S30) - M(R,S90)`

Values closer to zero indicate smaller degradation under stricter remote-generalization control.

### 26.2 Redundancy robustness

At fixed S:

`D = M(R30,S) - M(R90,S)`

Values closer to zero indicate less apparent dependence on redundant/closely related training examples.

### 26.3 Worst-case robustness

For each model family:

`W = min M(R,S)`

over the relevant 16 cells.

Worst-case performance is a useful complement to mean/rank summaries.

### 26.4 Error consistency

Using retained prediction tables, calculate per entity:

- number of models evaluating the entity;
- number/fraction misclassifying it;
- mean predicted positive probability;
- score variance across algorithms;
- score variance across representations;
- repeated error status.

These analyses identify sequences that are consistently difficult across modeling choices.

They are descriptive in the classification paper and can later support a dedicated determinants/patterns study.

---

# Part IX — Low-N analysis

## 27. VAL handling

VAL is retained as a biologically meaningful stress test.

Rules:

1. never hide low positive counts;
2. report TP/FN absolute counts alongside recall;
3. report train/validation/test positive support;
4. report bootstrap distributions cautiously;
5. avoid declaring small metric differences meaningful when they correspond to one or two entities;
6. emphasize rank/score behavior and consistency across R/S rather than only a single scalar value;
7. include multiclass VAL per-class metrics even when macro summaries appear stable.

Recommended figures:

- VAL positive support versus R/S;
- MCC/recall heatmaps annotated with positive counts;
- score distributions for VAL positives;
- model-agreement summaries.

---

# Part X — Positive–unlabeled learning

## 28. PU branch separation

PU learning remains a separate experiment registry, output tree, and analysis track.

Main PU experiment identifier:

`PU__<task>__Rxx__Syy__<representation>__<method>__<algorithm>`

No supervised metric table should mix ordinary binary and PU results without explicit labeling.

---

## 29. Naive P-vs-U diagnostic baseline

Run:

- LR;
- RF;
- KNN;
- XGBoost;
- MLP.

Use the same three predefined configurations as the supervised benchmark.

U is treated as temporary negative **only for optimization/diagnostic purposes**.

Metrics that assume U=negative must be prefixed or labeled as:

`observed_label_*`

Examples:

- `observed_label_mcc`;
- `observed_label_ap`;
- `observed_label_roc_auc`.

These are not biological specificity/precision estimates.

---

## 30. BaggingPU

Base learners:

- LR;
- RF;
- XGBoost.

Frozen method-level design:

- 30 internal bags;
- fixed random seed;
- sampled U:P ratio = **1:1**;
- each bag contains all fitting positives plus a resampled U subset of equal size;
- aggregate score = mean predicted positive probability across bags.

To avoid an uncontrolled PU hyperparameter explosion, use a single frozen middle-complexity base configuration per algorithm:

- LR02;
- RF02;
- XGB02.

This makes the PU-method comparison interpretable and avoids tuning against a validation target that is itself partially unlabeled.

Per test entity store:

- mean BaggingPU score;
- standard deviation across bags;
- optional 5th/50th/95th score quantiles;
- number of bags contributing.

Internal bag-level models are not stability repetitions.

---

## 31. Elkan–Noto

Base learners:

- LR;
- RF;
- XGBoost;
- MLP.

Use one frozen middle configuration:

- LR02;
- RF02;
- XGB02;
- MLP02.

### 31.1 Labeling-propensity estimation

To avoid estimating `c` on positives seen during fitting, use **group-aware cross-fitted estimation within outer train**.

Recommended contract:

- 5 cross-fitting folds;
- folds respect S-level homology clusters;
- each labeled positive receives an out-of-fold positive-label score;
- `c` is estimated from held-out positive predictions only;
- outer test is never used.

After estimating `c`:

1. fit the base classifier on the complete eligible outer train;
2. apply the frozen Elkan–Noto correction using the cross-fitted `c`;
3. score outer test.

Store:

- `c`;
- fold-level `c` estimates;
- number of positive entities per fold;
- corrected and uncorrected scores;
- clipping frequency if corrected probabilities exceed [0,1].

The exact cross-fitting assignment must be reused across representations for a given PU `task × R × S`.

---

## 32. nnPU

Base learner:

- MLP.

nnPU requires a positive-class prior:

`pi_P = P(Y=1)`

The prior estimator remains the only unresolved major PU methodological decision.

Rules already frozen:

- prior estimation uses outer-train information only;
- no outer-test prevalence;
- no outer-test score distribution;
- prior estimate and uncertainty must be persisted;
- if the prior cannot be estimated robustly, nnPU remains `PENDING_PRIOR` rather than being forced into the benchmark.

Until a dedicated prior-estimation protocol is frozen:

`nnPU status = PENDING_PRIOR`

---

## 33. PU evaluation

Primary PU evaluation focuses on ranking and recovery of labeled positives relative to U.

Store:

- held-out/cross-fitted positive recall;
- median P score;
- median U score;
- P and U score quantiles;
- fraction of positives above the U 90th/95th/99th percentiles;
- fraction of U above a frozen score threshold;
- score distributions;
- PU empirical risk when defined;
- BaggingPU score variability;
- Elkan–Noto `c`;
- nnPU prior if eventually enabled.

High-scoring U proteins are **candidates**, not asserted positives.

---

## 34. PU result bundle

Recommended path:

`results/training_pu/<task>/<R>_<S>/<representation>/<method>/<algorithm>/`

Required:

- `manifest.json`;
- `method_diagnostics.tsv`;
- `outer_test_scores.parquet`;
- `score_distribution_summary.tsv`;
- `runtime.json`;
- `warnings.log`.

Method-specific outputs:

### Naive

- configuration-level validation metrics/predictions;
- selected config;
- observed-label test diagnostics.

### BaggingPU

- bagging parameters;
- per-entity score mean/std/quantiles;
- bag fit failure count.

### Elkan–Noto

- cross-fit fold assignment;
- fold-level propensity estimates;
- final `c`;
- raw and corrected per-entity scores.

### nnPU

- prior estimate;
- prior-estimation manifest;
- PU empirical risk trajectory.

---

# Part XI — Challenge and discovery stages

## 35. Challenge sets

Challenge sets are never used for:

- parameter selection;
- threshold selection;
- representation selection;
- algorithm selection;
- PU class-prior estimation.

Frozen challenge groups include:

- Core multiple-nonexclusive;
- Core other lytic;
- UniProt lytic-related;
- UniProt ambiguous structural.

### 35.1 Global models

Report:

- predicted-positive fraction;
- score distribution;
- median/quantile score;
- agreement across representations/algorithms.

### 35.2 Subtype/multiclass

For ambiguous/multiple lytic sequences report:

- class probabilities;
- top class;
- prediction entropy;
- agreement/disagreement;
- rank consistency.

Do not force ambiguous sequences into clean single-label accuracy calculations.

### 35.3 Discovery

Discovery/mining is downstream of the completed benchmark.

Candidate selection must consume frozen final-model outputs and a separately defined selection/scoring contract.

The main classification paper may include a compact demonstration, but discovery must not alter benchmark definitions retrospectively.

---

# Part XII — Computational and reproducibility contract

## 36. Runtime/resource logging

Each model fit must record, where feasible:

- start/end timestamp;
- wall-clock fit time;
- validation/test prediction time;
- peak resident memory;
- CPU count;
- GPU use if any;
- representation dimension;
- input sparse/dense status;
- train/validation/test row counts;
- number of features;
- algorithm thread count;
- convergence state;
- warning count.

This enables explicit predictive-performance versus computational-cost analysis.

---

## 37. Failure handling

Every experiment family has one explicit state:

- `PENDING`
- `READY`
- `RUNNING`
- `DONE`
- `FAILED`
- `BLOCKED`
- `SKIP`
- `PENDING_PRIOR`

A failed job must preserve:

- exception/error text;
- factor identifiers;
- algorithm/config;
- representation;
- input hashes;
- runtime until failure;
- software versions;
- resource information when available.

No failed experiment may disappear from the registry.

---

## 38. Reproducibility freeze

Record:

- Python version;
- scikit-learn version;
- XGBoost version;
- pandas/numpy/scipy versions;
- Sylphy version;
- Roxy version;
- MMseqs2/BioSieve versions for upstream split provenance;
- random seed = 42;
- parameter-panel version;
- representation registry version/hash;
- experiment registry version/hash;
- outer-split hashes;
- eligibility-table hashes;
- inner-split hashes;
- representation file hashes;
- hardware/runtime metadata.

All stochastic algorithms use explicit random states wherever supported.

No repeated-seed stability benchmark is part of the main design.

---

# Part XIII — Global aggregation and analysis products

## 39. Supervised master table

After training, materialize one row per final supervised experiment family.

Recommended file:

`results/aggregated/supervised_master_results.tsv`

Minimum fields:

- `experiment_id`;
- `task`;
- `task_type`;
- `R`;
- `S`;
- `representation`;
- `algorithm`;
- `selected_config`;
- `n_outer_train`;
- `n_inner_train`;
- `n_validation`;
- `n_outer_test`;
- class counts;
- prevalence;
- MCC;
- AP;
- AP gain;
- ROC-AUC;
- F1;
- precision;
- recall;
- specificity;
- balanced accuracy;
- Brier;
- log loss;
- runtime;
- peak memory;
- status;
- warning count.

Multiclass uses corresponding multiclass fields.

---

## 40. Configuration-selection table

Materialize all configuration-level results:

`results/aggregated/supervised_configuration_results.tsv`

This table supports:

- selected-config frequency;
- hyperparameter sensitivity;
- selection margin;
- failure/convergence frequency;
- whether model family conclusions depend on one aggressive configuration.

---

## 41. Prediction index

Create a lightweight index:

`results/aggregated/prediction_file_index.tsv`

with:

- experiment ID;
- prediction file path;
- hash;
- row count;
- probability columns;
- status.

This prevents needing to scan thousands of folders during analysis.

---

## 42. Bootstrap/comparison tables

Recommended outputs:

- `metric_bootstrap_ci.tsv`;
- `paired_representation_deltas.tsv`;
- `paired_algorithm_deltas.tsv`;
- `esm2_capacity_deltas.tsv`;
- `factorial_rs_contrasts.tsv`;
- `common_test_rs_results.tsv`;
- `low_n_val_summary.tsv`;
- `compute_performance_tradeoff.tsv`;
- `error_consistency_summary.tsv`.

---

# Part XIV — Experiment scale

## 43. Supervised benchmark

Dimensions:

- 7 supervised tasks;
- 16 R/S cells;
- 9 representations;
- 5 algorithms.

Experiment families:

`7 × 16 × 9 × 5 = 5,040`

With three supervised configs:

- validation fits: **15,120**;
- final outer-train refits: **5,040**.

Total top-level supervised fits:

**20,160**

Validation predictions are retained for all three configs.

Outer-test predictions are retained for the selected config only.

---

## 44. PU benchmark

Per `PU task × R × S × representation`, the intended method rows are:

### Naive P-vs-U

5 algorithms × 3 configs = 15 validation evaluations.

### BaggingPU

3 fixed base learners = 3 method evaluations.

### Elkan–Noto

4 fixed base learners = 4 method evaluations.

### nnPU

1 MLP method = 1 evaluation, once prior estimation is frozen.

Total method/config evaluations per PU cell/representation:

`15 + 3 + 4 + 1 = 23`

Final PU experiment-family rows:

`5 + 3 + 4 + 1 = 13`

Overall PU family count:

`2 × 16 × 9 × 13 = 3,744`

nnPU rows remain registered as `PENDING_PRIOR` until enabled.

BaggingPU additionally performs 30 internal bags per method family.

Elkan–Noto additionally performs the frozen group-aware cross-fitting required for propensity estimation.

These internal fits must be tracked separately from top-level experiment counts.

---

# Part XV — Implementation stages and gates

## 45. Stage 10E — Predictive-universe freeze

**Input:** frozen Stage-08 cells + Stage-10A sequence universe  
**Output:** one eligible-entity table per cell + audit  
**Gate:** all 144 cells PASS

Required before any production training.

---

## 46. Stage 10F — Inner split freeze

**Input:** Stage-10E eligible outer-train entities + S-level clusters  
**Output:** group-aware 80/20 inner assignment per cell  
**Gate:** all 144 cells PASS

Key audit:

`shared_inner_train_validation_clusters = 0`

---

## 47. Stage 11A — Representation registry and completeness audit

**Input:** Roxy/one-hot/PLM artifacts  
**Output:** frozen registry + completeness matrix  
**Gate:** representation status explicitly READY/PENDING/INCOMPLETE/BLOCKED

No training script may infer file locations ad hoc.

---

## 48. Stage 11B — Experiment registry

Generate:

- complete supervised registry;
- complete PU registry;
- expected output locations;
- statuses.

No model training yet.

---

## 49. Stage 11C — Supervised smoke test

Before full execution, use a deliberately small but structurally representative subset:

- Global easy;
- one R/S cell;
- Roxy;
- one-hot;
- ESM2-t6;
- LR;
- RF;
- XGBoost.

Validate:

- UID joining;
- inner-split reuse;
- preprocessing isolation;
- configuration selection;
- metric calculation;
- probability export;
- manifest schema;
- runtime logging;
- restart/resume behavior.

Smoke-test outputs are discarded or clearly marked as test artifacts before production execution.

---

## 50. Stage 11D — Supervised production benchmark

Execution order is operational, not scientific.

Suggested:

1. global easy/hard/mixed;
2. subtype binary;
3. multiclass;
4. remaining representations as they become READY.

Because representation availability is registry-driven, currently available representations can be executed without waiting for all pending PLMs, provided the common biological eligibility universe has already been frozen.

---

## 51. Stage 11E — PU benchmark

Suggested order:

1. naive P-vs-U;
2. BaggingPU;
3. Elkan–Noto;
4. nnPU only after prior-estimation closure.

---

## 52. Stage 12 — Aggregation and statistical analysis

No new model training.

Substages:

- 12A: aggregate metrics/configuration results;
- 12B: bootstrap CIs;
- 12C: paired representation/algorithm comparisons;
- 12D: native R×S surfaces;
- 12E: common-test R/S sensitivity;
- 12F: Low-N analysis;
- 12G: compute-performance analysis;
- 12H: PU synthesis.

---

## 53. Stage 13 — Challenge/discovery evaluation

Only after Stage 12 benchmark freeze.

Challenge/discovery outputs are reported separately from canonical benchmark accuracy.

---

# Part XVI — Paper-level analysis map

## 54. Main questions and corresponding analyses

| Question | Primary analysis |
|---|---|
| Q1. Negative-evidence definition | easy vs hard vs mixed matched-factor surfaces |
| Q2. Redundancy | R90→R30 at fixed S |
| Q3. Remote generalization | S90→S30 at fixed R |
| Q4. R×S interaction | complete 4×4 surfaces + interaction contrast |
| Q5. Representation | paired within-cell Roxy/one-hot/PLM comparisons |
| Q6. PLM scale | paired ESM2 8M/35M/150M comparison |
| Q7. Algorithm | paired within-representation classifier comparison |
| Q8. Low-N | VAL support-aware analysis |
| Q9. PU | ranking/recovery diagnostics, not U-as-negative accuracy |
| Q10. Efficiency | performance vs runtime/memory/dimension |

---

## 55. Main-figure philosophy

The paper should avoid giant “best model” tables.

Likely high-value figure structures:

1. **Experimental-design figure:** evidence → task → R → S → representation → learner → analysis.
2. **R×S performance surfaces:** selected representations/algorithms, MCC/AP.
3. **Representation robustness:** paired deltas and rank distributions.
4. **Negative-definition effect:** easy/hard/mixed.
5. **Low-N subtype analysis:** VAL support + performance.
6. **PU ranking analysis.**
7. **Compute-performance trade-off.**
8. **Compact challenge/discovery demonstration.**

Full configuration tables, all 144-cell metrics, all runtime tables, and extensive bootstrap results belong in SI/repository.

---

# Part XVII — Frozen decisions

## 56. Frozen methodological decisions

The following are frozen unless a genuine blocking implementation issue is discovered:

- exact amino-acid sequence is the identity unit;
- stable SHA-256 `sequence_uid` is the representation join key;
- Stage-07 thresholds: R90/R70/R50/R30;
- Stage-08 thresholds: S90/S70/S50/S30;
- complete 4 × 4 R/S factorial;
- Stage-08 144-cell freeze is authoritative;
- outer test is immutable;
- biological predictive eligibility is canonical20 + length <=1024;
- representation incompleteness may not redefine eligibility;
- one fixed inner split per `task × R × S`;
- inner split respects S-level homology clusters;
- approximate inner split target 80/20;
- seed 42;
- nine main representations;
- five supervised algorithms;
- three supervised configurations per algorithm;
- no repeated-seed stability benchmark;
- preprocessing fitted only on fitting partitions;
- MCC is primary binary metric;
- AP is principal secondary binary ranking metric;
- multiclass MCC is primary multiclass metric;
- macro-F1 is secondary multiclass metric;
- decision threshold 0.5 for canonical threshold-dependent metrics;
- final per-sequence test predictions are mandatory;
- validation metrics/predictions for all configs are retained;
- paired bootstrap is restricted to identical evaluated entity sets;
- R/S native-test analyses are not mislabeled as paired;
- common-test sensitivity analysis is retained;
- runtime/resource metadata are scientific outputs;
- PU is a separate analytical track;
- U is never a confirmed negative class;
- BaggingPU uses 30 bags and 1:1 sampled U:P;
- PU-specific BaggingPU/Elkan–Noto methods use frozen middle base configurations;
- Elkan–Noto propensity estimation uses outer-train-only group-aware cross-fitting;
- challenge sets never train/select models;
- no production model serialization is required during benchmark training.

---

## 57. Remaining methodological decision

One major PU item remains open:

> **Positive-class prior estimation for nnPU**

Before nnPU is enabled, freeze:

- estimator;
- assumptions;
- parameters;
- group-aware validation/cross-fitting strategy if applicable;
- uncertainty reporting;
- acceptable prior range/failure rules;
- exact outer-train-only data use.

Until then:

`nnPU = PENDING_PRIOR`

This open decision does not block supervised training, naive PU, BaggingPU, or Elkan–Noto.

---

# Part XVIII — Final experimental architecture

```text
FROZEN BIOLOGICAL EVIDENCE
        ↓
FROZEN TASK DEFINITIONS
        ↓
INDEPENDENT REDUNDANCY CONTROL
R90 / R70 / R50 / R30
        ↓
FULL HOMOLOGY-AWARE OUTER FACTORIAL
S90 / S70 / S50 / S30
144/144 CELLS PASS
        ↓
COMMON BIOLOGICAL PREDICTIVE ELIGIBILITY
canonical20 + length ≤ 1024
        ↓
GROUP-AWARE INNER VALIDATION
complete S-level clusters, ~80/20, seed 42
        ↓
REPRESENTATION REGISTRY
Roxy / One-hot / 7 PLMs
        ↓
EXPERIMENT REGISTRY
Task × R × S × Representation × Algorithm
        ↓
CONFIGURATION EVALUATION
frozen three-config panels
        ↓
VALIDATION-ONLY SELECTION
        ↓
FULL OUTER-TRAIN REFIT
        ↓
FROZEN OUTER TEST
        ↓
PER-SEQUENCE PREDICTION FREEZE
        ↓
PAIRED + FACTORIAL + COMMON-TEST + LOW-N ANALYSES
        ↓
PU / CHALLENGE / DISCOVERY EXTENSIONS
```

The benchmark is therefore designed to separate the effects of biological evidence, sequence redundancy, homology-aware generalization, numerical representation, and classifier choice while preserving the raw prediction information required to investigate interactions, robustness, uncertainty, computational cost, and recurrent biological errors.

---

## 58. Implementation principle

The training engine must remain deliberately simple.

Scientific complexity belongs in the **frozen registries and contracts**, not in hidden trainer logic.

The implementation should follow:

```text
validate registry row
→ load frozen IDs
→ load frozen representation
→ fit frozen configurations
→ select on frozen validation
→ refit
→ predict frozen test
→ save everything required by the contract
→ stop
```

The trainer must not:

- create new biological labels;
- resplit data;
- redefine eligibility;
- repair missing embeddings by dropping entities;
- optimize thresholds on test;
- choose new hyperparameters;
- reinterpret U as negative;
- silently skip failures.

This separation between experimental design and model execution is essential for the intended methodological analysis and for reproducibility of the final PhageLysClass benchmark.

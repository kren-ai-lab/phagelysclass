# Scientific Workflow

**PhageLysClass | Source-code and data-flow guide**

The numbered scripts in `src/` preserve the original methodological progression. The workflow is designed to separate biological evidence, task construction, sequence-relatedness controls, numeric representations and predictive modeling. Each section below describes the purpose of the existing code, its canonical inputs and its original outputs.

## 1. Biological data and evidence states

The workflow begins with **PhageLysData v1.0 Core** as the curated positive resource and a taxonomy-defined UniProtKB phage-protein background. The UniProt background is retrieved broadly and classified *after* download, so functional search terms do not predetermine what counts as a negative.

| Stage | Existing script in `src/dataset_processing/` | Original inputs | Original outputs |
| --- | --- | --- | --- |
| **01** | `01_extract_phagelys_raw.py` | PhageLysData Core master (`core/phagelysdata_core.parquet`) | `data/raw/phagelysdata/` with Core, single-class, multiple-nonexclusive and other-lytic tables |
| **02** | `02_build_uniprot_background.py` | `config/queries/uniprot_background_queries_taxonomic.tsv`, UniProtKB REST API | Cursor-paginated phage-protein records and query manifests under `data/raw/uniprot_background/` |
| **03** | `03_classify_uniprot_background.py` | Raw UniProtKB background, Core positives and evidence rules | Accession-level evidence-state candidate pools under `data/processed/uniprot_background/` |
| **04** | `04_audit_uniprot_background.py` | Stage-03 classified records | Rule/annotation distributions, conflict summaries and review samples in `data/reports/` |
| **05** | `05_build_pnu_universes.py` | Stage-01 Core and Stage-03 evidence states | Exact-sequence P/N/U master under `data/processed/pnu_universes/` |

Biological evidence is retained as distinct states: `P_core` comes from PhageLysData; `N_easy` and `N_hard` are candidate negative strata with different types of alternative-function evidence; `U_explicit` and `U_other_unresolved` remain unlabeled. Target-like or structurally ambiguous annotations are segregated into challenge populations instead of being silently declared negative.

**Scientific rule:** *U means unlabeled, not confirmed non-lytic.* Exact-sequence decontamination of overlapping Core/UniProt entries occurs before downstream tasks are built.

## 2. Task definitions and redundancy controls

The Stage-06 scripts in [`../src/task_definitions/`](../src/task_definitions/) construct the task-specific datasets from the frozen P/N/U master. No feature computation, model selection, balancing or splitting is performed at this stage.

| Script | Scientific task |
| --- | --- |
| `06a_global_easy.py` | Core positives versus clearly alternative-function negative candidates |
| `06b_global_hard.py` | Core positives versus harder, lysis-adjacent negative candidates |
| `06c_global_mixed.py` | Global positives versus mixed easy/hard candidate negatives |
| `06d_subtype_binary.py` | Independent Endolysin, Depolymerase and VAL subtype recognition |
| `06e_multiclass_clean.py` | Conditional three-class functional discrimination |
| `06f_pu_global.py` | Positive–unlabeled conservative and broad settings |
| `06g_challenge_sets.py` | Ambiguous, non-exclusive and lytic-related diagnostic cohorts |

Outputs are materialized under `data/tasks/`, typically with `dataset.tsv.gz`, `summary.tsv` and `manifest.json`.

**Stage 07 — redundancy.** Scripts in [`../src/redundancy/`](../src/redundancy/) apply BioSieve/MMseqs2 separately within biological strata. Four independent identity thresholds (`H90`, `H70`, `H50`, `H30`) are each derived from the original Stage-06 dataset—not chained one after another. This preserves group composition while controlling within-stratum redundancy. Cross-stratum homology is handled at the splitting stage.

**Stage 08 — splitting.** [`08_build_factorial_homology_splits.py`](../src/splitting/08_build_factorial_homology_splits.py) builds the full redundancy × homology-split design: `R90/R70/R50/R30 × S90/S70/S50/S30`. Homology clustering is across the whole task, with the original target of approximately 70% outer train and 30% outer test, keeping cluster assignments together. The original files live under `data/splits/<task>/<R>_<S>/`.

**Stage 09 — reporting.** Scripts in [`../src/reporting/`](../src/reporting/) summarize the dataset flow and factorial split properties.

## 3. Numerical representations and fixed eligibility

The Stage-10 representation workflow uses a **shared exact-sequence universe** so that descriptors/embeddings are not redundantly computed for each task. The code lives in [`../src/representations/`](../src/representations/).

| Stage | Script | Role |
| --- | --- | --- |
| **10A** | `10a_build_unique_sequence_universe.py` | Consolidate unique sequences spanning supervised, PU and challenge universes |
| **10B** | `10b_characterize_unique_sequences_roxy.py` | Calculate the compact interpretable Roxy descriptor panel (48 features) |
| **10C** | `10c_encode_unique_sequences_onehot_sylphy.py` | Create a Sylphy one-hot representation with sparse CSR storage |
| **10D** | `10d_generate_plm_embeddings_sylphy.py` | Generate Sylphy protein-language-model embeddings with fixed pooling and precision |
| **10E** | `10e_freeze_common_predictive_universe.py` | Preserve a common, representation-independent predictive eligibility set |
| **10F** | `../training/10f_freeze_inner_group_splits.py` | Freeze deterministic, homology-group-aware inner train/validation partitions |

The predictive eligibility contract retained in Stage 10E is **canonical 20-amino-acid sequences of length at most 1,024 residues**. Ineligible records are reported; they are not reassigned between outer train/test. The Stage-10F inner split uses cluster-aware partitioning of the *outer-training* portion only. Representation availability must not silently redefine the eligible biological universe.

The representation registry is [`../config/representation_registry/representation_registry.tsv`](../config/representation_registry/representation_registry.tsv). See [Benchmarking and analysis](BENCHMARKING_AND_ANALYSIS.md) for the specific representations included in the paper's reported factorial comparison.

## 4. Training and result aggregation

| Stage | Existing source | Scientific responsibility |
| --- | --- | --- |
| **11 supervised** | [`11_train_supervised_factorial.py`](../src/training/11_train_supervised_factorial.py) | Fixed configuration panels, inner-validation selection, full outer-train refit and held-out test predictions |
| **11 PU** | [`11_train_pu_factorial.py`](../src/training/11_train_pu_factorial.py) | Separate PU regimes and methods while preserving the distinction between positive and unlabeled observations |
| **12 supervised QC** | [`12a_aggregate_and_qc_supervised.py`](../src/analysis/12a_aggregate_and_qc_supervised.py) | Aggregate predictions/metrics and audit experiment completion |
| **12 supervised tables** | [`12b_build_supervised_analysis_tables.py`](../src/analysis/12b_build_supervised_analysis_tables.py) | Construct factorial R/S contrasts and algorithm/representation comparisons |
| **12 PU QC** | [`12a2_audit_pu_training.py`](../src/analysis/12a2_audit_pu_training.py) | Audit PU experiment integrity |
| **12 PU tables** | [`12b2_build_pu_analysis_tables.py`](../src/analysis/12b2_build_pu_analysis_tables.py) | Aggregate PU diagnostics and contrasts without interpreting U as confirmed negative |

The original shell entry points are preserved in [`../runs/`](../runs/). Task-specific flags and data expectations remain defined by each script's own CLI.

## 5. Original directory convention

```text
data/raw/                   # Stage 01–02 source-derived records
data/processed/             # Evidence states and exact-sequence P/N/U master
data/tasks/                 # Frozen task-specific cohorts
data/redundancy/            # Independent R identity regimes
data/splits/                # Factorial R × S outer partitions
data/representations/       # Sequence universe, descriptors, one-hot and PLMs
data/predictive_universe/   # Frozen 10E eligible entities
data/training_splits/       # Frozen 10F inner partitions
results/training_supervised/
results/training_pu/
results/analysis_supervised/
results/analysis_pu/
```

These are **original run locations**, not claims that all generated datasets have been bundled into `release_repo`. The release contains source scripts and selected models/examples, while the large derived matrices and experiment trees remain outside it. Historical scripts may also expect the original working directory layout.

## 6. Detailed methodological records

- [Stages 01–04](PhageLysClass_Stages01-04_Data_Collection_Methods.md)
- [Task construction](PhageLysClass_Task_Definition_Methodology.md)
- [Stage 07 redundancy](PhageLysClass_Stage07_Redundancy_Methodology.md)
- [Training protocol v2](PhageLysClass_Training_Methodological_Design_v2.md)

**Documentation-only release:** this guide does not require, imply or report any new runs.

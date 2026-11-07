# Benchmarking and Analysis

**PhageLysClass | Experimental design and analytical notebooks**

The PhageLysClass benchmark evaluates the dependence of predictive performance on **biological evidence definition, training-task formulation, sequence redundancy, homology separation, representation and learning strategy**. It is not simply a leaderboard of classifier scores.

## 1. Different biological questions require different tasks

| Task family | Question addressed | Source of positive or target evidence |
| --- | --- | --- |
| `global_easy`, `global_hard`, `global_mixed` | Can known lytic proteins be distinguished from increasingly difficult, evidence-defined negative candidates? | Curated PhageLysData Core |
| `subtype_binary` | Is an individual lytic subtype compatible with a candidate sequence, relative to an explicit negative background? | Endolysin, Depolymerase or VAL subtype evidence |
| `multiclass_clean` | Among clean lytic subclasses, which functional class is most compatible? | Endolysin / Depolymerase / VAL |
| `pu_global` | What changes when uncertain background proteins are treated as unlabeled instead of confirmed negative? | Known P against conservative or broad U populations |
| `challenge_sets` | How does the system behave on ambiguous or non-exclusive evidence? | Segregated challenge cohorts; not ordinary training negatives |

The `subtype_binary` and `multiclass_clean` tasks answer different questions. A multiclass probability is conditional on consideration within the lytic-class space. The released inference interface only emits such subclass probabilities if its global core predicts `lytic`.

## 2. Factorial redundancy and homology design

The experiment independently varies two sequence-relatedness factors:

| Factor | Regimes | Meaning |
| --- | --- | --- |
| **R — within-task redundancy** | `R90`, `R70`, `R50`, `R30` | Maximum redundancy allowed when independently reducing task strata |
| **S — homology-aware separation** | `S90`, `S70`, `S50`, `S30` | Sequence-relatedness criterion used to form non-overlapping outer-split homology groups |

This creates a 4 × 4 **16-cell R × S surface** for each trainable task. The frozen outer test is kept separate from the homology-group-aware inner validation used for model configuration selection. Cross-cell comparisons may evaluate different test members and should not automatically be described as paired per-sequence tests. The analysis code retains both native-test contrasts and specialized sensitivity analyses.

## 3. Representations and predictive learners

The main factorial benchmark uses four registered representation families: `roxy_48`, `onehot_1024`, `esm2_t6_8m`, and `esm2_t30_150m`. The code also records broader representation-generation capabilities, but the existence of an encoder in the source does not imply it was included in every final benchmark cell.

The supervised trainer implements predefined panels of **Logistic Regression, Random Forest, K-Nearest Neighbors, XGBoost and Multilayer Perceptron**, with validation-only configuration selection. The positive–unlabeled branch is separate: its frozen base learners are LR and XGB and its methods include a naive P-versus-U diagnostic, BaggingPU and Elkan–Noto.

**Critical interpretation:** scores evaluated against observed P/U labels are *observed-label diagnostics*. The unlabeled set is not a certified biological negative set, so these scores must not be presented as unbiased biological discrimination metrics.

## 4. Where the analytical code lives

### Dataset exploration — `eda/`

| Notebook | Focus |
| --- | --- |
| `E0_build_phagelytic_EDA_universe_FIXED.ipynb` | Consolidated exploration universe |
| `E1_dataset_and_evidence_characterization.ipynb` | Class/evidence composition, sequence length and coverage |
| `E2_ROBUST_v2.ipynb` | Composition-level characterization |
| `E3_physicochemical_landscape_CORRECTED.ipynb` | Physicochemical descriptor landscapes |
| `E4_descriptor_space_geometry_ROBUST.ipynb` | Geometry and neighborhood structure in interpretable descriptors |
| `E5_v2_latent_landscape.ipynb` | Protein-language-model latent geometry |
| `E6_cross_representation_background_consistency.ipynb` | Cross-representation and background sensitivity |
| `PhageLysData_EDA_panel_builder_v9_F_scatter_wider.ipynb` | Final EDA visualization assembly |

### Model analysis — `notebooks/`

The retained research-question notebooks can be read in four broad groups. Several filenames contain final/fixed/v2 suffixes that reflect their original development history; they have not been renamed to avoid breaking references.

| Group | Typical research questions | Scientific focus |
| --- | --- | --- |
| **Task, evidence and R/S effects** | RQ1–RQ5 | Evidence definitions, redundancy, homology and interactions |
| **Representations, learners and PU** | RQ6–RQ15 | Capacity, model choice, task formulation, PU regime and method |
| **Robustness and strategy** | RQ16–RQ25 | Common-test sensitivity, uncertainty, runtime, transitions and cross-regime synthesis |
| **Model selection and deployment** | RQ26–RQ29 | Candidate portfolio, ensemble complementarity, frozen architecture, export and characterization |

A final benchmarking-panel builder is also retained in `notebooks/`. For deployment-specific interpretation, start with `RQ28_final_phagelysclass_architecture_READY.ipynb`, `29_RQ29A_production_model_export.ipynb`, and `29_RQ29B_model_characterization_SCHEMA_FIXED.ipynb`.

## 5. How the production release relates to the benchmark

The final production artifact is *not* one arbitrary winner from a single R/S cell. The released `models/phagelysclass_v1/manifest.json` records an architecture with two core functional decisions and three optional subtype-compatibility components. It also documents how selected configurations were consolidated using **validation-only evidence across benchmark cells**, without using outer-test results for the selection rule. The original final refit uses the frozen eligible task datasets; R/S are benchmark stress-test axes, not deployment parameters.

The model bundle contains its own manifest, checksums, task metadata and fitted estimators. See [Model and inference](MODEL_AND_INFERENCE.md) for the exact release semantics.

## 6. Availability of intermediate outputs

The notebooks were preserved as research code and analytical records. **Complete `data/` and `results/` trees were deliberately not copied into this compact repository folder.** Some notebooks read materialized tables and caches from those original locations. Do not assume a notebook can be restarted in a clean checkout without restoring the referenced inputs. Likewise, published figures stored under `figures/` are supplied artifacts, not evidence of a fresh rerun.

The extended methods are in [Predictive training and factorial analysis design](PhageLysClass_Training_Methodological_Design_v2.md).

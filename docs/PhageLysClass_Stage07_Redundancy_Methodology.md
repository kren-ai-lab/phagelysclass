# PhageLysClass — Redundancy Reduction Methodology

**Stage:** 07 — Task-specific sequence redundancy reduction  
**Strategy:** BioSieve + MMseqs2  
**Identity regimes:** H90, H70, H50, H30  
**Status:** Stage 07 — PASS  
**Input:** frozen task definitions from Stage 06  
**Output:** redundancy-controlled datasets for downstream homology-aware splitting

---

## 1. Purpose

The objective of Stage 07 is to quantify and reduce within-task sequence redundancy before constructing train/validation/test partitions.

This stage is intentionally separated from the subsequent splitting stage.

The workflow is:

```text
Frozen task definition
        ↓
Task-specific redundancy reduction
        ↓
H90 / H70 / H50 / H30 datasets
        ↓
Homology-aware splitting
```

The purpose is therefore **not** to create leakage-free train/test partitions yet. Instead, it generates progressively less redundant versions of each biological task while preserving its frozen class definition.

No balancing, representation calculation, train/test splitting, or model training is performed here.

---

# 2. Software strategy

Redundancy reduction is performed with:

```text
BioSieve
    └── MMseqs2 redundancy strategy
```

BioSieve provides the workflow layer, while MMseqs2 performs the sequence-similarity clustering/reduction.

Four sequence-identity regimes are generated:

```text
H90 = 90% identity
H70 = 70% identity
H50 = 50% identity
H30 = 30% identity
```

Each threshold represents an independent dataset version.

---

# 3. Independent threshold construction

The four redundancy levels are generated independently from the original Stage 06 dataset.

The workflow is:

```text
original task ──→ H90
              ├─→ H70
              ├─→ H50
              └─→ H30
```

The workflow is **not**:

```text
H90 → H70 → H50 → H30
```

This avoids making lower-identity datasets dependent on representatives selected at a previous threshold.

Each regime is therefore a direct reduction of the same frozen biological task.

---

# 4. Independent reduction by biological stratum

Redundancy is reduced independently within the biological groups composing each task.

This decision is central to the methodology.

For example, for the global mixed task:

```text
P_core  ──MMseqs2──→ representatives
N_easy  ──MMseqs2──→ representatives
N_hard  ─MMseqs2──→ representatives
                       ↓
                 recombined task
```

Similarly, the clean multiclass task is reduced independently for:

```text
endolysin
depolymerase
virion_associated_lytic_enzyme
```

and PU-learning datasets are reduced independently for their P and U evidence states.

### Why reduce within strata?

A global reduction across all classes could allow one biological group to eliminate representatives from another group.

For example:

- a negative sequence could become the representative of a cluster containing a positive sequence;
- the abundant endolysin class could displace sequences from a minority subtype;
- an unlabeled sequence could replace a known positive in a PU dataset.

To prevent this, representatives are selected within biological strata and subsequently recombined.

This preserves the biological task definition frozen in Stage 06.

---

# 5. Cross-stratum homology is not removed in Stage 07

A consequence of within-stratum reduction is that homologous proteins can still exist between different biological groups.

For example:

```text
P_core ↔ N_easy
P_core ↔ N_hard
P_core ↔ U
endolysin ↔ depolymerase
endolysin ↔ VAL
```

This is intentional.

Stage 07 addresses **within-group redundancy**, while the next methodological stage will address **cross-partition sequence relatedness and leakage** through homology-aware splitting.

Therefore:

> A dataset labeled H30 at Stage 07 is internally reduced at the 30% identity regime within each biological stratum, but this does not yet imply that future train and test sequences are separated below 30% identity.

---

# 6. Tasks processed

Redundancy reduction was applied independently to all trainable Stage 06 tasks:

```text
global_easy
global_hard
global_mixed

subtype_binary/
├── endolysin
├── depolymerase
└── virion_associated_lytic_enzyme

multiclass_clean

pu_global/
├── conservative
└── broad
```

The challenge sets were intentionally excluded from this stage.

They remain frozen external evaluation datasets and are not altered during ordinary task redundancy reduction.

---

# 7. Output organization

Each task receives an independent redundancy folder.

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

Subtype tasks retain a subtype-specific level:

```text
data/redundancy/subtype_binary/
├── endolysin/
├── depolymerase/
└── virion_associated_lytic_enzyme/
```

PU datasets similarly retain the learning regime:

```text
data/redundancy/pu_global/
├── conservative/
└── broad/
```

---

# 8. Output files

## `dataset.tsv.gz`

Contains the representatives retained by BioSieve/MMseqs2 for the corresponding threshold.

The original Stage 06 metadata are preserved, together with the redundancy threshold and biological redundancy group.

---

## `mapping.tsv.gz`

Stores the redundancy mapping produced during reduction, allowing removed sequences to be traced back to their representative sequence or cluster.

This file provides provenance for the reduction process.

---

## `summary.tsv`

Reports the number of sequences before and after redundancy reduction for every biological stratum.

Typical fields include:

```text
group
before
after
removed
retained_fraction
```

An overall task-level row is also generated.

---

## `manifest.json`

Stores the computational provenance of each redundancy run, including:

- task identifier,
- input dataset,
- input SHA-256,
- identity threshold,
- biological grouping strategy,
- BioSieve version,
- MMseqs2 version,
- counts before and after reduction,
- BioSieve group-level reports, and
- explicit indication that no balancing or splitting was performed.

---

## `redundancy_overview.tsv`

Summarizes the four identity regimes for a task:

```text
threshold
before
after
removed
retained_fraction
```

This provides a compact view of how sequence diversity changes from H90 to H30.

---

# 9. Implementation details

The Stage 07 scripts are organized independently by task:

```text
src/redundancy/
├── _redundancy_utils.py
├── 07a_reduce_global_easy.py
├── 07b_reduce_global_hard.py
├── 07c_reduce_global_mixed.py
├── 07d_reduce_subtype_binary.py
├── 07e_reduce_multiclass_clean.py
└── 07f_reduce_pu_global.py
```

The shared utility layer handles:

- BioSieve/MMseqs2 preflight checks,
- identity-threshold configuration,
- per-stratum execution,
- representative recovery,
- mapping integration,
- summary generation,
- output manifests, and
- reproducibility metadata.

The installed BioSieve MMseqs2 threshold parameter is probed before execution so that the pipeline does not silently assume an incompatible strategy configuration.

---

# 10. Global supervised tasks

## 10.1 Global easy

Starting dataset:

```text
P_core  = 11,867
N_easy  = 54,711
Total   = 66,578
```

After redundancy reduction:

| Threshold | P | N_easy | Total |
|---|---:|---:|---:|
| H90 | 5,013 | 25,734 | 30,747 |
| H70 | 3,654 | 18,354 | 22,008 |
| H50 | 2,555 | 12,583 | 15,138 |
| H30 | 1,705 | 7,511 | 9,216 |

The task remains viable under all four redundancy regimes.

---

## 10.2 Global hard

Starting dataset:

```text
P_core = 11,867
N_hard = 4,327
Total  = 16,194
```

After redundancy reduction:

| Threshold | P | N_hard | Total |
|---|---:|---:|---:|
| H90 | 5,013 | 1,894 | 6,907 |
| H70 | 3,654 | 1,429 | 5,083 |
| H50 | 2,555 | 1,106 | 3,661 |
| H30 | 1,705 | 860 | 2,565 |

The hard-negative universe retains substantial diversity even at H30.

---

## 10.3 Global mixed

Starting dataset:

```text
P_core              = 11,867
N_easy + N_hard     = 59,038
Total               = 70,905
```

After redundancy reduction:

| Threshold | P | N_easy | N_hard | Total |
|---|---:|---:|---:|---:|
| H90 | 5,013 | 25,734 | 1,894 | 32,641 |
| H70 | 3,654 | 18,354 | 1,429 | 23,437 |
| H50 | 2,555 | 12,583 | 1,106 | 16,244 |
| H30 | 1,705 | 7,511 | 860 | 10,076 |

The counts reproduce exactly the independently reduced strata from the easy and hard tasks, confirming consistent task construction.

---

# 11. Subtype-specific tasks

## 11.1 Endolysin

Positive counts:

```text
9,888 original
4,132 H90
2,965 H70
2,002 H50
1,266 H30
```

Endolysin sequences show substantial internal redundancy, but the task remains well populated even at H30.

---

## 11.2 Depolymerase

Positive counts:

```text
712 original
511 H90
453 H70
404 H50
350 H30
```

Depolymerases retain a comparatively high proportion of their sequences under increasingly stringent redundancy reduction.

At H30, approximately half of the original positive universe remains.

This indicates considerably greater sequence diversity than observed for endolysins.

---

## 11.3 Virion-associated lytic enzymes

Positive counts:

```text
241 original
91 H90
82 H70
65 H50
49 H30
```

VAL constitutes the main Low-N task after redundancy reduction.

The strongest reduction occurs already at H90, suggesting substantial high-identity redundancy in the original VAL positive universe.

H50 and especially H30 are retained as valid experimental regimes, but their downstream evaluation will require explicit Low-N consideration.

---

# 12. Clean multiclass task

The clean multiclass dataset contains only:

```text
endolysin
depolymerase
virion_associated_lytic_enzyme
```

After reduction:

| Threshold | Endolysin | Depolymerase | VAL | Total |
|---|---:|---:|---:|---:|
| H90 | 4,132 | 511 | 91 | 4,734 |
| H70 | 2,965 | 453 | 82 | 3,500 |
| H50 | 2,002 | 404 | 65 | 2,471 |
| H30 | 1,266 | 350 | 49 | 1,665 |

Because the classes are reduced independently, these counts exactly reproduce the corresponding subtype-specific positive universes.

The multiclass task becomes progressively more imbalanced as sequence redundancy is removed, with VAL constituting the limiting class.

---

# 13. PU-learning tasks

## 13.1 Conservative PU

Definition:

```text
P = P_core
U = U_explicit
```

After reduction:

| Threshold | P | U_explicit | Total |
|---|---:|---:|---:|
| H90 | 5,013 | 147,807 | 152,820 |
| H70 | 3,654 | 119,419 | 123,073 |
| H50 | 2,555 | 101,357 | 103,912 |
| H30 | 1,705 | 87,448 | 89,153 |

The unlabeled universe retains substantially more sequence diversity than the known-positive Core.

This causes the P:U imbalance to increase under more stringent identity reduction.

No balancing is performed at this stage.

---

## 13.2 Broad PU

Definition:

```text
P = P_core
U = U_explicit + U_other_unresolved
```

After reduction:

| Threshold | P | U_explicit | U_other_unresolved | Total |
|---|---:|---:|---:|---:|
| H90 | 5,013 | 147,807 | 69,415 | 222,235 |
| H70 | 3,654 | 119,419 | 52,176 | 175,249 |
| H50 | 2,555 | 101,357 | 39,824 | 143,736 |
| H30 | 1,705 | 87,448 | 28,328 | 117,481 |

Again, `U` remains unlabeled rather than negative.

The broad PU universe provides a large and sequence-diverse open-world background for later PU-learning experiments.

---

# 14. Main observations from Stage 07

The redundancy analysis reveals clear differences in the sequence diversity of the biological groups.

### Positive Core

The global positive Core is highly redundant:

```text
11,867 original
→ 5,013 H90
→ 3,654 H70
→ 2,555 H50
→ 1,705 H30
```

### Endolysins

Endolysins show particularly strong redundancy, retaining approximately 13% of their initial sequences at H30.

### Depolymerases

Depolymerases retain approximately half of their starting sequences at H30, suggesting a considerably more diverse sequence space.

### VAL

VAL contains substantial high-identity redundancy and becomes a Low-N classification problem after reduction:

```text
91 H90
82 H70
65 H50
49 H30
```

### Unlabeled proteins

The U universe is substantially more diverse than the known-positive Core, especially for `U_explicit`.

This difference becomes increasingly visible as identity thresholds become more stringent.

---

# 15. Interpretation of the H90–H30 regimes

The four redundancy regimes should be interpreted as progressively stronger controls on internal sequence redundancy:

```text
H90
    relatively permissive redundancy control

H70
    moderate family-level reduction

H50
    stringent redundancy reduction

H30
    highly stringent / remote sequence regime
```

They should **not yet** be interpreted as train/test homology guarantees.

That distinction belongs to the subsequent split stage.

---

# 16. Decisions frozen after Stage 07

The following methodological decisions are frozen:

1. Redundancy reduction is performed independently for every Stage 06 task.
2. Challenge sets are not reduced during Stage 07.
3. Four identity regimes are retained: H90, H70, H50, and H30.
4. Every threshold is generated independently from the original task.
5. Redundancy is reduced independently within biological strata.
6. Biological class definitions are never altered during redundancy reduction.
7. Cross-stratum homology is intentionally left unresolved at this stage.
8. BioSieve with the MMseqs2 strategy is used as the reduction workflow.
9. Removed sequences remain traceable to representatives through mapping files.
10. No balancing is performed.
11. No train/validation/test split is performed.
12. No representation or model information influences representative selection.
13. VAL H50/H30 are retained but flagged as Low-N regimes.
14. All H90–H30 datasets are propagated to the downstream split stage.

---

# 17. Transition to Stage 08

The outputs of Stage 07 form the direct inputs to the homology-aware splitting stage.

For each task and each redundancy regime:

```text
task
└── H90 / H70 / H50 / H30
        ↓
homology-aware split
        ↓
train / validation / test
```

The next stage must determine partitions while explicitly controlling sequence relatedness across train, validation, and test sets.

Only after those partitions are frozen should the project proceed to:

```text
dataset characterization
        ↓
protein representations
        ↓
one-hot baseline
        ↓
protein language model embeddings
        ↓
predictive modeling
```

This preserves the methodological hierarchy:

```text
biological definition
        ↓
task definition
        ↓
redundancy control
        ↓
homology-aware splitting
        ↓
data characterization
        ↓
representation
        ↓
modeling
```

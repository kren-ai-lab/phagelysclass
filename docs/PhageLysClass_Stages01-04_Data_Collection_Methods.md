# PhageLysClass — Data Collection and Evidence-State Construction

**Methods documentation for Stages 01–04**  
**Project:** PhageLysClass  
**Scope:** data acquisition, biological evidence-state construction, and pre-task audit  
**Status:** frozen for downstream task construction  
**Version:** methods freeze v1  
**Date:** 2026-08-22

---

## 1. Purpose and scope

This document records the complete methodology implemented in the first four data-construction stages of the PhageLysClass pipeline. These stages establish the biological universes used later for supervised classification, negative-definition experiments, homology-aware evaluation, and discovery analyses.

The workflow was designed around a central methodological principle:

> **absence of evidence for phage lytic activity is not equivalent to evidence of non-lytic function.**

Accordingly, the early pipeline does not force all non-positive phage proteins into a negative class. Instead, it separates:

- evidence-supported positive lytic proteins,
- proteins with explicit alternative functions,
- lysis-adjacent and structurally ambiguous proteins,
- proteins with direct target-like annotations outside the positive resource,
- explicitly uncharacterized proteins, and
- broader unresolved proteins.

The output of Stages 01–04 is therefore an **audited evidence-state universe**, not yet a final machine-learning dataset.

This separation is important because the central questions of the paper concern how negative definition, task formulation, protein representation, and sequence relatedness influence reliable recognition of phage lytic enzymes.

---

## 2. Workflow overview

The implemented data-collection workflow is:

```text
PhageLysData v1.0
        |
        | Stage 01
        v
Evidence-supported Core positives
        |
        +-----------------------------------------------+
                                                        |
UniProtKB phage protein universe                        |
        |                                               |
        | Stage 02                                      |
        v                                               |
Broad taxonomy-defined raw background                   |
        |                                               |
        | Stage 03                                      |
        v                                               |
Local evidence-state classification --------------------+
        |
        | Stage 04
        v
Audit of rules, conflicts, annotation quality,
target overlap, and random manual samples
```

Stages 01–04 intentionally do **not** perform:

- homology reduction,
- sequence-similarity filtering,
- class balancing,
- train/validation/test splitting,
- representation calculation,
- model training,
- final task-specific positive/negative labeling, or
- conversion of unlabeled proteins into negatives.

Those operations are handled downstream.

---

# 3. Stage 01 — Extraction of the PhageLysData evidence-supported Core

## 3.1 Objective

Stage 01 isolates the evidence-supported PhageLysData Core as the positive biological universe for the classification project.

The input is the frozen PhageLysData v1.0 Core table:

```text
core/phagelysdata_core.parquet
```

The project root is supplied either through:

```bash
PHAGELYSDATA_ROOT
```

or explicitly by command-line argument.

The corresponding implementation is:

```text
src/dataset_processing/01_extract_phagelys_raw.py
```

## 3.2 Positive-universe philosophy

The classification project does not reconstruct positive labels from external annotations. Instead, the positive universe is inherited from the curated evidence model of PhageLysData.

This provides three advantages:

1. positive labels are linked to a resource-level evidence model rather than a single annotation field;
2. exact-sequence entities are already consolidated;
3. the positive universe is frozen independently of downstream machine-learning performance.

This prevents circular task construction in which labels are modified after observing model behavior.

## 3.3 Extracted fields

Only fields required for downstream task definition are retained:

| Field | Role |
|---|---|
| `sequence_entity_id` | stable exact-sequence entity identifier |
| `sequence` | amino-acid sequence |
| `canonical_target_class` | resource-level target-class assignment |
| `canonical_class_basis` | basis of the class assignment |
| `assigned_tier` | evidence tier |

No descriptors, embeddings, similarity measures, or model-derived values are propagated into the raw classification project at this stage.

## 3.4 Biological classes

The evidence-supported Core contains four biological target assertions:

- `endolysin`
- `virion_associated_lytic_enzyme`
- `depolymerase`
- `other_phage_envelope_lytic_enzyme`

The resource additionally contains the derived label:

- `multiple_nonexclusive`

`multiple_nonexclusive` is not treated as a biological function. It denotes exact-sequence entities for which more than one non-predictive target assertion is supported.

## 3.5 Frozen counts

Stage 01 recovered the complete PhageLysData Core:

| Group | Exact-sequence entities |
|---|---:|
| Total Core | 11,867 |
| Single-class primary | 10,841 |
| `multiple_nonexclusive` | 999 |
| `other_phage_envelope_lytic_enzyme` | 27 |
| Endolysin | 9,888 |
| Depolymerase | 712 |
| Virion-associated lytic enzyme | 241 |

The count contract was verified directly against the frozen PhageLysData release.

## 3.6 Outputs

```text
data/raw/phagelysdata/
├── phagelysdata_core_raw.tsv.gz
├── phagelysdata_single_class_raw.tsv.gz
├── phagelysdata_multiple_nonexclusive_raw.tsv.gz
├── phagelysdata_other_lytic_raw.tsv.gz
└── summary TSV
```

## 3.7 Methodological guardrails

Stage 01 applies no sequence similarity logic.

In particular:

- no sequences are merged by homology,
- no sequences are removed for redundancy,
- no class balancing is performed,
- no multiclass exclusivity is imposed beyond the frozen resource labels.

Exact amino-acid identity remains the entity definition.

---

# 4. Stage 02 — Construction of the broad UniProtKB phage-protein background

## 4.1 Objective

Stage 02 constructs a broad raw universe of phage proteins from UniProtKB.

The design intentionally separates **retrieval** from **biological interpretation**.

The retrieval query therefore does not contain terms such as:

- hypothetical,
- uncharacterized,
- unknown,
- capsid,
- holin,
- lysin,
- depolymerase, or
- terminase.

Functional semantics are introduced only after the broad taxonomy-defined universe has been downloaded.

The implementation is:

```text
src/dataset_processing/02_build_uniprot_background.py
```

## 4.2 Why retrieval was taxonomy-defined

An initial exploratory strategy based on semantic UniProt queries produced highly uneven and vocabulary-dependent counts. Standardized UniProt protein naming meant that terms such as “unknown”, “protein of unknown function”, or “conserved hypothetical” had poor recall compared with broad taxonomic retrieval.

The workflow was therefore redesigned around the following principle:

> **retrieve broadly by phage taxonomy/context first; classify annotation evidence locally and reproducibly second.**

This avoids hard-coding the eventual U/N definitions into the acquisition query.

## 4.3 Query registry

The acquisition registry contains two enabled scopes:

### Caudoviricetes

```text
taxonomy_id:2731619
```

### Named phages outside Caudoviricetes

```text
(taxonomy_id:10239)
AND NOT (taxonomy_id:2731619)
AND ((taxonomy_name:phage) OR (taxonomy_name:bacteriophage))
```

The registry is stored separately from the code:

```text
config/queries/uniprot_background_queries_taxonomic.tsv
```

This makes the acquisition universe auditable and modifiable without changing Python logic.

## 4.4 Query counts before download

The UniProt count audit returned:

| Scope | UniProt records |
|---|---:|
| Caudoviricetes | 564,228 |
| Named non-Caudoviricetes phages | 5,967 |
| **Total raw background** | **570,195** |

The two query scopes are disjoint by construction.

## 4.5 UniProt release contract

The expected UniProtKB release was fixed to:

```text
2026_02
```

The script checks the release exposed by the UniProt REST API and can fail on mismatch unless explicitly overridden.

This temporal freeze is important because later discovery analyses can use newer UniProt releases as a distinct temporal validation universe.

## 4.6 Retrieved metadata

For each UniProt accession, the raw acquisition retained:

- accession,
- entry identifier,
- protein name,
- gene names,
- organism name,
- organism taxonomy identifier,
- taxonomic lineage,
- lineage identifiers,
- virus hosts,
- reviewed status,
- protein existence,
- annotation score,
- fragment status,
- sequence length,
- amino-acid sequence,
- function comment,
- catalytic-activity comment,
- EC number,
- GO identifiers,
- keywords,
- creation date,
- modification date,
- retrieval query identifier, and
- retrieval scope.

This metadata is retained because evidence-state construction requires more than protein names alone.

## 4.7 Native UniProt cursor pagination

The final acquisition implementation uses the native UniProtKB `/search` endpoint with cursor pagination rather than the long-lived `/stream` endpoint.

Default configuration:

```text
page size       = 500
retries/page    = 8
timeout/page    = 120 s
```

Each page is:

1. downloaded independently,
2. validated as TSV,
3. written atomically to a compressed page file,
4. recorded in a checkpoint, and
5. followed using UniProt's `Link: rel="next"` cursor.

This design was chosen because the full Caudoviricetes universe is too large to rely robustly on one long HTTP stream.

## 4.8 Checkpoint and resume behavior

Pages are stored under:

```text
data/raw/uniprot_background/pages/<query_id>/
```

with a `checkpoint.json`.

If the process is interrupted, the same command resumes from the last completed page rather than downloading the complete query again.

A query is not marked complete until:

- cursor traversal terminates,
- all pages are merged,
- the merged row count equals the UniProt `X-Total-Results` value.

## 4.9 Raw-output contract

Stage 02 produced:

```text
data/raw/uniprot_background/
├── pages/
├── queries/
├── uniprot_background_raw.tsv.gz
├── uniprot_background_query_manifest.tsv
└── uniprot_background_summary.tsv
```

The final raw table contains:

```text
570,195 accession-level records
```

At this point:

- accessions are not collapsed by sequence,
- no U/N labels are assigned,
- no PhageLysData decontamination is applied,
- no homology filtering is performed.

---

# 5. Stage 03 — Local evidence-state classification of the UniProt background

## 5.1 Objective

Stage 03 converts the broad UniProt raw universe into auditable evidence-state candidates.

The implementation is:

```text
src/dataset_processing/03_classify_uniprot_background.py
```

The rule registry is external to the code:

```text
config/rules/uniprot_background_evidence_rules.tsv
```

The script processes the complete UniProt background while retaining raw annotation provenance.

## 5.2 Conceptual distinction between labels and evidence states

Stage 03 does **not** assign final machine-learning labels.

Instead it assigns evidence states describing what can safely be concluded from the available UniProt annotation.

The final states are:

| Evidence state | Interpretation |
|---|---|
| `known_p_overlap` | exact sequence already present in PhageLysData Core |
| `lytic_related` | direct target-like lytic annotation outside exact Core overlap |
| `ambiguous_structural` | structural/adsorption proteins that can overlap biologically with VAL/DPO function |
| `n_candidate_hard` | lysis-adjacent, generally non-enzymatic proteins |
| `n_candidate_easy` | explicit alternative/non-target phage function |
| `u_uncharacterized` | explicitly hypothetical/uncharacterized/unknown protein |
| `u_domain_unknown` | DUF/UPF or unresolved-domain naming |
| `u_low_information` | generic ORF/predicted/putative naming |
| `u_other_unresolved` | remaining proteins not safely classified by the rules |

The naming `n_candidate_*` is intentional. These records remain **candidate negative pools**, not biological ground-truth negatives.

## 5.3 Exact positive overlap

Before local evidence-state use downstream, every UniProt sequence is checked against the PhageLysData Core by exact amino-acid identity.

If an exact match is found:

```text
evidence_state = known_p_overlap
```

This state has precedence over all local UniProt rules.

The local matches are nevertheless retained for audit.

## 5.4 Rule precedence

The effective precedence is:

```text
exact PhageLysData Core overlap
    >
direct lytic-related evidence
    >
ambiguous structural / adsorption evidence
    >
hard-negative candidate evidence
    >
easy-negative candidate evidence
    >
explicit unlabeled evidence
    >
other unresolved
```

The practical rule priority is encoded numerically in the external registry.

## 5.5 Rule registry

The final active rules are documented below.


| Rule | State | Group | Priority | Fields inspected | Biological interpretation |
|---|---|---|---:|---|---|
| `LYTIC_002` | `lytic_related` | `depolymerase` | 995 | protein_name, cc_catalytic_activity, keyword | Depolymerase terminology from protein identity/catalytic evidence; free-text Function [CC] is not used as a primary trigger. |
| `AMB_001` | `ambiguous_structural` | `tail_adsorption` | 900 | protein_name, cc_function, keyword | Tail/adsorption proteins can carry virion-associated lytic or depolymerase activity; do not treat as clean negatives. |
| `NHARD_001` | `n_candidate_hard` | `lysis_adjacent` | 800 | protein_name, cc_function, keyword | Lysis-adjacent phage proteins that are usually not lytic enzymes; retain as hard-negative candidates. |
| `NEASY_001` | `n_candidate_easy` | `capsid_head` | 700 | protein_name, cc_function, keyword | Explicit capsid/head structural function. |
| `NEASY_002` | `n_candidate_easy` | `packaging` | 700 | protein_name, cc_function, keyword | Explicit genome-packaging machinery. |
| `NEASY_003` | `n_candidate_easy` | `replication` | 700 | protein_name, cc_function, keyword | Explicit nucleic-acid replication/maintenance function. |
| `NEASY_004` | `n_candidate_easy` | `recombination_integration` | 700 | protein_name, cc_function, keyword | Explicit recombination/integration function. |
| `NEASY_005` | `n_candidate_easy` | `transcription_regulation` | 700 | protein_name, cc_function, keyword | Explicit transcription/regulation function. |
| `NEASY_006` | `n_candidate_easy` | `metabolism_modification` | 700 | protein_name, cc_function, keyword | Explicit non-lytic catalytic/metabolic functions. Generic protease/peptidase terms are excluded because phage lytic enzymes can themselves be peptidases. |
| `NEASY_007` | `n_candidate_easy` | `assembly_scaffold` | 700 | protein_name, cc_function, keyword | Explicit virion assembly/scaffolding role. |
| `U_001` | `u_uncharacterized` | `uncharacterized` | 500 | protein_name | Explicitly unresolved protein naming. |
| `U_002` | `u_domain_unknown` | `duf_upf` | 490 | protein_name, keyword | Domain/family assignment without resolved biological function. |
| `U_003` | `u_low_information` | `orf_generic` | 480 | protein_name | Generic low-information naming. |
| `LYTIC_001` | `lytic_related` | `lytic` | 1000 | protein_name, cc_catalytic_activity, ec, keyword | Direct or catalytic lytic-enzyme terminology. Free-text Function [CC] is excluded to avoid assigning holins/spanins merely because they mention endolysins. |


## 5.6 Final rule-design corrections

The Stage 04 audit identified two important failure modes in the first rule draft.

### 5.6.1 Free-text function descriptions were not used to infer lytic identity

The initial direct-lytic rule inspected `Function [CC]`.

This caused holins and spanins to be classified as lytic-related merely because their function descriptions mention release or action of endolysins.

The final lytic rules therefore use identity/catalysis-oriented fields:

```text
protein_name
cc_catalytic_activity
EC
keyword
```

depending on the rule.

Free-text function comments are not used as the primary trigger for direct lytic identity.

### 5.6.2 Generic protease/peptidase terms were removed from clean-negative evidence

Generic terms such as:

```text
protease
peptidase
```

were initially included in a clean alternative-function rule.

This was rejected because some phage lytic enzymes are themselves peptidases, including cell-wall peptidases.

The final `NEASY_006` rule therefore retains explicit alternative metabolic/catalytic functions but does not use generic `protease` or `peptidase` as sufficient evidence for a negative candidate.

## 5.7 Classification metadata retained for audit

For every UniProt row, Stage 03 records:

- normalized sequence,
- SHA-256 sequence hash,
- observed sequence length,
- canonical-20-amino-acid status,
- local evidence state,
- primary matched rule,
- all matched rule IDs,
- all matched rule groups,
- number of matched rules,
- number of matched rule groups,
- rule-conflict flag,
- exact PhageLysData Core overlap flag,
- PhageLysData Core match provenance, and
- final evidence state after precedence.

This makes every state assignment traceable to explicit rules.

## 5.8 Processing strategy

The 570,195-row background is processed in chunks of:

```text
50,000 rows
```

The master output is written incrementally to avoid unnecessary memory pressure.

Optional state-specific files can also be written.

## 5.9 Final accession-level evidence-state counts

After the Stage 04-informed rule correction, the final accession-level counts were:

| Evidence state | UniProt rows | Percent |
|---|---:|---:|
| `u_uncharacterized` | 276,099 | 48.42% |
| `u_other_unresolved` | 181,449 | 31.82% |
| `n_candidate_easy` | 68,061 | 11.94% |
| `ambiguous_structural` | 13,394 | 2.35% |
| `u_domain_unknown` | 12,803 | 2.25% |
| `known_p_overlap` | 9,087 | 1.59% |
| `n_candidate_hard` | 6,444 | 1.13% |
| `u_low_information` | 2,467 | 0.43% |
| `lytic_related` | 391 | 0.07% |

These counts are accession-level and therefore are not yet the exact-sequence universe used for modeling.

## 5.10 Interpretation of major states

### `u_uncharacterized`

This is the largest single state and captures explicitly uncharacterized, hypothetical, or unknown proteins.

It is treated as unlabeled, not negative.

### `u_other_unresolved`

This group contains records not safely assignable by the evidence rules.

It can include proteins with meaningful annotations that are simply outside the current evidence taxonomy.

Therefore:

> `u_other_unresolved` means unresolved by the classification rules, not necessarily biologically uncharacterized.

It is not converted to a negative class by default.

### `n_candidate_easy`

This state contains explicit alternative phage functions such as packaging, capsid/head, replication, recombination, transcription/regulation, assembly, and selected non-target metabolic functions.

It is the cleanest candidate negative pool.

### `n_candidate_hard`

This state contains lysis-adjacent but generally non-enzymatic proteins such as:

- holins,
- antiholins,
- spanins, and
- pinholins.

These are intentionally separated because they are biologically close to the target process and may constitute a more challenging negative definition.

### `ambiguous_structural`

This state includes proteins such as:

- tailspikes,
- tail fibers,
- receptor-binding proteins,
- baseplate proteins,
- adsorption proteins, and
- host-recognition proteins.

These are not used as clean negatives because some virion-associated lytic enzymes and depolymerases are functionally or architecturally associated with the same virion structures.

### `lytic_related`

This state contains proteins with direct lytic-enzyme-like identity or catalytic annotation that are not exact-sequence matches to the PhageLysData Core.

They are excluded from the negative and unlabeled pools.

They are not automatically promoted to positive ground truth because their evidence has not passed through the PhageLysData evidence model.

---

# 6. Stage 04 — Evidence-state audit and rule validation

## 6.1 Objective

Stage 04 is a descriptive audit designed to determine whether the Stage 03 evidence states behave as intended before downstream task construction.

The implementation is:

```text
src/dataset_processing/04_audit_uniprot_background.py
```

It does not modify any label or state.

## 6.2 Audit components

The audit evaluates:

1. global state composition,
2. dominant primary rules by state,
3. combinations of matched rule groups,
4. reviewed/unreviewed distribution,
5. protein-existence distribution,
6. annotation-score distribution,
7. fragment status,
8. retrieval scope,
9. PhageLysData class composition among exact overlaps,
10. composition of `lytic_related`,
11. rule conflicts, and
12. reproducible random samples from every state.

By default, the manual-inspection sample is:

```text
75 records per evidence state
random seed = 42
```

## 6.3 Audit outputs

```text
data/processed/uniprot_background/audit/
├── audit_state_overview.tsv
├── audit_primary_rules_by_state.tsv
├── audit_rule_group_combinations.tsv
├── audit_metadata_distributions.tsv
├── audit_annotation_score_by_state.tsv
├── audit_known_p_overlap_classes.tsv
├── audit_lytic_related_protein_names.tsv
├── audit_lytic_related_rules.tsv
├── audit_rule_conflicts_by_state.tsv
├── audit_random_samples.tsv
└── audit_report.md
```

## 6.4 Exact-overlap composition

At the UniProt accession level, 9,087 rows exactly matched PhageLysData Core sequences.

Their PhageLysData class composition was:

| Canonical target class | UniProt rows |
|---|---:|
| Endolysin | 7,016 |
| `multiple_nonexclusive` | 1,271 |
| Depolymerase | 606 |
| Virion-associated lytic enzyme | 167 |
| Other phage envelope lytic enzyme | 27 |

These counts are accession-level and must not be interpreted as unique PhageLysData entity coverage because several UniProt accessions can share one exact sequence.

Unique sequence-level overlap is calculated downstream.

## 6.5 Audit-driven rule correction

The first audit identified a severe conflict rate in `lytic_related`:

```text
907 / 1,376 = 65.9%
```

Manual inspection showed that holins and spanins were frequently classified as lytic-related because their `Function [CC]` descriptions mention endolysins.

After the rule correction described in Section 5.6:

```text
lytic_related = 391 records
conflicts     = 2 / 391 = 0.5%
```

This demonstrated that the corrected rules removed the major source of annotation leakage.

## 6.6 Final conflict profile

The final audit produced:

| Evidence state | Rule conflicts | Total | Fraction |
|---|---:|---:|---:|
| `known_p_overlap` | 1,352 | 9,087 | 14.9% |
| `lytic_related` | 2 | 391 | 0.5% |
| `ambiguous_structural` | 876 | 13,394 | 6.5% |
| `n_candidate_hard` | 2 | 6,444 | <0.1% |
| `n_candidate_easy` | 2,275 | 68,061 | 3.3% |
| `u_uncharacterized` | 0 | 276,099 | 0% |
| `u_domain_unknown` | 0 | 12,803 | 0% |
| `u_low_information` | 0 | 2,467 | 0% |
| `u_other_unresolved` | 0 | 181,449 | 0% |

The relatively high conflict rate among `known_p_overlap` does not invalidate those records because PhageLysData Core evidence has precedence over local UniProt text rules.

Instead, these cases illustrate why the evidence-supported positive resource is kept independent of the local background classifier.

## 6.7 Annotation-score profile

Median UniProt annotation scores were low across most states, which is expected for large phage-protein universes.

The final audit reported:

| State | Mean annotation score | Median |
|---|---:|---:|
| `known_p_overlap` | 1.49 | 1 |
| `lytic_related` | 1.00 | 1 |
| `ambiguous_structural` | 1.07 | 1 |
| `n_candidate_hard` | 1.16 | 1 |
| `n_candidate_easy` | 1.23 | 1 |
| `u_uncharacterized` | 1.00 | 1 |
| `u_domain_unknown` | 1.00 | 1 |
| `u_low_information` | 1.00 | 1 |
| `u_other_unresolved` | 1.05 | 1 |

This reinforces the decision not to equate missing target annotation with negative evidence.

---

# 7. Final Stage 01–04 data state

At the end of Stage 04, the project contains two distinct biological sources:

## 7.1 Evidence-supported positive source

```text
PhageLysData Core
11,867 exact-sequence entities
```

## 7.2 Broad background source

```text
UniProtKB phage universe
570,195 accession-level records
```

with locally audited evidence states:

```text
known_p_overlap                 9,087
lytic_related                     391
ambiguous_structural           13,394
n_candidate_hard                6,444
n_candidate_easy               68,061
u_uncharacterized             276,099
u_domain_unknown               12,803
u_low_information               2,467
u_other_unresolved            181,449
```

These are accession-level counts.

Exact-sequence consolidation is deliberately handled in the next stage and is not part of the Stage 01–04 freeze documented here.

---

# 8. Reproducibility and data lineage

## 8.1 Configuration outside code

Biological acquisition and classification decisions are stored independently from executable logic:

```text
config/
├── queries/
│   └── uniprot_background_queries_taxonomic.tsv
└── rules/
    └── uniprot_background_evidence_rules.tsv
```

This separation allows:

- transparent review of biological decisions,
- controlled updates without rewriting code,
- versioning of the evidence model, and
- straightforward reporting in supplementary methods.

## 8.2 Frozen rule provenance

Stage 03 writes the active rule registry beside the classified outputs:

```text
uniprot_background_evidence_rules_frozen.tsv
```

This ensures that the exact rules used for a given processed dataset are preserved even if the live configuration later changes.

## 8.3 Deterministic audit sampling

Manual-inspection samples use:

```text
random seed = 42
```

so the same records can be regenerated.

## 8.4 Sequence identity

Exact amino-acid sequence identity is used only for:

- PhageLysData entity identity,
- exact overlap detection between PhageLysData and UniProt.

No homology threshold is used during Stages 01–04.

This distinction is important because redundancy reduction and homology-aware splitting are experimental variables evaluated later.

---

# 9. Methodological rationale for the paper

The key methodological contribution of the data-construction workflow is that the negative universe is not defined by complementing the positive set.

A naïve binary formulation would effectively assume:

```text
not annotated as lytic = negative
```

The implemented workflow instead distinguishes at least four conceptually different situations:

```text
explicit alternative function
        -> candidate negative

lysis-adjacent / structurally ambiguous function
        -> challenge / hard candidate

explicitly unknown or hypothetical protein
        -> unlabeled

direct target-like annotation outside the curated positive resource
        -> exclude from negative/unlabeled construction
```

This allows the later experiments to ask whether apparent classifier performance depends on how strongly the negative class is biologically supported.

---

# 10. Paper-ready Methods draft

The following text can later be adapted directly into the manuscript.

## Data sources and positive reference universe

An evidence-supported positive universe was derived from the frozen PhageLysData v1.0 Core resource. Exact amino-acid sequence entities were retained together with their canonical target-class assignment, assignment basis, and evidence tier. The Core comprised 11,867 unique sequences, including 9,888 endolysins, 712 depolymerases, 241 virion-associated lytic enzymes, 27 other phage envelope lytic enzymes, and 999 entities carrying multiple non-exclusive target assertions. The `multiple_nonexclusive` designation was treated as a derived resource-level state rather than as an independent biological class. No similarity-based filtering or task-dependent relabeling was applied during positive-universe extraction.

## Broad phage background acquisition

A broad phage-protein background was retrieved from UniProtKB release 2026_02 using taxonomy-based queries rather than functional annotation terms. The primary query comprised all proteins assigned to Caudoviricetes (`taxonomy_id:2731619`), while a supplementary query retrieved named phage or bacteriophage proteins outside Caudoviricetes. Together, the queries returned 570,195 accession-level records. UniProt metadata included protein and gene names, taxonomy and host information, reviewed status, protein-existence evidence, annotation score, fragment status, sequence, function and catalytic-activity comments, EC numbers, GO terms, keywords, and record dates. Retrieval was implemented with the native UniProt REST search endpoint using cursor pagination with 500 records per page, per-page retries, checkpointing, and validation of the downloaded row count against `X-Total-Results`.

## Evidence-state assignment

UniProt records were classified locally using a version-controlled rule registry rather than function-specific retrieval queries. Exact sequence matches to the PhageLysData Core were marked as known positive overlaps and given precedence over local annotation rules. Remaining records were assigned to direct lytic-related, structurally ambiguous, hard-negative-candidate, easy-negative-candidate, explicitly unlabeled, or unresolved evidence states. Lysis-adjacent proteins such as holins, antiholins, and spanins were retained separately from clean alternative-function negatives, while tailspike, tail-fiber, receptor-binding, baseplate, adsorption, and host-recognition proteins were treated as structurally ambiguous because of their potential overlap with virion-associated lytic or depolymerase functions. Proteins described as hypothetical, uncharacterized, unknown, DUF/UPF-containing, or generic ORFs were retained as unlabeled rather than negative.

## Rule audit

Evidence-state assignments were audited using rule-level conflict counts, annotation metadata, exact-overlap composition, and reproducible random samples from each state. The initial direct-lytic rule was found to overclassify holins and spanins because free-text function descriptions commonly mention endolysins. Direct lytic identity was therefore restricted to protein-identity and catalytic-evidence fields, reducing the lytic-related conflict rate from 65.9% to 0.5%. Generic `protease` and `peptidase` terminology was also removed as sufficient clean-negative evidence because phage lytic enzymes can themselves exhibit peptidase activity. Following this audit, the accession-level background comprised 68,061 easy-negative candidates, 6,444 hard-negative candidates, 13,394 structurally ambiguous proteins, 391 lytic-related records, and 372,818 explicitly or otherwise unresolved proteins. These states were frozen before exact-sequence consolidation, task construction, homology filtering, balancing, or model evaluation.

---

# 11. What these stages do not establish

The following claims should not be made from Stages 01–04 alone:

- that every `n_candidate_easy` sequence is experimentally confirmed non-lytic;
- that every `n_candidate_hard` sequence is truly negative;
- that `u_uncharacterized` or `u_other_unresolved` are negative;
- that `lytic_related` records are validated positives;
- that annotation-state frequencies represent biological prevalence;
- that accession-level UniProt counts equal unique protein-sequence counts.

Instead, the stages establish an evidence-aware foundation from which controlled task definitions can be constructed.

---

# 12. Frozen decisions carried forward

The following decisions are considered frozen for downstream PhageLysClass experiments:

1. **PhageLysData Core is the primary positive universe.**
2. **UniProt retrieval is taxonomy-defined, not function-defined.**
3. **Absence of lytic annotation is not negative evidence.**
4. **Exact PhageLysData overlap has precedence over local annotation rules.**
5. **Lytic-related UniProt records are excluded from negative and unlabeled pools.**
6. **Tail/adsorption-associated proteins remain structurally ambiguous.**
7. **Holins/spanins and related lysis-adjacent proteins form a hard-negative candidate stratum.**
8. **Explicit alternative-function proteins form an easy-negative candidate stratum.**
9. **Explicitly unknown and unresolved proteins remain unlabeled.**
10. **No homology filtering occurs before biological evidence-state construction.**
11. **No task definition, balance, or split is allowed to alter these biological evidence states.**
12. **All later task definitions must be derived from these frozen universes rather than redefining them in place.**

---

## 13. Transition to downstream stages

The next stage operates on the audited evidence-state universe and performs exact-sequence consolidation into formal P/N/U candidate universes.

Only after that freeze should the project construct:

- global lytic-recognition tasks,
- subtype-specific binary tasks,
- clean multiclass tasks,
- multi-nonexclusive challenge sets,
- hard-negative challenge sets,
- PU-learning configurations,
- homology-aware generalization experiments, and
- later discovery analyses.

This separation preserves a clean data lineage:

```text
biological evidence
        ->
evidence-state universe
        ->
exact-sequence P/N/U freeze
        ->
task definition
        ->
homology/redundancy design
        ->
representation/model benchmarking
```


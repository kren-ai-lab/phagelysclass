# PhageLysClass — P/N/U Dataset Design and Evidence Contract
## Frozen methodological design v1.0

**Project:** PhageLysClass  
**Scope:** Detection and functional classification of phage lytic enzymes  
**Primary positive resource:** PhageLysData v1.0  
**External background resource:** UniProtKB, frozen release matched to the PhageLysData construction snapshot whenever possible  
**Primary purpose:** Define reproducible Positive (P), Negative (N), and Unlabeled (U) universes before any downstream representation, model-selection, hyperparameter-optimization, or performance analysis.

---

## 1. Purpose

PhageLysClass will evaluate how evidence definition, task formulation, protein representation, and sequence relatedness affect the reliable recognition of phage lytic enzyme function. The first methodological requirement is therefore not model selection, but a defensible definition of the biological universes from which prediction tasks will be constructed.

The central data-design principle is:

> **Absence of evidence is not equivalent to evidence of absence.**

Accordingly, proteins with explicit support for phage lytic function, proteins with explicit support for an alternative non-target function, and proteins lacking sufficient functional evidence must be treated as distinct evidence states.

The project therefore defines three primary modeling universes:

- **P — Positive:** proteins with non-predictive evidence supporting at least one target phage lytic function.
- **N — Negative:** phage proteins with explicit evidence supporting a biologically incompatible or clearly alternative non-target function.
- **U — Unlabeled:** phage proteins for which available evidence is insufficient to assign either P or N.

These definitions are frozen independently of downstream model performance.

---

## 2. Source basis and relationship to PhageLysData

PhageLysData v1.0 integrates 807,366 source observations into 759,105 exact-sequence entities and separates them into:

- **Core:** 11,867 evidence-supported entities;
- **Prediction Extension:** 745,092 prediction-only candidates;
- **Context:** 2,146 entities retained for provenance or biological context.

The Core contains four non-exclusive biological target classes:

1. `endolysin`
2. `virion_associated_lytic_enzyme`
3. `depolymerase`
4. `other_phage_envelope_lytic_enzyme`

The resource also preserves entities carrying more than one supported target assertion through the derived label `multiple_nonexclusive`.

PhageLysData uses exact amino-acid sequence identity to define entities. Sequence similarity, homology, domain architecture, predicted function, and clustering were not used to merge biological entities. This is critical for PhageLysClass because homology can subsequently be treated as an experimental variable rather than having been silently removed during resource construction.

The PhageLysData Core also contains evidence tiers distinguishing direct experimental support, high-confidence curated or independently supported evidence, annotation-based support, and ambiguous or incompletely resolved evidence. Enrichment layers such as InterProScan, Pfam, structural assets, physicochemical descriptors, and protein-language-model representations did not alter Core membership or biological assignments.

### 2.1 Consequence for PhageLysClass

PhageLysData Core will provide the positive biological universe.

The Prediction Extension will **not** be treated as negative evidence and will **not** be used as a primary unlabeled training set. Prediction-only membership indicates computational candidacy, not biological absence.

Context membership will likewise **not** be interpreted automatically as negative evidence.

External N and U universes will therefore be constructed independently from UniProtKB under a frozen query and evidence contract.

---

## 3. Non-negotiable methodological principles

### 3.1 Evidence-state separation

Every protein entering the classification framework must receive exactly one primary evidence state:

\[
P,\quad N,\quad U,\quad 	ext{or Excluded/Conflict}
\]

The assignment must be determined from biological evidence and frozen rules, never from model output.

### 3.2 Ground truth and prediction must remain separate

Prediction-derived annotations, PLM embeddings, InterProScan predictions, structural similarity, Foldseek hits, classifier scores, and other computational outputs cannot promote an entity from U to P or N during dataset construction.

They may later be used for:

- sensitivity analyses;
- orthogonal validation;
- interpretation;
- discovery demonstration.

### 3.3 Exact-sequence identity is resolved before homology

Exact duplicate sequences are consolidated before modeling.

Sequence similarity above or below any threshold is **not** an inclusion/exclusion criterion during P/N/U construction. Homology is preserved and subsequently measured explicitly during redundancy and generalization experiments.

### 3.4 Dataset definitions cannot be changed after model inspection

Once P/N/U are frozen, no downstream performance result may be used to:

- add or remove query families;
- change positive or negative terminology;
- redefine ambiguous proteins;
- exclude difficult protein families;
- change evidence thresholds solely because they reduce model performance.

Any post-freeze biological correction must be recorded as a versioned amendment with a documented rationale independent of predictive performance.

### 3.5 Modeling eligibility is distinct from biological membership

A protein can belong biologically to P, N, or U but fail a particular representation contract.

For example, the current PhageLysData representation contract requires canonical amino acids and sequence length ≤1,024 residues. Such restrictions define a **modeling subset**, not the biological universe itself.

---

# 4. Positive universe — P

## 4.1 Primary source

The positive universe is derived exclusively from the evidence-supported **PhageLysData Core v1.0**.

Current Core composition is:

| Resource-level class | Entities |
|---|---:|
| Endolysin | 9,888 |
| Depolymerase | 712 |
| Virion-associated lytic enzyme | 241 |
| Other phage-envelope lytic enzyme | 27 |
| Multiple non-exclusive | 999 |
| **Total Core** | **11,867** |

These counts describe the frozen resource and are not automatically the final modeling counts.

## 4.2 Positive biological definition

An entity belongs to P when PhageLysData contains at least one valid **non-predictive target assertion** supporting a target class.

The modeling layer should use the underlying target assertions rather than relying only on the derived canonical class.

For each positive entity, the following binary target vector should be reconstructed:

\[
Y_i=
[y_{\mathrm{endo}},
y_{\mathrm{VAL}},
y_{\mathrm{depol}},
y_{\mathrm{other}}]
\]

where each element indicates whether the entity has valid support for that target assertion.

This avoids treating `multiple_nonexclusive` as a biological class.

## 4.3 Positive subuniverses

### P_single

Entities with exactly one supported primary target among:

- endolysin;
- virion-associated lytic enzyme;
- depolymerase.

These proteins form the clean primary universe for conventional subtype classification.

### P_multi

Entities carrying support for two or more target classes.

The current resource-level `multiple_nonexclusive` group contains 999 entities, but membership for PhageLysClass must be reconstructed from assertion-level labels.

**Primary role:** challenge set for comparing independent binary and multiclass task formulations.

P_multi must initially remain outside the primary training set used to compare binary one-vs-rest and multiclass formulations.

### P_other

Entities supported as `other_phage_envelope_lytic_enzyme` without sufficient evidence for assignment to endolysin, VAL, or depolymerase.

Because the current class contains only 27 entities, it will not be treated as a primary trainable subtype.

Primary interpretation:

\[
	ext{lytic-positive, subtype unresolved}
\]

**Primary role:** challenge set for Stage-1 lytic detection and subtype uncertainty.

## 4.4 Evidence-tier handling

Evidence tier is retained as a first-class variable for every P entity.

At dataset-construction time:

- no Tier is automatically discarded solely because it is weaker;
- Tier 4 is not silently converted to U;
- Tier information does not become a feature for the sequence classifier;
- evidence-tier sensitivity will be evaluated after task construction.

The primary benchmark may later define a stricter evidence subset, but that decision must be based on the existing evidence contract and class viability, not predictive performance.

Recommended fields:

```text
sequence_entity_id
sequence
endolysin
virion_associated_lytic_enzyme
depolymerase
other_phage_envelope_lytic_enzyme
label_cardinality
evidence_tier
evidence_summary
source_count
source_provenance
modeling_eligible
```

---

# 5. External universe construction from UniProtKB

## 5.1 Temporal consistency

The initial N/U construction should preferentially use **UniProtKB release 2026_02**, matching the frozen UniProtKB release used by PhageLysData v1.0.

This prevents P from being defined using one annotation state while N/U are defined using a later state containing information that was unavailable during resource construction.

Later UniProtKB releases are reserved primarily for temporal/external discovery experiments.

## 5.2 Phage-protein eligibility

Every external candidate for N or U must first satisfy a phage-association criterion.

The criterion must be implemented through a frozen combination of UniProtKB taxonomic and/or viral-host metadata appropriate to bacteriophage proteins.

The precise query syntax will be stored in the query registry.

The biological requirement is:

> the model should distinguish target lytic proteins from other proteins in a biologically relevant phage-protein universe, not from arbitrary proteins from unrelated organisms.

---

# 6. Negative universe — N

## 6.1 Conceptual definition

N contains phage proteins for which available evidence supports an alternative function sufficiently strongly that the protein can reasonably be treated as non-target for the prediction task.

Formally:

\[
N = \{x:
E(	ext{alternative function}|x)	ext{ is sufficient}
\land
E(	ext{target lytic function}|x)=0
\}
\]

The defining requirement is **positive evidence for another function**, not absence of a lytic annotation.

## 6.2 N_easy — biologically clear negatives

N_easy contains proteins whose known functions are clearly distinct from the target lytic enzyme space.

Candidate query families include:

### Virion assembly and genome packaging
- major capsid proteins;
- scaffold proteins where clearly annotated;
- portal proteins;
- terminase large and small subunits;
- genome-packaging ATPases with unambiguous annotation.

### Genome replication and maintenance
- DNA polymerases;
- RNA polymerases where phage encoded;
- helicases;
- primases;
- ligases;
- nucleases with clearly non-envelope roles;
- replication initiation proteins.

### Integration and recombination
- integrases;
- excisionases;
- recombinases;
- resolvases with clear annotation.

### Transcriptional regulation
- transcription factors;
- repressors;
- antirepressors;
- sigma-factor-like proteins;
- transcriptional regulators with explicit function.

### Metabolic/accessory functions
Only functions that are clearly incompatible with the target lytic enzyme definitions and whose annotation is sufficiently specific.

### Intended use

N_easy provides a controlled evidence-backed negative universe and allows evaluation of the relatively broad question:

> Can phage lytic enzymes be distinguished from phage proteins with clearly unrelated functions?

---

# 7. Hard-negative universe — N_hard

## 7.1 Motivation

A classifier can obtain deceptively high performance if positives are compared only with biologically distant proteins.

N_hard therefore contains non-target proteins from functional or genomic neighborhoods closer to lytic systems.

Potential families include:

- holins;
- antiholins;
- spanins;
- selected non-enzymatic virion proteins;
- selected structural proteins involved in infection;
- selected cell-envelope-associated proteins with explicit non-target roles.

### Critical caution

Tail proteins, tail fibers, receptor-binding proteins, baseplate proteins, tailspikes, and related structural categories cannot be treated as automatically negative because:

- virion-associated lytic enzymes can be physically associated with the virion;
- depolymerases can occur in tail-associated and receptor-binding architectures.

Therefore, these families can enter N_hard only when the available evidence supports a clearly non-target function and all target-evidence exclusion rules are satisfied.

### Intended use

N_hard asks the more difficult question:

> Can the model recognize lytic enzyme function after trivial phage-versus-nonphage and gross structural-context signals are minimized?

N_hard is a negative subset, not an independent evidence state.

---

# 8. Unlabeled universe — U

## 8.1 Conceptual definition

U contains phage proteins for which available annotation is insufficient to support either target lytic function or a sufficiently specific alternative function.

Formally:

\[
U = \{x:
E(	ext{target}|x)	ext{ insufficient}
\land
E(	ext{alternative non-target}|x)	ext{ insufficient}
\}
\]

Typical candidates include:

- hypothetical proteins;
- uncharacterized proteins;
- proteins of unknown function;
- generic phage proteins;
- weak or non-specific functional annotations.

## 8.2 U is intentionally allowed to contain hidden positives

The U universe is not expected to be biologically pure.

Indeed, the possibility that U contains undiscovered lytic enzymes is central to the methodological question.

Therefore:

- U is not interpreted as negative ground truth;
- U cannot be used as a negative-labelled evaluation set;
- computational evidence discovered after U construction must not retroactively purge all plausible positives from U.

This intentional uncertainty permits a principled comparison between:

1. evidence-backed supervised learning;
2. conventional absence-as-negative learning;
3. positive-unlabeled learning.

## 8.3 What should be excluded from U

A UniProt record should not enter U if it already contains explicit target evidence identifiable by the frozen target-evidence registry.

Examples include records directly annotated with terminology or controlled annotations supporting:

- endolysin;
- virion-associated lytic enzyme;
- peptidoglycan-degrading phage enzyme;
- phage depolymerase;
- explicitly target-compatible catalytic or substrate functions according to the PhageLysData biological scope.

The exclusion registry should reuse and expand, where necessary, the same families of evidence used to build the PhageLysData UniProt query registry:

- explicit protein identities;
- functional annotations;
- catalytic annotations;
- Gene Ontology terms;
- Enzyme Commission terms;
- substrate terminology;
- structural terminology;
- phage-associated context.

However, the mere presence of a computationally predicted domain should not automatically convert an otherwise unlabeled protein into P.

---

# 9. Evidence assignment precedence

External proteins are assigned using the following precedence:

\[
oxed{P > Conflict/Exclude > N > U}
\]

### Rule 1 — Positive overlap

If an external sequence is an exact sequence match to a PhageLysData positive entity, it cannot be N or U.

It remains linked to P provenance or is removed from the external pool.

### Rule 2 — Conflicting explicit evidence

If a record contains both credible target evidence and credible alternative-function evidence that cannot be resolved under frozen rules, it enters:

\[
	ext{Excluded/Conflict}
\]

It must not be forced into N or U.

### Rule 3 — Negative evidence

If target evidence is absent and a qualifying explicit alternative function exists, assign N.

### Rule 4 — Unlabeled

If neither qualifying positive nor qualifying negative evidence is present, assign U.

---

# 10. Target-evidence exclusion registry

A frozen exclusion registry must be created before N/U materialization.

Each rule should contain:

```text
rule_id
evidence_source
field
term_or_identifier
target_class
match_type
rationale
action
version
```

Possible evidence sources include:

- UniProt protein name;
- function comments;
- catalytic-activity annotations;
- keywords;
- GO;
- EC;
- domain/family cross-references when they are part of the frozen annotation record;
- source-specific phage terminology.

The registry is intended primarily to prevent obvious target proteins from contaminating N.

It is not intended to prove that every protein lacking these signals is negative.

---

# 11. Negative-function query registry

A separate registry defines which alternative-function families qualify for N.

Recommended schema:

```text
query_id
negative_family
difficulty_level
uniprot_release
query_string
required_fields
inclusion_logic
target_exclusion_rules
biological_rationale
expected_risk
retrieval_date
status
```

`difficulty_level` should distinguish at least:

- `easy`
- `hard`

Every retrieved record remains traceable to the query or query family that admitted it.

Overlapping retrieval by several negative queries does not increase biological confidence automatically.

---

# 12. Exact-sequence consolidation

External UniProt records must be consolidated by normalized exact amino-acid sequence before modeling.

Multiple accessions mapping to an identical sequence become one external sequence entity while retaining:

- all accessions;
- all annotations;
- all query memberships;
- all relevant taxa/phage associations;
- source provenance.

This mirrors the central entity principle used by PhageLysData and prevents duplicated sequences from acting as independent observations.

---

# 13. Cross-universe decontamination

After sequence consolidation:

\[
P \cap N = arnothing
\]

\[
P \cap U = arnothing
\]

\[
N \cap U = arnothing
\]

under exact sequence identity.

All exact-sequence overlaps must be resolved before representation generation or model splitting.

## 13.1 Prediction Extension handling

An external N/U sequence matching a PhageLysData Prediction Extension entity should be flagged.

Recommended default:

- do not automatically promote it to P;
- do not automatically treat prediction-only membership as biological evidence;
- record the overlap as `prediction_extension_overlap = TRUE`.

For the strictest primary analysis, such sequences may be excluded from N to avoid using previously predicted lytic candidates as asserted negatives.

Whether they remain in U can be evaluated under a predeclared sensitivity analysis.

## 13.2 Context handling

Context overlaps should similarly be flagged rather than interpreted automatically.

---

# 14. Homology is preserved, not pre-cleaned

No sequence is removed from P/N/U merely because it is homologous to another sequence.

After P/N/U freeze, a global homology map will be constructed with MMseqs2 or an equivalent reproducible sequence-clustering procedure.

This map will support:

- redundancy analysis;
- cluster-disjoint splitting;
- maximum identity-to-training calculations;
- detection of label-conflicting clusters.

This separation is essential:

> **dataset definition determines biological evidence; homology analysis determines evaluation difficulty.**

---

# 15. Label-conflicting homology clusters

Homology clustering may reveal clusters containing:

- P and N entities;
- different P functional classes;
- single-label and multi-label positives;
- U mixed with known functional classes.

These clusters are scientifically informative and must not be silently corrected.

For the primary gold-standard evaluation:

- P/N conflicts may be quarantined from strict benchmark construction if the biological assertions cannot be reconciled;
- cross-class positive clusters remain valid unless they reveal an evidence error;
- all conflict decisions must be logged.

A dedicated `cluster_conflicts.tsv` should preserve them for later analysis and potential use in the patterns/interpretability paper.

---

# 16. Representation eligibility

The biological P/N/U universes should be constructed before applying representation-specific constraints.

For comparisons involving the existing PhageLysData representation contract, a common modeling subset should use identical eligibility rules across P/N/U.

The current PhageLysData contract includes:

- sequences containing only the 20 canonical amino acids;
- maximum sequence length of 1,024 residues.

Each entity should therefore receive:

```text
biological_member = TRUE/FALSE
representation_eligible = TRUE/FALSE
eligibility_reason
```

No biological label is changed because a sequence is representation-ineligible.

---

# 17. Class prevalence and sampling

P/N/U will be materialized at their naturally retrieved prevalence.

No permanent downsampling, oversampling, or balancing is allowed during universe construction.

Class balancing belongs to training only.

Primary evaluation sets should retain a frozen biologically meaningful prevalence or use explicitly defined balanced challenge sets alongside prevalence-preserving evaluation.

This avoids confusing a sampling decision with a biological evidence decision.

---

# 18. P/N/U-derived prediction tasks

The P/N/U design supports several deliberately different learning assumptions.

## 18.1 Evidence-backed supervised detection

\[
P \;	ext{vs}\; N
\]

This is the primary gold-standard binary detection task.

## 18.2 Easy-negative detection

\[
P \;	ext{vs}\; N_{easy}
\]

Measures broad discrimination.

## 18.3 Hard-negative detection

\[
P \;	ext{vs}\; N_{hard}
\]

Measures biologically difficult discrimination.

## 18.4 Absence-as-negative formulation

\[
P \;	ext{vs}\; U
\]

This intentionally treats unlabeled proteins as negatives during training and evaluates the consequences of that common assumption.

U is **not** used as negative ground truth during final evaluation.

## 18.5 Mixed conventional formulation

\[
P \;	ext{vs}\; (N+U)
\]

Optional analysis representing a conventional background-pool construction.

## 18.6 Positive–Unlabeled learning

\[
P + U
\]

U remains explicitly unlabeled.

All RQ1 strategies should ultimately be assessed against evidence-backed held-out P/N data.

---

# 19. P-derived subtype task design

Subtype tasks are independent of external N/U once Stage-1 lytic detection has been defined.

## 19.1 Independent binary subtype models

Train using P_single only:

\[
Endolysin \;vs\; other\ single\ lytic\ classes
\]

\[
VAL \;vs\; other\ single\ lytic\ classes
\]

\[
Depolymerase \;vs\; other\ single\ lytic\ classes
\]

P_multi remains unseen during training and is subsequently used as a biological challenge set.

Primary question:

> Can independent binary detectors recover multiple supported functions in proteins never shown as multi-functional during training?

## 19.2 Multiclass subtype model

Train on the same P_single universe:

\[
Endolysin / VAL / Depolymerase
\]

P_multi is then evaluated using:

- predictive entropy;
- top-1 probability;
- top-2 probability;
- top-1/top-2 margin;
- top-2 label concordance with supported assertions.

Primary question:

> Does non-exclusive biological evidence manifest as increased predictive ambiguity when a mutually exclusive task formulation is imposed?

## 19.3 P_other challenge

P_other is evaluated primarily as:

- Stage-1 lytic-positive recovery;
- Stage-2 unresolved or uncertain subtype behavior.

It is not used as a fourth primary subtype unless future data support a statistically and biologically justified class.

---

# 20. Dataset freezing and temporal discovery separation

The P/N/U universe belongs to the **development time point**.

A later UniProtKB release must not be used to modify P/N/U after freeze except through a documented dataset version update.

Later releases are reserved for a separate discovery universe:

\[
D = 	ext{new, previously unseen phage proteins}
\]

This makes temporal discovery possible without contaminating model development.

A particularly valuable retrospective temporal validation may arise if proteins that were unlabeled at the development snapshot later receive supported lytic annotations.

---

# 21. Quality-control gates

P/N/U construction cannot proceed to representation generation until all mandatory gates pass.

## Gate 1 — Source closure
- frozen UniProt release recorded;
- query registry versioned;
- retrieval date recorded;
- raw response files/checksums available.

## Gate 2 — Sequence integrity
- sequences normalized;
- exact duplicates consolidated;
- invalid/missing sequences explicitly flagged.

## Gate 3 — Evidence closure
- every retained entity assigned P, N, or U;
- every exclusion has a reason;
- every N entity linked to qualifying alternative-function evidence;
- every U entity lacks qualifying positive and negative evidence by the frozen rules.

## Gate 4 — Cross-universe exclusivity
- no exact sequence occurs simultaneously in P, N, and/or U;
- overlaps with PhageLysData Prediction Extension and Context are quantified.

## Gate 5 — Positive assertion closure
- P_single, P_multi, and P_other are materialized from assertion-level evidence;
- `multiple_nonexclusive` is not used as a biological target class.

## Gate 6 — Task feasibility
- counts available by P subtype;
- counts available for N_easy, N_hard, and U;
- sequence-length and representation-eligibility distributions evaluated;
- class viability documented before any split threshold is chosen.

## Gate 7 — No model-dependent decisions
A signed/frozen manifest confirms that P/N/U definitions were completed before predictive performance was inspected.

---

# 22. Required outputs

Recommended files:

```text
data_design/
├── p_positive_universe.tsv
├── p_single.tsv
├── p_multi_challenge.tsv
├── p_other_challenge.tsv
├── n_negative_universe.tsv
├── n_easy.tsv
├── n_hard.tsv
├── u_unlabeled_universe.tsv
├── excluded_conflicts.tsv
├── exact_sequence_overlap_report.tsv
├── prediction_extension_overlap.tsv
├── context_overlap.tsv
├── target_exclusion_registry.tsv
├── negative_query_registry.tsv
├── uniprot_retrieval_manifest.tsv
├── master_pnu_universe.tsv
├── pnu_summary.tsv
├── pnu_freeze_manifest.json
└── checksums.sha256
```

The master table should contain, at minimum:

```text
entity_id
sequence
sequence_hash
primary_state
state_subtype
functional_labels
label_cardinality
evidence_tier
negative_family
negative_difficulty
uniprot_accessions
phage_taxonomy
host_context
source_provenance
query_ids
phagelysdata_overlap
prediction_extension_overlap
context_overlap
representation_eligible
exclusion_status
exclusion_reason
```

---

# 23. Planned methodological analyses enabled by the frozen P/N/U design

The P/N/U architecture is designed specifically to support the following scientific questions without redefining datasets later.

### RQ1 — Negative evidence
Compare evidence-backed negatives, hard negatives, absence-as-negative formulations, and PU learning.

### RQ2 — Task formulation
Compare independent binary subtype models and multiclass models using P_single, then challenge them with P_multi and P_other.

### RQ3 — Representation/model strategies
Generate the same representation contract for external N/U and evaluate representation × learner combinations.

### RQ4 — Generalization
After P/N/U freeze, build homology clusters and quantify performance degradation as train–test similarity decreases.

### RQ5 — Demonstration
Apply the frozen model to a temporally/external discovery universe and use orthogonal annotation, structure, catalytic plausibility, and targeted docking where biologically justified.

---

# 24. Decisions explicitly deferred

The following decisions must not be prematurely embedded into P/N/U construction:

- final learner;
- final PLM;
- hyperparameter-search strategy;
- class balancing method;
- final homology thresholds;
- final PU algorithm;
- probability threshold;
- candidate-ranking function;
- structural validation pipeline;
- docking protocol.

These belong to later milestones.

---

# 25. Frozen methodological statements for future manuscript reuse

### Positive definition
“Positive entities were derived from the evidence-supported PhageLysData Core and required at least one valid non-predictive assertion supporting a target phage lytic function. Functional labels were reconstructed from assertion-level evidence rather than the derived resource-level canonical class.”

### Negative definition
“Negative proteins were retrieved independently from a frozen UniProtKB phage-protein universe and required explicit annotation supporting an alternative non-target function. Lack of a lytic annotation alone was not considered negative evidence.”

### Unlabeled definition
“Proteins lacking sufficient evidence for either target lytic function or a qualifying alternative function were retained as unlabeled rather than converted to negatives.”

### Homology principle
“Sequence homology was not used to define biological evidence states. Exact-sequence consolidation preceded modeling, whereas redundancy and train–test relatedness were quantified separately and treated as experimental variables.”

### Freeze principle
“All biological inclusion, exclusion, and evidence-assignment rules were frozen before model training and were not modified in response to downstream predictive performance.”

---

# 26. Immediate implementation plan

Before any model training, the project should complete the following sequence:

1. Extract assertion-level P from PhageLysData Core.
2. Materialize P_single, P_multi, and P_other.
3. Build the frozen target-evidence exclusion registry.
4. Design N_easy query families.
5. Design N_hard query families.
6. Design U retrieval rules.
7. Freeze UniProt query registry.
8. Retrieve UniProtKB 2026_02 candidates.
9. Consolidate external records by exact sequence.
10. Apply P/N/U precedence and conflict rules.
11. Audit overlaps with all PhageLysData universes.
12. Apply representation-eligibility annotations without changing biological membership.
13. Produce master P/N/U universe and QC summary.
14. Freeze checksums and dataset-design manifest.
15. Only then proceed to homology analysis, representations, task construction, and modeling.

---

## Final design statement

PhageLysClass will not define the negative class as “everything that is not known to be a phage lytic enzyme.” Instead, the project separates evidence-supported positives, evidence-supported alternative functions, and genuinely unresolved proteins into P, N, and U universes under frozen provenance-aware rules. This separation allows the effect of negative-class definition itself to become an experimental question and provides the foundation for subsequent binary, multiclass, positive-unlabeled, representation, and homology-aware analyses without allowing model performance to redefine the underlying biology.

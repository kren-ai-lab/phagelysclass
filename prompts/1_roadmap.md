Sí, papu. **Ahora sí estamos en el punto correcto para congelar el roadmap antes de escribir una sola línea seria de modelado.** La estructura del repo está perfecta para partir; no la inflaría todavía. Primero definiría el contrato científico y después cada módulo del repo aparece naturalmente.

Además, tenemos una base bastante privilegiada: PhageLysData ya separó evidencia de predicción, consolidó por identidad exacta —no por homología— y nos entrega 11,867 entidades Core, incluyendo 9,888 endolisinas, 712 depolimerasas, 241 VAL, 999 `multiple_nonexclusive` y 27 `other_phage_envelope_lytic_enzyme`; 11,259 Core tienen además el contrato completo de representaciones. 

# PhageLysClass — Full Roadmap v1.0

## Objetivo científico

Yo congelaría desde ya esta como **pregunta paraguas**:

> **How do evidence definition, task formulation, protein representation, and sequence relatedness determine reliable recognition of phage lytic enzyme function, and can these models identify previously unannotated but functionally plausible proteins?**

Y debajo quedan **cinco RQs**, cada una produciendo una pieza necesaria para la siguiente:

| RQ      | Pregunta                                                         | Qué queremos aprender                                                         |
| ------- | ---------------------------------------------------------------- | ----------------------------------------------------------------------------- |
| **RQ1** | What constitutes a negative?                                     | efecto de evidencia negativa vs ausencia de evidencia                         |
| **RQ2** | How should the functional task be formulated?                    | binarios independientes vs multiclass y comportamiento de casos no exclusivos |
| **RQ3** | Which representation/model strategy captures lytic function?     | interacción representación × learner                                          |
| **RQ4** | How does prediction generalize as sequence similarity decreases? | redundancia, leakage y extrapolación                                          |
| **RQ5** | Can the framework identify previously unannotated candidates?    | demostración prospectiva + validación in silico                               |

Y **no** incluimos como RQ: HPO, “qué algoritmo gana”, benchmark contra SUBLYME/DposFinder, motifs profundos ni ranking sofisticado.

---

# M0 — Scientific freeze

Antes de datasets.

Aquí congelamos formalmente:

**Scope biológico**

[
{\text{Endolysin},\text{VAL},\text{Depolymerase}}
]

`other_phage_envelope_lytic_enzyme` queda como **lytic-positive unresolved subtype**, no como cuarta clase entrenable.

`multiple_nonexclusive` deja de ser una clase artificial y vuelve a las assertions originales.

**Regla fundamental**

[
\text{absence of evidence}\neq\text{evidence of absence}
]

**Regla de evaluación**

Ningún threshold, representación, algoritmo o definición puede modificarse después de mirar el final test.

### Deliverable

`prompts/1_definitions.md`

Este debería convertirse realmente en el **frozen scientific contract** del proyecto.

---

# M1 — Construct the modeling universes

Esta es probablemente la etapa más crítica del paper.

Vamos a construir:

[
P,;N,;U,;C,;D
]

## `P` — Positive evidence

Fuente:

**PhageLysData Core v1.0**.

Pero recuperamos las **assertions**, no simplemente `canonical_class`.

Generaremos:

```text
sequence_entity_id
endolysin
virion_associated_lytic_enzyme
depolymerase
other_phage_envelope_lytic_enzyme
evidence_tier
```

### `P_single`

Proteínas con exactamente una de:

* Endolysin
* VAL
* Depolymerase

Este será el universo principal para RQ2.

### `P_multi`

Las ~999 `multiple_nonexclusive`.

**No entran al entrenamiento principal.**

Son un **biological challenge set**.

### `P_other`

Las 27 `other_phage_envelope_lytic_enzyme`.

Otro challenge set:

[
Lytic=yes,\quad subtype=unresolved
]

### Evidence sensitivity

Mantendría Tier como metadata.

No congelaría todavía “Tier 4 fuera del training” hasta ver exactamente las assertions y tamaños por clase; sí congelaría que **el efecto de evidence tier será auditable y nunca ignorado**.

---

# M2 — Build `N` and `U` from UniProt

Acá necesitamos construir el otro lado del problema de manera deliberada.

Usaría inicialmente **UniProtKB 2026_02**, porque es la release congelada en PhageLysData y evita inconsistencias temporales con el universo positivo. 

## `N_easy`

Phage proteins con **evidencia explícita de otra función**, por ejemplo:

capsid, portal, terminase, replication, DNA/RNA polymerase, helicase, primase, integrase, recombinase, packaging, transcriptional machinery, etc.

Esto responde:

> Can we separate lytic enzymes from clearly unrelated phage proteins?

## `N_hard`

Proteínas de fagos biológicamente próximas:

* holins;
* antiholins;
* spanins;
* selected structural proteins;
* selected tail/baseplate proteins;
* other envelope-associated proteins,

**pero solamente cuando exista evidencia suficiente para afirmar función alternativa**.

Aquí tendremos una exclusión muy estricta de cosas como `tailspike`, `tail fiber`, `RBP`, etc. si existe riesgo de mezclar depolimerasas o VAL.

## `U`

Phage proteins cuya función no permite clasificarlas ni como P ni como N:

* hypothetical;
* uncharacterized;
* unknown function;
* weak/generic annotations.

Y aquí **no hacemos InterPro/Pfam filtering agresivo posterior**, porque queremos conservar la posibilidad de que `U` contenga positivos desconocidos.

### Necesitamos un query registry congelado

Cada consulta deberá tener:

```text
query_id
uniprot_query
intended_universe
biological_rationale
inclusion_rule
exclusion_rule
release
retrieval_date
```

Esto será material Supplementary precioso.

---

# M3 — Contamination and provenance control

Antes de modelar:

[
P\cap N=\varnothing
]

[
P\cap U=\varnothing
]

por exact sequence.

También excluiremos secuencias externas que aparezcan en cualquier universo relevante de PhageLysData cuando eso pueda introducir circularidad.

Pero **no eliminaremos homólogos simplemente por ser homólogos**.

La homología no se “limpia” silenciosamente.

Se mide y después se convierte en variable experimental.

### Output

`master_modeling_universe.tsv`

con al menos:

```text
sequence_id
sequence
origin
evidence_state
functional_labels
uniprot_accession
source
evidence_tier
eligibility
```

---

# M4 — RQ1: What constitutes a negative?

Ahora sí construimos tres formulaciones principales.

### E1 — Evidence-backed supervised

[
P;vs;N
]

### E2 — Absence-as-negative

[
P;vs;U
]

Representa la práctica común:

> no está anotada como lítica → negativo.

### E3 — Positive–Unlabeled

[
P+U
]

tratando explícitamente U como unlabeled.

Y probablemente un análisis complementario:

[
P;vs;(N+U)
]

### Evaluación justa

El gold test debe contener **P y N con evidencia**.

No evaluamos contra U fingiendo conocer su verdad.

### Resultado esperado

No buscamos demostrar que PU siempre gana.

Buscamos medir:

[
\Delta Performance =
f(\text{negative definition})
]

y posteriormente:

[
\Delta Performance =
f(\text{negative definition},\text{sequence distance})
]

Esto puede ser uno de los resultados principales del paper.

---

# M5 — RQ2: How should functional classification be formulated?

Esta etapa utiliza únicamente `P_single` para training.

## T1 — Independent binary models

Tres modelos:

[
Endo;vs;non-Endo
]

[
VAL;vs;non-VAL
]

[
DPO;vs;non-DPO
]

donde los negativos vienen de las **otras clases líticas single-label**.

### Después:

Pasamos `P_multi`, jamás visto como multi durante entrenamiento.

Pregunta:

> **Can independent binary detectors recover multiple supported functions?**

Evaluaremos:

* number of detected labels;
* label recall;
* exact-set agreement;
* partial-set agreement.

Este experimento me parece **muy potente**.

---

## T2 — Multiclass

[
Endo/VAL/DPO
]

utilizando exactamente el mismo `P_single`.

Después aplicamos las ~999 `P_multi`.

Pero ahora el modelo está obligado a:

[
\sum_kP(y_k)=1.
]

Por tanto analizamos:

* entropy;
* top-1 probability;
* top-2 probabilities;
* top-1/top-2 margin;
* class ambiguity.

### Hipótesis

Los `P_multi` deberían mostrar más:

[
predictive\ uncertainty
]

que proteínas single-label equivalentes.

Esto convierte la comparación binary vs multiclass en una **pregunta biológica**, no en una competición de métricas.

---

## T3 — Explicit multilabel

Lo dejaría **secundario/condicional**.

Podemos entrenarlo con las assertion matrices completas, pero no quiero perder la belleza de utilizar `P_multi` como challenge completamente externo a T1/T2.

Si aporta algo nuevo → entra.

Si solamente replica los binarios → Supplementary o fuera.

---

# M6 — Representation contract

Aquí tenemos una ventaja enorme.

PhageLysData ya distribuye 11 PLMs sobre el mismo conjunto y contrato de pooling/precision. 

Para los nuevos N/U generamos exactamente el mismo contrato.

Yo **sí exploraría las 11 representaciones inicialmente**, porque ya es barato conceptualmente y nos permite aprovechar el resource.

Pero después agrupamos/congelamos finalistas.

Además:

### Baseline A

Physicochemical/sequence descriptors.

### Baseline B

Sequence similarity / nearest homolog.

Eso es obligatorio conceptualmente.

Queremos distinguir:

[
\text{functional learning}
]

de

[
\text{nearest homolog retrieval}.
]

---

# M7 — Modeling strategies

Definimos formalmente:

[
Strategy =
Representation\times Algorithm
]

No “modelo”.

## Algorithms

Yo congelaría cuatro familias:

* **Logistic Regression**
* **SVM**
* **RF or XGBoost**
* **MLP**

Después de preflight decidimos RF vs XGB, o conservamos ambos si computationally cheap.

### Screening

Primero hyperparameters moderados/fijos.

### HPO

Solo después de identificar configuraciones razonables.

Mismo budget para todos los finalistas.

Y HPO **no es RQ**.

---

# M8 — Global homology map

Aquí hacemos MMseqs2 sobre el universo completo relevante.

Generamos al menos:

[
C90
]

[
C70
]

[
C50
]

[
C30.
]

También almacenamos pairwise/local nearest similarity necesaria para:

[
I_{max}(x,Train).
]

Importante:

### Redundancy y split son dos experimentos distintos.

---

# M9 — RQ4A: Effect of redundancy

Partimos del universo exact-sequence unique.

Construimos distintas versiones del training set:

[
NR100
]

[
NR90
]

[
NR70
]

[
NR50.
]

No prometo NR30 porque VAL puede colapsar.

Preguntamos:

> **How does internal dataset redundancy affect apparent model performance and strategy ranking?**

Aquí será interesante ver si cambia incluso:

[
best\ representation
]

o:

[
best\ algorithm.
]

---

# M10 — RQ4B: Generalization

Este sí es uno de los corazones del paper.

Cada **strategy** tendrá:

[
Random
]

[
H90
]

[
H70
]

[
H50
]

[
H30
]

cuando el tamaño de clase lo permita.

### Regla congelada

No eliminamos un threshold porque el resultado quede feo.

Solo se elimina si existe un **predefined viability failure**, por ejemplo insuficientes positivos en una clase.

---

## Generalization curves

Para cada strategy:

[
MCC_S(h)
]

[
AUPRC_S(h)
]

en función de:

[
h.
]

Pero además tendremos la versión que me parece todavía más importante:

[
I_{max}(x)=
\max_{t\in Train}Identity(x,t)
]

y estimaremos:

[
Performance_S=f(I_{max})
]

en bins, posiblemente:

[

> 90
> ]

[
70-90
]

[
50-70
]

[
30-50
]

[
<30.
]

### Esto debe hacerse para **cada strategy relevante**.

Exactamente como dijiste.

No queremos simplemente:

> “ESM2+SVM ganó”.

Podemos encontrar:

> ESM2+SVM obtiene el mayor random score pero ProtT5+LR pierde menos performance bajo remote homology.

Eso es científicamente mucho más interesante.

---

# M11 — Locked model selection

Solo después de RQ1–RQ4.

No necesariamente existe un único winner.

Podemos tener:

### Best conventional classifier

máxima performance global.

### Best generalizer

mejor comportamiento bajo baja identidad.

### Best calibrated / stable strategy

si observamos diferencias relevantes.

Pero para la **demo final** congelamos una estrategia o un ensemble simple antes de tocar `D`.

---

# M12 — Discovery universe

Acá hacemos la demo, no Paper 4.

Preferencia:

[
\boxed{\text{UniProt release posterior a 2026_02}}
]

Tomamos proteínas de fagos:

* nuevas respecto a PhageLysData;
* no presentes por exact sequence;
* nunca utilizadas en P/N/U;
* sin anotación lítica explícita.

Eso genera:

[
D.
]

Idealmente también buscamos una pequeña validación temporal:

> proteins unannotated at time (t_0) that received a lytic annotation at (t_1).

Si existen, es casi una mini **time-split prospective validation**.

Pero no la hacemos requisito.

---

# M13 — Discovery demonstration

Aplicamos el classifier congelado.

No hacemos todavía:

[
Novelty+Uncertainty+Patterns+Developability
]

como función de ranking formal.

Eso es Paper 4.

Aquí usamos reglas simples y declaradas antes:

* positive lytic prediction;
* reproducibility across repeats/strategies;
* subtype signal;
* sequence novelty;
* manageable length/quality.

Elegimos algo como:

[
6-12\ candidates.
]

Idealmente con representación de Endo/VAL/DPO.

---

# M14 — Orthogonal in-silico validation

Esta debe responder:

> **Does evidence not used by the classifier support the predicted function?**

## Layer 1 — Functional annotation

Después de seleccionar:

* InterProScan;
* Pfam;
* HMMER.

## Layer 2 — Structure

* existing AFDB/PDB when possible;
* prediction otherwise;
* Foldseek structural neighbors.

PhageLysData ya tiene InterProScan para 97.3% del Core y estructuras existentes para 745 entidades, así que tenemos un excelente **positive reference universe**. 

## Layer 3 — Catalytic plausibility

Cuando corresponda:

* expected catalytic fold;
* catalytic residues;
* spatial arrangement;
* pocket.

## Layer 4 — Docking

Solo si existe una hipótesis biológica defendible.

### Endo / VAL

fragmentos de peptidoglicano compatibles con la catalytic family.

### DPO

oligosacárido específico únicamente cuando host/serotype/substrate lo justifiquen.

**No docking genérico decorativo.**

Y siempre:

[
Candidate
]

versus

[
KnownPositive
]

idealmente más:

[
NegativeControl.
]

---

# M15 — Final robustness / sensitivity

Antes del freeze final:

* evidence tier sensitivity;
* class balance sensitivity;
* cluster threshold sensitivity;
* repeated seeds;
* calibration;
* failure analysis;
* `P_multi`;
* `P_other`;
* conflict clusters.

Todo lo secundario puede ir al SI.

---

# M16 — Results freeze and journal gate

Recién ahí decidimos journal.

### Scenario A

Resultados correctos pero relativamente esperables:

→ **Protein Science**.

### Scenario B

Vemos efectos conceptuales fuertes:

* definición de negativo cambia sustancialmente las conclusiones;
* binary models recuperan overlaps reales;
* multiclass uncertainty refleja `multiple_nonexclusive`;
* estrategias intercambian ranking conforme baja homología;
* remote candidates reciben fuerte soporte estructural/catalítico;

→ **paramos y reevaluamos target hacia arriba**.

No hay razón para decidir hoy el techo.

---

# Métricas congeladas

### Binary detection / OVR

Primary:

[
MCC,\quad AUPRC
]

Secondary:

F1, sensitivity, specificity, balanced accuracy, AUROC.

### Multiclass

Primary:

[
Macro\text{-}F1,\quad MCC,\quad Macro\text{-}AUPRC
]

más per-class metrics.

### Multi-nonexclusive challenge

* true-label recovery;
* predicted label cardinality;
* exact/partial agreement;
* entropy;
* top-2 agreement;
* probability margin.

### Generalization

No solamente score absoluto.

También:

[
\Delta M_{Random\rightarrow Hx}
]

y pendiente/degradación de cada strategy.

---

# Qué queda explícitamente fuera del Paper 2

Esto es importante congelarlo ya:

| Tema                                    | Destino |
| --------------------------------------- | ------- |
| systematic motif discovery              | Paper 3 |
| physicochemical design rules            | Paper 3 |
| structural pattern atlas                | Paper 3 |
| deep XAI / residue pattern discovery    | Paper 3 |
| formal discovery ranking                | Paper 4 |
| sophisticated uncertainty-aware ranking | Paper 4 |
| mining masivo                           | Paper 4 |
| Prediction Extension prioritization     | Paper 4 |
| de novo generation                      | Paper 5 |
| experimental validation                 | Paper 5 |

Paper 2 únicamente **demuestra que discovery es posible**.

---

# Cómo lo mapearía a tu repo

Tu estructura actual está bien. No metería veinte directorios nuevos.

```text
phagelysclass/
├── LICENSE
├── README.md
├── prompts/
│   ├── 1_definitions.md
│   └── 2_full_roadmap.md
│
└── src/
    ├── dataset_processing/
    │   ├── 01_build_positive_universe.py
    │   ├── 02_query_uniprot_negative_unlabeled.py
    │   ├── 03_build_master_universe.py
    │   ├── 04_build_homology_clusters.py
    │   └── 05_build_discovery_universe.py
    │
    ├── task_definitions/
    │   ├── 01_detection_tasks.py
    │   ├── 02_binary_subtype_tasks.py
    │   ├── 03_multiclass_tasks.py
    │   └── 04_challenge_sets.py
    │
    ├── representation_strategies/
    │   ├── 01_generate_descriptors.py
    │   ├── 02_generate_embeddings.py
    │   └── 03_validate_representation_contract.py
    │
    └── training_models/
        ├── 01_screen_strategies.py
        ├── 02_optimize_finalists.py
        ├── 03_evaluate_generalization.py
        ├── 04_evaluate_challenge_sets.py
        └── 05_fit_frozen_models.py
```

Posteriormente podemos agregar `src/discovery_demo/`, pero **no lo añadiría hasta llegar realmente a M12**.

---

# Orden de trabajo real desde hoy

Y para mí este es el orden que debemos respetar:

1. **Freeze RQs and definitions.**
2. **Design P/N/U biologically.**
3. **Freeze UniProt queries.**
4. **Materialize P/N/U.**
5. **Audit overlaps and class sizes.**
6. **Freeze T1/T2/T3 tasks.**
7. **Build homology map.**
8. **Only then start representations/models.**
9. **Screen strategies.**
10. **Run generalization landscape for every relevant strategy.**
11. **Freeze final classifier.**
12. **Touch discovery universe for the first time.**
13. **Perform orthogonal structural/functional validation.**
14. **Freeze paper results.**

Y yo pondría una regla de oro en `1_definitions.md`:

> **No downstream modeling result is allowed to redefine the biological composition of P, N, U, challenge, or discovery universes.**

Eso nos protege de adaptar la definición del problema para conseguir mejores scores.

Con esto, papu, **PhageLysClass ya tiene una hoja de ruta completa desde evidencia hasta demostración**. El siguiente paso lógico no es programar embeddings ni modelos: es bajar **M1–M3 a especificaciones exactas de datasets y queries UniProt**, porque ahí se juega buena parte de la calidad científica del trabajo.

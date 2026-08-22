Sí, papu. **Ahora sí congelaría el diseño científico v1.0.** Y con tu observación de binarios vs multiclass queda bastante más interesante que lo que veníamos armando.

El Paper 2 no será simplemente *“classify phage lytic enzymes”*. La pregunta grande queda:

> **How do evidence definition, task formulation, protein representation and sequence relatedness determine our ability to recognize phage lytic enzyme function and extrapolate toward previously unannotated proteins?**

PhageLysData nos da una base especialmente buena porque el Core ya separa soporte no predictivo de la Prediction Extension, conserva las afirmaciones funcionales originales y sus evidencias, consolida únicamente por identidad exacta y distribuye 11 representaciones PLM sobre un mismo contrato. Además, el propio Data Report deja explícito que su clasificación actual es solo ilustrativa, con split estratificado random, sin control de homología ni HPO. 

# Paper 2 — Frozen scientific design v1.0

## Las cinco preguntas centrales

| RQ      | Pregunta                                                           | Datos principales                             | Experimento central                                      |
| ------- | ------------------------------------------------------------------ | --------------------------------------------- | -------------------------------------------------------- |
| **RQ1** | **What constitutes a negative?**                                   | PhageLysData `P` + UniProt `N/U`              | evidence-negative vs absence-as-negative vs PU           |
| **RQ2** | **How should function be formulated?**                             | Core single-class + 999 multiple-nonexclusive | binary one-vs-rest vs multiclass + ambiguity challenge   |
| **RQ3** | **What representations/modeling strategies encode function?**      | PLMs + descriptors                            | representation × algorithm                               |
| **RQ4** | **How much performance survives sequence divergence?**             | P/N/U + MMseqs2 clusters                      | redundancy + random/H90/H70/H50/H30                      |
| **RQ5** | **Can the framework recognize previously unannotated candidates?** | later/external UniProt universe               | prospective prediction + orthogonal in-silico validation |

Estas cinco cuentan **una sola historia**:

[
Evidence
\rightarrow
Task
\rightarrow
Representation
\rightarrow
Generalization
\rightarrow
Application
]

---

# RQ1 — What constitutes a negative?

Esta queda congelada como uno de los dos experimentos metodológicos principales.

## `P` — positive evidence

Fuente:

[
\boxed{\text{PhageLysData Core}}
]

Pero no utilizaría todo indiscriminadamente.

### Development positives

Primariamente:

[
P_{strict}=Tier,1+Tier,2+Tier,3
]

porque Tier 4 representa la evidencia ambigua/incompletamente resuelta.

Eso nos deja, **antes de aplicar eligibility de representación**:

* 588 depolymerases;
* 9,875 endolysins;
* 220 VAL;
* 27 other;

aproximadamente, derivados de la tabla Tier × class del recurso. La cifra final para modelado se recalculará tras aplicar el mismo contrato de secuencia a todos los universos. 

### Challenge positives

[
P_{challenge}=Tier,4
]

incluyendo las **999 `multiple_nonexclusive`**.

Esto es mucho mejor que meter automáticamente Tier 4 al training.

---

## `N` — explicit negative evidence

Fuente:

[
\boxed{\text{UniProtKB}}
]

Idealmente la misma release **2026_02** utilizada por PhageLysData para que la definición inicial sea temporalmente consistente. 

Pero `N` significa:

> **evidence of another phage-protein function**

y no simplemente `NOT lysin`.

### `N_easy`

Funciones inequívocamente diferentes:

capsid, portal, terminase, polymerase, helicase, primase, integrase, recombinase, packaging proteins, etc.

### `N_hard`

Proteínas biológicamente más próximas:

holins, antiholins, spanins y determinadas proteínas estructurales/de interacción de fago **solo cuando exista evidencia explícita de una función no enzimática**.

No pondría alegremente `tailspike` dentro de negativos porque podemos meter depolimerasas reales.

---

## `U` — unlabeled / absence of evidence

También UniProtKB:

[
U=
\text{phage proteins without sufficient evidence for P or N}
]

Principalmente:

* hypothetical;
* uncharacterized;
* weakly annotated;
* proteins without sufficient functional assignment.

Y aquí está la gracia:

[
U\neq N.
]

---

## Las tres estrategias de RQ1

Todas se evalúan contra **el mismo gold test compuesto por P + N**, nunca tratamos `U_test` como ground truth negativo.

### Strategy E1 — evidence-supervised

[
P;vs;N
]

### Strategy E2 — absence-as-negative

[
P;vs;U
]

imitando la práctica común de considerar no anotado = negativo.

### Strategy E3 — PU

[
P+U
]

con Positive–Unlabeled learning.

Podemos tener también:

[
P;vs;(N+U)
]

como análisis convencional complementario.

### Pregunta que responderemos

> **Does the way negative evidence is defined materially alter functional recognition and generalization?**

Esto me parece muy fuerte porque es extrapolable a muchos problemas de function prediction.

---

# RQ2 — Binary or multiclass?

Acá tu observación mejora mucho el diseño.

## No escogería uno antes de probar.

**Compararía las dos formulaciones porque representan hipótesis biológicas diferentes.**

Y no utilizaría las 999 `multiple_nonexclusive` para entrenar ninguna de las dos inicialmente.

Esas 999 se convierten en un **challenge set excepcional**.

---

## Strategy T1 — Independent binary models

Una vez que Stage 1 dice:

[
Lytic=yes
]

entrenamos tres detectores independientes:

[
f_E(x)=P(Endolysin)
]

[
f_V(x)=P(VAL)
]

[
f_D(x)=P(Depolymerase)
]

Training:

* single-class Endo;
* single-class VAL;
* single-class DPO.

Negativos de cada modelo:

> las otras clases líticas limpias.

Los 999 `multiple_nonexclusive` **no participan en entrenamiento ni selección**.

### Después hacemos algo precioso

Los pasamos por los tres binarios.

Una proteína puede producir:

[
Endo=1,\quad VAL=1,\quad DPO=0.
]

Y comparamos eso contra las **assertions originales de PhageLysData**.

Entonces preguntamos:

> **Can independent functional detectors recover biologically non-exclusive assignments that were never shown during training?**

Eso me encanta.

---

# Strategy T2 — Multiclass

Entrenamos:

[
Endo;/;VAL;/;DPO
]

solo sobre proteínas single-class.

Otra vez:

**999 multiple_nonexclusive fuera del training.**

Después las aplicamos al modelo multiclass.

Pero como softmax obliga:

[
\sum P_k=1,
]

nos interesa estudiar:

* predictive entropy;
* top-1/top-2 probabilities;
* probability margin;
* top-2 label concordance;
* calibration.

### Hipótesis

Las proteínas `multiple_nonexclusive` deberían mostrar:

[
H(P(y|x))_{multiple}

>

H(P(y|x))_{single}.
]

Es decir:

> **la ambigüedad biológica debería manifestarse como incertidumbre predictiva cuando obligamos al modelo a elegir una única clase.**

Esta pregunta es bastante sexy.

---

# Y eso nos permite comparar binary vs multiclass de una manera NO trivial

No preguntamos solamente:

> cuál da mayor F1.

Preguntamos:

### Binary

¿Puede reconocer **coexistencia funcional**?

### Multiclass

¿La biología no exclusiva emerge como **uncertainty**?

Eso es muchísimo más interesante.

---

# ¿Multilabel explícito?

### Lo dejaría como análisis secundario, no como eje principal.

Después podemos entrenar:

[
[Endo,VAL,DPO]
]

usando las assertion matrices completas.

Pero eso le revela durante training la existencia de múltiples labels y pierde el experimento precioso de las 999 secuencias nunca vistas como overlap.

Así que:

**primary comparison**

[
Binary\ OVR \quad vs\quad Multiclass
]

**secondary**

[
Multilabel
]

si aporta.

---

# Las 27 `other_phage_envelope_lytic_enzyme`

No haría una clase con n=27.

Las mantendría como un segundo **out-of-definition challenge set**.

Preguntamos:

> ¿Stage 1 las reconoce como líticas?

y:

> ¿Stage 2 muestra baja confianza / múltiples señales / ninguna asignación clara?

Eso es biológicamente mucho más razonable.

---

# RQ3 — Representation × algorithm

Acá sí definiría explícitamente **strategy** como:

[
S_{ij}=
Representation_i+Algorithm_j.
]

Y desde este momento todos los análisis posteriores operan sobre **estrategias**, no sobre “el mejor modelo”.

Eso responde exactamente a lo que dices.

---

## Representaciones

PhageLysData ya tiene, sobre 11,259 Core elegibles:

* 11 PLMs;
* one-hot;
* physicochemical properties. 

Yo haría screening de **todas las representaciones PLM disponibles** porque ya tenemos el asset.

Eso explota realmente el valor del resource paper.

Añadiría:

### Baseline

physicochemical / composition descriptors.

### Homology baseline

MMseqs2 nearest-neighbor label transfer.

One-hot podemos incorporarlo si encontramos un learner razonable sin introducir una arquitectura completamente distinta; no lo haría obligatorio.

---

# Algoritmos

Congelaría cuatro familias:

[
\boxed{LR}
]

linear.

[
\boxed{SVM}
]

kernel/maximum-margin.

[
\boxed{RF;or;XGBoost}
]

tree ensemble.

[
\boxed{MLP}
]

nonlinear embedding learner.

Yo probablemente usaría **XGBoost**, pero esa es una decisión técnica que podemos cerrar después de preflight.

No necesitamos diez algoritmos.

---

# HPO

Sí.

Pero:

> **HPO es parte del protocolo, no una research question.**

Pipeline:

[
train
\rightarrow
internal\ validation
\rightarrow
HPO
\rightarrow
freeze
\rightarrow
test.
]

Mismo budget de optimización por strategy.

Nada toca el test.

Y guardamos default-vs-optimized como Supplementary porque puede revelar cuánto importa realmente HPO.

---

# RQ4 — Redundancy and generalization

Acá haría dos experimentos explícitamente separados.

## A. Redundancy

PhageLysData usa **exact-sequence identity** para formar entidades; no usa homología para fusionarlas. Eso es ideal para empezar este experimento. 

Sobre:

[
P+N+U
]

generamos clusters MMseqs2:

[
C90,;C70,;C50,;C30.
]

Y construimos versiones:

[
NR100,;NR90,;NR70,;NR50.
]

No necesariamente NR30 porque podemos destruir VAL.

Esto responde:

> **How does training-set redundancy influence estimated performance?**

---

# B. Generalization

Y esto es distinto.

Para **cada strategy representación + algoritmo** tendremos:

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

donde `Hxx` significa cluster-disjoint train/test al threshold correspondiente.

### Esto queda congelado:

> **Every modeling strategy gets its own generalization curve.**

Por ejemplo:

[
MCC_S(h)
]

o

[
AUPRC_S(h)
]

con:

[
h={Random,90,70,50,30}.
]

Eso nos deja comparar no solamente quién gana:

> ¿qué estrategia **degrada más lentamente** cuando la homología desaparece?

Ahí está la pregunta interesante.

---

# Y además calculamos la versión continua

Para cada test protein:

[
I_{max}(x)=
\max_{t\in Train}Identity(x,t).
]

Eso nos permite obtener:

[
Performance_S=f(I_{max})
]

para cada strategy.

Los bins pueden ser:

* > 90%
* 70–90%
* 50–70%
* 30–50%
* <30%.

### Esta debería ser una de las figuras principales.

Porque podemos encontrar casos como:

**Strategy A**

[
0.95\rightarrow0.65
]

**Strategy B**

[
0.92\rightarrow0.81.
]

Entonces A es el mejor classifier convencional, pero B es el mejor **functional generalizer**.

Eso es mucho más interesante para protein science.

---

# Split viability gate

No voy a congelar que obligatoriamente habrá H30 si VAL no lo permite.

Congelo el procedimiento:

> usamos H90/H70/H50/H30 siempre que exista un mínimo estadísticamente útil de positivos por clase en validation/test.

Si un threshold colapsa una clase:

* se reporta;
* no se fuerza;
* el último threshold viable se convierte en el límite de evaluación.

Incluso eso es información sobre el sequence space.

---

# Algo importante: cluster conflicts

MMseqs2 puede producir clusters donde aparezcan:

[
P+N
]

o diferentes functional labels.

Esos casos **no se resuelven automáticamente**.

Los marcamos:

[
C_{conflict}
]

y los:

* excluimos del gold benchmark primario;
* analizamos posteriormente como challenge/error set.

Podrían ser biológicamente muy interesantes para Paper 3.

---

# Métricas congeladas

## Stage 1 binary detection

Primarias:

[
\boxed{MCC}
]

[
\boxed{AUPRC}
]

Secundarias:

* F1;
* balanced accuracy;
* sensitivity;
* specificity;
* AUROC.

---

## Stage 2 binary OVR

Por clase:

* MCC;
* AUPRC;
* F1;
* sensitivity.

Y para multiple-nonexclusive:

### `multi-hit recall`

> cuántos labels verdaderos son recuperados.

### exact-set / partial-set agreement

cuando las assertions lo permitan.

---

## Stage 2 multiclass

Primarias:

* macro-F1;
* multiclass MCC;
* macro-AUPRC.

Además:

* per-class recall;
* confusion matrix.

Para las 999:

* entropy;
* probability margin;
* top-2 agreement;
* calibration/uncertainty distributions.

---

# Todas las comparaciones llevan repetición

No seed 42 y listo.

Congelaría:

[
5\ independent\ repeats
]

como mínimo para screening.

Y:

[
10
]

para finalistas cuando sea computacionalmente razonable.

Reportamos:

[
95%,CI
]

y efectos pareados entre estrategias.

---

# RQ5 — Discovery demonstration

Esta no será todavía la metodología de mining.

Será simplemente:

> **Can the final classifier nominate plausible unseen candidates?**

---

## `D` — discovery space

No usaría Prediction Extension como demo primaria.

El recurso contiene **745,092 prediction-only candidates**, por lo que encontrar candidatos ahí no significaría necesariamente que nosotros los descubrimos; otra herramienta ya los propuso. 

Preferiría:

[
\boxed{\text{later UniProtKB release}}
]

posterior a 2026_02.

Definimos:

[
D=
\text{new phage proteins}
]

que:

* no aparecen por exact sequence en PhageLysData v1.0;
* no participaron en P/N/U;
* no tienen anotación lítica explícita;
* pasan QC.

---

# Selección de candidatos

Sin inventar todavía el ranking de Paper 4.

Utilizamos criterios simples y predefinidos:

1. predicted as lytic;
2. subtype prediction reproducible;
3. agreement entre varias strategies fuertes;
4. low similarity to PhageLysData positives.

Idealmente:

[
<30% identity
]

con cobertura razonable.

Seleccionamos quizá:

[
2-3
]

por clase.

Total:

[
6-9
]

proteínas.

Perfecto para una demo.

---

# Validación in silico

Acá congelaría una jerarquía de evidencia.

## V1 — domain evidence

Después de la predicción:

* InterProScan;
* Pfam/HMMER.

Es información que **no participó en el predictor**, por lo que puede ser evidencia ortogonal.

---

## V2 — structural plausibility

Obtenemos estructura disponible o la predecimos.

Después:

[
Foldseek
]

contra estructuras funcionalmente relevantes.

Buscamos:

* fold compatible;
* domain architecture;
* catalytic-domain geometry.

---

## V3 — catalytic plausibility

Cuando conocemos la familia:

* catalytic residues;
* pocket;
* spatial arrangement.

---

## V4 — ligand/substrate interaction

### Endolysin / VAL

Docking con fragmentos de peptidoglicano compatibles con el mecanismo esperado.

### Depolymerase

Docking **solo si conocemos un polisacárido/host/serotype razonablemente justificable**.

Si no:

> structural/domain evidence is the validation.

No docking decorativo.

---

# Controles de validación

Para cada candidato:

[
Candidate
]

versus

[
Known\ Core\ positive
]

versus, cuando sea razonable,

[
N_{hard}.
]

Esto hace que docking/estructura tengan una referencia y no sean simplemente un score aislado.

---

# Qué NO hacemos todavía

Esto queda congelado para que los papers no se pisen:

| Tema                                              | Paper |
| ------------------------------------------------- | ----- |
| deep motif / pattern discovery                    | **3** |
| mechanistic feature interpretation                | **3** |
| structural landscape systematic analysis          | **3** |
| formal candidate scoring/ranking                  | **4** |
| uncertainty-aware mining                          | **4** |
| mining masivo de Prediction Extension/metagenomes | **4** |
| de novo generation/design                         | **5** |
| experimental validation                           | **5** |

Paper 2 **abre la puerta** a discovery. No la cruza completa.

---

# Las figuras ya salen casi solas

No congelaría todavía el número final, pero veo **4 figuras principales + 1 tabla**.

### **Figure 1 — Evidence and task formulation**

P/N/U + Stage 1 + binary vs multiclass + challenge sets.

### **Figure 2 — Representation/model strategies**

representation × algorithm × task.

Aquí podemos meter heatmaps y efectos de negative definition.

### **Figure 3 — Generalization landscape**

La gran figura:

[
Performance
\ vs
homology
]

para **cada strategy**.

Más binary/multiclass challenge de las 999 si visualmente cabe.

### **Figure 4 — Discovery and orthogonal validation**

unseen proteins → predictions → structural/domain support → docking/interaction.

### **Table 1**

Datasets, evidence states, task definitions, sizes y experimental roles.

Y el resto:

Supplementary.

---

# Workflow congelado

La secuencia de trabajo debería ser:

**M0 — Freeze biological questions**

✅ ahora.

**M1 — Construct `P/N/U`**

PhageLysData + frozen UniProt query registry.

**M2 — Freeze task definitions**

Stage 1 + binary/multiclass Stage 2 + challenge sets.

**M3 — Build global homology map**

P/N/U + MMseqs2 thresholds.

**M4 — Generate representation contract**

PhageLysData assets + same embeddings for N/U.

**M5 — Representation × algorithm screening**

Todas las strategies.

**M6 — HPO and final strategy freeze**

Sin tocar final test.

**M7 — Generalization experiments**

Random/H90/H70/H50/H30 + continuous identity curves para **cada strategy**.

**M8 — Multiple/nonexclusive challenge**

Binary multi-hit vs multiclass uncertainty.

**M9 — External/temporal discovery**

UniProt post-2026_02.

**M10 — Orthogonal in-silico validation**

domains → structure → catalytic plausibility → docking where justified.

**M11 — Results freeze + journal decision**

Ahí miramos lo que tenemos.

---

## Y no congelaría todavía el journal

**Protein Science es el mínimo target que estamos diseñando.**

Pero si obtenemos, por ejemplo:

* efecto fuerte y generalizable de P/N/U;
* diferencias claras binary vs multiclass en las 999 secuencias ambiguas;
* algunas strategies que mantengan rendimiento sorprendentemente bien bajo H30;
* candidatos realmente remotos;
* soporte estructural/catalítico convincente;

entonces **paramos antes de enviar y reevaluamos target**.

Porque ahí puede haber un paper bastante más fuerte que simplemente un Methods/Application para Protein Science.

### Para mí, esta es la versión que debemos congelar.

Especialmente me gustan tres resultados potenciales que le dan personalidad propia:

[
\boxed{\text{negative evidence matters}}
]

[
\boxed{\text{non-exclusive biology emerges differently in binary vs multiclass models}}
]

[
\boxed{\text{model quality is better defined by degradation with homology than by one random-split score}}
]

Si conseguimos demostrar esas tres cosas y luego cerrar con una pequeña demo de discovery + validación ortogonal, **Paper 2 queda muy redondo** y nos deja exactamente donde queremos estar para los Papers 3 y 4.

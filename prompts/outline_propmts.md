Sí, papu. Para **Paper 2 → Protein Science**, yo congelaría desde ya una filosofía: **el objeto científico es la función proteica; ML es el instrumento**. La revista acepta estudios computacionales de función proteica, pero para un Full-Length Article quiere hallazgos que avancen la comprensión de las proteínas, no solo un predictor técnicamente correcto. El formato además favorece un paper compacto: ~5,000 palabras y ≤10 elementos principales. ([Wiley Online Library][1])

Y protegería explícitamente los papers siguientes: **Paper 2 reconoce**, Paper 3 explica patrones/determinantes, Paper 4 mina/rankea, Paper 5 diseña.

# Tesis central del Paper 2

No:

> *We developed a machine-learning classifier for phage lytic enzymes.*

Sí:

> **Reliable recognition of phage lytic enzyme function depends strongly on how functional evidence, negative examples, sequence redundancy, representation, and evolutionary generalization are handled.**

La pregunta general sería:

> **How reliably can phage lytic enzyme function be recognized from protein sequence, particularly beyond the sequence space represented during training?**

Eso tiene mucho más sabor a *Protein Science*.

---

# Las 4 preguntas científicas que deberían gobernar todo

## RQ1 — ¿Cómo debe definirse el problema biológico?

No podemos asumir que:

[
\text{not annotated as lytic}=\text{non-lytic}.
]

Queremos determinar cómo cambia el reconocimiento cuando distinguimos:

* evidencia positiva de función lítica;
* evidencia negativa real cuando exista;
* proteínas de fagos con otra función conocida;
* proteínas simplemente **unlabeled**.

Aquí entra PU learning, pero como una **pregunta sobre la definición biológica de la tarea**, no como gimmick de ML.

### Hipótesis

> La definición de negativos tendrá un impacto comparable o superior al cambio de algoritmo.

Esta me parece una de las hipótesis más fuertes del paper.

---

# RQ2 — ¿Qué representación captura mejor función bajo divergencia?

Comparar bajo **exactamente las mismas particiones**:

* sequence/homology baseline;
* simple sequence composition/descriptors;
* physicochemical descriptors;
* PLM embeddings.

No necesitamos 15 PLMs.

Yo haría **3 familias suficientemente diferentes**, por ejemplo:

[
\text{ESM2}
\qquad
\text{ProtT5}
\qquad
\text{ESM-C}
]

más los baselines.

Podemos incluso incluir dos tamaños de ESM2 si queremos preguntar por model scale, pero no lo convertiría en zoológico de embeddings.

### Hipótesis

> Protein language representations should improve recognition particularly as sequence similarity to the training set decreases.

Esta frase vende muchísimo mejor el paper que:

> ESM2 had the highest F1.

---

# RQ3 — ¿Cuánto de la performance depende de redundancia y homología?

Este debería ser probablemente **el experimento central del artículo**.

Separaría dos cosas que la literatura frecuentemente mezcla:

### Redundancy

[
D_{100},D_{90},D_{70},D_{50}
]

### Partition strategy

[
Random
\quad vs\quad
Cluster/Homology-aware
]

No son lo mismo.

Y evaluaría el deterioro de performance cuando nos alejamos del training.

### Hipótesis

> Conventional random partitions overestimate recognition performance because homologous proteins occupy both training and test sets.

Pero iría todavía un paso más allá:

[
Performance=f(\max Identity_{train})
]

Eso nos permite responder la pregunta protein-centered:

> **¿hasta qué distancia evolutiva sigue siendo reconocible la función?**

Ahí veo una figura central del paper.

---

# RQ4 — ¿Podemos reconocer correctamente las distintas funciones líticas?

Aquí entra finalmente la arquitectura jerárquica.

Yo no usaría:

[
phage/nonphage
\rightarrow
class
]

porque nuestra aplicación real será trabajar dentro de espacios de proteínas de fagos.

Haría:

[
\boxed{\text{Lytic/envelope-degrading vs other phage protein}}
]

y después:

[
\boxed{
\text{Endolysin}
\qquad
\text{VAL}
\qquad
\text{Depolymerase}
}
]

con formulación **multilabel cuando la evidencia lo requiera**.

Y tendría como baseline:

[
\text{flat multiclass}
]

para poder demostrar si la jerarquía realmente aporta algo.

No asumamos que lo hará.

Eso es importante.

---

# La arquitectura experimental completa

Yo la diseñaría en etapas, **no como un factorial gigante**.

## Stage A — Frozen classification universe

Partir exclusivamente del universo congelado de PhageLysData.

Importantísimo:

### Ground truth ≠ prediction

No debería entrar como truth ninguna etiqueta producida por nuestros propios modelos o por la prediction extension.

Mantendría:

* sequence;
* class/evidence;
* provenance;
* evidence tier;
* exclusivity/non-exclusivity;
* soporte ortogonal cuando exista.

Y generaría un **classification manifest congelado**.

Cada proteína tendría algo como:

```text
entity_id
sequence
functional_labels
evidence_level
evidence_sources
label_origin
redundancy_cluster_90
redundancy_cluster_70
redundancy_cluster_50
redundancy_cluster_30
```

Este manifest debería ser completamente reproducible.

---

# Stage B — Definir positivos, negativos y unlabeled

Aquí haría tres escenarios.

### B1. Conventional supervised

Positivos:

> supported PhageLysData lytic enzymes.

Negativos:

> phage proteins con función conocida diferente de la función objetivo.

### B2. Broad-background supervised

Positivos iguales.

Background:

> otras proteínas de fagos sin evidencia de función lítica.

Aquí deliberadamente asumimos unlabeled = negative para medir qué ocurre.

### B3. PU-learning

[
P=\text{known lytic}
]

[
U=\text{phage proteins without positive lytic evidence}
]

No necesitamos probar veinte métodos PU.

Uno o dos enfoques robustos bastan.

### Resultado importante

No solamente:

> PU performed better.

Sino:

[
\Delta Performance =
f(\text{negative definition})
]

---

# Un detalle que considero fundamental: hard negatives

Construiría un subconjunto de negativos difíciles:

* otras proteínas del mismo proteoma de fago;
* proteínas estructurales;
* tail-associated proteins;
* otras enzimas;
* proteínas próximas funcionalmente cuando sea posible.

Porque distinguir una endolisina de una proteína humana random sería absurdo.

Nuestro problema real es:

> **dada una proteína de fago, ¿es una de las proteínas funcionales que buscamos?**

Esto hace al paper muchísimo más creíble.

---

# Stage C — Comparar formulaciones de la tarea

Usaría las mismas secuencias para probar:

### C1 — Flat multiclass

[
Endo/VAL/Depol/Other
]

como baseline comparable con approaches existentes.

### C2 — Hierarchical

[
Lytic
\rightarrow
Endo/VAL/Depol
]

### C3 — Hierarchical multilabel

[
Lytic
\rightarrow
[y_{endo},y_{VAL},y_{depol}]
]

permitiendo:

[
\sum y_k > 1.
]

Esto aprovecha directamente una característica que ya existe en el recurso: **la evidencia funcional no siempre es mutuamente exclusiva**.

Pero no vendería la jerarquía como novedad. SUBLYME ya tiene una jerarquía.

La pregunta es:

> **Which formulation best reflects the biological evidence and generalizes most reliably?**

---

# Stage D — Representation benchmark

Primero congelamos un split estándar y después tocamos representaciones.

### Baseline 1 — Sequence similarity

Obligatorio.

Por ejemplo:

[
\text{nearest homolog}
\rightarrow
\text{label transfer}.
]

MMseqs2/BLAST.

Esto es importantísimo para Protein Science.

Queremos saber:

> ¿cuándo ML aporta información por encima de simple homología?

### Baseline 2 — sequence composition

AAC / dipeptide o algo equivalente.

### Baseline 3 — physicochemical descriptors

Aprovechar la capa que ya tenemos en PhageLysData.

### PLMs

Yo elegiría máximo **3–4 representaciones**.

No más.

---

# Stage E — Model benchmark

Tampoco haría zoológico.

## Linear probe

**Logistic regression**

Es científicamente útil porque mide cuánto está ya organizado el espacio de representación.

## Kernel

**SVM**

Además nos conecta directamente con SUBLYME.

## Tree ensemble

Yo elegiría:

**XGBoost o RF**, no necesariamente ambos.

## Nonlinear embedding model

**small MLP**

Entonces tenemos:

[
linear
/
kernel
/
tree
/
neural
]

Perfectamente suficiente.

---

# Hyperparameter optimization

Aquí reduciría muchísimo la importancia que inicialmente le estábamos dando.

Sí debemos hacerlo.

Pero no debe convertirse en un capítulo.

## Diseño

Primero hacemos screening.

Después seleccionamos:

* top representations;
* top learners.

Y recién ahí:

[
Nested\ CV + Optuna
]

dentro de training.

Nunca tocar test.

Incluso podemos preguntar en Supplementary:

[
\Delta_{HPO}
============

Metric_{optimized}-Metric_{default}
]

para mostrar cuánto realmente aporta optimizar.

Si resulta que la ganancia es 0.01, maravilloso:

> representation/data decisions mattered more than extensive tuning.

Eso es un hallazgo mucho más interesante.

---

# Stage F — Redundancy experiment

Aquí sí haría un experimento explícito.

Por ejemplo:

[
NR100
\rightarrow
NR90
\rightarrow
NR70
\rightarrow
NR50
]

manteniendo metodología idéntica.

Medimos:

* número de secuencias;
* clases retenidas;
* prevalence;
* performance;
* variability.

Eso nos muestra:

[
Performance=f(redundancy)
]

Pero **sin confundirlo con generalización**.

---

# Stage G — Split experiment

Compararía como mínimo:

### Random stratified

El escenario convencional.

### Homology-aware

clusters completos asignados a:

[
train/validation/test.
]

Y me gustaría tener por lo menos dos niveles:

[
H50
]

y

[
H30.
]

Si la cobertura queda suficiente.

No necesariamente necesitamos H90/H70/H50/H30 todos como experimentos completos.

Podemos generar todos los clusters y mostrar la curva, pero reservar:

**Random / H50 / H30**

como escenarios principales.

Eso mantiene el paper manejable.

---

# Stage H — El experimento estrella: generalization landscape

Esto sí lo quiero.

Para cada secuencia de test calculamos:

[
I(x)=
\max_{t\in Train}
Identity(x,t).
]

Agrupamos:

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
<30%.
]

Y mostramos:

[
MCC
]

[
AUPRC
]

[
Recall
]

versus máximo similarity-to-training.

### Esa figura responde:

> **How far from known sequence space can phage lytic enzyme function still be recognized?**

Para mí esta podría perfectamente ser **Figure 4 o Figure 5 y la figura más citada del trabajo**.

---

# Stage I — Independent evaluation

Esto podría subir bastante nuestras posibilidades en Protein Science.

No confiaría exclusivamente en CV.

Construiría un **final locked independent test set** que no participa en:

* representation selection;
* model selection;
* HPO;
* threshold optimization.

Idealmente de una de estas formas:

### Opción A — source-disjoint

Retener completamente determinadas fuentes/publicaciones.

### Opción B — evidence-disjoint

Proteínas experimentalmente/literature-supported mantenidas completamente fuera del training.

### Opción C — temporal holdout

Secuencias incorporadas después de la congelación del training resource.

La opción **B** me parece especialmente bonita si tenemos suficientes.

Porque tendríamos:

[
Train =
\text{resource evidence}
]

versus

[
Test =
\text{orthogonally supported proteins}.
]

Eso evita circularidad.

---

# Importantísimo: evaluar “known-domain” vs “domain-remote”

Esto le da mucho Protein Science sin comernos Paper 3.

No necesito descubrir patrones.

Simplemente estratificamos test:

### Canonical-domain supported

proteínas con dominios líticos reconocibles mediante Pfam/InterPro/HMM.

### No canonical-domain support

proteínas cuya evidencia funcional existe pero donde esos métodos no reconocen claramente el dominio esperado.

Entonces preguntamos:

[
Performance_{known-domain}
\quad vs\quad
Performance_{domain-remote}.
]

Si el PLM mantiene capacidad en el segundo grupo, tenemos un resultado precioso:

> **learned sequence representations recover functional signal beyond readily detectable canonical domain similarity.**

Si cae brutalmente, también es científicamente valioso:

> estos predictores siguen dependiendo mucho más del espacio conocido de lo que sugieren sus scores globales.

Win-win.

---

# Stage J — Calibration, pero sin gastar Paper 4

Haría **calibración básica**, porque estamos hablando de predicción fiable.

Por ejemplo:

* reliability curves;
* Brier score;
* ECE;
* temperature/isotonic calibration cuando corresponda.

Pero **NO haría todavía**:

* conformal discovery;
* epistemic/aleatoric decomposition profunda;
* candidate ranking;
* novelty score;
* uncertainty-driven mining.

Eso pertenece al Paper 4.

Acá solamente demostramos:

> predicted probabilities can/cannot be interpreted reliably.

---

# Métricas

No usaría accuracy como headline.

## Binary detection

Primarias:

[
MCC
]

[
AUPRC
]

Secondary:

* balanced accuracy;
* sensitivity;
* specificity;
* F1;
* AUROC.

Especialmente AUPRC porque el problema real será desbalanceado.

## Functional multilabel

Primarias:

* macro-AUPRC;
* macro-F1;
* per-class AUPRC.

Secondary:

* micro-F1;
* precision/recall por clase;
* Hamming loss.

### Para hierarchy

Reportaría dos niveles:

1. conditional functional classification entre verdaderos lytic enzymes;
2. **end-to-end performance**, contando errores producidos en Stage 1.

Porque si Stage 1 pierde una proteína, Stage 2 jamás puede rescatarla.

---

# Variabilidad y estadística

Nada de:

> seed 42, listo.

Como mínimo:

[
n=5
]

repeticiones independientes por configuración importante.

Idealmente 10 para las finales.

Reportaría:

[
mean/median + 95% CI.
]

Y para comparaciones pareadas usaría bootstrap/permutation o modelos estadísticos sencillos sobre los splits.

Más importante:

> **effect size > p-value.**

Ejemplo:

[
\Delta MCC_{random\rightarrow H30}=-0.23
]

es mucho más informativo que `p < 0.001`.

---

# Benchmark con herramientas existentes

Para Protein Science, **ya no lo considero opcional**.

Pero debe ser justo.

## Lysin / Endo / VAL

Comparar cuando la tarea coincida con:

**SUBLYME/PhaLP 2.0.**

## Depolymerase

Cuando corresponda:

* DposFinder;
* DePP;
* DepoScope/PhageDPO si es ejecutable de manera razonable.

Pero no haría una tabla grotesca diciendo:

> ours 0.94 > theirs 0.89.

Las definiciones de tareas cambian.

El benchmark serio debe usar:

> **the same independent evaluation set whenever technically possible.**

Y donde no sea posible:

> explicitly report the mismatch.

Eso da mucha credibilidad.

---

# Qué NO pondría en Paper 2

Esto es igual de importante.

### ❌ Pattern discovery profundo

Motifs, posiciones, arquitecturas y reglas de diseño:

**Paper 3.**

### ❌ Structure-function atlas

Paper 3.

### ❌ Ranking de cientos de candidatos

Paper 4.

### ❌ Novelty-aware scoring

Paper 4.

### ❌ Conformal/ensemble uncertainty discovery

Paper 4.

### ❌ Large-scale environmental mining

Paper 4.

### ❌ Generative modeling

Paper 5.

### ❌ Experimental validation de de novo sequences

Paper 5.

---

# Sí permitiría un pequeño error analysis

Porque necesitamos demostrar que entendemos qué está haciendo el clasificador.

Por ejemplo:

* performance por familia;
* performance por identity-to-train;
* performance con/sin canonical domain support;
* performance por evidence tier;
* confusion endolysin ↔ VAL;
* comportamiento de `multiple_nonexclusive`.

Eso no invade Paper 3.

Simplemente caracteriza **dónde funciona y dónde falla el reconocimiento**.

---

# Las figuras principales

Yo intentaría **6 figuras**, no más.

## Figure 1 — From PhageLysData evidence to prediction tasks

Una figura conceptual preciosa:

[
PhageLysData
\rightarrow
evidence filtering
\rightarrow
P/U/N
\rightarrow
Stage,1
\rightarrow
Stage,2.
]

Mostrar clases, overlaps y tamaños.

---

## Figure 2 — Task definition shapes apparent performance

Comparar:

* strict negatives;
* broad background;
* PU;
* flat vs hierarchical/multilabel.

Mensaje:

> **What we call a negative materially changes the problem.**

---

## Figure 3 — Representation and learner benchmark

Heatmap:

[
representation\times model
]

para:

* MCC;
* macro-AUPRC.

Nada más.

---

## Figure 4 — Redundancy and partitioning

Probablemente:

[
NR100/90/70/50
]

vs.

[
Random/H50/H30.
]

Mensaje:

> conventional evaluation can substantially overestimate functional recognition.

---

## Figure 5 — Recognition across sequence divergence

Mi favorita.

X:

[
maximum\ identity\ to\ training.
]

Y:

[
performance.
]

Curvas por representación.

Podemos agregar known-domain / domain-remote.

**Esta es la figura Protein Science.**

---

## Figure 6 — Independent test and biological limits

Final locked test.

Mostrar:

* binary detection;
* Endo;
* VAL;
* DPO;
* specialized baselines;
* calibration;
* failure categories.

Cerrar con un modelo conceptual de:

[
known-like
\rightarrow
remote-but-recognizable
\rightarrow
unresolved.
]

Sin hacer todavía mining.

---

# Supplementary fuerte

Aquí mandaría el volumen:

* todas las hyperparameter grids;
* todos los PLMs;
* todos los algorithms secundarios;
* threshold analyses;
* full confusion matrices;
* per-seed metrics;
* per-evidence-tier results;
* detailed class prevalence;
* cluster statistics;
* alternative redundancy definitions;
* complete classification reports;
* source-specific performance;
* calibration curves adicionales.

El main paper debe ser **muy limpio**.

---

# Outline del manuscrito

## **Introduction**

### P1 — Biological problem

Diversidad y relevancia de:

* endolysins;
* VALs;
* depolymerases.

### P2 — Why prediction remains difficult

* divergent sequences;
* modular evolution;
* inconsistent annotations;
* imbalance;
* unknown proteins.

### P3 — Limitations of current predictors

No atacar papers individuales.

El mensaje:

> reported predictive performance is difficult to compare because task definitions, negative sets, redundancy and homology control differ substantially.

### P4 — Opportunity

PhageLysData permite controlar esas variables en un universo evidence-aware.

### P5 — Objective

> determine how task formulation, representation, redundancy and sequence divergence affect reliable recognition of phage lytic enzyme function.

---

# Results

### **1. Evidence-aware formulation of phage lytic enzyme recognition tasks**

### **2. Functional task and negative-set definition strongly influence apparent predictive performance**

### **3. Protein representations differ in their ability to encode lytic enzyme function**

### **4. Sequence redundancy inflates conventional estimates of predictive performance**

### **5. Homology-aware evaluation reveals the limits of functional generalization**

### **6. Functional recognition persists—or deteriorates—across divergent and domain-remote protein space**

### **7. Independent evaluation establishes a reference framework for reliable phage lytic enzyme prediction**

Dependiendo de resultados podemos fusionar 6–7.

---

# Discussion

Yo la estructuraría alrededor de cuatro mensajes, no alrededor de algoritmos.

### 1. Functional prediction starts with task definition

“No evidence” no significa “negative”.

### 2. Representation quality matters most outside close homology

Lo importante no es que ESM gane 0.02 en random split.

Es **dónde empieza a separarse de homology transfer**.

### 3. Random performance is a poor proxy for discovery capability

Este probablemente será un mensaje fuerte.

### 4. Reliable classification has explicit boundaries

Definir dónde el modelo deja de ser confiable prepara perfectamente el Paper 4.

Y ahí podemos terminar diciendo:

> these boundaries provide the foundation for subsequent uncertainty- and novelty-aware candidate prioritization.

Sin hacer todavía ese paper.

---

# El “minimum publishable strong paper”

Si queremos mantenerlo rápido, para mí hay **7 cosas no negociables**:

1. **Evidence-clean ground truth** sin predicciones circulares.
2. **Biologically meaningful phage-protein negatives/unlabeled data.**
3. **Flat vs hierarchical/multilabel task formulation.**
4. **Homology baseline + descriptors + 2–3 PLMs.**
5. **Random vs homology-aware evaluation.**
6. **Performance as a function of sequence distance to training.**
7. **One truly independent locked test set.**

Con eso, **yo sí lo veo competitivo para Protein Science**.

PU learning, full HPO analysis, muchas representaciones, calibration exhaustiva, etc. pueden mejorar el trabajo, pero no son el núcleo.

---

## Y creo que hay una formulación todavía más potente para venderlo

El paper no debería concluir:

> “Our best model obtained an MCC of 0.XX.”

Debería concluir algo parecido a:

> **Phage lytic enzyme function is highly recognizable within known sequence space, but apparent performance depends strongly on evidence definition and sequence relatedness. Homology-aware evaluation identifies the conditions under which learned protein representations retain functional signal beyond close sequence similarity.**

Eso **sí es un resultado sobre proteínas**.

Y eso encaja muy bien con lo que Protein Science dice querer: avances en comprensión de **protein function**, incluyendo estudios computacionales, y no necesariamente una nueva arquitectura de ML. ([Wiley Online Library][1])

Para mí, papu, **ese es el Paper 2 que debemos congelar antes de escribir código nuevo**. Y la mejor parte es que no necesita convertirse en un monstruo: la dificultad está en diseñar bien los conjuntos y la evaluación, no en entrenar modelos enormes.

[1]: https://onlinelibrary.wiley.com/page/journal/1469896x/homepage/forauthors.html?utm_source=chatgpt.com "Protein Science"

--- 

Sí, papu. **Con el texto exacto de Protein Science delante, yo también me movería a *Methods and Applications***, pero eso cambia la pregunta del paper.

No debería ser “hagamos todos los análisis posibles con PhageLysData”. Tiene que existir **un método concreto, reutilizable y suficientemente general**, demostrado con phage lytic enzymes. Y coincido contigo: **3 figuras multipanel + 1 tabla** me parece mucho mejor que 6 figuras.

Además, hay una limitación que debemos reconocer desde el principio: **PhageLysData por sí solo no necesariamente contiene el universo negativo necesario para entrenar un detector “lytic vs non-lytic”**. Tenemos un recurso extraordinario de positivos/evidencia funcional, pero para detección probablemente necesitaremos construir *un único background externo bien definido* de proteínas de fagos no líticas/unlabeled. No veinte datasets nuevos.

## Yo reformularía el paper así

El producto metodológico sería algo del estilo:

> **An evidence- and homology-aware framework for hierarchical protein function classification under incomplete functional annotation**

y **phage lytic enzymes serían la aplicación principal**.

Eso le da la “broad applicability” que pide *Methods and Applications*: la metodología no depende de que la proteína sea una endolisina. En principio sirve para cualquier problema en que tengamos:

* positivos de diferente calidad;
* negativos imperfectos o ausencia de verdaderos negativos;
* clases funcionales desbalanceadas;
* funciones potencialmente no mutuamente exclusivas;
* redundancia/homología que puede inflar la evaluación.

Ahí ya tenemos algo más general que “otro predictor de lysins”.

---

# Primero congelaría solo **tres preguntas**

Nada más por ahora.

## **Q1. ¿Cómo formular correctamente el problema de clasificación?**

Esta sería la pregunta metodológica principal.

Tenemos que comparar una formulación convencional con una formulación que respete la evidencia del recurso.

Por ejemplo:

### Formulación convencional

[
\text{Endolysin / VAL / Depolymerase / Other}
]

single-label.

versus

### Formulación propuesta

[
\text{Lytic enzyme detection}
]

seguido por:

[
\text{Endolysin / VAL / Depolymerase}
]

permitiendo multilabel cuando corresponda.

Pero aquí aparece el primer **data requirement real**:

### Para Stage 1 necesitamos un background externo.

PhageLysData puede proporcionarnos muy bien:

[
P=\text{supported lytic enzymes}
]

pero necesitamos obtener:

[
N/U=\text{other phage proteins}
]

Yo haría **una sola adquisición externa congelada**, probablemente desde UniProtKB y/o INPHARED, con reglas transparentes.

No buscaría decenas de recursos.

Idealmente:

* proteínas de fagos con anotaciones claramente no líticas → **N**;
* proteínas de fagos sin evidencia suficiente para decidir → **U**.

Entonces podemos incluso comparar:

[
P\text{ vs }N
]

con

[
P\text{ vs }U
]

o implementar una estrategia PU sencilla.

Pero **PU sería una extensión de esta pregunta**, no un paper dentro del paper.

### Lo que queremos aprender

No:

> PU beats SVM by 0.03.

Sino:

> **Treating lack of annotation as evidence of absence changes the inferred protein-function boundary.**

Eso es metodológicamente interesante y generalizable.

---

# **Q2. ¿Qué información de secuencia necesitamos para reconocer la función?**

Aquí sí utilizaría muy bien lo que ya tenemos.

PhageLysData ya contiene distintas capas que permiten comparar familias de representación sin salir a fabricar demasiadas cosas adicionales:

### Sequence-derived descriptors

Nuestros descriptores fisicoquímicos.

### PLM representations

Un conjunto **muy pequeño y representativo**, no todo Sylphy.

Algo como:

* ESM2;
* ProtT5;
* ESM-C.

Incluso podría reducirlo a dos PLMs si los resultados preliminares muestran redundancia.

Y usar exactamente **el mismo clasificador downstream** para comparar representaciones.

Esto es importante.

No haría:

[
8\ PLMs \times 8\ classifiers \times 200\ HPOs.
]

Porque eso convierte el paper en benchmarking.

Haría algo tipo:

[
Representation
\rightarrow
\text{fixed/simple learner}
]

y después, una vez elegida la representación, optimizamos solamente el modelo final.

### Pregunta científica

> **Does learned protein representation improve functional discrimination relative to conventional sequence-derived properties, particularly for divergent proteins?**

Eso conecta directamente con protein science.

---

# **Q3. ¿Cuánto depende nuestra estimación de performance de la similitud entre proteínas?**

Esta para mí sigue siendo **obligatoria** y probablemente la pregunta más potente.

Pero también la simplificaría.

No necesitamos:

NR100 × NR90 × NR70 × NR50 × random × H90 × H70 × H50 × H30 × 10 modelos.

Eso sería infernal.

Haría algo mucho más limpio.

### Dataset redundancy

Un dataset principal con un criterio razonable de reducción de redundancia.

Y quizá comparar:

[
NR90
]

versus

[
NR50
]

como análisis de sensibilidad.

### Evaluation strategy

Dos condiciones principales:

[
\boxed{Random/stratified}
]

versus

[
\boxed{Homology-aware}
]

Por ejemplo cluster-disjoint a 30–50%, según lo que soporte el número de muestras por clase.

Eso basta para responder:

> **How much does conventional random partitioning overestimate functional generalization?**

Y adicionalmente podemos calcular para cada proteína del test:

[
I_{\max}(x)=
\max_{t\in Train}Identity(x,t)
]

para observar cómo cambia la predicción con distancia al training.

Eso nos da muchísimo sin multiplicar experimentos.

---

# Con esas tres preguntas ya tenemos el paper

Y fíjate que son metodológicamente coherentes:

[
\boxed{
\text{What is the task?}
}
]

[
\downarrow
]

[
\boxed{
\text{How should proteins be represented?}
}
]

[
\downarrow
]

[
\boxed{
\text{How should generalization be evaluated?}
}
]

Eso puede transformarse en un **framework reproducible**.

No necesitamos responder todavía:

* qué patrones determinan actividad;
* cuáles son los candidatos nuevos;
* cómo rankearlos;
* qué dominios son novedosos;
* cómo diseñar nuevas enzimas.

Todo eso queda limpio para #3–#5.

---

# ¿Qué necesitamos obtener externamente?

Aquí sería muy conservador.

## Necesario

### **1. Background de proteínas de fagos**

Este es el único bloque importante que hoy probablemente falta.

Idealmente tres grupos:

**P — Positive**

PhageLysData supported core.

**N — Negative**

Proteínas de fagos con evidencia funcional incompatible con lytic enzyme.

Por ejemplo proteínas inequívocamente:

* capsid;
* terminase;
* DNA polymerase;
* helicase;
* transcriptional regulators;
* etc.

**U — Unlabeled**

Proteínas de fagos sin anotación lítica pero tampoco evidencia negativa suficiente.

Eso nos permitiría estudiar la definición de tarea sin inventarnos negativos.

### Fuentes posibles

Podemos resolverlo probablemente con:

* UniProtKB phage proteins;
* INPHARED;
* quizá PHROG annotation si hace falta.

Pero escogería **una fuente primaria**, no las tres por defecto.

---

## No necesario para este paper

No iría a buscar:

* estructuras adicionales;
* AlphaFold masivo;
* metagenomas;
* nuevos dominios;
* experimental annotations manuales a gran escala;
* taxonomía compleja;
* host ranges;
* kinetic information.

Eso pertenece a otros trabajos.

---

# Tampoco haría benchmark de herramientas

Estoy contigo.

SUBLYME responde:

> lysin vs non-lysin → Endo vs VAL.

DposFinder responde otra tarea.

DePP responde DPO.

DeepMineLys tiene otra construcción de dataset.

Comparar:

[
F1_{\text{ours}} > F1_{\text{theirs}}
]

sería metodológicamente bastante dudoso porque **los ground truths no son iguales**.

Y para un Methods paper me parece incluso más elegante hacer **ablations internas controladas**:

* conventional labels vs evidence-aware;
* random vs homology-aware;
* descriptors vs PLM;
* flat vs hierarchical.

Ahí cada comparación tiene exactamente el mismo universo de datos.

Eso vale más que una tabla con seis servidores.

Podemos mencionarlos en Discussion como contexto.

---

# ¿Qué sería entonces la “novedad metodológica”?

Esta parte hay que cuidarla, porque *Methods and Applications* exige **novelty + broad applicability**.

No diría que inventamos:

* hierarchical learning;
* PU learning;
* PLMs;
* homology-aware splitting.

Individualmente no son nuevos.

La metodología sería la **integración formal de tres decisiones que normalmente se tratan por separado**:

[
\boxed{
\text{evidence-aware task construction}
}
]

*

[
\boxed{
\text{representation-controlled functional modelling}
}
]

*

[
\boxed{
\text{homology-aware evaluation}
}
]

dentro de un protocolo reproducible para problemas de **protein functional classification with incomplete annotation**.

Eso sí puede tener broad applicability.

El software podría recibir una tabla:

```text
protein_id
sequence
label
evidence_state
```

y generar:

```text
task definitions
positive / negative / unlabeled sets
homology clusters
train / val / test
representations
models
evaluation
```

Es decir, PhageLysData sería **el showcase**, no una dependencia del framework.

Eso me gusta bastante más para *Methods and Applications*.

---

# Y haría exactamente **3 figuras + 1 tabla**

## **Figure 1 — Method**

Una multipanel.

**A. Evidence states**

[
P/N/U
]

**B. Hierarchical task**

[
Lytic?
\rightarrow
Endo/VAL/DPO
]

**C. Representation**

descriptors → PLM.

**D. Homology-aware evaluation**

clusters → train/val/test.

Una sola figura deja claro **qué aportamos metodológicamente**.

---

# **Figure 2 — What changes the classifier?**

La figura cuantitativa principal.

### Panel A

Task definition.

Conventional negatives vs strict negatives vs PU si funciona.

### Panel B

Representation.

Descriptors vs ESM2 vs ProtT5/ESM-C.

### Panel C

Flat vs hierarchical/multilabel.

### Panel D

Per-class results.

Todo puede entrar en una figura compacta.

El mensaje sería:

> **task definition and representation determine the apparent separability of functional classes.**

---

# **Figure 3 — Does it generalize?**

Esta es la figura fuerte para Protein Science.

### Panel A

Random vs homology-aware performance.

### Panel B

[
Performance
\text{ vs }
maximum\ identity\ to\ training
]

### Panel C

Confusion/recall por:

* Endo;
* VAL;
* DPO.

### Panel D

Quizá calibration/error landscape, si aporta.

O simplemente distribución de correct/incorrect predictions en el espacio de similitud.

El mensaje:

> **sequence-function recognition progressively deteriorates outside homologous sequence space, and the extent depends on representation/task formulation.**

Eso ya es un resultado bastante protein-centered.

---

# **Table 1 — Dataset and experimental design**

Aquí descargaría muchísimo contenido.

| Component       | Definition                            |  N | Notes        |
| --------------- | ------------------------------------- | -: | ------------ |
| Positive        | supported lytic                       |  … | PhageLysData |
| Negative        | curated non-lytic phage               |  … | external     |
| Unlabeled       | phage proteins without lytic evidence |  … | external     |
| Endolysin       | …                                     |  … |              |
| VAL             | …                                     |  … |              |
| DPO             | …                                     |  … |              |
| Representations | descriptors / PLMs                    |    |              |
| Splits          | random / homology-aware               |    |              |

Así no gastamos figuras explicando números.

---

# El paper además puede ser bastante corto

Yo lo imagino con resultados en apenas cuatro secciones:

### **1. An evidence-aware framework for phage lytic enzyme classification**

Data + method.

### **2. Task formulation and representation determine functional separability**

Q1 + Q2.

### **3. Homology-aware evaluation reveals the limits of sequence-based functional recognition**

Q3.

### **4. A reusable workflow for protein classification under incomplete annotation**

Generalización metodológica + software.

Y Discussion corta.

Esto cabe cómodamente alrededor de **4,000–4,500 palabras**.

---

# Algo que reconsideraría: HPO

Lo haría, pero **no lo estudiaría como pregunta científica**.

Pipeline:

1. seleccionar task;
2. seleccionar representations;
3. seleccionar 2–3 learners razonables;
4. HPO del finalista;
5. locked test.

No necesitamos una sección:

> “Effect of hyperparameter optimization”.

No aporta a la historia.

---

# Learners

También simplificaría muchísimo.

Yo tendría:

* Logistic Regression;
* SVM;
* XGBoost o RF.

Quizá MLP únicamente para PLMs.

Con eso basta.

Si el SVM con ProtT5 gana, maravilloso.

No necesitamos demostrar que probamos todo sklearn.

---

# Y PU learning: no lo congelaría todavía como obligatorio

Primero hay que mirar cómo podemos construir `N` y `U`.

Si conseguimos:

* **N fuerte y suficiente**, PU puede ser innecesario.
* `N` muy pequeño pero `U` enorme → PU se vuelve muy atractivo.

Por eso, **la primera tarea antes de programar modelos debería ser auditar qué background podemos realmente construir**.

Ese resultado determina el diseño.

---

# Así que congelaría el scope de Paper 2 de esta manera

### **Core obligatorio**

**1. Evidence-aware task definition.**

**2. Hierarchical/multilabel classification.**

**3. Descriptors vs a small number of PLM representations.**

**4. Random vs homology-aware evaluation.**

**5. Generalization as a function of sequence similarity.**

**6. Reproducible implementation/software.**

### **Condicional**

PU learning.

### **Fuera**

Benchmark externo, discovery, ranking, patrones profundos, estructuras, design, experimental validation.

---

Y para mí hay una decisión que debemos tomar **antes de hacer cualquier otra cosa**:

> **¿PhageLysData nos permite construir únicamente el universo positivo y funcional, o también contiene evidencia suficiente para declarar proteínas realmente no-líticas?**

Mi lectura actual es que **para un verdadero detector Stage 1 tendremos que agregar un background externo de proteínas de fagos**, mientras que para Stage 2 PhageLysData prácticamente nos da todo.

Entonces el próximo diseño que yo haría no es todavía el modelo: es definir exactamente **P / N / U y las tres tareas predictivas**. Una vez congelemos eso, el resto del paper cae casi solo.

---

Sí, papu. **Ahí sí: me había ido demasiado al “minimum viable paper” y le quité justamente la parte que hacía interesante la pregunta.** Si estamos dispuestos a usar PhageLysData como universo positivo y construir inteligentemente desde UniProt los universos `N/U`, podemos recuperar el diseño fuerte **sin convertirlo en un monstruo**.

Para *Methods and Applications*, yo haría que el método responda **cuatro preguntas centrales**. No seis, no diez.

## La columna vertebral que me convence

### **RQ1 — What is a negative for phage lytic enzyme prediction?**

Esta me parece probablemente **la pregunta más sexy del trabajo**.

Partimos de:

[
P=\text{PhageLysData supported lytic enzymes}
]

y desde UniProt construimos dos universos conceptualmente distintos:

[
N=\text{phage proteins with evidence for another function}
]

[
U=\text{phage proteins without evidence of lytic function}
]

Entonces podemos comparar tres formulaciones:

[
P;vs;N
]

**strict supervised**

versus

[
P;vs;(N+U)
]

**conventional “everything else is negative”**

versus

[
P;vs;U
]

mediante **PU learning**.

La pregunta no es “qué algoritmo gana”.

Es:

> **How does the biological definition of the negative class alter what a model learns as phage lytic function?**

Y la hipótesis:

> **Absence of lytic annotation is not equivalent to evidence of non-lytic function, and treating both as negatives should change apparent performance and decision boundaries.**

Esto es generalizable mucho más allá de fagos.

---

# **RQ2 — How should the functional problem itself be formulated?**

Acá recuperaría tu idea original de jerarquía.

Comparamos:

### Flat

[
Endolysin/VAL/Depolymerase/Other
]

versus

### Hierarchical

[
\text{Lytic?}
\rightarrow
\text{functional class}
]

versus, cuando la evidencia lo permita:

### Hierarchical multilabel

[
\text{Lytic?}
\rightarrow
[Endo,\ VAL,\ Depol]
]

La pregunta es:

> **Does a biologically structured task formulation improve recognition relative to forcing all proteins into mutually exclusive functional classes?**

No estamos diciendo que inventamos hierarchical learning —SUBLYME ya existe—.

Lo novedoso es que **la definición deriva explícitamente de estados de evidencia**, incluyendo proteínas con funciones no necesariamente exclusivas.

Además, esto conecta RQ1 y RQ2:

[
\boxed{
\text{Evidence}
\rightarrow
\text{Task definition}
\rightarrow
\text{Prediction}
}
]

Ese puede ser el núcleo del método.

---

# **RQ3 — Which representation preserves functional information beyond close sequence similarity?**

Acá sí representación, pero limpia.

No haría 10 PLMs.

Compararía **familias conceptualmente distintas**:

[
\text{Sequence similarity}
]

[
\text{Physicochemical / sequence descriptors}
]

[
\text{PLM representation A}
]

[
\text{PLM representation B}
]

Quizás:

* nearest-neighbor / MMseqs2 baseline;
* descriptors;
* ESM2;
* ProtT5 o ESM-C.

Y usamos inicialmente **el mismo learner** para todos.

La pregunta:

> **Do learned protein representations encode phage lytic function beyond information recoverable from conventional sequence similarity or physicochemical descriptors?**

Pero lo interesante viene cuando esto se cruza con RQ4.

---

# **RQ4 — How much apparent predictive performance survives sequence divergence?**

Esta sigue siendo para mí **la figura fuerte del paper**.

Separamos:

### Redundancy control

[
NR100,\ NR90,\ NR70,\ NR50
]

No necesariamente todas en el main; podemos seleccionar 2–3 después.

versus

### Evaluation strategy

[
Random
]

[
Homology\text{-}aware
]

Y después calculamos explícitamente:

[
S(x)=\max_{t\in Train}Identity(x,t)
]

para cada proteína test.

La verdadera pregunta pasa a ser:

> **How does functional recognition decay as proteins move away from the sequence space represented during training?**

Eso es muchísimo más interesante que:

> H30 gives lower F1 than random.

Y podemos estudiar:

[
Performance =
f(
Task,
Representation,
Sequence\ similarity
)
]

No hace falta hacer el factorial completo de todas las combinaciones. Diseñamos etapas.

---

# Las cuatro preguntas cuentan una historia perfecta

Mira cómo queda:

### **1. What constitutes evidence?**

[
P/N/U
]

↓

### **2. How should that evidence define the prediction task?**

[
flat/hierarchical/multilabel
]

↓

### **3. How should proteins be represented?**

[
sequence/descriptors/PLM
]

↓

### **4. Does the resulting model generalize beyond known protein space?**

[
random/homology-aware/distance-to-train
]

Eso **sí me parece Methods and Applications**.

No es un benchmark de modelos.

Es un **framework para diseñar protein-function classification problems when functional annotation is incomplete and sequence redundancy is substantial.**

Phage lytic enzymes son el caso de estudio.

---

# Y UniProt nos permite hacerlo mucho mejor

Estoy completamente contigo en esto.

No necesitamos que PhageLysData contenga todos los negativos.

De hecho, metodológicamente puede ser incluso mejor **no construirlos desde PhageLysData**, porque queremos una definición independiente.

Yo imaginaría tres queries/conjuntos externos.

## **N — Evidence-backed non-lytic**

Proteínas de fagos con función anotada inequívocamente diferente.

Por ejemplo familias de:

* capsid proteins;
* terminases;
* DNA/RNA polymerases;
* helicases;
* replication proteins;
* transcriptional proteins;
* packaging proteins;
* integrases/recombinases;

etc.

Preferiblemente reviewed cuando el tamaño lo permita.

Aquí tenemos:

[
\text{evidence of another function}
]

no simplemente ausencia de lysin.

---

## **U — No evidence**

Proteínas de fagos donde:

* no aparece evidencia lítica;
* no podemos afirmar que sean negativas;
* podrían ser hypothetical/uncharacterized o anotaciones débiles.

Y excluimos cualquier cosa con:

* lysin;
* endolysin;
* peptidoglycan hydrolase;
* muralytic;
* depolymerase;
* tailspike/depolymerase-related terminology según la tarea;
* dominios inequívocamente asociados cuando corresponda.

Así:

[
U \neq N.
]

Eso es fundamental.

---

## **Hard-N**

Incluso crearía un tercer subconjunto:

[
N_{hard}
]

con proteínas funcionalmente/contextualmente próximas:

* tail proteins;
* virion structural proteins;
* cell-envelope-interacting proteins;
* enzymes de fagos;
* proteins from the **same genera/proteomes** cuando podamos.

Entonces podemos preguntar:

> Does a model still recognize lytic function when trivial phage/background signals are removed?

Ese experimento me encanta.

---

# Y ahora sí PU learning se vuelve científicamente justificado

No:

> “metimos PU porque es fancy”.

Sino porque tenemos:

[
P=\text{evidence positive}
]

[
N=\text{evidence negative}
]

[
U=\text{absence of evidence}.
]

Eso es literalmente la situación para la que PU learning existe.

Incluso podríamos tener una conclusión muy interesante aunque PU **no gane**.

Por ejemplo:

> Strict supervised learning with evidence-backed negatives generalized better than PU learning.

También sería un hallazgo.

Lo importante es que por primera vez podemos **comparar estados de evidencia de manera controlada**.

---

# ¿Dónde entra la optimización de hiperparámetros?

No como pregunta.

Después de responder preliminarmente RQ1–RQ3:

[
best\ task
+
best\ representation
+
few\ learners
]

hacemos HPO únicamente sobre las configuraciones finalistas.

Por ejemplo:

* LR;
* SVM;
* RF/XGB;
* small MLP.

Y listo.

HPO va a Methods/Supplementary.

---

# ¿Benchmark contra SUBLYME/DposFinder?

Sigo estando contigo:

### **No como experimento principal.**

Las tareas son demasiado diferentes.

En Discussion:

> SUBLYME solves X, DposFinder solves Y, whereas our study addresses how evidence and sequence relatedness affect the formulation and evaluation of the prediction problem.

Perfecto.

No necesitamos fabricar una competencia artificial.

---

# Y seguimos perfectamente dentro de **3 figuras + 1 tabla**

## **Figure 1 — Evidence defines the prediction problem**

Multipanel:

**A.** PhageLysData `P` + UniProt `N/U`.

**B.** strict N / broad N / PU.

**C.** flat vs hierarchical/multilabel.

**D.** homology-aware methodological workflow.

Esta figura cuenta **todo el método**.

---

## **Figure 2 — What determines functional recognition?**

Puede ser una figura matriz muy fuerte.

**A.** P/N vs P/(N+U) vs PU.

**B.** flat vs hierarchy.

**C.** descriptors vs PLMs.

**D.** quizás interaction:

[
Task\times Representation
]

La pregunta visual:

> **What matters most: how we define the biological task or which representation/model we use?**

Eso es muy bonito.

---

## **Figure 3 — From apparent performance to real generalization**

**A.**

Random vs homology-aware.

**B.**

[
Performance
\text{ vs }
max\ sequence\ identity\ to\ train
]

**C.**

por clase:

Endo / VAL / Depol.

**D.**

quizá strict-N vs PU en remote sequence space.

Este último podría ser muy potente:

> **Does evidence-aware task formulation matter more as proteins become more remote?**

Eso conecta todo el paper.

---

## **Table 1 — Experimental universe**

Una sola tabla grande:

| Universe/task | Definition                      | Source       | Evidence requirement         |  N |
| ------------- | ------------------------------- | ------------ | ---------------------------- | -: |
| P             | supported phage lytic           | PhageLysData | positive evidence            |    |
| N             | confirmed alternative function  | UniProtKB    | explicit functional evidence |    |
| Hard-N        | functionally/contextually close | UniProtKB    | explicit                     |    |
| U             | no lytic evidence               | UniProtKB    | unlabeled                    |    |
| Endolysin     | …                               | PhageLysData | …                            |    |
| VAL           | …                               |              |                              |    |
| DPO           | …                               |              |                              |    |

Y listo.

---

# Incluso veo una pregunta integradora muy buena

Podría ser casi la frase de Abstract:

> **Are advances in phage lytic enzyme prediction driven more by model choice, or by how functional evidence and sequence relatedness are handled?**

Esa pregunta **sí me parece sexy**.

Y probablemente la respuesta termine siendo algo como:

[
\boxed{
\text{Data/task design}

>

\text{algorithmic sophistication}
}
]

que además está muy alineado con buena metodología moderna.

---

# Yo congelaría ahora estas cuatro RQs

**RQ1. How does the definition of negative and unlabeled evidence affect phage lytic enzyme detection?**

**RQ2. Does biologically structured hierarchical/multilabel classification improve functional assignment compared with conventional flat classification?**

**RQ3. Which protein representations retain functional information as sequence similarity decreases?**

**RQ4. How strongly do redundancy and train–test homology influence apparent versus out-of-homology generalization?**

Y una **meta-question** que une todo:

> **Which methodological decisions are required for reliable protein-function classification when positive evidence is strong but negative evidence is incomplete?**

Ahí, papu, **volvemos al paper que originalmente me entusiasmó**, pero ahora con una ventaja: ya sabemos exactamente qué nos da PhageLysData y qué único bloque tenemos que traer inteligentemente desde UniProt. No necesitamos bajar la ambición científica; necesitamos **diseñar bien los universos P/N/U**.

---

Sí, papu. **Con Full-Length Article yo recuperaría la versión ambiciosa**, pero con una frontera clara: Paper 2 debe llegar hasta **“podemos reconocer de forma confiable y demostrar utilidad para discovery”**, sin resolver todavía *por qué* funcionan las proteínas, cómo rankear todo el universo, ni cómo diseñarlas.

Creo que la pregunta global podría ser:

> **Can phage lytic enzyme function be learned reliably from protein sequence when functional annotation is incomplete and sequence homology is explicitly controlled, and can such models identify previously unrecognized, functionally plausible candidates?**

Eso ya da un paper completo.

---

# Las 5 preguntas que yo congelaría

## **RQ1 — What constitutes a negative for phage lytic enzyme prediction?**

Esta sigue siendo para mí **la pregunta más diferenciadora**.

Construimos tres estados:

[
P=\text{positive evidence}
]

desde el core soportado de PhageLysData.

[
N=\text{evidence of an alternative/non-lytic phage function}
]

desde UniProt mediante queries controladas.

[
U=\text{phage protein with no evidence of lytic function}
]

también desde UniProt, pero **sin convertir ausencia de evidencia en evidencia negativa**.

Entonces estudiamos:

[
P;vs;N
]

strict supervised,

[
P;vs;(N+U)
]

conventional negative assignment,

y

[
P;vs;U
]

con una formulación PU.

### Pregunta real

> **Does treating absence of annotation as negative evidence change what the model learns as phage lytic function?**

### Hipótesis

Sí, y probablemente el efecto sea especialmente importante en proteínas alejadas del training set.

Esto además tiene **broad applicability** para protein-function prediction en general.

---

# **RQ2 — How should phage lytic enzyme function be formulated computationally?**

Aquí va tu idea jerárquica, pero bien planteada.

### Flat

[
Endolysin/VAL/Depolymerase/Other
]

versus:

### Hierarchical

[
\text{Phage lytic/envelope-degrading enzyme?}
]

↓

[
Endolysin / VAL / Depolymerase
]

versus:

### Hierarchical multilabel

[
[y_{Endo},y_{VAL},y_{DPO}]
]

cuando la evidencia permita funciones no exclusivas.

Esto último aprovecha algo importante de PhageLysData: **no tenemos por qué forzar toda la biología a categorías mutuamente excluyentes**.

### Pregunta

> **Does an evidence-consistent hierarchical/multilabel formulation provide more reliable functional assignment than conventional flat classification?**

No vendemos la jerarquía como invento. Vendemos **la formulación de la tarea desde el estado de evidencia**.

Y reportamos dos resultados:

[
\text{conditional Stage-2 performance}
]

y

[
\text{true end-to-end performance}.
]

---

# **RQ3 — What representation of protein sequence is required to recognize lytic function?**

Acá sí volvería a algo un poco más completo de lo que discutimos recién.

No 20 PLMs, pero sí comparar **familias de información**:

[
\text{sequence similarity}
]

vs.

[
\text{sequence/physicochemical descriptors}
]

vs.

[
\text{protein language representations}.
]

Dentro de PLM sí podemos poner 3–4 representantes:

**ESM2 / ProtT5 / ESM-C / quizá Ankh**, dependiendo de disponibilidad y coste.

Lo importante es evaluar con las **mismas particiones y learners**.

### Pregunta

> **Do learned representations encode functional information that remains recoverable when conventional sequence similarity becomes weak?**

Esto conecta representación con biología y evita el típico:

> ESM-C = 0.918, ProtT5 = 0.912.

Lo interesante será:

[
\Delta_{\mathrm{PLM-homology}}
==============================

f(\text{sequence divergence}).
]

---

# **RQ4 — How much predictive performance survives removal of sequence redundancy and close homology?**

Este es el segundo eje grande.

Y separaría explícitamente:

[
\textbf{dataset redundancy}
\neq
\textbf{train-test homology}.
]

Podemos investigar diferentes reducciones de redundancia:

[
NR100,;NR90,;NR70,;NR50
]

y diferentes regímenes de partición:

[
Random
]

versus

[
Cluster/Homology-aware.
]

No necesariamente entrenaremos **todas las combinaciones factorialmente**. Primero caracterizamos el efecto y después elegimos escenarios representativos.

Pero el resultado estrella sería continuo.

Para cada proteína test:

[
S(x)=\max_{t\in Train}
Identity(x,t)
]

y entonces:

[
Performance=f(S(x)).
]

Por ejemplo:

[

> 90,;70-90,;50-70,;30-50,;<30%.
> ]

### Pregunta

> **At what distance from known sequence space does phage lytic enzyme function cease to be reliably recognizable?**

Esa pregunta me parece extremadamente *Protein Science*.

Porque ya no hablamos de split.

Hablamos de **persistencia de señal funcional bajo divergencia de secuencia**.

---

# **RQ5 — Can the resulting framework identify previously unrecognized phage lytic enzymes?**

Y acá cerramos el paper de manera aplicada.

Pero **sin comernos Paper 4**.

Creamos un universo de discovery independiente:

[
D=\text{phage proteins without current lytic annotation}
]

de UniProt, que **no haya participado en entrenamiento ni como U de PU**.

Idealmente podemos hacer incluso un holdout muy limpio por:

* proteomas completos;
* taxa;
* release/time;
* fuente.

Aplicamos **una única configuración congelada** después de RQ1–RQ4.

Y elegimos un pequeño conjunto, por ejemplo:

[
6-12\ candidates
]

no 5,000.

Los elegimos mediante criterios simples:

**high predicted probability + low similarity to PhageLysData positives + representation of the different functional classes.**

Nada de score sofisticado todavía.

Porque:

[
Probability+Novelty+Uncertainty+Patterns+Ranking
]

es **Paper 4**.

Aquí solamente queremos demostrar:

> **The classification framework can nominate biologically plausible proteins outside the annotated space.**

---

# Y acá entra la validación in silico

Yo no haría simplemente:

> candidate → docking → yay, −7.3 kcal/mol.

Eso es débil.

Haría una **validación ortogonal escalonada**.

Para cada candidato seleccionado:

[
\text{sequence prediction}
]

↓

[
\text{domain/function evidence}
]

↓

[
\text{structure}
]

↓

[
\text{catalytic plausibility}
]

↓

[
\text{substrate interaction}
]

Es decir, intentar demostrar que **información que el classifier no utilizó directamente** es consistente con la función predicha.

---

## 1. Domain-level validation

Buscar posteriormente:

* InterPro;
* Pfam;
* HMMER;
* conserved catalytic domains.

Pero **estos resultados no entran al classifier**.

Son evidencia ortogonal.

Un candidato especialmente interesante sería:

[
\text{high ML confidence}
+
\text{low sequence identity}
+
\text{weak/no canonical HMM evidence}
]

pero con evidencia estructural posterior.

---

## 2. Structural validation

Predecimos estructura de los candidatos externos y buscamos:

* Foldseek;
* DALI si corresponde;
* structural neighbors;
* catalytic fold similarity;
* domain organization.

Esto puede ser incluso más convincente que docking.

Un candidato con:

[
<30%\ sequence\ identity
]

pero una arquitectura estructural compatible con una familia lítica sería una demo preciosa.

---

# 3. Catalytic-site plausibility

Para candidatos donde conocemos el tipo catalítico esperado:

* localizar pocket;
* residuos catalíticos;
* geometría;
* conservación espacial.

Esto es importantísimo porque dos proteínas pueden compartir fold pero no necesariamente función.

---

# 4. Docking

Sí lo haría, **pero únicamente donde tenga sentido biológico**.

Para endolisinas/VALs ya existe precedente reciente para docking contra fragmentos de peptidoglicano usando AutoDock Vina, 3D-RISM, HADDOCK y posteriormente MD. ([PubMed Central (PMC)][1])

Se han utilizado, entre otros:

* NAG;
* NAM;
* NAG–NAM;
* muramyl dipeptide;
* peptidoglycan-derived peptides;
* fragmentos mayores de PG. ([Nature][2])

Incluso un trabajo de 2026 combina docking contra muramyl dipeptide, MD y MM/PBSA con ensayos experimentales de endolisinas. ([Nature][3])

### Pero hay una sutileza importante

No usaría **el mismo ligando para todas las endolisinas**.

Primero inferimos mecanismo/domain family:

[
glycosidase
]

[
amidase
]

[
endopeptidase
]

etc.

Y seleccionamos un fragmento compatible con el tipo de enlace que debería reconocer.

Entonces docking tiene una hipótesis funcional.

---

# Y con depolymerases sería todavía más cuidadoso

Para una DPO, decir:

> “docking contra un polisacárido cualquiera”

no significa mucho.

La actividad depende mucho de:

* cápsula;
* EPS;
* LPS;
* composición del polisacárido;
* host/serotype.

Por eso, para DPO haría:

### si conocemos host/substrate context

→ docking contra oligosacárido relevante.

### si no conocemos substrate

→ **structural/domain validation**, no forzaría docking.

Eso hace el estudio más serio.

---

# Incluso podríamos hacer MD solo para 2–3 candidatos estrella

No hace falta para todos.

Por ejemplo:

[
Candidate_{Endo}
]

[
Candidate_{VAL}
]

[
Candidate_{DPO}
]

si DPO tiene un substrate conocido.

Docking genera pose.

Después:

[
50-100,ns
]

de MD —o lo que definamos metodológicamente después— y una estimación comparativa tipo MM/PBSA.

No vendería esto como prueba de actividad.

Lo llamaría:

> **orthogonal in silico functional support.**

Eso es científicamente más correcto.

Docking reciente en endolisinas se ha utilizado precisamente como apoyo para mecanismos de reconocimiento de peptidoglicano, y cuando se ha comparado con actividad experimental puede ser consistente con diferencias funcionales, pero no sustituye esa validación. ([PubMed][4])

---

# Entonces la arquitectura conceptual completa queda preciosa

[
\boxed{
\textbf{Evidence}
}
]

**What is positive, negative and unlabeled?**

↓

[
\boxed{
\textbf{Task}
}
]

**How should function be formulated?**

↓

[
\boxed{
\textbf{Representation}
}
]

**What information encodes function?**

↓

[
\boxed{
\textbf{Generalization}
}
]

**Does it survive sequence divergence?**

↓

[
\boxed{
\textbf{Application}
}
]

**Can it find something previously unannotated?**

↓

[
\boxed{
\textbf{Orthogonal validation}
}
]

**Is the prediction structurally and biochemically plausible?**

Eso sí es un **Full-Length Article**.

---

# Y mantiene perfectamente separados los siguientes papers

| Paper                  | Pregunta                                                                    |
| ---------------------- | --------------------------------------------------------------------------- |
| **2 — Classification** | **Can we recognize them reliably?**                                         |
| **3 — Patterns**       | **What molecular features distinguish and define them?**                    |
| **4 — Mining**         | **How do we systematically rank and discover the best unknown candidates?** |
| **5 — Design**         | **Can we intentionally create functional new enzymes?**                     |

En Paper 2 hacemos una **small discovery demonstration**, no un discovery platform.

En Paper 3 sí diseccionamos:

* motifs;
* residues;
* physicochemical signatures;
* architectures;
* learned patterns;
* design rules.

En Paper 4 utilizamos esos patrones junto con classifier, uncertainty, novelty, structural support, etc., para un **ranking formal y mining masivo**.

Y Paper 5 los convierte en diseño.

---

# Qué haría con HPO y algoritmo

Ahora que tenemos esta historia, **dejan de ser research questions**.

Son decisiones técnicas necesarias.

Podemos explorar:

[
LR,\ SVM,\ RF/XGB,\ MLP
]

hacer tuning correctamente sobre train/validation, seleccionar y congelar.

Eso va en Methods/Supplementary.

La pregunta científica no es:

> which hyperparameter gives better MCC?

---

# Y tampoco haría benchmark contra SUBLYME como centerpiece

Seguimos comparando cosas distintas.

Sí incluiría **Table 1 — State of the art and methodological scope**, donde mostremos conceptualmente:

| Method              | Negatives explicit | PU | Multilabel | Homology-aware | Representation comparison | Discovery | Orthogonal validation |
| ------------------- | -----------------: | -: | ---------: | -------------: | ------------------------: | --------: | --------------------: |
| existing approaches |                  … |  … |          … |              … |                         … |         … |                     … |
| **This work**       |                  ✓ |  ✓ |          ✓ |              ✓ |                         ✓ |         ✓ |                     ✓ |

Eso posiciona la contribución sin hacer un benchmark metodológicamente dudoso.

---

# Y ahora no me limitaría artificialmente a tres figuras

Tampoco necesitamos diez.

Yo veo **4 figuras grandes + 1 tabla** como punto dulce:

**Figure 1 — Evidence-aware task formulation.**

P/N/U + flat/hierarchical/multilabel + workflow.

**Figure 2 — Task definition and representation.**

RQ1–RQ3.

**Figure 3 — Redundancy, homology and functional generalization.**

RQ4, con `performance vs identity-to-train` como centerpiece.

**Figure 4 — Prospective discovery demonstration and orthogonal validation.**

Candidatos → sequence novelty → structural/domain support → docking/MD para candidatos representativos.

**Table 1 — existing approaches vs methodological questions addressed here.**

Con eso tenemos un paper **mucho más completo sin transformarlo en un monstruo**.

---

## Para mí, el mensaje final del Paper 2 debería ser este

No:

> *We achieved 0.94 F1 for classification of phage lytic enzymes.*

Sino:

> **Reliable recognition of phage lytic enzyme function requires explicit treatment of functional evidence and sequence relatedness. When these factors are controlled, learned protein representations retain predictive signal beyond close homologs and can nominate previously unannotated proteins with independent structural and substrate-level support.**

Ahí sí, papu. **Eso tiene cuerpo de Full-Length Article para Protein Science**, deja una demo concreta de discovery y, al mismo tiempo, nos deja servido el pase para que Paper 3 explique *qué patrones lo hacen posible* y Paper 4 convierta esa capacidad en un verdadero sistema de mining y priorización.

[1]: https://pmc.ncbi.nlm.nih.gov/articles/PMC11597070/?utm_source=chatgpt.com "Molecular Docking of Endolysins for Studying Peptidoglycan Binding Mechanism - PMC"
[2]: https://www.nature.com/articles/s41598-025-27272-7?utm_source=chatgpt.com "Combined experimental and in silico elucidation of KP27 endolysin reveals a phage derived antibacterial with pH and thermal robustness | Scientific Reports"
[3]: https://www.nature.com/articles/s41598-026-57871-x?utm_source=chatgpt.com "Integrated computational and experimental benchmarking of Bacillus phage endolysins reveals the relationship between peptidoglycan-fragment recognition descriptors and antibacterial performance | Scientific Reports"
[4]: https://pubmed.ncbi.nlm.nih.gov/40680964/?utm_source=chatgpt.com "Improved lytic action of engineered phage-encoded endolysin and docking insights into its action on bacterial peptidoglycan - PubMed"

---

Sí, papu. **Ahora que tengo delante el Data Report, podemos dejar de diseñar esto en abstracto y construir el Paper 2 encima de la arquitectura real de PhageLysData.** Y hay varias cosas que cambian para bien.

La ventaja enorme es que PhageLysData ya separa explícitamente **Core, Prediction Extension y Context**, conserva las afirmaciones originales y sus linajes de evidencia, y además no utilizó similitud/homología para definir entidades. El Core tiene 11,867 secuencias con soporte no predictivo, mientras la Prediction Extension contiene 745,092 candidatos sostenidos exclusivamente por predicciones. Eso nos permite separar **ground truth, unlabeled/background y discovery space** de una manera bastante más limpia que si partiéramos desde cero. 

Pero hay una regla que congelaría inmediatamente:

> **Prediction Extension ≠ negative. Context ≠ negative.**

El propio Data Report deja explícito que ausencia de target support no implica ausencia biológica y que la Prediction Extension son candidatos computacionales. Usarlos como negativos sería justamente cometer el error que queremos estudiar. 

---

# Yo construiría cinco universos de datos

Esto nos ordena prácticamente todo el paper:

[
\boxed{P}
\quad
\boxed{N}
\quad
\boxed{U}
\quad
\boxed{D}
\quad
\boxed{V}
]

donde:

* **P = positive evidence**
* **N = explicit non-target evidence**
* **U = unlabeled / no evidence**
* **D = discovery universe**
* **V = orthogonal validation universe**

Y cada uno viene de una fuente diferente por una razón deliberada.

---

## 1. `P` — Positive universe

### Fuente: **PhageLysData Core v1.0**

No necesitamos reconstruir positivos.

El Core ya contiene exactamente lo que queremos: secuencias con al menos una afirmación objetivo **no exclusivamente predictiva**, conectadas a provenance y evidence lineage. 

Tenemos:

| canonical class                   |     n |
| --------------------------------- | ----: |
| Endolysin                         | 9,888 |
| Depolymerase                      |   712 |
| VAL                               |   241 |
| Other phage-envelope lytic enzyme |    27 |
| Multiple non-exclusive            |   999 |

Esto inmediatamente nos dice algo importante.

### No usaría `canonical_class` directamente como ground truth de Stage 2

Porque `multiple_nonexclusive` **no es una función biológica**.

Es una etiqueta derivada del resource porque una misma secuencia tiene soporte no predictivo para más de una clase. El Data Report preserva justamente esas afirmaciones en vez de obligarlas a una única clase. 

Entonces para Paper 2 deberíamos regresar a la **assertion layer**:

```text
sequence_entity_id
    ├── endolysin = 0/1
    ├── VAL = 0/1
    ├── depolymerase = 0/1
    └── other_lytic = 0/1
```

Eso nos permite construir la tarea correctamente.

---

# También usaría los tiers, pero con cuidado

El Core tiene:

* Tier 1 = 90
* Tier 2 = 784
* Tier 3 = 9,836
* Tier 4 = 1,157

y la distribución no es simétrica: Tier 1 son exclusivamente depolimerasas y la mayoría de `multiple_nonexclusive` cae en Tier 4. 

Por eso **no podemos decir simplemente “entrenamos solo con Tier 1+2”** porque destruiríamos VAL y gran parte del universo endolysin.

Yo definiría:

### `P_main`

[
Tier1+Tier2+Tier3
]

como positivos principales.

### `P_ambiguous`

[
Tier4
]

como **challenge/sensitivity set**, no necesariamente training primario.

Esto añade incluso una pregunta secundaria interesante:

> Does classification reliability deteriorate for proteins supported by incomplete or conflicting functional evidence?

Pero no hace falta convertirla en RQ principal.

---

# 2. `N` — True negative evidence

### Fuente: **UniProtKB**, mediante queries explícitamente negativas

Acá sí debemos construir algo nuevo.

Y creo que debemos reutilizar una fortaleza del propio PhageLysData: ya existe un **frozen query registry** para recuperar positivos mediante:

* protein identity;
* function;
* catalytic terminology;
* GO;
* EC;
* substrate terminology;
* structural terminology;
* phage context. 

Perfecto.

Podemos construir **la operación complementaria**, pero no simplemente:

```text
NOT lysin
```

Eso sería horrible.

Queremos:

[
N =
\text{phage protein}
+
\text{explicit evidence for another function}
]

---

## Yo dividiría `N` en dos familias

### `N_easy`

Funciones inequívocamente diferentes:

* major capsid proteins;
* portal proteins;
* terminase subunits;
* DNA polymerases;
* helicases;
* primases;
* integrases;
* recombinases;
* transcriptional regulators;
* nucleotide-metabolism proteins;
* packaging proteins.

Estas proteínas nos dan un universo **evidence-backed non-lytic**.

---

### `N_hard`

Y esta me parece muchísimo más interesante.

Proteínas biológicamente próximas al sistema lítico pero que **no son las enzimas objetivo**:

* holins;
* antiholins;
* spanins;
* some baseplate proteins;
* virion structural proteins;
* non-enzymatic tail proteins;
* receptor-binding proteins sin evidencia depolimerasa.

Con muchísimo cuidado con las últimas categorías porque tail proteins/tailspikes pueden contener VALs o DPOs.

### La idea

[
P\ vs\ N_{easy}
]

pregunta:

> ¿reconocemos una proteína lítica frente a una proteína de fago cualquiera?

Mientras:

[
P\ vs\ N_{hard}
]

pregunta:

> **¿reconocemos realmente función lítica frente a proteínas que comparten contexto biológico?**

Esta segunda tarea es bastante más interesante.

---

# Regla de exclusión muy importante para `N`

Todo candidato negativo debería excluirse si su secuencia exacta aparece en:

[
Core \cup PredictionExtension
]

Yo probablemente excluiría incluso **todo PhageLysData**, incluyendo Context, para mantener los universos conceptualmente independientes.

Porque el recurso ya contiene 759,105 entidades exactas. 

Después controlaremos homología mediante clustering.

No excluiría proteínas simplemente porque tengan 30–50% de identidad a un positivo.

Eso artificialmente haría la tarea fácil.

---

# 3. `U` — Unlabeled universe

También desde **UniProtKB**, pero con una definición radicalmente diferente.

[
U =
\text{phage protein}
--------------------

## \text{positive evidence}

\text{explicit alternative-function evidence}
]

Es decir:

> **no sabemos qué es.**

Aquí entrarían principalmente:

* hypothetical protein;
* uncharacterized protein;
* poorly annotated proteins;
* generic phage protein annotations;
* proteínas sin GO/EC/function suficientes para asignarlas.

Y aquí hay una sutileza fundamental:

### NO excluiría automáticamente proteínas porque InterPro/Pfam sugieran algo lítico.

Porque eso sería usar información que justamente podría revelar que dentro de `U` existen positivos ocultos.

Para construir `U` usamos principalmente **evidencia anotacional disponible en el momento de definir la tarea**.

Después podemos utilizar:

* InterPro;
* Pfam;
* structure;
* Foldseek;

como **evidencia independiente posterior**.

Eso deja muchísimo más limpia la historia.

---

# Con P / N / U respondemos RQ1

Entonces:

### Task A

[
P\ vs\ N
]

**evidence-backed supervised**

### Task B

[
P\ vs\ (N+U)
]

**absence-is-negative conventional formulation**

### Task C

[
P + U
]

**PU formulation**

Y quizás:

### Task D

[
P\ vs\ N_{hard}
]

**biologically difficult discrimination**

Entonces la pregunta:

> **What constitutes a negative for phage lytic enzyme prediction?**

ya tiene un diseño experimental concreto.

Y me gusta muchísimo.

---

# 4. RQ2 — ¿De dónde salen las clases?

Aquí no necesitamos absolutamente nada externo.

### Todo viene de PhageLysData Core.

Pero desde las **biological assertions**, no desde `canonical_class`.

Tendríamos:

## Level 1

[
P(\text{phage lytic enzyme})
]

usando:

**P vs N/U.**

## Level 2

Tres heads:

[
P(\text{endolysin})
]

[
P(\text{VAL})
]

[
P(\text{depolymerase})
]

y permitimos:

[
\sum y_k > 1.
]

`other_phage_envelope_lytic_enzyme`, con solamente **27 secuencias**, no debería ser una clase ML independiente. 

Yo haría:

> **lytic-positive but unresolved subtype**

para esas 27.

Eso incluso es conceptualmente bonito porque el hierarchical model puede decir:

```text
Lytic = yes
Subtype = unresolved
```

en vez de inventar una clase con n=27.

---

# El `multiple_nonexclusive` del ejemplo actual desaparece como clase

Esto me parece importante.

El Data Report actualmente demuestra cuatro clases y junta:

```text
other_phage_envelope_lytic_enzyme
+
multiple_nonexclusive
```

como `other`. Pero el propio manuscrito deja claro que es solamente una demostración de interoperabilidad, con un split random y sin control de homología ni HPO. 

Para Paper 2:

**NO hacemos eso.**

Los 999 `multiple_nonexclusive` vuelven a sus labels originales.

Eso por sí solo ya es una evolución metodológica muy natural desde Data Report → classification paper.

---

# 5. RQ3 — Representaciones

Aquí PhageLysData prácticamente nos hizo el trabajo.

Tenemos **11 PLMs + one-hot sobre exactamente las mismas 11,259 secuencias elegibles**, cubriendo:

* Ankh;
* ESM-C;
* ESM-2;
* Mistral-Prot;
* ProtBERT;
* ProtT5. 

Además tenemos:

* amino-acid composition;
* physicochemical descriptors;
* InterPro;
* Pfam.

Así que para **P** ya está prácticamente todo.

---

## Lo único que debemos hacer

Para `N` y `U`:

[
sequence
\rightarrow
\text{same representation contract}
]

mediante **Sylphy**.

Mismos:

* checkpoints;
* pooling;
* precision;
* sequence-length constraints.

Eso es importantísimo porque convierte a PhageLysData en el **reference representation contract** del nuevo trabajo.

Yo incluso consideraría explorar las **11 representaciones inicialmente**.

No necesariamente todas terminarán en las figuras.

Computacionalmente ya sabemos cómo producirlas y científicamente nos permite preguntar:

> Is generalization dependent on PLM family or scale?

Después las figuras pueden mostrar solamente familias o finalists.

---

# 6. RQ4 — Redundancia y generalización

### Fuente: las secuencias combinadas

[
P+N+U
]

No necesitamos ninguna base nueva.

Pero PhageLysData usa **exact sequence identity** para construir entidades y explícitamente **no utilizó homología para fusionarlas**. 

Eso es ideal.

Partimos de una base sin exact duplicates y hacemos nosotros:

[
MMseqs2
]

sobre el universo de clasificación completo.

Podemos obtener:

[
C_{90},C_{70},C_{50},C_{30}
]

o los thresholds que finalmente decidamos.

Y de la misma ejecución derivamos dos experimentos diferentes.

### A — Redundancy

¿Qué ocurre cuando representamos cada cluster mediante uno/varios representatives?

### B — Generalization

¿Qué ocurre cuando **ningún cluster se distribuye entre train/test**?

Eso responde limpiamente:

[
Performance=f(redundancy)
]

versus:

[
Performance=f(train-test\ homology)
]

Y de los alignments guardamos también:

[
S_i=\max Identity(x_i,Train)
]

para construir la curva continua de generalización.

---

# 7. `D` — Discovery universe

Aquí cambiaría algo importante respecto a lo que veníamos hablando.

### Yo NO usaría inicialmente Prediction Extension como “discovery”.

Es tremendamente tentador porque contiene **745,092 candidatos**, pero son explícitamente **prediction-only candidates** provenientes de sistemas que forman parte del ecosistema que PhageLysData integró. 

Entonces, si nuestro modelo encuentra uno:

> “nuevo candidato”

en realidad podría ser simplemente:

> “otro predictor ya lo había encontrado”.

Prediction Extension nos va a servir muchísimo más adelante para **Paper 4 / mining**.

---

# Para la demo del Paper 2 construiría `D` desde UniProt

Pero **independiente de N y U usados para training**.

Idealmente:

### opción más bonita

Un **UniProt release posterior al frozen release 2026_02** que utilizó PhageLysData.

El Data Report congela UniProtKB 2026_02. 

Entonces:

[
D =
\text{new phage proteins in a later UniProt snapshot}
]

que:

1. no existen por exact sequence en PhageLysData v1.0;
2. no tienen anotación lítica explícita;
3. no participaron en N/U;
4. pasan las reglas básicas de calidad.

Eso es una **temporal discovery demonstration**.

Conceptualmente es mucho más limpio que tomar random proteins de la misma snapshot.

---

# Incluso podemos tener una pequeña joya adicional

Si entre el freeze 2026_02 y la release posterior existen secuencias que:

### en la versión antigua

eran:

> uncharacterized / no evidence

y posteriormente pasan a ser anotadas como:

> endolysin / VAL / depolymerase,

tendríamos una especie de **retrospective temporal validation**.

Eso sería espectacular:

> Could the model identify proteins before their functional annotation became available?

No dependería de que haya miles.

Con unas pocas ya sería una validación preciosa.

Pero la considero **opportunistic**, no requisito del paper.

---

# 8. `V` — Orthogonal validation

Aquí podemos aprovechar PhageLysData otra vez.

El Core ya contiene:

* InterPro matches para 11,549;
* Pfam para 9,300;
* PDB para 77;
* AlphaFold DB para 672;
* 745 con alguna estructura existente. 

Eso nos permite construir un **reference library** de positivos conocidos.

---

## Para los nuevos candidatos D

Después de hacer la predicción —y solamente después—:

### V1 — functional/domain support

Correr:

* InterProScan;
* Pfam/HMMER.

Comparar arquitectura con Core.

---

### V2 — structural support

Generar estructura si no existe.

Después:

[
CandidateStructure
\xrightarrow{Foldseek}
KnownCoreStructures/PDB/reference
]

y evaluar:

* fold similarity;
* domain arrangement;
* catalytic-domain similarity.

Esto es mucho más fuerte como evidencia ortogonal que solamente docking.

---

### V3 — catalytic plausibility

Para los candidatos estrella:

* active-site position;
* catalytic residues;
* pocket geometry;
* residue arrangement.

Ahí empezamos a preguntarnos:

> ¿la estructura no solo parece una proteína de esa familia, sino que preserva una maquinaria catalítica plausible?

---

### V4 — docking / interaction

Y recién ahí docking.

Con **substrato específico para la función**.

Para endolysin/VAL:

[
\text{appropriate peptidoglycan fragment}
]

según catalytic family.

Para DPO:

[
\text{host-compatible oligosaccharide}
]

solamente cuando tengamos host/substrate suficientemente definido.

El Data Report conserva precisamente **phage and host information** dentro de las relaciones de entidad, aunque la disponibilidad deberá auditarse candidato por candidato. 

Si no sabemos el substrate de una DPO:

**no docking inventado**.

Fold/domain/catalytic-site support es mejor.

---

# Los controles del docking también pueden salir de PhageLysData

Eso me gusta bastante.

Para cada candidato:

### Positive control

Core protein:

* misma clase;
* idealmente mismo catalytic family;
* estructura PDB cuando exista.

### Negative control

Una proteína `N_hard` estructuralmente razonable pero no lítica.

Entonces no reportamos simplemente:

> docking score = −8.2 kcal/mol.

Reportamos:

[
Candidate
\quad vs\quad
KnownPositive
\quad vs\quad
NegativeControl.
]

Muchísimo más defendible.

---

# Entonces el mapa completo queda así

| Universe        | Fuente                                       | Rol                                         |
| --------------- | -------------------------------------------- | ------------------------------------------- |
| **P**           | PhageLysData Core                            | positive ground truth                       |
| **P ambiguous** | Core Tier 4                                  | evidence sensitivity/challenge              |
| **N_easy**      | UniProt phage proteins                       | explicit alternative function               |
| **N_hard**      | UniProt phage proteins                       | biologically adjacent negatives             |
| **U**           | UniProt phage proteins                       | no sufficient positive or negative evidence |
| **D**           | later/external UniProt snapshot              | unseen discovery space                      |
| **V_known**     | PhageLysData InterPro/Pfam/PDB/AFDB          | orthogonal positive references              |
| **V_candidate** | newly generated annotation/structure/docking | candidate validation                        |

**Prediction Extension queda fuera del training ground truth.**

Y yo la reservaría principalmente para Paper 4.

---

# Ahora las cinco preguntas tienen fuentes reales

## **RQ1**

### What constitutes a negative?

[
P_{PhageLysData}
+
N_{UniProt}
+
U_{UniProt}
]

---

## **RQ2**

### How should function be formulated?

[
Assertions_{PhageLysData}
\rightarrow
hierarchical/multilabel
]

Nada externo.

---

## **RQ3**

### What representation carries functional information?

[
11,PLMs + physicochemical
]

ya existentes para Core,

*

mismas representaciones generadas para UniProt N/U.

---

## **RQ4**

### Does recognition survive sequence divergence?

[
MMseqs2(P+N+U)
]

→ redundancy / homology-controlled splits / distance-to-training.

Nada externo adicional.

---

## **RQ5**

### Can the framework identify previously unrecognized proteins?

[
D_{temporal/external\ UniProt}
]

→ model

→ novel candidates

→ InterPro/Pfam

→ structure/Foldseek

→ catalytic site

→ targeted docking where biologically justified.

---

# Y hay una sexta pregunta que yo no convertiría necesariamente en RQ, pero que puede emerger

La arquitectura del Data Report nos permite preguntar:

> **Does evidence quality interact with sequence novelty?**

Porque tenemos:

[
Tier
\times
IdentityToTrain
\times
PredictionAccuracy.
]

Podría ocurrir que los Tier 4 sean difíciles simplemente porque están peor anotados, no porque el modelo falle biológicamente.

Eso sería un resultado interesante para Discussion o Supplementary.

---

## Lo primero que haría antes de desarrollar modelos

No entrenaría nada todavía.

Construiría un documento/manifest de **dataset design** con seis bloques:

1. `P_definition.tsv`
2. `N_query_registry.tsv`
3. `U_query_registry.tsv`
4. `exclusion_registry.tsv`
5. `D_discovery_definition.tsv`
6. `task_definitions.tsv`

Y especialmente congelaría las **queries de UniProt antes de mirar resultados de ML**.

Eso es importantísimo metodológicamente: no queremos redefinir N/U después porque una definición “funciona mejor”.

Primero definimos biológicamente:

[
P/N/U/D
]

después modelamos.

### Y ahora sí siento que Paper 2 está realmente conectado con PhageLysData

El Data Report ya termina diciendo explícitamente que sus workflows **no constituyen benchmarks controlados por homología ni validación biológica independiente**. 

Paper 2 empieza exactamente ahí:

> **PhageLysData provides the evidence architecture; Paper 2 asks what happens when that architecture is used to formulate and rigorously evaluate protein-function prediction.**

Eso hace que ambos papers se conecten muy naturalmente sin duplicarse.

Y para mí el siguiente paso debería ser **diseñar concretamente P/N/U/D**, empezando por las familias de queries que utilizaremos en UniProt y las reglas de exclusión/inclusión. Ahí ya podemos congelar la primera fase experimental del proyecto.

---

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


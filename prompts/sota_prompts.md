Sí, papu. Hice una revisión bastante amplia, incluyendo trabajos históricos, herramientas específicas de **endolisinas/VALs**, predictores de **depolimerasas**, pipelines de *mining*, PLMs, estructura y métodos generales de anotación de proteínas de fagos. La literatura reciente cambió muchísimo entre 2024 y 2026.

La conclusión importante es esta:

> **El campo ya pasó de “¿puedo clasificar una lysin?” a “¿puedo encontrar proteínas funcionales muy alejadas del espacio conocido y priorizarlas para validación?”**

Pero sigue teniendo huecos metodológicos grandes en **definición de tarea, negativos, redundancia, leakage, comparación de representaciones, incertidumbre y open-world discovery**. Ahí veo claramente espacio para nosotros.

---

# 1. Primero: no existe realmente un único problema “phage lytic enzyme prediction”

Hay que separar por lo menos tres familias funcionales.

### Endolysins

Son enzimas que degradan **peptidoglicano**, típicamente al final del ciclo lítico. Suelen presentar arquitecturas modulares con:

* enzymatically active domains, EAD;
* cell-wall binding domains, CBD;
* linkers;
* en algunos casos regiones de interacción con membrana.

La modularidad y *domain shuffling* son muy fuertes, por lo que identidad global de secuencia puede ser baja incluso conservando función. Una revisión crítica de 2024 mostró precisamente que conjuntos relativamente pequeños de HMMs de dominios bien curados pueden explicar buena parte de las endolisinas conocidas. ([ASM Journals][1])

### Virion-associated lysins, VALs/VAPGHs

También atacan peptidoglicano, pero durante **entrada del fago**, y frecuentemente forman parte o están asociadas a estructuras del virión.

Comparten algunos EADs con endolisinas, pero el contexto estructural y arquitectónico puede ser completamente diferente.

### Depolymerases

Acá hay que tener cuidado semánticamente.

Muchas depolimerasas no degradan peptidoglicano: degradan **cápsula, EPS, LPS u otros polisacáridos superficiales**. Frecuentemente son parte de tailspikes/receptor-binding proteins.

Por eso:

[
\text{endolysin/VAL biology}
\neq
\text{depolymerase biology}
]

aunque todas pertenezcan al universo funcional que estamos capturando en PhageLysData.

Esto importa mucho para el modelado: **una única multiclass plana puede estar imponiendo una ontología artificial al problema**.

---

# 2. Primera generación: homología, dominios y descriptores clásicos

Antes del boom de deep learning, el problema se abordaba fundamentalmente mediante:

[
sequence
\rightarrow
homology/domain/descriptor
\rightarrow
classification
]

## Homología y HMMs

El enfoque más biológicamente directo ha sido:

* BLAST/DIAMOND;
* HMMER;
* Pfam;
* InterPro;
* perfiles de dominios catalíticos;
* combinación de EAD/CBD;
* reglas basadas en SignalP/TMHMM, longitud y arquitectura.

Esto sigue siendo importante porque las lysins son proteínas modulares.

La revisión crítica de Bałdysz et al. mostró que **13 modelos HMM representativos** ya capturan una fracción importante de las lysins bien anotadas y comparó directamente esta estrategia con herramientas como phiBiScan. ([ASM Journals][1])

El problema es evidente:

> los métodos basados en homología son excelentes dentro del espacio conocido, pero deterioran cuando entramos en el *remote-homology / protein-dark-matter space*.

Y esto es particularmente grave en fagos: una revisión de 2025 estima que alrededor del **65% de las proteínas de fagos no puede recibir una anotación funcional confiable**. ([ASM Journals][2])

---

# 3. Los primeros modelos ML para lytic enzymes

## Lypred — 2016

Uno de los predictores históricos importantes.

Representaba proteínas mediante una versión de **pseudo-amino-acid composition, PseAAC**, incluyendo:

* g-gap dipeptides;
* propiedades fisicoquímicas;
* correlaciones de orden de secuencia.

Después:

[
PseAAC
\rightarrow
ANOVA\ feature\ selection
\rightarrow
SVM
]

La base se redujo mediante CD-HIT al **40% de identidad**. Los positivos incluían lyases bacterianas y de fagos; los negativos provenían también de fagos. Reportó accuracy ≈84.8% y AUROC ≈0.926 mediante jackknife validation. ([PubMed][3])

Es históricamente interesante porque ya reconocía dos cuestiones que siguen siendo relevantes:

1. redundancia;
2. representación de la secuencia.

Pero no era específicamente una formulación moderna de *phage lytic enzyme discovery*.

---

## CWLy-pred — 2020

Representaba inicialmente cada proteína con **473 features** y después realizaba una selección agresiva de variables mediante MRMD, quedándose con apenas seis.

El clasificador era un **SVM**.

Reportó:

* ACC 93.1%;
* sensitivity 85.3%;
* specificity 94.8%;
* MCC 0.775;
* AUC 0.900.

Los propios autores reconocieron como limitación **la ausencia de un independent test set**. ([ScienceDirect][4])

Esta generación tenía la filosofía:

> *better handcrafted representation + feature selection + conventional classifier.*

Eso prácticamente dejó de ser SOTA una vez aparecieron los PLMs.

---

# 4. PhaLP 1.0: aparece por primera vez una clasificación funcional muy parecida a la nuestra

PhaLP, en 2021, es un precedente muy importante.

Partió de 7,957 secuencias únicas y realizó una clasificación semimanual usando:

* UniProt annotations;
* GO;
* conserved domains;
* protein names.

Obtuvieron:

* 1,604 VALs;
* 2,492 endolysins.

Después:

[
sequence
\overset{SeqVec}{\longrightarrow}
1024D
\rightarrow
RandomForest
]

con 100 árboles y balanceo de clases.

En stratified 10-fold CV reportaron:

[
F1 \approx 0.984
]

([PubMed Central (PMC)][5])

### Lo importante metodológicamente

Esto representa el primer salto claro desde:

**handcrafted descriptors**

hacia:

**learned protein representations.**

Pero el paper no estaba diseñado para estudiar sistemáticamente:

* efecto de representación;
* efecto de homología;
* distintos esquemas de negativos;
* OOD generalization;
* uncertainty.

---

# 5. Depolymerases: un desarrollo computacional casi paralelo

Aquí la literatura es mucho más activa de lo que parece.

---

## DePP — 2023

Uno de los primeros predictores específicos de depolimerasas con ML.

Lo interesante es que los autores fueron deliberadamente conservadores:

**solamente ~50 depolimerasas experimentalmente verificadas** para los positivos.

Usaron features derivadas de la secuencia, incluyendo di/tripeptides, y compararon:

* SVM;
* Random Forest.

Hubo nested validation para hyperparameter selection. El modelo final fue RF con unas 1,500 trees y polynomial feature transformation.

Accuracy final alrededor de:

[
0.90
]

([Springer Link][6])

### Pero DePP introdujo algo que me gusta mucho

No evaluó solamente:

> DP / non-DP.

También preguntó:

> **¿en qué posición rankea la verdadera depolimerasa dentro del proteoma completo del fago?**

En 155 genomas, la verdadera DPO aparecía aproximadamente:

* top-1: 40.6%;
* top-3: ~70%;
* top-10: ~87%;
* top-20: ~95%.

([PubMed Central (PMC)][7])

**Eso es una evaluación orientada a discovery.**

Y creo que nosotros deberíamos recuperar esa idea.

---

# 6. DeepMineLys — 2024: cambia el estándar de “discovery”

Este paper es probablemente uno de los precedentes más importantes para nosotros.

DeepMineLys desarrolló una arquitectura **two-track CNN con una tarea de cuatro clases**, diseñada específicamente para encontrar lysins en microbiomas. ([ScienceDirect][8])

El mensaje metodológico central fue explícitamente:

> **interpolative performance is not enough; extrapolative ability matters.**

En un independent dataset:

[
F1 = 0.84
]

superando los métodos comparados por ~20.8 puntos porcentuales. ([ScienceDirect][8])

Luego aplicaron el modelo a **370,849 proteínas** de microbiomas humanos e identificaron aproximadamente:

* 18,514 endolysins;
* 20,340 VALs;
* 6,784 holins. ([PubMed Central (PMC)][9])

Pero lo que transforma el paper es la validación prospectiva:

> 16 candidatos expresados → **11 activos**.

Y el mejor mostró actividad 6.2× superior a hen egg-white lysozyme. ([ScienceDirect][8])

### Lección

DeepMineLys estableció que un buen paper de discovery no termina en AUROC.

Termina en:

[
model
\rightarrow
unseen\ sequence\ space
\rightarrow
candidate\ prioritization
\rightarrow
validation
]

---

# 7. DepoScope — 2024: PLM + residue-level prediction

Este es conceptualmente muy bueno.

DepoScope usa:

[
ESM2
+
CNN
]

pero hace simultáneamente:

### Sequence-level prediction

¿Es depolymerase o no?

### Amino-acid-level prediction

¿Dónde está el dominio funcional?

Es decir:

[
\text{protein classification}
+
\text{token classification}
]

([PubMed Central (PMC)][10])

La construcción del dataset fue bastante sofisticada:

* > 500k proteínas de INPHARED;
* base propia de dominios degradadores de polisacáridos;
* InterPro;
* estructura;
* ESMFold;
* Foldseek;
* folds tipo β-helix / β-propeller, etc.

Además evaluaron la generalización construyendo splits con CD-HIT a múltiples identidades:

[
25%,30%,35%,40%,45%,65%,70%,75%,80%,85%
]

y usaron Bayesian optimization/Optuna para tuning.

### Esta idea nos interesa mucho

No porque debamos copiar DepoScope, sino porque demuestra que:

> **sequence identity threshold puede convertirse en una variable experimental, no solamente en un preprocessing step.**

Eso coincide exactamente con lo que veníamos pensando para nuestro paper.

---

# 8. SpikeHunter — 2024: PLM y control de homología muy serio

SpikeHunter no predice depolimerasas directamente; detecta **tailspike proteins**, muchas de las cuales funcionan como receptor-binding/depolymerase proteins.

Arquitectura:

[
ESM2_{650M}
\rightarrow
1280D\ embedding
\rightarrow
3-layer\ MLP
]

Dataset:

* 1,912 tailspikes;
* 200,732 otras proteínas de fagos.

Pero lo excelente es el split:

> **ninguna proteína de un set podía tener >30% identity con alguna proteína de otro set.**

([PubMed Central (PMC)][11])

Obtuvo aproximadamente:

* F1 0.990;
* MCC 0.990;
* AUPRC 0.997.

Y después se aplicó masivamente a genomas bacterianos/profagos para estudiar especificidad de hospedero. ([PubMed Central (PMC)][11])

### Este paper es muy relevante para nuestro diseño

Porque demuestra que **30%-identity-disjoint evaluation ya existe como estándar fuerte en este subcampo**.

---

# 9. DeepLysin — 2024: descubrimiento como pipeline, no como simple classifier

DeepLysin aborda un problema ligeramente diferente:

> dada una proteína que parece una lysin, ¿es probablemente antibacteriana?

El pipeline tiene dos módulos.

### Candidate mining

Usa:

* prophage identification;
* gene calling;
* CD-HIT;
* HMMsearch;
* lytic-domain filtering;
* DeepTMHMM.

Después:

### Activity prediction

Representaciones handcrafted:

* AAC;
* BPNC;
* GAAC;
* CTD;
* AAE;
* AAI,

y un **stacking ensemble** compuesto por:

* ExtraTrees;
* ANN;
* KNN;
* XGBoost;
* Logistic Regression;

con otra LR como meta-classifier. ([Wiley Online Library][12])

Después minaron >4,000 posibles lysins de fagos de *S. aureus*. ([Wiley Online Library][12])

Experimentalmente:

[
17\ candidates
\rightarrow
7\ active
]

incluyendo LLysSA9, posteriormente probado in vitro e in vivo. ([PubMed][13])

### Limitación metodológica interesante

Este trabajo mezcla:

[
\text{domain-based candidate definition}
+
\text{activity prediction}.
]

Por lo tanto, no responde exactamente al problema:

> “Given arbitrary phage protein X, is X a previously unknown lytic enzyme?”

---

# 10. PhageDPO — 2025

PhageDPO vuelve a un enfoque más clásico:

[
features
\rightarrow
SVM
]

pero construye su positive universe usando secuencias vinculadas a **siete dominios relacionados con depolimerasas**, complementadas con DPOs validadas experimentalmente.

Reporta:

* accuracy 96%;
* recall 97%;
* precision 94%;
* F1 96%.

Además valida casos de literatura y datos nuevos. ([ScienceDirect][14])

Tiene una limitación conceptual importante para cualquier benchmark que hagamos:

> si el dataset positivo se construye en parte mediante dominios conocidos, el modelo puede estar evaluando qué tan bien recupera precisamente ese universo definido por dominios.

No significa que esté mal, pero sí que **la definición de ground truth condiciona la dificultad del problema**.

---

# 11. PDP-Miner — 2025: narrowing the search space

PDP-Miner adopta una estrategia diferente:

[
bacterial\ genome
\rightarrow
prophage
\rightarrow
tail\ proteins
\rightarrow
DePP
\rightarrow
domain\ validation
]

Es esencialmente un wrapper orientado a discovery que restringe primero el universo biológico.

En 1,294 genomas de *Pseudomonas* encontró diez candidatos de alta confianza. ([PubMed Central (PMC)][15])

Esto ilustra otro paradigma:

> **no intentar distinguir la proteína objetivo de todo UniProt; reducir previamente el espacio al contexto biológicamente plausible.**

Eso será importante cuando discutamos nuestros negativos.

---

# 12. Estructura: SPAED y Phold

Entre 2025 y 2026 aparece otro cambio fuerte.

## SPAED

SPAED está específicamente diseñado para segmentar dominios de endolisinas usando la **PAE matrix de AlphaFold**.

No es un predictor de lysin/non-lysin.

Hace:

[
AlphaFold\ PAE
\rightarrow
hierarchical\ clustering
\rightarrow
domain\ boundaries
]

En 376 endolisinas manualmente delimitadas obtiene:

* mean IoU 96%;
* domain-boundary score 89%;

frente a 94%/70% de Chainsaw. ([OUP Academic][16])

Esto evidencia que **la estructura aporta información ortogonal a la secuencia para estas proteínas modulares**.

---

## Phold — 2026

No es específico de lytic enzymes, pero es muy relevante.

Phold representa estructura mediante:

[
sequence
\overset{ProstT5}{\longrightarrow}
3Di
\overset{Foldseek}{\longrightarrow}
structural\ homologs
]

y busca contra >1.36 millones de estructuras de proteínas de fagos.

En su benchmark:

* Pharokka/MMseqs2 anotó ~34.7%;
* Pharokka/HMM ~37.7%;
* Phold/ProstT5 ~49.4%;
* Phold/structures ~51.5%.

Y también aumentó las anotaciones dentro de la categoría **lysis**. ([PubMed Central (PMC)][17])

Esto demuestra que para proteínas extremadamente divergentes:

[
sequence\ homology
<
structural\ homology
]

puede ser muy real.

---

# 13. Genome context también está emergiendo

Y esto puede ser muy relevante a futuro para nosotros.

## GOPhage

GOPhage combina:

[
ESM2\ protein\ embeddings
+
genomic\ context
+
Transformer
]

para predecir función de proteínas de fagos.

La idea es aprovechar que los genomas de fagos tienen **arquitectura modular y orden génico conservado**.

Superó métodos previos en:

* proteínas divergentes: +6.78%;
* funciones infrecuentes: +13.05%.

Además pudo identificar cientos de posibles holins sin homología convencional. ([PubMed Central (PMC)][18])

---

## Phynteny / Phynteny Transformer

El enfoque reciente combina:

* PLM embeddings;
* positional encoding;
* BiLSTM;
* Transformer;
* synteny/genome organization.

Fue entrenado en >280k genomas de fagos y predice las nueve grandes categorías PHROG, incluyendo **lysis**. ([BioRxiv][19])

No resuelve nuestro problema fino, pero demuestra algo importante:

> la secuencia de una proteína no es necesariamente toda la información disponible.

Podemos tener:

[
P(y|sequence, genome\ context)
]

en lugar de únicamente

[
P(y|sequence).
]

---

# 14. DposFinder — 2026: probablemente el SOTA actual de depolymerase prediction

Este es especialmente importante.

Publicado en **Genome Medicine en abril de 2026**.

Dataset inicial:

* 384 depolimerasas bien caracterizadas/experimentales;
* 3,373 negativos.

Los negativos incluían:

* lysins;
* bacteriocins;
* lysozymes;
* lipases;
* adhesins;
* otras proteínas de fagos.

Después eliminaron secuencias con ≥60% identidad mediante CD-HIT, quedando:

* 244 DPO;
* 817 non-DPO.

([Springer Link][20])

Arquitectura:

[
ESM2
\rightarrow
additional\ Transformer
\rightarrow
classification
]

y atención para interpretar regiones importantes.

Independent test:

[
61\ DP + 205\ non-DP
]

con:

* AUROC 0.991;
* AUPRC 0.955;
* sensitivity DPO ≈90.2%;
* specificity/non-DPO ≈98%.

([Springer Link][20])

Además hicieron algo metodológicamente muy bueno:

reentrenaron en **su propio dataset**:

* DePP;
* PhageDPO;
* SpikeHunter;
* DepoScope;

para una comparación controlada.

DposFinder obtuvo:

[
F1\approx0.952.
]

([Springer Link][20])

Y después viene lo realmente fuerte:

### prospective discovery

Secuenciaron fagos de *K. pneumoniae* y experimentalmente confirmaron **seis depolimerasas nuevas**, todas con <50% identidad respecto de las conocidas. ([Springer Link][20])

Además predicen el **capsular serotype target**.

Y minaron >440k genomas de fagos obteniendo una base con >100k candidatos.

Este es probablemente el competidor más serio del lado DPO.

---

# 15. SUBLYME / PhaLP 2.0 — 2026: el competidor directo de nuestra idea jerárquica

Este es el paper que debemos estudiar con más cuidado.

Publicado en **Database en 2026**.

SUBLYME tiene exactamente:

[
\boxed{\text{Stage 1: lysin vs non-lysin}}
]

seguido de:

[
\boxed{\text{Stage 2: endolysin vs VAL}}
]

([pmc.ncbi.nlm.nih.gov][21])

### Dataset Stage 1

* 15,387 lysins;
* 47,288 non-lysin phage proteins.

Los negativos vienen de INPHARED y excluyen anotaciones potencialmente vinculadas con funciones líticas.

Representation:

[
ProtT5
]

Model:

[
RBF-SVM
]

incluso con **default hyperparameters**.

### Lo bueno: homology-aware splitting

Clustering con MMseqs2:

[
30%\ identity,;70%\ coverage
]

y **un cluster completo pertenece solamente a train o test**.

Eso es excelente diseño.

Evaluación:

[
10\times10\ fold
]

para Stage 1.

Para endolysin/VAL:

* ~4.4k VAL;
* ~11k endolysin;

con 50×5-fold repeated CV.

Ambos modelos:

[
F1\approx97.8%
]

([PubMed Central (PMC)][21])

### Discovery

Lo aplicaron a EnVhog:

**46 millones de proteínas**.

Resultado:

[
743,000
]

nuevas putative lysins, llevando PhaLP 2.0 a ~754k entradas y aumentando el número de clusters conocidos aproximadamente 40×. ([pmc.ncbi.nlm.nih.gov][21])

---

# 16. Pero SUBLYME deja huecos clarísimos

Y acá aparece nuestra oportunidad.

### 1. Solo una representación

ProtT5.

No preguntan:

> ¿cuánto depende el resultado de la representación?

### 2. Solo un algoritmo

RBF-SVM.

No es malo; de hecho funciona excelente.

Pero no sabemos interacción:

[
representation\times learner
]

### 3. Un único threshold de homología

30%.

No estudian sistemáticamente:

[
100\rightarrow90\rightarrow70\rightarrow50\rightarrow30
]

ni la degradación continua de performance con distancia al training.

### 4. Negativos tratados como negativos reales

Una proteína de fago que no esté anotada como lysin se utiliza como non-lysin después de filtros.

Pero:

[
\text{absence of lysin annotation}
\neq
\text{evidence of non-lysin}
]

Este es exactamente el argumento para **PU learning**.

### 5. La segunda tarea es mutuamente exclusiva

La propia implementación de SUBLYME señala que:

[
P(endo)+P(VAL)=1
]

y siempre asigna una de ellas. ([GitHub][22])

Eso es problemático para nuestro universo porque sabemos que existen casos:

* múltiples anotaciones;
* arquitecturas mixtas;
* proteínas funcionalmente ambiguas;
* evidencia discordante.

### 6. No incluye depolymerases

SUBLYME/PhaLP están fundamentalmente centrados en peptidoglycan lysins:

* endolysin;
* VAL.

Nuestro universo es mucho más amplio.

### 7. No hay uncertainty formal

Hay `predict_proba` / confidence.

Pero:

[
classifier\ probability
\neq
epistemic\ uncertainty.
]

No hay, por ejemplo:

* ensemble uncertainty;
* conformal prediction;
* OOD score;
* distance-aware confidence;
* calibrated selective prediction.

### 8. No distingue annotation de discovery

Una predicción con 99% de confianza puede ser simplemente:

> “esto se parece mucho a lo que ya conocía”.

No necesariamente:

> “esto constituye una nueva familia funcional convincente”.

Ese distinction es clave.

---

# 17. LysinFusion — 2026, por ahora preprint

También hay que tenerlo en radar.

LysinFusion combina cuatro descriptores:

* CKSAAP;
* CTDD;
* APAAC;
* CTDC;

con:

[
CNN + Transformer.
]

([GitHub][23])

Algo interesante es cómo construyeron los negativos:

> proteínas no-líticas provenientes **de los mismos virus que los positivos**, excluyendo unknowns.

Eso es bastante más difícil que coger proteínas random de UniProt. ([ResearchGate][24])

En un independent benchmark de 148 proteínas experimentalmente validadas:

* ACC ≈0.811;
* AUC ≈0.892;
* MCC ≈0.623.

Reportan superar DeepMineLys y reducir false positives. ([Life Science Network][25])

Además hacen occlusion/LIME para interpretación.

Pero, **al 21 de agosto de 2026, sigue siendo preprint**, así que lo separaría del SOTA peer-reviewed. ([Life Science Network][25])

---

# 18. Entonces la evolución del campo se puede resumir así

| Era        | Representation                                | Model                | Problema principal          |
| ---------- | --------------------------------------------- | -------------------- | --------------------------- |
| ~2009–2018 | PseAAC / AAC / physicochemical                | Fisher/SVM           | lysin recognition           |
| 2019–2022  | richer sequence descriptors / SeqVec          | SVM/RF               | functional classification   |
| 2023       | handcrafted + optimized ML                    | RF/SVM               | specialized DPO prediction  |
| 2024       | CNN / ESM2                                    | CNN/MLP              | extrapolative discovery     |
| 2024–25    | PLM + domain/structure                        | CNN/Transformer      | remote functional detection |
| 2025–26    | PLM + structure/context                       | Transformer/Foldseek | viral dark matter           |
| 2026       | PLM hierarchical / interpretable transformers | SVM/Transformer      | large-scale discovery       |

Y la tendencia es claramente:

[
\boxed{
homology
\rightarrow
descriptors
\rightarrow
embeddings
\rightarrow
PLMs
\rightarrow
structure/context
\rightarrow
discovery
}
]

---

# 19. Ahora viene lo más importante: ¿qué NO está resuelto?

Después de revisar todo esto, para mí hay **ocho vacíos bastante claros**.

## A. Task definition

Cada paper parece definir una tarea distinta:

* “cell wall lyase”;
* lysin;
* endolysin;
* VAL;
* DPO;
* tailspike;
* lysis protein;
* antibacterial lysin.

Por tanto, comparar F1 entre papers casi no tiene sentido.

No existe todavía una formulación realmente fuerte y explícita de:

[
\text{phage-associated envelope-degrading enzyme}
\rightarrow
{\text{endolysin, VAL, DPO, ...}}
]

con evidencia y relaciones no mutuamente excluyentes.

**PhageLysData nos da exactamente eso.**

---

# 20. B. Negative definition

Este me parece uno de los huecos más importantes.

Se han utilizado:

### Easy negatives

proteínas generales de UniProt.

### Phage negatives

otras proteínas de fagos.

### Hard negatives

otras proteínas del **mismo genoma/fago**.

### Functionally close negatives

lysozymes, adhesins, bacteriocins, tail proteins, etc.

Pero casi nadie pregunta sistemáticamente:

[
\text{How does negative definition change the model?}
]

Y tampoco encontré en esta búsqueda dirigida un predictor específico de lysins/depolimerasas que formule explícitamente la tarea mediante **positive-unlabeled learning**.

Ahí veo un hueco real.

---

# 21. C. Redundancy y split se suelen confundir

Esto es importantísimo.

Dos operaciones diferentes:

### Dataset redundancy control

[
D \rightarrow NR90,\ NR70,\ NR50...
]

### Evaluation partitioning

[
clusters
\rightarrow
train/val/test
]

Muchos trabajos hacen una y asumen que resolvieron la otra.

Los mejores trabajos actuales —SpikeHunter y SUBLYME— sí hacen separación fuerte por homología. ([PubMed Central (PMC)][11])

Pero no existe una caracterización sistemática de:

[
Performance =
f(
redundancy,\ split,\ representation
)
]

Y eso es exactamente nuestra idea.

---

# 22. D. Los scores publicados son casi incomparables

Tenemos:

* 98% F1 SUBLYME;
* 99% SpikeHunter;
* 95% DposFinder;
* 84% DeepMineLys;
* 81% LysinFusion;
* 90% DePP.

Eso **no significa que SUBLYME sea mejor que DeepMineLys**.

Porque cambian simultáneamente:

* labels;
* negatives;
* sequence space;
* identity thresholds;
* prevalence;
* split;
* evaluation dataset;
* biological question.

Esto nos permite construir un mensaje muy potente:

> **model performance cannot be interpreted independently of the experimental definition of the prediction task.**

---

# 23. E. Representation benchmarking casi no existe

Los papers suelen escoger:

* PseAAC;
* SeqVec;
* ProtT5;
* ESM2;

y siguen adelante.

Pero no preguntan sistemáticamente:

[
AAC
\quad vs\quad
descriptors
\quad vs\quad
ESM2
\quad vs\quad
ESMC
\quad vs\quad
ProtT5
\quad vs\quad
structure
]

bajo **exactamente los mismos splits**.

Ese experimento sigue siendo valioso.

Y nosotros tenemos una ventaja enorme porque ya hemos construido el recurso precisamente con múltiples capas de representación.

---

# 24. F. Sequence novelty no está cuantificada adecuadamente

Los buenos trabajos intentan demostrar remote discovery.

Pero faltaría una evaluación del tipo:

[
P =
f\left(
\max_{t\in Train}
Identity(x,t)
\right)
]

o:

[
P =
f(
embedding\ distance
)
]

Es decir:

> performance vs. distance from training knowledge.

Eso sería muchísimo más informativo que un único F1.

---

# 25. G. Incertidumbre prácticamente no está explotada

En esta revisión dirigida **no encontré un predictor dedicado de phage lytic enzymes que use formalmente**:

* conformal prediction;
* uncertainty decomposition;
* calibrated abstention;
* deep ensembles for epistemic uncertainty;
* OOD-aware selective prediction;

como pieza central del discovery.

Hay scores y probabilities, sí.

Pero no:

[
\text{I predict class X}
+
\text{I know how little I know}.
]

Para nuestro paper esto sigue siendo un hueco excelente.

---

# 26. H. Open-world classification

Todos los clasificadores especializados esencialmente asumen:

[
x\in{known\ classes}.
]

SUBLYME incluso fuerza:

[
Endolysin \oplus VAL.
]

Pero en el mundo real tendremos:

* proteins outside known labels;
* multiple functional classes;
* ambiguous annotations;
* new domain combinations;
* potentially new lytic mechanisms.

Por eso me atrae mucho para nosotros algo como:

[
\text{known}
/
\text{uncertain}
/
\text{novel}
]

en vez de solamente:

[
argmax(P_k).
]

---

# 27. La literatura general de proteínas de fagos nos da otra pista

GOPhage, Phynteny y Phold están convergiendo hacia tres fuentes complementarias:

[
\boxed{
Sequence +
Structure +
Genome\ context
}
]

([PubMed Central (PMC)][17])

Yo **no metería todo eso en nuestro primer paper**, porque se volvería monstruoso.

Pero sí dejaría claro que el problema puede verse como:

[
f(sequence)
]

para esta primera versión y posteriormente:

[
f(sequence,\ structure,\ context).
]

---

# 28. Mapa competitivo directo para nuestro paper

Si reduzco todo el estado del arte a quienes realmente competirían contra nosotros:

| Tool             | Target           | Representation        | Model           | Homology-aware               | Discovery           | Experimental     |
| ---------------- | ---------------- | --------------------- | --------------- | ---------------------------- | ------------------- | ---------------- |
| **PhaLP 1**      | Endo vs VAL      | SeqVec                | RF              | limitado                     | database annotation | No               |
| **DePP**         | DPO              | handcrafted           | RF              | parcialmente                 | genome ranking      | No               |
| **DeepMineLys**  | Lysins           | DL features           | CNN             | extrapolation-oriented       | metagenome          | **Yes**          |
| **DepoScope**    | DPO + domain     | ESM2                  | CNN             | **Yes, multiple thresholds** | large-scale         | No               |
| **SpikeHunter**  | Tailspike        | ESM2                  | MLP             | **30% disjoint**             | genomic mining      | No               |
| **DeepLysin**    | active lysins    | handcrafted           | stacking        | moderate                     | candidate mining    | **Yes**          |
| **PhageDPO**     | DPO              | engineered            | SVM             | —                            | genome prediction   | validation       |
| **DposFinder**   | DPO + serotype   | ESM2                  | Transformer     | NR60 + independent           | **Yes**             | **Yes**          |
| **SUBLYME**      | lysin → Endo/VAL | ProtT5                | SVM             | **30% cluster split**        | **huge-scale**      | domain-based     |
| **LysinFusion*** | lysin            | handcrafted multiview | CNN-Transformer | curated                      | prediction          | independent test |

*preprint as of August 2026. ([PubMed Central (PMC)][5])

---

# 29. Y entonces reformulo bastante nuestra idea inicial

Después de este estado del arte, **yo ya no vendería nuestra contribución como un nuevo predictor jerárquico**.

SUBLYME nos mata esa novedad directamente.

Tampoco vendería:

> *“first PLM predictor for phage lytic enzymes”*

porque evidentemente ya no lo somos.

Ni:

> *“first uncertainty-assisted biological discovery”*

sin precisar muchísimo.

## La oportunidad real es otra

Para mí, el paper debería preguntar:

> **What makes computational prediction of phage lytic enzymes reliable beyond the known sequence space?**

Y PhageLysData nos permite responderlo porque podemos manipular experimentalmente:

[
\boxed{
Task\ definition
}
]

[
\times
\boxed{
Negative\ evidence
}
]

[
\times
\boxed{
Redundancy
}
]

[
\times
\boxed{
Split\ strategy
}
]

[
\times
\boxed{
Representation
}
]

[
\rightarrow
\boxed{
Generalization
}
]

y después:

[
\boxed{
Uncertainty + Novelty
\rightarrow
Discovery
}
]

**Eso no lo veo resuelto por ninguno de los trabajos actuales.**

---

# 30. Y hay una oportunidad especialmente bonita con PhageLysData

La mayoría de estos papers tuvo que construir un dataset específicamente para su predictor.

Nosotros hicimos el camino inverso:

> **primero construimos un recurso evidence-aware y armonizado; después utilizamos ese recurso para estudiar cómo la definición de los datos determina la capacidad predictiva.**

Esa diferencia es enorme conceptualmente.

No sería:

> “here is another classifier.”

Sería:

> **“Here is what happens to apparent and real discovery performance when we systematically change the evidence, redundancy, representation and generalization regime.”**

Eso sí lo veo diferenciable incluso frente a **SUBLYME y DposFinder**.

---

## Mi lectura final del estado del arte

Hay tres generaciones actuales coexistiendo:

**annotation methods**

[
BLAST/HMM/domain/structure
]

**classification methods**

[
features/PLMs
\rightarrow
SVM/RF/CNN/Transformer
]

y

**discovery systems**

[
classification
+
remote\ search
+
candidate\ ranking
+
experimental/orthogonal\ validation.
]

Lo que prácticamente todavía no veo es una cuarta:

[
\boxed{
\textbf{evidence-aware + generalization-aware + uncertainty-aware discovery}
}
]

con un análisis sistemático del efecto de **cómo construyes la tarea**.

Y creo, papu, que **ese debería ser nuestro espacio**, no competir por meter un Transformer más grande que DposFinder ni una arquitectura más fancy que DeepMineLys. Tenemos algo bastante más defendible: un recurso que permite experimentar con las decisiones que todos los demás papers normalmente fijan una sola vez y después dan por correctas.

Si quieres, el siguiente paso puede ser **convertir este estado del arte en una matriz exacta “qué hicieron / qué no hicieron / qué haremos nosotros”**, y de ahí congelamos las **research questions e hipótesis del paper** antes de tocar un solo modelo.

[1]: https://journals.asm.org/doi/10.1128/aem.02361-23?utm_source=chatgpt.com "“Tear down that wall”—a critical evaluation of bioinformatic resources available for lysin researchers | Applied and Environmental Microbiology"
[2]: https://journals.asm.org/doi/10.1128/mmbr.00022-25?utm_source=chatgpt.com "Computational function prediction of bacteria and phage proteins | Microbiology and Molecular Biology Reviews"
[3]: https://pubmed.ncbi.nlm.nih.gov/27437396/?utm_source=chatgpt.com "Identification of Bacterial Cell Wall Lyases via Pseudo Amino Acid Composition - PubMed"
[4]: https://www.sciencedirect.com/science/article/pii/S0888754320309915?utm_source=chatgpt.com "CWLy-pred: A novel cell wall lytic enzyme identifier based on an improved MRMD feature selection method - ScienceDirect"
[5]: https://pmc.ncbi.nlm.nih.gov/articles/PMC8310338/?utm_source=chatgpt.com "PhaLP: A Database for the Study of Phage Lytic Proteins and Their Evolution - PMC"
[6]: https://link.springer.com/article/10.1186/s12859-023-05341-w?utm_source=chatgpt.com "DePolymerase Predictor (DePP): a machine learning tool for the targeted identification of phage depolymerases | BMC Bioinformatics | Springer Nature Link"
[7]: https://pmc.ncbi.nlm.nih.gov/articles/PMC10199479/?utm_source=chatgpt.com "DePolymerase Predictor (DePP): a machine learning tool for the targeted identification of phage depolymerases - PMC"
[8]: https://www.sciencedirect.com/science/article/pii/S2211124724009124?utm_source=chatgpt.com "DeepMineLys: Deep mining of phage lysins from human microbiome - ScienceDirect"
[9]: https://pmc.ncbi.nlm.nih.gov/articles/PMC12031009/?utm_source=chatgpt.com "Phage Endolysins as an Alternative Biocontrol Strategy for Pathogenic and Spoilage Microorganisms in the Food Industry - PMC"
[10]: https://pmc.ncbi.nlm.nih.gov/articles/PMC11326577/?utm_source=chatgpt.com "DepoScope: Accurate phage depolymerase annotation and domain delineation using large language models - PMC"
[11]: https://pmc.ncbi.nlm.nih.gov/articles/PMC11034027/?utm_source=chatgpt.com "Large-scale genomic survey with deep learning-based method reveals strain-level phage specificity determinants - PMC"
[12]: https://advanced.onlinelibrary.wiley.com/doi/10.1002/advs.202404049?utm_source=chatgpt.com "Discovery of Antimicrobial Lysins from the “Dark Matter” of Uncharacterized Phages Using Artificial Intelligence - Zhang - 2024 - Advanced Science - Wiley Online Library"
[13]: https://pubmed.ncbi.nlm.nih.gov/38899839/?utm_source=chatgpt.com "Discovery of Antimicrobial Lysins from the \"Dark Matter\" of Uncharacterized Phages Using Artificial Intelligence - PubMed"
[14]: https://www.sciencedirect.com/science/article/pii/S0010482525001866?utm_source=chatgpt.com "PhageDPO: A machine-learning based computational framework for identifying phage depolymerases - ScienceDirect"
[15]: https://pmc.ncbi.nlm.nih.gov/articles/PMC12579547/?utm_source=chatgpt.com "PDP-Miner: an AI/ML tool to detect prophage tail proteins with depolymerase domains across thousands of bacterial genomes - PMC"
[16]: https://academic.oup.com/bioinformatics/article/41/10/btaf531/8262958?utm_source=chatgpt.com "SPAED: harnessing AlphaFold output for accurate segmentation of phage endolysin domains | Bioinformatics | Oxford Academic"
[17]: https://pmc.ncbi.nlm.nih.gov/articles/PMC12774648/?utm_source=chatgpt.com "Protein structure-informed bacteriophage genome annotation with Phold - PMC"
[18]: https://pmc.ncbi.nlm.nih.gov/articles/PMC11751364/?utm_source=chatgpt.com "GOPhage: protein function annotation for bacteriophages by integrating the genomic context - PMC"
[19]: https://www.biorxiv.org/content/10.1101/2025.07.28.667340v2?utm_source=chatgpt.com "Synteny-aware functional annotation of bacteriophage genomes with Phynteny | bioRxiv"
[20]: https://link.springer.com/article/10.1186/s13073-026-01657-3?utm_source=chatgpt.com "DposFinder: an interpretable transformer model for predicting phage-derived polysaccharide depolymerases and their host capsular serotypes | Genome Medicine | Springer Nature Link"
[21]: https://pmc.ncbi.nlm.nih.gov/articles/PMC13237592/?utm_source=chatgpt.com "PhaLP 2.0: extending the community-oriented phage lysin database with a SUBLYME pipeline for metagenomic discovery - PMC"
[22]: https://github.com/Rousseau-Team/sublyme?utm_source=chatgpt.com "GitHub - Rousseau-Team/sublyme: Prediction of bacteriophage lysins in metagenomic datasets. · GitHub"
[23]: https://github.com/sinuo560/LysinFusion?utm_source=chatgpt.com "GitHub - sinuo560/LysinFusion: A hybrid CNN–Transformer framework for lysin prediction · GitHub"
[24]: https://www.researchgate.net/publication/402525053_LysinFusion_Integrating_Multi-Feature_Encoding_and_Hybrid_CNN-Transformer_Architecture_for_Phage_Lysin_Prediction?utm_source=chatgpt.com "(PDF) LysinFusion: Integrating Multi-Feature Encoding and Hybrid CNN-Transformer Architecture for Phage Lysin Prediction"
[25]: https://www.lifescience.net/preprints/11609/lysinfusion-integrating-multi-feature-encoding-and/?utm_source=chatgpt.com "LysinFusion: Integrating Multi-Feature Encoding and Hybrid CNN-Transformer Architecture for Phage Lysin Prediction"


---
Sí, papu. **Ahora que veo la línea completa de cinco papers, yo dejaría de escoger journals paper por paper y construiría una escalera editorial.** Hay contribuciones que son naturalmente más metodológicas, otras más de protein science y otras con potencial de salto grande si aparece validación experimental.

Y en ese esquema **CSBJ prácticamente sale de mi estrategia**: full OA caro y no nos da una ventaja suficiente frente a alternativas híbridas/subscription.

## Cómo repartiría los cinco trabajos

| #     | Paper          | Contribución principal                               | Target que reservaría           |
| ----- | -------------- | ---------------------------------------------------- | ------------------------------- |
| **1** | PhageLysData   | recurso / evidence-aware dataset                     | **Frontiers in Bioinformatics** |
| **2** | Classification | detección + clasificación funcional                  | **PROTEINS / Protein Science**  |
| **3** | Patterns       | reglas, diferencias, determinantes de función/diseño | **Protein Science / JCIM**      |
| **4** | Mining         | ranking + novelty + uncertainty + discovery          | **Briefings in Bioinformatics** |
| **5** | De novo        | diseño + validación experimental                     | **ACS Catalysis** / ACS SynBio  |

Y acá hay una decisión estratégica interesante entre #2 y #3.

---

# 2. Clasificación: aquí NO gastaría Briefings

Totalmente de acuerdo contigo.

Después de ver SUBLYME, DposFinder, DeepMineLys, etc., el paper de clasificación puede ser **muy bueno**, pero sigue siendo:

> construcción rigurosa de modelos de clasificación para phage lytic enzymes.

Yo ya no intentaría venderlo como nuestro gran paper metodológico.

Haría:

* task definition;
* negative definition;
* representación;
* redundancia;
* random vs homology-aware split;
* unos pocos modelos;
* benchmark;
* quizás PU como análisis secundario;
* interpretación sencilla;
* modelo final.

Y se acabó.

### Para este paper ahora me gusta muchísimo **PROTEINS**

Porque tiene literalmente en scope:

> structure, function, **computation**, genetics and design

y acepta explícitamente nuevos approaches computacionales y hasta **prediction reports**. ([Wiley Online Library][1])

Además:

* IF 3.3;
* **8 días** medianos a primera decisión;
* no hay cargos obligatorios;
* OA es opcional. ([Wiley Online Library][2])

Esto es casi exactamente lo que buscamos para el #2:

> **rápido + respectable + protein-oriented + no APC duro + no necesidad de inventar una red neuronal.**

### ¿Y Protein Science?

**También entra perfectamente.**

Protein Science tiene IF 5.6, acceptance rate publicada de 53% y 33 días medianos a primera decisión. Además acepta explícitamente estudios computacionales sobre función, estructura, diseño y evolución proteica. ([Wiley Online Library][3])

Y la opción tradicional cuesta actualmente **$750 para non-members**, no miles de euros de APC. ([Wiley Online Library][4])

Así que si me dices:

> *“me tinca mucho Protein Science y quiero probar ahí el classifier”*

te diría **sí, absolutamente**.

Pero hay una razón por la que estratégicamente quizá la guardaría.

---

# 3. El paper de patrones puede ser incluso MÁS Protein Science

Este:

> **patrones y diferencias entre phage lytic enzymes + reconocimiento de características relevantes para diseño y clasificación**

es extremadamente *Protein Science*.

Porque ya no estamos preguntando:

> ¿qué modelo predice mejor?

Estamos preguntando:

> **¿qué distingue funcional, fisicoquímica, secuencial y estructuralmente a estas familias proteicas?**

Podemos integrar:

* sequence patterns;
* domain architectures;
* physicochemical properties;
* PLM space;
* structural space;
* catalytic/domain organization;
* EAD/CBD organization;
* endolysin vs VAL vs depolymerase;
* conserved vs divergent regions;
* family-specific signatures;
* important positions/regions;
* relationships between pattern and function;
* potentially design constraints.

Eso toca directamente **protein function, structure, design and evolution**, que son precisamente las áreas que Protein Science declara explícitamente. ([Wiley Online Library][5])

Por eso, **si optimizo el portafolio**, haría:

### Paper 2 → **PROTEINS**

### Paper 3 → **Protein Science**

Y creo que los dos quedarían perfectamente diferenciados.

---

# 4. Aquí sí reservaría Briefings in Bioinformatics

Este es el trabajo que yo no “quemaría” antes.

Porque el paper #4 puede ser bastante más interesante metodológicamente:

[
\text{PhageLysData}
\rightarrow
\text{patterns}
\rightarrow
\text{classifier}
\rightarrow
\text{candidate universe}
]

y entonces:

[
\text{activity likelihood}
+
\text{novelty}
+
\text{uncertainty}
+
\text{pattern consistency}
+
\text{distance from known space}
]

[
\downarrow
]

[
\boxed{\text{candidate ranking}}
]

Aquí sí aparece algo diferente de SUBLYME/DposFinder.

No queremos simplemente:

> probability > 0.95.

Queremos un **discovery framework**.

Por ejemplo:

[
Score(x)=
w_1P_{\mathrm{lytic}}
+
w_2P_{\mathrm{class}}
+
w_3Novelty
+
w_4PatternSupport
-----------------

w_5Uncertainty
]

o una formulación no necesariamente lineal, pero con esa filosofía.

Y podemos estudiar:

* ranking stability;
* calibration;
* uncertainty;
* conformal/selective prediction;
* OOD;
* sequence novelty;
* structural novelty;
* candidate diversity;
* prioritization.

Esto sí es perfectamente defendible como **Problem Solving Protocol / Case Study / methodological application** en BiB, categorías que actualmente aceptan nuevas aplicaciones metodológicas sobre problemas biológicos reales. ([OUP Academic][6])

### Acá sí gastaría BiB.

No en el classifier.

---

# ¿Y ACS Catalysis para el paper 4?

Aquí pondría una condición fuerte.

### Si el #4 es solamente computacional:

**No mandaría a ACS Catalysis.**

ACS Catalysis dice explícitamente que en computational enzymology los resultados computacionales deberían generalmente ir respaldados por datos experimentales. Además, para enzyme discovery/engineering quieren evidencia de **catalytic performance**, incluyendo aspectos como eficiencia, selectividad, estabilidad y especificidad de sustrato. ([American Chemical Society Publications][7])

Entonces:

> ranking → 20 candidatos computacionales

no me parece ACS Catalysis.

Pero:

> ranking → candidatos novedosos → expresión → actividad → cinética/caracterización

**ahí cambia completamente.**

---

# 5. El paper de diseño de novo es donde realmente veo ACS Catalysis

Este es potencialmente **el paper grande de la serie**.

Imagínate:

[
PhageLysData
\rightarrow
learned\ design\ rules
]

[
\downarrow
]

[
de\ novo\ generation/design
]

[
\downarrow
]

[
ranking + developability + structural filtering
]

[
\downarrow
]

[
experimental\ expression
]

[
\downarrow
]

[
lysis / catalytic\ activity
]

Si tenemos diseños **verdaderamente novedosos**, y experimentalmente caracterizamos:

* enzymatic activity;
* substrate specificity;
* catalytic efficiency si podemos;
* stability;
* dose-response;
* bacteriolytic activity;
* Gram+/Gram− context cuando corresponda;
* quizá comparación con naturales;

entonces **ACS Catalysis se vuelve un target completamente legítimo**.

ACS Catalysis tiene actualmente IF **13.6** y explícitamente busca descubrimiento, caracterización e ingeniería de enzimas y nuevas metodologías computacionales para generar enzimas con propiedades nuevas o mejoradas. ([American Chemical Society Publications][8])

Y están publicando ML para problemas enzimáticos en 2026; por ejemplo, un trabajo reciente de predicción de parámetros cinéticos mediante ML fue recibido el 19 de mayo, revisado el 10 de junio, aceptado el 15 de junio y publicado online el 23 de junio de 2026. ([American Chemical Society Publications][9])

Así que tampoco es necesariamente una revista insoportablemente lenta cuando el paper funciona.

---

# Si el diseño #5 queda más “engineering” que “catalysis”

Entonces tenemos **ACS Synthetic Biology**.

Su scope incluye explícitamente:

* design;
* engineering;
* optimization of biological systems;
* computational methods for biological design. 

IF actual 4.5. ([American Chemical Society Publications][10])

Pero si tenemos **enzimas de novo catalíticamente activas y bien caracterizadas**, yo intentaría ACS Catalysis antes.

---

# Y JCIM, ¿dónde queda?

Para mí es nuestro **comodín premium**.

JCIM acepta explícitamente:

* ML sobre datos biológicos;
* pattern recognition;
* molecular modeling;
* computer-aided molecular design;
* biopharmaceutical applications. ([American Chemical Society Publications][11])

Entonces podría recibir:

### Paper #3

si lo hacemos muy computational:

> *sequence/structure patterns governing functional specialization of phage lytic enzymes*

con ML interpretability + structural/representation analysis.

O:

### Paper #4

si BiB dice que no:

> *uncertainty- and novelty-aware computational mining of phage lytic enzymes.*

Yo **no lo gastaría en el classifier #2** a menos que el paper quede inesperadamente metodológico.

---

# Por eso mi mapa editorial quedaría así

### **Paper 1 — Resource**

**PhageLysData**

→ **Frontiers in Bioinformatics**

Ya está.

---

### **Paper 2 — Prediction**

**Reliable classification of phage lytic enzymes**

Target:

**PROTEINS → BMC Bioinformatics / similar**

o, si quieres intentar algo más alto:

**Protein Science → PROTEINS.**

Este debe salir **rápido**.

---

### **Paper 3 — Understanding**

**Sequence, structural and physicochemical determinants of phage lytic enzyme functional diversity**

Target:

### **Protein Science**

Backup:

**JCIM → PROTEINS**

Para mí este es *mucho más Protein Science* que el #2.

---

### **Paper 4 — Discovery**

**Evidence-, novelty- and uncertainty-aware mining of phage lytic enzymes**

Target:

### **Briefings in Bioinformatics**

Backup:

**JCIM**

Aquí ponemos la metodología más interesante.

---

### **Paper 5 — Design**

**De novo design of functional phage lytic enzymes**

Si hay validación experimental fuerte:

### **ACS Catalysis**

y dependiendo de los resultados:

**ACS Synthetic Biology / JMB**

JMB explícitamente favorece trabajos multidisciplinarios que combinen computación y experimento para obtener insight mecanístico/funcional. ([ScienceDirect][12])

---

# Y hay algo todavía más importante: evitar que parezca salami slicing

Los cinco papers tienen que responder **preguntas distintas**, no ser “PhageLysData parts 1–5”.

Yo congelaría desde ya esta separación:

**1 — What data exist?**

→ resource.

**2 — Can we reliably recognize them?**

→ classification.

**3 — What makes them different?**

→ patterns / determinants.

**4 — Where are the unknown ones?**

→ mining / discovery.

**5 — Can we create new ones?**

→ design / experiment.

Eso es una progresión preciosa:

[
\boxed{
DATA
\rightarrow
PREDICT
\rightarrow
UNDERSTAND
\rightarrow
DISCOVER
\rightarrow
DESIGN
}
]

Y además cada paper **produce insumos para el siguiente sin ser dependiente de él**.

Eso, para mí, es la línea PhageLysData/PhageLytic completa.

### Entonces, para el paper que estamos diseñando ahora mismo —el #2—

Yo ya **no perseguiría BiB ni ACS Catalysis**.

Mi decisión práctica sería:

[
\boxed{\textbf{Protein Science o PROTEINS}}
]

Si quieres **mejor journal y te tinca**: **Protein Science**.

Si quieres **máxima velocidad, cero APC obligatorio y sacar la wea del escritorio**: **PROTEINS**, que actualmente reporta 8 días medianos a primera decisión. ([Wiley Online Library][2])

Y personalmente, viendo el portafolio completo, **me inclino ligeramente por PROTEINS para #2 y guardaría Protein Science para #3**, porque ahí el encaje científico es casi perfecto.

[1]: https://onlinelibrary.wiley.com/page/journal/10970134/homepage/forauthors.html?utm_source=chatgpt.com "Author Guidelines - Proteins: Structure, Function, and Bioinformatics - Wiley Online Library"
[2]: https://onlinelibrary.wiley.com/journal/10970134?utm_source=chatgpt.com "PROTEINS: Structure, Function, and Bioinformatics | Protein Science Journal | Wiley Online Journal"
[3]: https://onlinelibrary.wiley.com/journal/1469896x?utm_source=chatgpt.com "Protein Science - Wiley Online Library"
[4]: https://onlinelibrary.wiley.com/page/journal/1469896x/homepage/forauthors.html?utm_source=chatgpt.com "Protein Science"
[5]: https://onlinelibrary.wiley.com/journal/1469896X?utm_source=chatgpt.com "Protein Science - Wiley Online Library"
[6]: https://academic.oup.com/bib/pages/author-guidelines?utm_source=chatgpt.com "Author guidelines | Briefings in Bioinformatics | Oxford Academic"
[7]: https://pubs.acs.org/doi/10.1021/acscatal.6b01749?utm_source=chatgpt.com "ACS Catalysis and the Scope of Papers Sought in Three Catalysis Subdisciplines: Biocatalysis and Enzymology, Molecular Catalysis for Organic Synthesis, and Heterogeneous Photocatalysis | ACS Catalysis | ACS Publications"
[8]: https://pubs.acs.org/journal/accacs?utm_source=chatgpt.com "ACS Catalysis | ACS Publications"
[9]: https://pubs.acs.org/doi/10.1021/acscatal.6c03874?utm_source=chatgpt.com "Enzyme Kinetic Parameter Prediction via Catalytic Pocket-Augmented Machine Learning | ACS Catalysis | ACS Publications"
[10]: https://pubs.acs.org/journal/asbcd6?utm_source=chatgpt.com "ACS Synthetic Biology | ACS Publications"
[11]: https://pubs.acs.org/jcisd8/pages/About?utm_source=chatgpt.com "About | Journal of Chemical Information and Modeling | ACS Publications"
[12]: https://www.sciencedirect.com/journal/journal-of-molecular-biology/publish/guide-for-authors?utm_source=chatgpt.com "Guide for authors - Journal of Molecular Biology - ISSN 0022-2836 | ScienceDirect.com by Elsevier"

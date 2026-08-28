Sí bebé. Ahora sí conviene **congelar el roadmap científico antes de tocar un notebook**, porque si diseñamos bien esta capa, el paper deja de ser “otra clasificación de proteínas de fagos” y pasa a ser un estudio de **generalización en protein-function prediction bajo cambios de evidencia, redundancia, homología, representación y label uncertainty**.

La regla que usaría para todo el roadmap es:

> **Exploramos exhaustivamente el tensor completo, pero el paper reporta los principios generales que emergen de él.**

El objetivo para *Briefings in Bioinformatics* no sería encontrar solamente el mejor modelo. Sería responder **cuándo, por qué y bajo qué condiciones la predicción funcional generaliza**.

---

# Roadmap analítico completo — estrategia *Briefings*

## 0. Universo experimental congelado

Tenemos dos espacios que deben mantenerse conceptualmente separados hasta que hagamos comparaciones válidas.

### Supervised

$$
7\ \text{tasks}\times
16\ (R,S)\times
4\ \text{representations}\times
5\ \text{algorithms}
=
2240
$$

con:

$$
R,S\in\{30,50,70,90\}.
$$

Cada `task × representation × algorithm` contiene el factorial completo:

$$
4R\times4S=16\text{ celdas}.
$$

Por superficie disponemos de:

$$
{16\choose2}=120
$$

comparaciones entre celdas, **24 contrastes de R a S fijo**, **24 de S a R fijo** y **36 interacciones R×S**.

### PU

$$
2\ \text{U definitions}\times
16\ (R,S)\times
4\ \text{representations}\times
3\ \text{learning paradigms}\times
2\ \text{base learners}
=
768.
$$

Paradigmas:

```text
naive P-vs-U
BaggingPU
Elkan–Noto
```

Aquí `naive P-vs-U` es precisamente nuestro **control tradicional**, porque convierte operacionalmente U en clase 0 mientras los otros métodos reconocen explícitamente label uncertainty.

---

# RQ1 — ¿Cómo determina la definición del problema lo que parece “aprendible”?

Esta debería abrir Results.

Queremos demostrar que **la dificultad aparente no es una propiedad intrínseca de “la función”**, sino también de cómo definimos evidencia positiva, negativa y desconocida.

| Comparación                   | Qué responde                                  |
| ----------------------------- | --------------------------------------------- |
| `global_easy vs global_hard`  | Efecto de la naturaleza del negativo          |
| `global_easy vs global_mixed` | Efecto de ampliar el background               |
| `global_hard vs global_mixed` | Efecto de combinar backgrounds                |
| subtype binary vs multiclass  | Efecto de cambiar la formulación del problema |
| `PU conservative vs broad`    | Efecto de ampliar el universo desconocido     |

### Análisis

Para `global_*`:

$$
\Delta M=
M_{\mathrm{hard}}-M_{\mathrm{easy}}
$$

etc., para **cada**:

```text
R × S × representation × algorithm
```

No solo medias globales.

Debemos mirar:

```text
MCC
AP
AP gain over prevalence
recall
specificity
```

porque ya sabemos que las prevalencias son diferentes.

### Análisis especialmente fuerte: backgrounds anidados

`global_mixed` contiene componentes de easy/hard, así que podemos hacer comparaciones sobre entidades compartidas.

Por ejemplo:

```text
easy vs mixed
→ common positives + shared N_easy
```

y:

```text
hard vs mixed
→ common positives + shared N_hard
```

Eso nos permite separar mejor:

> **task-definition effect**

de:

> **test-composition effect**.

### Binary vs multiclass

Para:

```text
endolysin
depolymerase
VAL
```

compararemos:

```text
subtype binary
vs
one-vs-rest behavior inside multiclass
```

y, mejor aún, la recuperación de **los mismos positivos cuando estén presentes en ambos test sets**.

### Resultado buscado para BIB

No:

> Hard task performs better than easy.

Sino:

> **The apparent learnability of protein function changes substantially with the definition of the comparison class.**

---

# RQ2 — ¿Qué hace realmente la reducción de redundancia R?

Aquí las **curvas serán centrales**.

Para cada S fijo:

$$
M(R\mid S=s)
$$

con:

```text
R = 90 → 70 → 50 → 30
```

y cuatro curvas correspondientes a:

```text
S90
S70
S50
S30
```

### No miraremos solo R90→R30

Analizaremos los seis contrastes:

```text
90→70
90→50
90→30
70→50
70→30
50→30
```

a cada S.

Eso permite descubrir:

* thresholds;
* comportamiento no monotónico;
* regiones de saturación;
* puntos donde empieza el colapso;
* diferencias entre representaciones;
* diferencias entre tareas.

### Pero hay una confusión importante que tenemos que resolver

Reducir R también reduce:

```text
n_train
n_positive
number of clusters
effective diversity
```

Por lo tanto añadiremos:

$$
\Delta M_R
\quad\text{vs}\quad
\Delta n_{\mathrm{train}}
$$

y, cuando esté disponible:

$$
\Delta M_R
\quad\text{vs}\quad
\Delta n_{\mathrm{clusters}}.
$$

La pregunta importante será:

> **¿El efecto de R es simplemente data quantity, o queda un efecto asociado a redundancy structure después de considerar cuánto training data eliminamos?**

Para BIB esto es muchísimo más interesante que “lower R gives lower MCC”.

---

# RQ3 — ¿Qué hace realmente la separación por homología S?

Exactamente la construcción inversa.

Para cada R fijo:

$$
M(S\mid R=r)
$$

con:

```text
S = 90 → 70 → 50 → 30
```

y cuatro curvas:

```text
R90
R70
R50
R30
```

Analizaremos nuevamente todos:

```text
90→70
90→50
90→30
70→50
70→30
50→30
```

### Pregunta principal

> **How rapidly does functional recognition deteriorate as evaluation becomes increasingly remote in sequence space?**

Y queremos saber si esa degradación depende de:

```text
task
representation
algorithm
```

---

# RQ4 — ¿R y S son realmente factores independientes?

Acá está probablemente una de las contribuciones más fuertes.

Tenemos las **36 interacciones 2×2** de cada superficie.

Conceptualmente:

$$
I=
\Delta_R(S_{\mathrm{strict}})
-
\Delta_R(S_{\mathrm{loose}})
$$

y equivalentemente desde el efecto S.

La pregunta no es solamente si `interaction ≠ 0`.

Queremos descubrir patrones:

> Cuando R ya es restrictivo, ¿S deja de importar?

> Cuando S es muy restrictivo, ¿la pérdida de redundancia se vuelve más dañina?

> ¿Las representaciones aprendidas reducen esa interacción?

> ¿Las funciones biológicas presentan superficies diferentes?

---

# RQ5 — ¿Cuál factor pesa más: R o S?

Para cada:

```text
task × representation × algorithm
```

calcularemos resúmenes de:

$$
E_R=\mathrm{mean}(|\Delta_R|)
$$

y:

$$
E_S=\mathrm{mean}(|\Delta_S|).
$$

Después:

$$
E_R/E_S.
$$

Pero además lo haremos por intervalo:

```text
R90→70 vs S90→70
R70→50 vs S70→50
R50→30 vs S50→30
```

Esto nos permitirá obtener mapas tipo:

```text
R-dominated regime
S-dominated regime
mixed regime
```

Eso ya empieza a construir nuestro concepto de **generalization regimes**.

---

# RQ6 — ¿Cuál es la geometría completa de la superficie R×S?

No basta con curvas separadas.

Cada modelo tiene una superficie:

$$
M(R,S)
$$

representada por las 16 celdas.

Exploraremos tres vistas complementarias.

| Vista                 | Pregunta                            |
| --------------------- | ----------------------------------- |
| Curvas de R a S fijo  | ¿cómo responde a redundancy?        |
| Curvas de S a R fijo  | ¿cómo responde a remote evaluation? |
| Superficie/matriz R×S | ¿cómo interactúan ambas?            |

Además definiremos la degradación respecto a la condición permisiva:

$$
D(R,S)=M(R,S)-M(90,90).
$$

Así podremos comparar modelos que parten de performances distintas según:

> **cuánto rendimiento consiguen conservar al desplazarnos por el espacio de generalización.**

---

# RQ7 — ¿La representación cambia la forma de la generalización?

Este es uno de los análisis más importantes para BIB.

Tenemos una secuencia conceptual perfecta:

```text
Roxy
→ physicochemical information

One-hot
→ explicit sequence information

ESM2-T6
→ small pretrained representation

ESM2-T30
→ larger pretrained representation
```

Para cada:

```text
task × R × S × algorithm
```

calcularemos todos los contrasts:

```text
OneHot − Roxy
T6 − Roxy
T30 − Roxy

T6 − OneHot
T30 − OneHot

T30 − T6
```

### Pero la pregunta fuerte será

No:

> ¿Cuál representación da más MCC?

Sino:

> **Which representations retain functional information as redundancy and train–test relatedness are progressively removed?**

Por tanto estudiaremos:

$$
G_{\mathrm{rep}}(R,S)
$$

como **superficie de ganancia**.

---

# RQ8 — ¿Escalar ESM2 de T6 a T30 sirve realmente?

Esta merece análisis independiente.

Definimos:

$$
G_{\mathrm{scale}}(R,S)
=
M_{T30}(R,S)-M_{T6}(R,S).
$$

Y lo estudiamos en función de:

```text
task
algorithm
R
S
R×S
```

Las tres conclusiones posibles son igualmente interesantes:

```text
T30 > T6 everywhere
→ model capacity broadly matters
```

```text
T30 > T6 mainly at strict R/S
→ capacity matters under difficult generalization
```

```text
T30 ≈ T6
→ a small PLM already captures most useful functional signal
```

También añadiremos:

```text
worst-case gain
mean gain
win rate
loss rate
rank stability
```

---

# RQ9 — ¿Cuánto importa realmente el algoritmo downstream?

Compararemos los cinco learners a representación fija.

No queremos otro leaderboard LR/RF/KNN/XGB/MLP.

Queremos contrastar:

$$
\text{algorithm effect}
$$

contra:

$$
\text{representation effect}.
$$

Calcularemos por celda:

```text
range across representations
range across algorithms
```

y después su relación.

Además:

```text
winner count
mean rank
worst-case rank
regret
performance range
```

### Pregunta Briefings

> **Once proteins are represented appropriately, how much predictive behavior is still determined by classifier complexity?**

Si KNN sobre ESM queda cerca de XGB/MLP, eso es conceptualmente muy potente.

---

# RQ10 — ¿Existen generalization regimes específicos de función?

Separaremos:

```text
global lytic recognition
endolysin
depolymerase
VAL
multiclass
```

Para cada función analizaremos:

```text
R sensitivity
S sensitivity
R×S interaction
representation gain
T30-T6 gain
algorithm dependence
```

Esto puede llevarnos a identificar algo tipo:

```text
DPO → robust functional regime
Endolysin → intermediate
VAL → evidence-limited regime
```

sin imponer esos resultados de antemano.

---

# RQ11 — VAL como caso explícito de Low-N

VAL debe ir separado de los resúmenes globales cuando la varianza dependa fuertemente de support.

Para cada R/S tendremos:

```text
n_train positives
n_test positives
MCC
recall
precision
F1
```

y analizaremos:

$$
M\quad\text{vs}\quad n_{\mathrm{positive}}.
$$

También:

```text
representation gain vs support
R effect vs support
S effect vs support
```

Así podremos distinguir:

> generalization failure

de:

> sampling volatility.

---

# RQ12 — ¿Qué cambia al reconocer explícitamente label uncertainty?

Esta es la gran extensión que separa el paper actual del data-centric.

En el mismo universo P/U tenemos:

```text
naive P-vs-U
BaggingPU
Elkan–Noto
```

manteniendo fijos:

```text
task
R
S
representation
base learner
```

Por tanto es una comparación **controlada de learning paradigm**.

Definiremos:

$$
G_{\mathrm{BaggingPU}}
=
M_{\mathrm{BaggingPU}}-M_{\mathrm{naive}}
$$

y:

$$
G_{\mathrm{EN}}
=
M_{\mathrm{EN}}-M_{\mathrm{naive}}.
$$

Aquí `M` será un diagnóstico PU válido, no una falsa interpretación de U como negativo.

Actualmente tenemos:

```text
p_recall_above_u_q90
p_recall_above_u_q95
p_recall_above_u_q99
```

y desde los scores podemos añadir ranking diagnostics.

---

# RQ13 — Curvas PU en función de R y S

Tal como dijiste, **necesitamos curvas**, no solo tablas.

### A R fijo/S fijo

Para cada S:

```text
x = R

naive
BaggingPU
Elkan–Noto
```

y viceversa:

```text
x = S

naive
BaggingPU
Elkan–Noto
```

para:

```text
conservative
broad

Roxy
One-hot
T6
T30

LR
XGB
```

### Preguntas

> ¿PU ayuda más cuando baja R?

> ¿PU ayuda más cuando baja S?

> ¿existe un crossover donde naive funciona bien en condiciones fáciles pero PU gana en condiciones difíciles?

Eso puede ser uno de los resultados más fuertes del paper.

---

# RQ14 — ¿Cómo depende la ganancia PU de R×S?

Así como tenemos representación-gain surfaces, tendremos:

$$
G_{\mathrm{PU}}(R,S).
$$

Eso permite estudiar:

```text
paradigm × redundancy
paradigm × homology
paradigm × representation
paradigm × learner
```

y especialmente:

$$
R\times S\times\text{paradigm}.
$$

No necesitamos necesariamente modelar una interacción estadística de tres vías formal; primero la estudiaremos como contrastes factoriales y superficies de gain.

---

# RQ15 — ¿Influye la definición del unlabeled universe?

Compararemos:

```text
PU conservative
vs
PU broad
```

en términos de:

```text
positive enrichment
ranking
score distributions
PU gain
method agreement
```

Pero existe un análisis especialmente limpio:

`U_explicit` está contenido en el universo broad.

Podemos comparar **los mismos U_explicit** puntuados bajo modelos entrenados con:

```text
conservative U
vs
broad U
```

y analizar:

```text
Spearman correlation
score shift
top 1% overlap
top 5% overlap
rank displacement
```

Eso responde:

> **Does expanding the uncertain background alter what the model considers functionally plausible?**

---

# RQ16 — Representation × learning paradigm

Esto puede ser otra contribución muy buena.

Queremos estudiar si:

$$
G_{\mathrm{PU}}(T30)
>
G_{\mathrm{PU}}(Roxy)
$$

o al revés.

Preguntas:

> ¿Una representación fuerte reduce la necesidad de PU?

> ¿O una representación rica permite que PU aproveche mejor la estructura de U?

> ¿T30 estabiliza los rankings más que T6 bajo PU?

Esto es **effect-of-effect**, no simple leaderboard.

---

# RQ17 — Robustez: ¿qué significa realmente “ganar”?

No habrá un único ganador.

Para cada familia calcularemos:

| Concepto                | Definición                         |
| ----------------------- | ---------------------------------- |
| Peak                    | mejor celda                        |
| Average                 | promedio 16 R×S                    |
| Median                  | mediana 16 R×S                     |
| Worst case              | mínima performance                 |
| Range                   | max − min                          |
| Win count               | número de celdas ganadas           |
| Mean rank               | ranking promedio                   |
| Worst rank              | peor ranking                       |
| Regret                  | distancia al ganador de cada celda |
| R sensitivity           | magnitud de efecto de R            |
| S sensitivity           | magnitud de efecto de S            |
| Interaction sensitivity | magnitud R×S                       |

Entonces podremos distinguir:

```text
best peak
best average
best robust
best remote-generalization
best redundancy-resistant
best efficient
```

Eso responde mucho mejor “¿qué modelo ganó?”.

---

# RQ18 — Common-test sensitivity

Este análisis es **obligatorio** antes de hacer afirmaciones fuertes sobre R/S.

Los test sets cambian con R/S.

Por tanto los native-test deltas:

$$
M(R_1,S_1)-M(R_2,S_2)
$$

mezclan:

```text
model/generalization effect
+
test-composition effect
```

### Solución

Para comparaciones de S a R fijo:

```text
intersection of test entities
across S90/S70/S50/S30
```

y recalculamos performance únicamente sobre entidades compartidas.

Cuando la intersección de cuatro niveles sea pequeña:

```text
pairwise intersections
```

Lo mismo como sensitivity secundaria para R.

Después comparamos:

```text
native-test effect
vs
common-test effect
```

Si ambos cuentan la misma historia, nuestra conclusión queda mucho más blindada.

---

# RQ19 — Paired uncertainty / bootstrap

No necesitamos hacer bootstrap absurdo sobre las 48 mil comparaciones.

Lo aplicaremos donde el contraste es científicamente principal.

### Supervised

Paired bootstrap cuando existen las mismas entidades:

```text
T30 vs T6
T6 vs Roxy
key algorithm comparisons
common-test R/S contrasts
```

### PU

Resampling estratificado de:

```text
P
U
```

para:

```text
BaggingPU vs naive
EN vs naive
T30 vs T6
```

### Multiclass

Bootstrap estratificado por clase.

Reportaremos:

```text
Δ metric
95% CI
directional consistency
```

en lugar de perseguir p-values masivamente.

---

# RQ20 — ¿Las curvas son consistentes o hay cambios de régimen?

Para cada superficie estudiaremos no solo endpoints sino **stepwise changes**:

```text
90→70
70→50
50→30
```

Eso permite identificar:

```text
early degradation
mid-range collapse
late collapse
stable regime
non-monotonic regime
```

Podemos después clasificar descriptivamente las superficies.

Si los patrones son muy claros, una fase exploratoria secundaria podría agrupar **response profiles de 16 celdas** para ver si ciertas tareas/representaciones caen en los mismos archetypes.

Esto sería exploratorio, no una afirmación central salvo que resulte muy interpretable.

---

# RQ21 — Mixed R×S trajectories

Además de las curvas ortogonales, exploraremos la dificultad conjunta.

La más intuitiva:

```text
R90/S90
R70/S70
R50/S50
R30/S30
```

pero **nunca será usada para inferir aisladamente R o S**, porque ambos cambian.

También podemos usar:

$$
D(R,S)
$$

respecto de `R90/S90` para obtener curvas de **performance retention under joint strictness**.

Otra vista útil será agrupar celdas por niveles comparables de joint strictness y mostrar cuánto varían según cómo se distribuye esa strictness entre R y S.

Eso ayuda a demostrar que:

> `R50/S90` y `R90/S50` no necesariamente representan la misma dificultad.

---

# RQ22 — Performance vs cantidad y estructura de información

Esta merece notebook propio.

Queremos relacionar performance con:

```text
n_train
n_test
positive support
prevalence
number of sequence clusters
R
S
```

Esto es crucial para interpretar R.

Podemos construir modelos descriptivos bloqueados por:

```text
task
representation
algorithm
```

y preguntar cuánto del efecto de R permanece después de considerar training size.

No trataría R/S automáticamente como variables continuas en la inferencia principal; las mantendríamos como factores ordenados y usaríamos tendencias continuas solo como diagnóstico.

---

# RQ23 — Compute/performance

Secundario, pero útil.

Para cada learner:

```text
performance
vs
training runtime
```

con:

```text
representation
task
R/S
```

Eso podría reforzar conclusiones tipo:

> a simple learner operating on a strong representation approaches the performance of more complex models.

Caveat importante:

**el runtime del embedding debe separarse del runtime downstream**.

Si no tenemos tiempos completos de embedding, no venderemos esto como end-to-end computational efficiency.

---

# RQ24 — Candidate prioritization como demostración, no como otro paper

Nada de docking por obligación.

Después de congelar qué modelos/paradigmas son robustos, seleccionamos candidatos de U usando:

```text
high predictive/ranking support
+
consensus across robust models
+
persistence across R/S
+
agreement between PU methods
+
agreement T6/T30
+
physicochemical support
+
latent-space support
+
nontrivial sequence relatedness
```

### Análisis de espacio latente

Queremos saber si las candidatas:

* no son OOD extremas;
* ocupan regiones soportadas por positivos;
* no son simplemente duplicados próximos de positivos conocidos;
* permanecen separadas de regiones claramente negativas cuando corresponda.

### Fisicoquímica

Compararemos sus Roxy descriptors con la distribución positiva:

```text
percentile profiles
multivariate distance
range compatibility
```

Sin afirmar mecanismo.

### Conclusión correcta

> prioritized for experimental characterization.

No:

> computationally validated lytic enzymes.

Docking, estructura mecanística, motifs profundos, etc. se reservan para **mining/discovery**.

---

# RQ25 — Síntesis de “generalization regimes”

Esta debe ser la culminación de todo el análisis.

A partir de RQ1–24 queremos identificar regímenes como:

| Regime                       | Característica                          |
| ---------------------------- | --------------------------------------- |
| **Evidence-rich**            | muchos modelos funcionan                |
| **Redundancy-sensitive**     | R domina                                |
| **Homology-sensitive**       | S domina                                |
| **Joint-stringency**         | interacción R×S fuerte                  |
| **Representation-dependent** | PLM gain aumenta                        |
| **Capacity-dependent**       | T30 gana a T6                           |
| **Label-uncertain**          | PU supera naive                         |
| **Evidence-limited**         | support domina, e.g. potencialmente VAL |

No debemos imponer que todos existan.

**Los datos deben decirnos cuáles aparecen realmente.**

Esta sería la principal diferenciación respecto del paper data-centric.

---

# Contrato estadístico del análisis

Antes de notebooks, congelaría estas reglas:

| Comparación                                       | Naturaleza                                         |
| ------------------------------------------------- | -------------------------------------------------- |
| Representación dentro de mismo task/R/S/algorithm | **Paired**                                         |
| Algoritmo dentro de mismo task/R/S/representation | **Paired**                                         |
| PU method dentro de mismo task/R/S/rep/learner    | **Paired**                                         |
| R/S native tests                                  | **No asumir paired**                               |
| Common-test R/S                                   | **Paired**                                         |
| Different task definitions                        | Generalmente **not paired**, salvo subsets comunes |
| Conservative vs broad PU                          | Comparación global + subset U común                |
| VAL                                               | Siempre support-aware                              |

Y para AP:

> siempre reportar prevalence/AP baseline o AP gain.

---

# Cómo dividiría finalmente los notebooks

No haría uno por cada plot. Haría **workbenches científicos**:

| Notebook | Contenido                             |
| -------- | ------------------------------------- |
| **00**   | QC + global experimental map          |
| **01**   | Task/evidence formulation             |
| **02**   | Complete R curves                     |
| **03**   | Complete S curves                     |
| **04**   | R×S surfaces and interactions         |
| **05**   | R-vs-S dominance + data-size effects  |
| **06**   | Representation effects                |
| **07**   | T30 vs T6 scaling                     |
| **08**   | Algorithm dependence + robustness     |
| **09**   | Functional subtype regimes            |
| **10**   | VAL / Low-N                           |
| **11**   | Binary vs multiclass formulation      |
| **12**   | PU paradigm                           |
| **13**   | PU R/S curves + PU gain surfaces      |
| **14**   | Conservative vs broad U               |
| **15**   | Representation × PU paradigm          |
| **16**   | Common-test sensitivity               |
| **17**   | Paired bootstrap / uncertainty        |
| **18**   | Compute/performance                   |
| **19**   | Candidate prioritization              |
| **20**   | Final generalization-regime synthesis |

Son muchos notebooks, sí, pero acá **cada uno responde una pregunta distinta**. No serían 20 notebooks haciendo el mismo heatmap.

---

# Qué probablemente iría al main paper

Todavía no seleccionaría figuras finales, pero la arquitectura probable sería:

| Main figure | Historia                                                        |
| ----------- | --------------------------------------------------------------- |
| **Fig. 1**  | Experimental design: task × R × S × representation × paradigm   |
| **Fig. 2**  | Task formulation / evidence definition                          |
| **Fig. 3**  | Complete R/S generalization landscape                           |
| **Fig. 4**  | Representation + PLM scaling across R×S                         |
| **Fig. 5**  | Functional regimes / subtype behavior                           |
| **Fig. 6**  | Traditional vs PU and PU-gain landscape                         |
| **Fig. 7**  | Generalization-regime synthesis + small candidate demonstration |

El resto puede ir brutalmente bien al SI porque tenemos un factorial enorme.

---

# La prueba de fuego para *Briefings*

Al terminar el análisis deberíamos poder responder **sí** a varias de estas afirmaciones:

> **R and S affect generalization through distinguishable mechanisms.**

> **Their effects are not universally additive.**

> **Protein representation changes the shape of the generalization landscape.**

> **PLM scaling has a regime-dependent—or demonstrably negligible—effect.**

> **Classifier choice explains less variation than upstream problem formulation/representation, si eso efectivamente ocurre.**

> **Functional classes occupy distinct generalization regimes.**

> **Treating unlabeled proteins as negatives changes predictive behavior, and PU provides identifiable advantages under specific regimes.**

> **These effects cannot be inferred from a single redundancy threshold or a single homology split.**

Si encontramos **varias de esas conclusiones con evidencia consistente**, yo iría a *Briefings* sin demasiada duda.

Si al final resulta que todas las superficies son prácticamente paralelas, T30≈T6, PU≈naive y las interacciones son mínimas, el trabajo sigue siendo bueno, pero ahí **Protein Science vuelve a ser el target natural**.

Ese sería para mí el roadmap completo. El siguiente paso ya no sería programar todavía: sería convertir **RQ1–RQ25 en una matriz operacional** con `input exacto → cálculo → estadístico → salida esperada → notebook`, y recién ahí implementar los workbenches uno por uno.

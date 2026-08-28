## RQ1 — How does evidence/background definition reshape apparent functional learnability?

La definición del background tuvo un efecto sustancial sobre la dificultad aparente de reconocer proteínas líticas. A través del factorial completo `R×S`, `global_hard` mostró una separabilidad consistentemente mayor que `global_easy` y `global_mixed`, mientras que `global_easy` y `global_mixed` se comportaron prácticamente como el mismo problema. En MCC, `hard−easy` presentó una mediana de aproximadamente **+0.17**, mientras que `mixed−easy` fue esencialmente nulo (**≈ −0.004**). El uso de AP normalizado por prevalencia confirmó esta conclusión: `hard−easy` alcanzó una mediana de **ΔAP-skill ≈ +0.153**, `hard−mixed ≈ +0.158`, y `mixed−easy ≈ −0.005`. Esta normalización fue necesaria porque el AP gain bruto está limitado por `1−prevalence` y, por tanto, penalizaba artificialmente la tarea `global_hard`, cuya prevalencia positiva es considerablemente mayor. 

| Hallazgo                      |                                              Resultado principal | Interpretación                                                                                      |
| ----------------------------- | ---------------------------------------------------------------: | --------------------------------------------------------------------------------------------------- |
| **Hard vs easy**              |                         ΔMCC ≈ **+0.17**; ΔAP-skill ≈ **+0.153** | La naturaleza del background cambia fuertemente la separabilidad                                    |
| **Mixed vs easy**             |                        ΔMCC ≈ **−0.004**; ΔAP-skill ≈ **−0.005** | Añadir `N_hard` al background amplio prácticamente no cambia el problema                            |
| **Hard vs mixed**             |                                           ΔAP-skill ≈ **+0.158** | `N_hard` aislado define un régimen predictivo distinto                                              |
| **Shared easy→mixed**         |                                                ΔMCC ≈ **−0.004** | El null result no se explica por diferencias de test composition                                    |
| **Shared hard→mixed**         | ΔMCC ≈ **−0.229**, ΔROC-AUC ≈ **−0.035**, ΔAP-skill ≈ **−0.043** | Hay pérdida real de discrimination, pero el efecto sobre decisiones a threshold fijo es mucho mayor |
| **Representation dependence** |     One-hot > Roxy > T6 > T30 en sensibilidad al task definition | Las representaciones aprendidas son progresivamente más estables frente a cambios del background    |

El análisis sobre **entidades test compartidas** fue especialmente informativo. `easy` y `mixed` permanecieron prácticamente indistinguibles cuando se evaluaron exactamente las mismas proteínas, descartando que su similitud fuera consecuencia de cambios en composición del outer test. En contraste, la diferencia `hard→mixed` persistió e incluso se amplificó para MCC, demostrando que la formulación de la tarea modifica genuinamente el comportamiento predictivo. Sin embargo, las pérdidas en ROC-AUC y AP-skill fueron bastante menores que las pérdidas en MCC y balanced accuracy, indicando que task formulation afecta **tanto la discriminación como, de forma más pronunciada, el operating point y la distribución de scores**. El análisis extendido fue diseñado precisamente para distinguir estas dos componentes. 

El efecto más fuerte de RQ1 fue su dependencia con la representación. La mediana de `hard−easy` en AP-skill disminuyó aproximadamente de **+0.476 para one-hot**, a **+0.238 para Roxy**, **+0.090 para ESM2-T6** y solamente **+0.024 para ESM2-T30**. De forma concordante, al pasar de `hard` a `mixed`, los scores de positivos conocidos disminuyeron mucho más con one-hot y Roxy que con los PLMs; T30 mostró el menor desplazamiento de scores y la mayor conservación del ranking. Esto sugiere que **la sensibilidad del reconocimiento funcional a la construcción del background disminuye progresivamente con representaciones más informativas**, aunque todavía debemos separar este comportamiento de un posible ceiling effect en análisis posteriores de representación y PLM scaling.

También apareció una señal de interacción con redundancy: la diferencia entre task formulations tendió a aumentar a medida que `R` se volvió más restrictivo, mientras que la dependencia con `S` fue más débil y menos monotónica. No conviene desarrollar esa conclusión dentro de RQ1 porque será analizada formalmente en las preguntas específicas de R, S y `R×S`.

### Conclusión de RQ1

> **The apparent learnability of phage lytic function is strongly dependent on how the comparison class is defined. A lysis-adjacent background constitutes a substantially more separable learning regime than the broader background, whereas adding this subset to the broader negative universe has negligible impact. Importantly, sensitivity to task formulation decreases markedly with increasingly informative protein representations, reaching near-invariance with ESM2-T30. Common-entity analyses further show that these effects reflect both genuine changes in discrimination and stronger shifts in score/decision behavior, rather than differences in test-set composition alone.**

## RQ2 — Effect of progressive redundancy reduction

La reducción progresiva de redundancia tuvo un efecto global negativo pero **mucho menor que la magnitud de la reducción del dataset**. Al pasar de `R90` a `R30`, la mediana global fue de aproximadamente **ΔMCC = −0.070**, con deterioro en ~**82%** de las comparaciones. Sin embargo, a `R30` permanece solo ~**31% de los training entities** originales. Esto indica que una gran fracción de las secuencias eliminadas aporta información redundante y que el rendimiento no escala proporcionalmente con el tamaño bruto del training set.

La estructura de la evidencia restante fue más informativa que el simple número de secuencias. La reducción de R eliminó entidades más rápido que riqueza relativa de sequence space: aproximadamente, `R70` conserva ~**72% de entidades pero casi toda la riqueza relativa de clusters**, `R50` ~**50% de entidades y ~86% de cluster richness**, y `R30` ~**31% de entidades pero ~56% de cluster richness**. Por tanto, las primeras etapas de redundancy reduction eliminan principalmente repetición dentro de regiones ya representadas, mientras que en regímenes muy estrictos comienza también una pérdida sustancial de diversidad. Los clusters se interpretan aquí como **relative S-level cluster richness**, no como conjuntos literalmente preservados entre R, porque cada universo R fue construido y reclusterizado independientemente. 

El efecto de R fue fuertemente **task-dependent**. Endolysin mostró la mayor dependencia de redundancia, con `R90→R30` de aproximadamente **−0.143 MCC** y deterioro prácticamente universal. `global_easy` y `global_mixed` mostraron pérdidas cercanas a **−0.10 MCC**, mientras que `global_hard` fue considerablemente más robusto, con ~**−0.046 MCC**. Depolymerase constituyó una excepción clara: el rendimiento fue esencialmente estable e incluso mejoró ligeramente bajo redundancy reduction, con una mediana cercana a **+0.018 MCC**. VAL y multiclass mostraron respuestas predominantemente no monotónicas, compatibles con regímenes de evidencia limitada o heterogeneidad entre clases.

No apareció un **breakpoint universal** de redundancia. Los pasos `R90→R70`, `R70→R50` y `R50→R30` tuvieron pérdidas medianas similares, y el intervalo con mayor deterioro cambió entre tasks y modelos. Aproximadamente la mitad de las curvas fueron monotonic decreasing y casi la otra mitad non-monotonic, mostrando que el efecto de redundancy reduction depende del contexto funcional y del modelo, no de un único threshold de identidad.

La **representación** moduló fuertemente la robustez frente a pérdida de redundancia. La sensibilidad siguió aproximadamente la jerarquía:

> **one-hot > Roxy > ESM2-T6 > ESM2-T30**

en magnitud de pérdida. Para `R90→R30`, one-hot mostró una mediana cercana a **−0.095 MCC**, frente a ~**−0.057** para T30. Más importante aún, al normalizar por la cantidad de training entities eliminados, la pérdida mediana por cada 10% de datos removidos disminuyó de aproximadamente **0.010 MCC para one-hot** a **0.006 para T30**. Esto indica que las representaciones pretrained explotan de manera más eficiente la evidencia no redundante restante y que la mayor sensibilidad de representaciones simples no puede explicarse solo por reducción de sample size. La lógica de este análisis está alineada con el diseño metodológico de comparar robustness frente a R dentro de condiciones controladas. 

La cantidad de datos eliminada, por sí sola, explicó muy poco de la pérdida predictiva global: la asociación entre fracción de entidades removidas y MCC loss fue cercana a cero en promedio. En cambio, **performance retention mostró una asociación descriptiva mayor con relative cluster richness que con entity retention**, sugiriendo que conservar diversidad de sequence space puede ser más importante que conservar un gran número de ejemplos redundantes.

Finalmente, el efecto de R estuvo claramente condicionado por S. La pérdida `R90→R30` fue mayor bajo condiciones de split más permisivas y disminuyó hacia `S30`, lo que sugiere que la importancia marginal de redundancy reduction cambia cuando train y test ya están fuertemente separados por homología. Esto constituye una señal clara de una posible **R×S interaction**, que debe analizarse formalmente en la RQ correspondiente y no cerrarse dentro de RQ2.

### Take-home message

> **Redundancy reduction does not simply make the dataset smaller. It changes the structure of the available evidence. Large numbers of sequences can be removed while preserving a much larger fraction of sequence-space richness and predictive performance, but the cost of doing so depends strongly on biological task and protein representation. Pretrained representations, particularly ESM2-T30, are substantially more robust to the loss of redundant sequence support than simpler sequence or physicochemical encodings.**

## RQ3 — Effect of train–test homology separation

La separación progresiva por homología tuvo un efecto global **real pero moderado** sobre la performance. Para `S90→S30`, la mediana global fue de aproximadamente **ΔMCC = −0.017**, con deterioro en ~**71%** de las comparaciones. Este efecto fue considerablemente menor que el observado para redundancy reduction en RQ2, lo que indica que `R` y `S` no capturan exactamente el mismo mecanismo de generalización.

No apareció un **breakpoint universal de homología**. Los pasos `S90→S70`, `S70→S50` y `S50→S30` mostraron pérdidas progresivamente mayores, pero el threshold dominante cambió entre tasks y modelos. Además, aproximadamente **63% de las curvas fueron non-monotonic**, lo que sugiere que la generalización remota depende fuertemente de la topología particular de los clusters y de qué grupos completos terminan en train y test, no de un único cutoff de identidad.

El efecto de S fue claramente **task-dependent**. `global_mixed`, multiclass y endolysin mostraron pérdidas moderadas bajo `S90→S30`, mientras que `global_hard` fue nuevamente uno de los regímenes más robustos, con una mediana cercana a **−0.007 MCC** y una pérdida mínima también en AP-skill. VAL mostró la mayor volatilidad y debe seguir tratándose como un **evidence-limited regime**. Depolymerase fue comparativamente poco sensible a homology separation.

La **representación** volvió a modular fuertemente la robustez. La sensibilidad absoluta a `S90→S30` siguió aproximadamente:

> **one-hot > Roxy > ESM2-T6 > ESM2-T30**

con T30 mostrando la menor degradación. Más importante aún, la **ventaja relativa de T30 aumentó a medida que S se volvió más estricto**. Por ejemplo, el gap mediano `T30−T6` creció desde aproximadamente **+0.069 MCC en S90** hasta **+0.091 en S30**, mientras que la ventaja frente a one-hot aumentó de ~**+0.37** a ~**+0.46 MCC**. Esto sugiere que la capacidad de las representaciones pretrained se vuelve progresivamente más valiosa cuando disminuye la similitud local entre train y test.

El tamaño relativo de train y test permaneció prácticamente constante entre S, con test fractions cercanas a **0.30** en todos los niveles. Por tanto, el deterioro observado no puede explicarse por cambios gruesos en sample size. Lo que sí cambió de forma marcada fue la **granularidad/topología de los clusters**: hacia `S30` disminuyó sustancialmente el número de clusters y aumentó el número de entidades por cluster. En consecuencia, S debe interpretarse como un cambio en la estructura de la partición homológica y no simplemente como un “split más difícil”.

Los cambios de prevalence entre S fueron pequeños y no explicaron la variación de MCC; la asociación entre `Δprevalence` y `ΔMCC` fue prácticamente nula. Esto apoya que el efecto de S refleja diferencias reales en el régimen de generalización y no simples cambios de balance de clases.

### Hallazgo central de RQ3

El efecto de S estuvo **fuertemente condicionado por R**. Para `S90→S30`, la mediana de ΔMCC fue aproximadamente:

```text
R90   -0.047
R70   -0.031
R50   -0.006
R30    ~0.000
```

Así, cuando el universo conserva alta redundancia (`R90`), imponer una separación train–test más estricta tiene un costo claro y consistente. En cambio, cuando la redundancia ya ha sido reducida agresivamente (`R30`), el efecto marginal adicional de S prácticamente desaparece.

### Take-home message

> **Remote homology separation imposes a genuine but globally modest generalization penalty. This penalty is highly task- and representation-dependent, and pretrained representations—especially ESM2-T30—retain substantially more performance as train–test similarity decreases. Most importantly, the cost of stricter homology separation collapses as redundancy control becomes more aggressive, strongly indicating that R and S do not act independently and motivating a formal R×S interaction analysis.**

## RQ4 — Interaction between redundancy reduction and homology separation

RQ4 confirmó que **R y S no actúan de forma aditiva**. A través de los **5,040 rectángulos 2×2** posibles del factorial, la interacción MCC mostró una mediana positiva de aproximadamente **+0.015**, con cerca de **66% de interacciones positivas**. En el contraste extremo `R90/R30 × S90/S30`, la señal fue mucho más fuerte: la mediana de interacción alcanzó aproximadamente **+0.051 MCC** y **~89% de las 140 superficies** mostraron interacción positiva.

Dado que tanto redundancy reduction como homology separation producen en general pérdidas de performance, una interacción positiva implica que:

> **the combined strict R/S regime performs better than expected from simply adding the isolated penalties of R and S.**

En otras palabras, **los efectos de R y S se atenúan parcialmente entre sí**.

Esta atenuación fue cuantitativamente importante. Entre las superficies donde el modelo aditivo predecía una pérdida conjunta de al menos `0.05 MCC`, la interacción recuperó en mediana aproximadamente **33% del penalty esperado**. Cerca de **73%** de esas superficies recuperaron más del 25% de la pérdida aditiva, y ~**24%** recuperaron más de la mitad.

La descomposición geométrica de las superficies 4×4 mostró que R sigue siendo el principal determinante del landscape, explicando en mediana aproximadamente **70.5% de la variación de MCC**, mientras que S explicó ~**9.1%**. Sin embargo, la interacción `R×S` explicó aproximadamente **18.3%**, es decir, alrededor del doble de la contribución independiente de S. Por tanto, una descripción basada únicamente en main effects habría perdido una parte sustancial de la estructura de generalización.

El peso de la interacción fue nuevamente **task-dependent**. Depolymerase, multiclass y VAL mostraron las mayores fracciones relativas de estructura no aditiva, mientras que global easy/mixed y endolysin estuvieron más dominados por el efecto principal de R. Aun así, excluyendo el régimen inestable de VAL, las medianas de interacción extrema fueron positivas para prácticamente todas las tareas, indicando que la atenuación R×S no depende de una sola función.

La **representación** también moduló fuertemente la interacción. La magnitud extrema siguió aproximadamente:

> **one-hot > Roxy > ESM2-T6 > ESM2-T30**

con medianas cercanas a **+0.122**, **+0.058**, **+0.046** y **+0.027 MCC**, respectivamente. De forma consistente, la fracción de variación de la superficie atribuible a interacción disminuyó desde ~25% para one-hot hasta ~9% para T30. Esto refuerza el patrón observado desde RQ1: **representaciones pretrained más informativas producen generalization landscapes más estables frente a cambios en la definición de la tarea y del protocolo de evaluación**.

La geometría local mostró que **no existe un único hotspot universal de interacción**. El rectángulo `R70→R50 × S70→S50` presentó la mayor interacción mediana local, pero solo fue dominante en una minoría de superficies. La no-aditividad está distribuida por diferentes regiones del factorial dependiendo de task, representation y learner, por lo que no existe un único threshold de identidad capaz de resumir el fenómeno.

El resultado más importante conecta directamente con RQ3. El penalty `S90→S30` cambia aproximadamente de:

```text
R90   ΔMCC ≈ -0.047
R30   ΔMCC ≈  0.000
```

mientras que el efecto `R90→R30` también disminuye cuando S se vuelve estricto:

```text
S90   ΔMCC ≈ -0.084
S30   ΔMCC ≈ -0.040
```

Por tanto, **hacer más estricto uno de los factores reduce el efecto marginal del otro**. R y S capturan componentes diferentes del protocolo, pero sus penalties se comportan como **partially substitutive constraints** sobre la generalización.

La interacción tampoco fue un artefacto exclusivo de MCC o del threshold de clasificación. El contraste extremo se mantuvo positivo con **AP skill (~+0.040)** y **ROC-AUC (~+0.008)**, con fuerte correlación entre las magnitudes de interacción de las distintas métricas. Para multiclass, macro-F1 también mostró una interacción positiva fuerte. Esto demuestra que la no-aditividad afecta tanto las decisiones a threshold fijo como la calidad del ranking/discriminación.

Existe cierta contribución de **ceiling/headroom**, especialmente en algunas representaciones, pero no explica por completo el fenómeno. La interacción permanece positiva a través de tareas, métricas y representaciones, y la desaparición progresiva del S effect al reducir R está directamente observada.

### Take-home message

> **Redundancy control and homology-separated evaluation are not independent constraints. Their joint effect is systematically less detrimental than predicted by an additive model, because aggressive redundancy reduction progressively attenuates the additional cost of remote train–test separation. The magnitude of this non-additivity depends strongly on biological task and representation, and more informative pretrained representations produce substantially more stable R×S generalization landscapes.**

Y la consecuencia metodológica más fuerte para el paper sería:

> **A single redundancy threshold combined with a single homology split cannot characterize model generalization, because the apparent effect of one factor depends on where the other factor is fixed.**

## RQ5 — R versus S dominance and the role of data quantity

RQ5 mostró que **R es globalmente el principal driver de cambio predictivo**, pero esa dominancia no es universal ni constante a través del landscape. La magnitud mediana del efecto absoluto fue aproximadamente **0.085 MCC para R** frente a **0.026 MCC para S**, con un ratio mediano `R/S ≈ 2.55`. El dominance index global fue cercano a **0.45**, y usando el cutoff principal cerca de **67% de las superficies fueron R-dominant**, ~28% balanced y solo ~5% S-dominant.

Esta conclusión fue robusta a una segunda definición independiente basada en la descomposición factorial `SS_R` versus `SS_S`, con alta concordancia entre ambas métricas de dominancia. Por tanto, la conclusión de que R domina globalmente no depende de una única formulación descriptiva.

El grado de dominancia fue fuertemente **task-dependent**. Endolysin y las tareas globales fueron marcadamente R-driven, con dominance indices altos y una gran mayoría de superficies clasificadas como R-dominant. En cambio, depolymerase mostró un régimen más intermedio, mientras que multiclass y VAL se comportaron como landscapes mucho más balanced y fuertemente interaction-dependent.

Un hallazgo central fue que **la identidad del factor dominante puede cambiar dentro de una misma superficie R×S**. El dominance index mediano aumentó progresivamente al avanzar hacia regímenes más estrictos, desde aproximadamente **0.28 en el nivel 90** hasta ~**0.67 en el nivel 30**, y la fracción de superficies R-dominant aumentó desde ~54% hasta ~78%. Aproximadamente **41% de las superficies cambiaron de signo de dominance en algún punto**, mostrando que no es correcto resumir todo el benchmark con una afirmación fija del tipo “R siempre importa más que S”.

La representación también modificó la estructura de dominancia. ESM2-T30 mostró el dominance index más alto, pero no porque fuera más sensible a R. Al contrario, T30 presentó los menores efectos absolutos tanto de R como de S; su alta dominancia relativa apareció porque **S se vuelve particularmente poco influyente con representaciones más fuertes**, dejando R como la principal fuente residual de variación. Esto refuerza la idea de que representaciones pretrained producen landscapes más robustos y menos sensibles a cambios en el protocolo.

RQ5 también clarificó que **R y S modifican componentes diferentes del problema**. El cambio `R90→R30` redujo el training universe en aproximadamente **69% de las entidades**, mientras que `S90→S30` modificó el tamaño de train/test en prácticamente **0%**. Así, R altera la cantidad y estructura de la evidencia disponible, mientras que S altera principalmente la **topología de la partición homológica**.

La pérdida asociada a R no escaló de forma simple con la cantidad de datos eliminados. La relación entre fracción de training removido y MCC loss fue prácticamente nula en promedio, confirmando que **R dominance no puede reducirse a un mero sample-size effect**. Además, la pérdida por unidad de evidencia eliminada fue menor para representaciones pretrained que para one-hot o Roxy.

En contraste, S cambió fuertemente la estructura de clusters manteniendo casi constante el tamaño del dataset. Hacia `S30` disminuyó el número relativo de clusters y aumentó el número de entidades por cluster, consistente con una partición más coarse-grained. Esto apoya la interpretación de S como un factor de **remote-evaluation topology** más que de cantidad de datos.

La dominancia también se preservó al usar métricas alternativas. MCC y AP-skill mostraron dominance indices muy similares y una concordancia considerable, mientras que multiclass con macro-F1 permaneció mucho más balanced. Por tanto, la historia no depende únicamente del threshold de clasificación.

Finalmente, el **driver map** mostró que los regímenes S-dominant fueron raros y aparecieron casi exclusivamente en superficies con alta interacción R×S. Esto indica que S rara vez actúa como un main effect limpio y dominante; cuando parece dominar, suele hacerlo dentro de landscapes fuertemente no aditivos.

### Take-home message

> **Redundancy reduction is the dominant global driver of predictive change, but dominance is itself regime-dependent. R primarily modifies the amount and redundancy structure of the available evidence, whereas S modifies train–test homology topology. The apparent dominant factor can change within the same model surface, especially in subtype and multiclass regimes, while globally defined tasks and endolysin remain consistently R-driven. S-dominant regimes are rare and largely confined to strongly interaction-driven landscapes.**

Y para mí la frase conceptual más importante es:

> **Generalization difficulty is not controlled by a single protocol factor; it emerges from the balance between evidence redundancy, homology-separated evaluation, and their interaction.**

## RQ6 — Effect of protein representation on predictive generalization

RQ6 mostró que la **representación de proteínas es uno de los principales determinantes tanto de performance como de robustez de generalización**. En comparaciones estrictamente pareadas dentro de la misma `task × R × S × algorithm`, la jerarquía fue muy consistente:

> **ESM2-T30 > ESM2-T6 > Roxy > one-hot**

La mediana global de MCC fue aproximadamente **0.902 para T30**, **0.829 para T6**, **0.675 para Roxy** y **0.496 para one-hot**. Más importante aún, los paired deltas mostraron ventajas muy consistentes: `T30−one-hot ≈ +0.41 MCC`, `T30−Roxy ≈ +0.24`, y `T30−T6 ≈ +0.08`, con T30 superando a one-hot y Roxy en prácticamente todas las comparaciones y a T6 en ~93% de los matched cells.

La superioridad de T30 no estuvo restringida a condiciones permisivas. En el régimen más estricto `R30_S30`, T30 mantuvo la mayor performance, con una mediana de MCC cercana a **0.88**, frente a ~0.78 para T6, ~0.62 para Roxy y ~0.39 para one-hot. De forma consistente, T30 también mostró el mejor **worst-case performance** a través de las superficies completas R×S.

Además de mayor performance, T30 produjo **landscapes más estables**. El rango mediano de MCC dentro de cada superficie fue menor para T30 que para las demás representaciones, indicando menor sensibilidad a cambios en redundancy reduction y homology separation.

La ventaja de las representaciones pretrained aumentó a medida que la evaluación se volvió más estricta, especialmente a lo largo de S. Por ejemplo, la diferencia mediana `T30−T6` aumentó aproximadamente de **+0.069 MCC en S90** a **+0.091 en S30**, mientras que la ventaja frente a one-hot aumentó de ~**+0.37** a ~**+0.46 MCC**. Esto indica que **representation quality becomes progressively more important for remote generalization**.

Frente a redundancy reduction, T30 y T6 mostraron pérdidas direccionales medianas similares, pero T30 tuvo menor sensibilidad absoluta y menos extremos. Frente a homology separation, la jerarquía de robustez fue más clara:

> **T30 > T6 > Roxy > one-hot**

con T30 mostrando la menor degradación absoluta bajo `S90→S30`.

La magnitud del beneficio de T30 fue claramente **task-dependent**. Las mayores ganancias sobre T6 aparecieron en global easy, global mixed y endolysin, mientras que `global_hard` mostró una ganancia mucho menor, consistente con un régimen que ya es altamente separable y cercano a saturación con T6. VAL mostró diferencias grandes pero debe seguir interpretándose con cautela por su bajo soporte.

La misma jerarquía se mantuvo con métricas alternativas. En tareas binarias, AP-skill reprodujo claramente `T30 > T6 > Roxy > one-hot`, y en multiclass macro-F1 mostró el mismo orden. Por tanto, la ventaja de representación no depende únicamente de MCC o del threshold de clasificación.

Finalmente, T30 obtuvo el mejor rank de representación en aproximadamente **93% de los matched cells**, estuvo dentro del top-2 prácticamente siempre y mostró el mejor comportamiento tanto promedio como worst-case. Esto confirma que su ventaja no depende de unas pocas configuraciones favorables, sino que es ampliamente consistente a través del benchmark.

### Take-home message

> **Protein representation strongly determines both predictive accuracy and generalization robustness. Pretrained ESM2 representations consistently outperform one-hot and physicochemical encodings, while ESM2-T30 achieves the strongest overall, worst-case, and remote-generalization performance. Importantly, the advantage of stronger representations increases as evaluation becomes more stringent, indicating that representation quality becomes increasingly consequential as redundant and closely homologous support is removed.**

Y la frase conceptual que yo guardaría para el paper es:

> **Better protein representations do not merely improve predictive performance; they flatten the generalization landscape.**

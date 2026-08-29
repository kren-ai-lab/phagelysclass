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

## RQ7 — Effect of increased ESM2 capacity on generalization

RQ7 mostró que aumentar la capacidad dentro de la familia ESM2, desde **T6 (~8M)** hasta **T30 (~150M)**, produce una mejora clara pero **fuertemente dependiente del régimen**. En las 560 comparaciones pareadas, la ganancia mediana fue de aproximadamente:

$$
\boxed{\Delta MCC_{T30-T6}\approx +0.077}
$$

T30 superó numéricamente a T6 en ~**93%** de los benchmark cells, fue mejor por más de `0.01 MCC` en ~**83%**, y mostró una mejora superior a `0.05 MCC` en ~**62%**. Los casos donde T6 fue mejor por más de `0.01 MCC` fueron muy raros (~**2%**) y no formaron un régimen coherente.

La ganancia de capacidad fue claramente **task-dependent**. Los mayores beneficios aparecieron en:

```text
Endolysin      ~+0.134 MCC
Global easy    ~+0.113
Global mixed   ~+0.104
```

mientras que multiclass y depolymerase mostraron mejoras moderadas y `global_hard` quedó prácticamente saturado:

```text
Global hard    ~+0.008 MCC
```

Así, RQ7 separó claramente **capacity-sensitive regimes** de tareas donde **T6 ya es prácticamente suficiente**.

La ventaja de T30 aumentó al hacer la evaluación más estricta. Con homology separation:

```text
S90   ~+0.069
S70   ~+0.075
S50   ~+0.078
S30   ~+0.091 MCC
```

y entre los corners extremos:

```text
R90_S90   ~+0.055 MCC
R30_S30   ~+0.098 MCC
```

Por tanto, **el valor adicional de la representación de mayor capacidad crece cuando la generalización requerida se vuelve más remota y menos apoyada por redundancia local**.

Con redundancy reduction también aumentó la ventaja de T30 desde R90 hacia R50, antes de alcanzar aproximadamente un plateau en R50–R30. Esto indica que la capacidad adicional se vuelve más valiosa a medida que se elimina soporte redundante, aunque el efecto no sigue una progresión estrictamente monotónica.

T30 también mejoró de forma importante el **worst-case generalization**. La elevación mediana del peor MCC de cada superficie 4×4 fue aproximadamente:

$$
\boxed{+0.109\ MCC}
$$

y las superficies T30 fueron además ligeramente más planas que las de T6. Esto indica que la capacidad adicional no solo eleva la performance promedio, sino que **sube el piso de generalización en condiciones adversas**.

T30 tendió también a ser más robusto frente a R y S, aunque este efecto fue menos universal que la ganancia directa de performance. T30 mostró menor sensibilidad absoluta a R en ~**64%** de los casos y menor sensibilidad a S en ~**65%**. Por tanto, la evidencia principal de RQ7 es el aumento de performance y worst-case robustness, no una reducción uniforme de sensibilidad.

Un resultado clave fue la fuerte relación entre performance basal de T6 y beneficio de T30:

$$
\boxed{\rho\approx -0.82}
$$

Cuanto más cerca estaba T6 de saturar una tarea, menor era la mejora producida por T30. Esto fue especialmente evidente en `global_hard`, mientras que endolysin, global easy y global mixed mostraron todavía gran headroom y grandes ganancias de capacidad. La interpretación prudente es que **capacity gains exhibit strong diminishing returns as the smaller PLM approaches task saturation**.

En las condiciones con headroom suficiente, T30 recuperó en mediana aproximadamente **42% del margen de MCC disponible**, alcanzando valores cercanos a 57–59% en endolysin y las tareas global easy/mixed.

Las métricas alternativas confirmaron la misma historia. AP-skill mostró una ganancia mediana cercana a **+0.063**, multiclass macro-F1 aproximadamente **+0.043**, y las ganancias entre métricas estuvieron fuertemente correlacionadas. Por tanto, el efecto de capacidad no depende exclusivamente de MCC.

### Regímenes identificados

**Capacity-sensitive:** global easy, global mixed, endolysin.
T30 aporta mejoras grandes y consistentes.

**Moderately capacity-sensitive:** multiclass y depolymerase.
T30 generalmente gana, pero el beneficio es menor.

**Capacity-saturated:** global hard.
T6 ya alcanza performance muy alta y T30 aporta poco.

**Evidence-limited:** VAL.
Las ganancias pueden ser grandes, pero deben interpretarse con cautela por el bajo soporte.

### Take-home message

> **Increasing ESM2 capacity from T6 to T30 provides substantial additional functional-generalization value across most benchmark regimes, but the benefit is not uniform. Gains are largest where the smaller representation retains predictive headroom and where evaluation becomes more stringent, while highly separable tasks approach saturation with T6 and show sharply diminishing returns.**

Y la frase conceptual que guardaría para el paper es:

> **Larger pretrained representations are most valuable where generalization remains representation-limited; once the task approaches saturation, additional PLM capacity yields sharply diminishing returns.**

## RQ8 — Effect of downstream algorithm choice and robustness

RQ8 mostró que el **clasificador downstream sí afecta la performance**, pero su importancia es claramente menor que la de la representación y además depende fuertemente del feature space utilizado.

A nivel global, MLP presentó la mayor MCC mediana (~**0.817**), seguido de XGB (~0.796), KNN (~0.784), RF (~0.742) y LR (~0.717). Sin embargo, las diferencias entre los tres mejores algoritmos fueron pequeñas en comparaciones pareadas: `MLP−KNN ≈ +0.007 MCC`, `MLP−XGB ≈ +0.006` y `XGB−KNN ≈ 0`. Por tanto, **MLP, XGB y KNN forman un upper tier competitivo**, sin evidencia de un algoritmo universalmente superior.

El comportamiento del algoritmo fue fuertemente **representation-dependent**. Con Roxy, XGB fue claramente competitivo; con one-hot, KNN funcionó mucho mejor que RF; con T6, MLP/XGB/KNN quedaron próximos; y con T30 prácticamente todos los clasificadores convergieron hacia performances altas. Por tanto:

> **there is no representation-independent learner hierarchy.**

Uno de los hallazgos más fuertes fue que la **importancia del algoritmo disminuye drásticamente a medida que mejora la representación**. El rango mediano entre el mejor y peor algoritmo fue aproximadamente:

```text
One-hot    ~0.306 MCC
Roxy       ~0.236
ESM2-T6    ~0.115
ESM2-T30   ~0.051
```

Así, una representación fuerte reduce de forma muy marcada la dependencia del classifier downstream. Esto indica que **better representations absorb part of the predictive burden that otherwise falls on the learning algorithm**.

LR ilustró especialmente bien este fenómeno. Frente a Roxy, modelos no lineales podían mejorar más de `+0.20 MCC` respecto de LR; con T30, esas diferencias se redujeron típicamente a solo ~`+0.02–0.04 MCC`. Esto muestra que un decision boundary simple puede volverse bastante competitivo cuando el embedding ya separa bien el problema.

En el régimen más estricto `R30_S30`, MLP mantuvo la mejor performance agregada (~**0.786 MCC**), seguido de XGB y KNN. MLP también mostró el mejor **worst-case MCC** a través de las superficies completas R×S, por lo que aparece como una opción fuerte en términos de performance absoluta y robustez.

Sin embargo, el algoritmo con menor degradación frente a R o S no fue necesariamente el de mejor performance. LR mostró las menores sensibilidades absolutas a `R90→R30` y `S90→S30`, pero partiendo desde un nivel de performance inferior. Por tanto:

> **absolute predictive quality and robustness to evaluation perturbations are distinct properties.**

Dentro de T30, MLP combinó performance alta con muy buena robustez a R y S, lo que lo vuelve una opción particularmente fuerte dentro del embedding más informativo.

El ganador retrospectivo también cambió con frecuencia a través del landscape. En ~**64%** de las superficies `task × representation`, el algoritmo ganador cambió entre `R90_S90` y `R30_S30`. Sin embargo, T30 mostró bastante más estabilidad entre corners extremos que Roxy o T6. Esto indica que la identidad del “mejor algoritmo” puede ser sensible al evaluation regime, especialmente cuando las representaciones son más débiles.

El resultado conceptual más fuerte de RQ8 fue la comparación directa entre **representation spread** y **algorithm spread**. A través de todos los `task × R × S` conditions:

```text
median algorithm spread        ~0.176 MCC
median representation spread   ~0.423 MCC
```

con un ratio mediano cercano a:

$$
\boxed{2.8\times}
$$

y la representación produjo un spread mayor que el algoritmo en **100% de las condiciones**.

Por tanto, la elección de representación fue consistentemente mucho más importante que la elección del clasificador.

AP-skill y multiclass macro-F1 reprodujeron la misma historia: las diferencias entre MLP, XGB y KNN fueron pequeñas, mientras RF y LR mostraron mayor dependencia del tipo de representación.

### Take-home message

> **Downstream learner choice matters, but its effect is strongly representation-dependent and substantially smaller than the effect of representation itself. MLP, XGB, and KNN form a broadly competitive upper tier, while stronger protein representations sharply reduce performance differences among classifiers. Across all benchmark regimes, choosing how proteins are represented is more consequential than choosing the downstream learner.**

Y la frase conceptual que guardaría para el paper:

> **A strong representation reduces the burden on the downstream learner; what the model sees matters substantially more than which classifier partitions the resulting feature space.**

## RQ9 — Functional generalization regimes across biological tasks

RQ9 mostró que **no existe un único régimen de generalización para protein-function prediction**. Cada formulación biológica responde de forma distinta a redundancy control, homology separation, representation quality, PLM capacity, classifier choice y evidencia disponible.

El hallazgo más claro fue que **global easy y global mixed ocupan prácticamente el mismo régimen**. Sus fingerprints fueron muy similares en MCC, sensibilidad a R y S, dependencia de representación, T30−T6, algorithm spread y worst-case performance. Esto confirma que agregar el background lysis-adjacent al broad background no cambia sustancialmente la estructura del problema predictivo.

En contraste, **global hard forma un régimen completamente distinto**. Presenta MCC alto (~0.90), worst-case elevado (~0.85), baja sensibilidad a R y S, poca dependencia de representación adicional y un gain T30−T6 casi nulo (~+0.008 MCC). Por tanto, representa un régimen:

> **high-floor, protocol-robust and largely representation-saturated.**

Esto demuestra que la definición del background negativo no solo altera la performance absoluta, sino que puede cambiar completamente **qué limita la generalización**.

### Endolysin: redundancy- and representation-limited

Endolysin mostró uno de los landscapes más dependientes de decisiones upstream:

* `|R| ≈ 0.143 MCC`;
* `T30−T6 ≈ +0.134`;
* `T30−one-hot ≈ +0.51`;
* representation range ~0.51;
* worst-case MCC ~0.56.

Por tanto:

> **Endolysin prediction is strongly redundancy-sensitive and representation-limited.**

La tarea depende sustancialmente tanto del soporte redundante como de la calidad de la representación.

### Depolymerase: high-floor but interaction-dependent

Depolymerase presentó una firma muy diferente:

* MCC global ~0.84;
* strict MCC ~0.85;
* worst-case ~0.79;
* poca degradación neta con R o S;
* `T30−T6 ≈ +0.02`;
* interacción R×S relativamente importante;
* mayor dependencia del classifier que otras tareas de alto desempeño.

Por tanto, DPO parece un régimen:

> **high-floor, comparatively capacity-saturated, but interaction- and learner-dependent.**

No hay evidencia para interpretar su leve mejora nativa bajo R más estricto como un beneficio biológico de eliminar redundancia.

### VAL: evidence-limited and unstable

VAL fue el régimen más inestable del benchmark:

* worst-case MCC ~**0.19**;
* surface range ~**0.45**;
* `|S| ≈ 0.078`;
* extreme R×S interaction ~**0.34**;
* T30−T6 ~+0.12;
* representation range ~**0.71**;
* algorithm range ~**0.32**.

Sin embargo, esta aparente sensibilidad extrema debe interpretarse principalmente como:

> **an evidence-limited generalization regime.**

Su bajo soporte positivo puede amplificar los efectos aparentes de R, S, representation y algorithm. Por eso RQ10 será necesaria antes de hacer claims fuertes sobre VAL.

### Multiclass: stable in magnitude, complex in geometry

Multiclass mostró:

* MCC global ~0.81;
* strict MCC ~0.82;
* worst-case ~0.74;
* efectos netos pequeños de R y S;
* T30−T6 ~+0.03;
* representation effects todavía mayores que algorithm effects.

Aunque RQ2–RQ4 mostraron mucha non-monotonicity e interacción local, el efecto extremo neto fue relativamente pequeño.

Así, multiclass parece:

> **a comparatively stable but locally non-monotonic discrimination regime.**

Esto demuestra que una superficie puede tener geometría compleja sin presentar una gran degradación extrema.

### Los subtipos no constituyen un único régimen

RQ9 dejó claro que “subtype prediction” no es una categoría metodológicamente homogénea:

* **Endolysin:** redundancy- and representation-limited.
* **Depolymerase:** high-floor and interaction-dependent.
* **VAL:** evidence-limited and highly unstable.

Por tanto, la función biológica concreta determina qué factor limita la generalización.

### Representation sigue dominando dentro de cada task

En todas las tareas, el spread inducido por cambiar representación fue mayor que el inducido por cambiar classifier.

El ratio `representation range / algorithm range` fue aproximadamente:

* global hard: **5.2×**
* global easy: **3.6×**
* endolysin: **3.0×**
* multiclass: **2.9×**
* global mixed: **2.9×**
* VAL: **2.2×**
* depolymerase: **1.6×**

Incluso en los regímenes donde el learner importa, **representation choice sigue siendo el factor model-side dominante**.

### AP skill confirma los mismos regímenes

Las sensibilidades con AP skill reprodujeron esencialmente los fingerprints obtenidos con MCC:

* endolysin y global easy/mixed siguen siendo R-sensitive;
* global hard sigue siendo muy estable;
* VAL sigue mostrando la mayor sensibilidad a S;
* DPO sigue mostrando poca degradación neta.

Por tanto, los functional regimes no son artefactos de una sola métrica.

## Take-home message

> **Protein-function prediction does not exhibit a single generalization regime. Different biological task formulations fail for different reasons: some are limited primarily by redundant evidence, others by representation capacity, others by non-additive evaluation effects, and low-support tasks additionally by evidence scarcity.**

Y la frase conceptual más fuerte para guardar:

> **Generalization difficulty is a property of the biological learning problem itself, not merely of the model used to solve it.**

## RQ10 — VAL under Low-N: evidence scarcity and metric fragility

RQ10 confirmó que **VAL constituye un régimen genuinamente evidence-limited**. Tras el control de redundancia y la separación homológica, el número de positivos disponibles queda extremadamente reducido:

* **Train positives:** ~31–75
* **Test positives:** ~14–32
* **Positive test homology clusters:** ~12–26

Por tanto, el problema no es solo que haya pocas secuencias, sino que también existen **pocos componentes homólogos positivos independientes**.

### La evaluación queda pobremente resuelta

Con solo 14–32 positivos en test, un único VAL correctamente o incorrectamente clasificado puede cambiar recall en aproximadamente:

$$
\boxed{3.1\%-7.1\%}
$$

Esto significa que diferencias aparentemente grandes entre modelos pueden depender de apenas uno o dos proteins.

La conclusión conceptual más importante es:

> **Low-N constrains not only what can be learned, but also how precisely generalization can be measured.**

### Learning scarcity y evaluation scarcity son problemas distintos

La reducción de R disminuye directamente la cantidad de evidencia VAL disponible para aprendizaje. Por otro lado, S puede redistribuir pequeños grupos homólogos entre train y test, produciendo cambios importantes en el número de positivos evaluados aun cuando el tamaño global del dataset cambie muy poco.

El soporte positivo de test mostró una asociación relativamente fuerte con la performance observada, por ejemplo:

$$
\rho(n_{test+}, MCC_{median})\approx0.67
$$

y asociaciones también claras con recall, balanced accuracy, AP y ROC-AUC.

Estas relaciones son descriptivas, no causales, porque al cambiar S también cambia la composición del test.

### La escasez afecta la estabilidad de las métricas de forma metric-dependent

La volatilidad no puede reducirse únicamente al threshold de clasificación. Ranking metrics también mostraron sensibilidad al bajo soporte, especialmente ROC-AUC.

Por tanto:

> **Evidence scarcity affects different evaluation metrics differently, rather than producing a uniform increase in variance.**

### Existe prediction fragility, no solo metric fragility

El análisis por secuencia mostró que los modelos frecuentemente **no coinciden en qué VAL proteins son reconocibles**.

La mediana del fraction-of-models vote para positivos fue cercana a **0.40**, y aproximadamente **65% de las exposiciones positivas** fueron `model-dependent`.

En algunas condiciones estrictas prácticamente no hubo positivos reconocidos de manera estable por ≥80% de los modelos, mientras una fracción importante aparecía como `stable miss`.

Esto demuestra dos fenómenos distintos:

> **Metric fragility:** pocos positivos hacen que cada protein tenga un peso enorme en la métrica.

> **Prediction fragility:** diferentes representaciones y classifiers reconocen diferentes subsets de esos pocos positivos.

### Algunos VAL proteins son persistentemente difíciles

Las secuencias positivas que reaparecen en varios test regimes muestran considerable desacuerdo entre modelos, con algunos proteins sistemáticamente mal reconocidos.

Esto indica que la incertidumbre no está distribuida uniformemente entre todos los positivos.

### Representaciones fuertes mantienen señal incluso en Low-N

A pesar de la escasez extrema, ESM2-T30 mantuvo una ventaja muy fuerte:

$$
\boxed{T30-T6\approx+0.12\ MCC}
$$

$$
\boxed{T30-onehot\approx+0.69\ MCC}
$$

En algunas celdas con apenas ~14 positivos, one-hot podía fallar prácticamente todos los VAL, mientras T30 recuperaba una fracción importante de ellos.

Por tanto:

> **Evidence scarcity increases uncertainty without eliminating learnability.**

Sin embargo, T30 **no resuelve el problema Low-N**. Su ventaja no aumentó sistemáticamente a medida que disminuyó el soporte, por lo que no hay evidencia para afirmar que el PLM compensa proporcionalmente la escasez.

### Low-N tampoco explica por sí solo la dependencia del classifier

No se observó una relación fuerte entre positive support y algorithm spread. Así, la variabilidad entre modelos VAL parece surgir de una combinación de:

* evidence scarcity;
* representation;
* learner;
* test composition;
* homology topology.

No puede atribuirse únicamente a `n`.

### VAL es mucho más escaso que los otros subtypes

En mediana:

| Task         | Train + |  Test + |
| ------------ | ------: | ------: |
| Endolysin    |   ~1668 |    ~680 |
| Depolymerase |    ~277 |    ~117 |
| **VAL**      | **~47** | **~18** |

Así, VAL posee aproximadamente **6–7× menos positivos test que depolymerase** y casi **40× menos que endolysin**.

Esto respalda claramente su clasificación como **evidence-limited functional regime**.

### El problema también aparece dentro de multiclass

VAL mantiene soporte muy bajo dentro de `multiclass_clean`, aproximadamente **14–36 positivos test**.

En los regímenes más estrictos, VAL recall y F1 disminuyen claramente aun cuando el MCC multiclass global sigue siendo relativamente bueno.

Esto demuestra que:

> **Aggregate multiclass performance can conceal severe uncertainty in a low-support functional class.**

## Take-home message

> **VAL prediction occupies a distinct evidence-limited generalization regime. Stringent sequence control leaves only a few dozen positive proteins and homology components for both learning and evaluation, causing individual proteins to exert disproportionate influence on reported metrics and producing substantial disagreement across model families. Nevertheless, strong pretrained representations retain meaningful predictive signal, showing that evidence scarcity increases uncertainty and measurement fragility without eliminating learnability.**

Y la frase conceptual que yo guardaría para el paper es:

> **Sparse biological evidence constrains not only what a model can learn, but how precisely its generalization can be measured.**

## RQ11 — Binary subtype recognition vs multiclass functional discrimination

RQ11 mostró que **la formulación de tarea cambia de manera profunda el problema de generalización**, incluso cuando las clases biológicas nominales son las mismas. El efecto no tiene una dirección universal: multiclass favorece fuertemente a endolysin, mejora moderadamente DPO y perjudica claramente a VAL.

### Principales hallazgos

* **Endolysin**

  * Binary F1 ≈ **0.719**
  * Multiclass F1 ≈ **0.969**
  * Δ ≈ **+0.252**
  * Multiclass fue mejor en prácticamente el **100%** de las condiciones.
  * La mejora se explica por aumentos tanto de recall como de precision.

* **Depolymerase**

  * Binary F1 ≈ **0.843**
  * Multiclass F1 ≈ **0.900**
  * Δ ≈ **+0.058**
  * También se beneficia del multiclass, pero de manera mucho más moderada.

* **VAL**

  * Binary F1 ≈ **0.400**
  * Multiclass F1 ≈ **0.299**
  * Δ ≈ **−0.112**
  * El deterioro está impulsado especialmente por una caída fuerte en precision.
  * Esto muestra que VAL puede separarse del background no lítico, pero resulta mucho más difícil distinguirlo de otras funciones líticas.

### El comparison universe explica gran parte del efecto

El número de positivos por subtipo es similar entre formulaciones, pero cambia radicalmente **contra qué proteínas se realiza la decisión**.

Binary pregunta:

> “¿Es este protein compatible con este subtype frente a un background no lítico?”

Multiclass pregunta:

> “¿Cuál de estas funciones líticas mutuamente exclusivas describe mejor este protein?”

Por tanto, el cambio de formulación modifica el decision problem, no simplemente el número de clases.

### VAL se confunde principalmente con depolymerase

La matriz de confusión multiclass mostró aproximadamente:

* True DPO → DPO: **~98%**
* True Endolysin → Endolysin: **~89%**
* True VAL → VAL: **~25%**
* True VAL → DPO: **~64%**

Así, el principal cuello de botella funcional es:

> **VAL vs depolymerase discrimination.**

VAL prácticamente no se confunde con endolysin.

### T30 mejora la discriminación, pero no elimina el problema

Para VAL multiclass:

* one-hot: reconocimiento prácticamente nulo;
* Roxy: ~22%;
* T6: ~36%;
* T30: ~47%.

La confusión VAL→DPO disminuye con representaciones más fuertes, pero incluso T30 deja aproximadamente la mitad de los VAL asignados como DPO.

Por tanto:

> **better representations improve functional resolution, but do not remove the intrinsic formulation difficulty of VAL.**

### El resultado persiste incluso con T30 y bajo R30/S30

Con T30:

* Endolysin: Δ multiclass−binary ≈ **+0.060**
* DPO: ≈ **+0.047**
* VAL: ≈ **−0.231**

Y en `R30_S30`:

* Endolysin: ≈ **+0.086**
* DPO: ≈ **+0.046**
* VAL: ≈ **−0.280**

Por tanto, la penalización multiclass de VAL no puede explicarse simplemente por una representación débil o por redundancy/homology leakage.

### Task formulation cambia también el generalization fingerprint

El cambio binary→multiclass no solo mueve F1; cambia qué factores limitan cada tarea.

El caso más claro es **endolysin**:

* Binary:

  * |R| ≈ 0.136
  * representation range ≈ 0.467
  * algorithm range ≈ 0.151
  * T30−T6 ≈ 0.119

* Multiclass:

  * |R| ≈ 0.019
  * representation range ≈ 0.036
  * algorithm range ≈ 0.013
  * T30−T6 ≈ 0.004

Así, el mismo subtipo pasa de ser **redundancy- and representation-sensitive** a un régimen prácticamente **saturated and robust** simplemente al cambiar el decision problem.

### Los rankings de modelos sí se transfieren bastante bien

Aunque performance y sensibilidad cambian, las configuraciones fuertes suelen seguir siendo fuertes:

* Endolysin: Spearman ≈ **0.89**
* DPO: ≈ **0.87**
* VAL: ≈ **0.86**

Dentro de cada celda R×S, las correlaciones son incluso mayores.

Esto indica:

> **task formulation changes absolute difficulty and sensitivity, but does not completely reshuffle which model configurations are strong.**

### Un buen multiclass MCC puede esconder una clase mal resuelta

A través de las familias multiclass:

* median global MCC ≈ **0.81**
* median worst-class F1 ≈ **0.30**
* median best-class F1 ≈ **0.97**
* class-F1 range ≈ **0.65**

Por tanto, un buen metric global puede coexistir con una clase funcional pobremente resuelta, especialmente VAL.

## Take-home message

> **Binary subtype recognition and multiclass functional discrimination are not interchangeable formulations of the same problem. The formulation changes both absolute performance and the entire generalization landscape, with strongly class-dependent effects.**

Y la frase conceptual más fuerte:

> **A protein may be easy to recognize as compatible with a functional category yet difficult to resolve against biologically related alternatives.**

O, más general:

> **A protein-function benchmark is defined not only by its labels, but by the biological alternatives those labels are required to discriminate.**

## RQ12 — Positive–unlabeled learning paradigm

RQ12 mostró que **el principal efecto del paradigma PU no es necesariamente mejorar la recuperación de positivos ya conocidos, sino cambiar cómo se organiza y prioriza el espacio de proteínas no anotadas**.

### Principales hallazgos

* **Naive P-vs-U y BaggingPU recuperan positivos conocidos de forma casi equivalente.**

  * P ≥ U90: ~0.868 vs ~0.867
  * P ≥ U95: ~0.805 vs ~0.798
  * P ≥ U99: ~0.584 vs ~0.576
  * Los deltas medianos son mínimos y gran parte de las comparaciones caen dentro de una zona de practical tie.

* Por tanto:

> **BaggingPU does not provide a systematic benchmark-wide improvement in known-positive recovery relative to naive P-vs-U.**

### BaggingPU sí cambia de manera real la priorización de U

Aunque el recovery agregado sea parecido, BaggingPU reorganiza sustancialmente el ranking de proteínas unlabeled:

* U-score Spearman vs naive: **~0.875**
* tie-aware top-1% Jaccard: **~0.63**
* tie-aware top-5% Jaccard: **~0.60**

Así, aproximadamente **40% de los candidatos del extremo superior cambian** entre naive y BaggingPU.

Esto demuestra:

> **Similar aggregate recovery does not imply similar discovery behavior.**

La elección del paradigma PU puede cambiar materialmente qué proteínas se seleccionarían para validación experimental.

### Elkan–Noto actúa principalmente como score correction

Elkan–Noto conserva prácticamente por completo el ordenamiento naive:

* U-score Spearman: **~0.999995**
* top-5% tie-aware overlap: **1.0**
* top-1% overlap coefficient: **1.0**

Esto indica que EN **no produce una repriorización sustancial del proteoma U**, sino que principalmente transforma la escala de scores mediante propensity correction.

### El aparente gain U99 de Elkan–Noto está dominado por clipping y ties

Elkan–Noto muestra aparentemente mayor P ≥ U99 (~0.683), pero la auditoría tie-aware reveló:

* median U99 threshold = **1.0**
* ~**2.2%** de los U quedan exactamente en score 1.0
* ~**66%** de los P quedan exactamente en 1.0
* strict P > U99 = **0**
* ~**82%** del recovery U99 mediano depende de ties en el boundary

Por tanto:

> **Elkan–Noto's apparent U99 advantage reflects ceiling compression and boundary ties rather than improved resolved ranking above the extreme U tail.**

### Elkan–Noto pierde resolución en el extremo superior

En top-1%, el pool tie-inclusive de EN contiene aproximadamente **2.2%** de los U, porque muchos quedan empatados en 1.0.

El top-1% naive queda prácticamente contenido dentro de ese plateau de EN.

Así, EN no selecciona necesariamente candidatos distintos; **reduce la resolución entre los candidatos más altos**.

### BaggingPU es internamente bastante estable

Across its 30 internal bags:

* median per-protein score SD ≈ **0.020**
* mean score SD ≈ **0.051**

Pero mayor variabilidad entre bags se asocia con menor positive recovery:

$$
\rho \approx -0.70
$$

Esto sugiere que cuando la predicción depende fuertemente de qué U proteins son temporalmente tratados como negativos, la recuperación de positivos también se vuelve menos consistente.

### No existe un PU method winner universal

Los tres paradigmas hacen cosas distintas:

* **Naive:** baseline simple y estable.
* **BaggingPU:** recovery similar, pero candidate ranking realmente distinto.
* **Elkan–Noto:** ranking casi idéntico al naive, pero score scale corregida y comprimida en el extremo superior.

Por tanto, no hay soporte para una jerarquía simple tipo:

$$
EN > BaggingPU > Naive
$$

La elección depende del objetivo del sistema.

## Take-home message

> **Explicit PU learning does not necessarily improve the recovery of already annotated positives, but it can substantially alter how unlabeled proteins are prioritized for discovery. BaggingPU genuinely reorders the candidate space, whereas Elkan–Noto largely preserves the naive ranking while changing score scale and creating upper-tail ties through propensity correction.**

Y la frase conceptual que guardaría:

> **When annotation is incomplete, the choice of learning paradigm may matter less for recovering what is already known than for deciding what should be investigated next.**

Claro bebé. Para tus notas, dejaría **RQ13** así:

## RQ13 — PU generalization across redundancy and homology constraints

RQ13 mostró que, incluso bajo aprendizaje positive–unlabeled, **la reducción de redundancia (`R`) domina claramente el generalization landscape**, mientras que la separación homológica (`S`) tiene un efecto bastante menor y más irregular.

### Principales hallazgos

* **Efecto de R90→R30 sobre U95 recovery**

  * Δ mediana ≈ **−0.073**
  * ~**96%** de los contrastes son negativos
  * efecto absoluto mediano ≈ **0.073**

* **Efecto de S90→S30**

  * Δ mediana ≈ **−0.008**
  * ~**72%** de los contrastes son negativos
  * efecto absoluto mediano ≈ **0.012**

Así, el efecto típico de `R` es aproximadamente **6 veces mayor** que el de `S`.

> **PU recovery is primarily redundancy-sensitive rather than homology-separation-sensitive.**

### La pérdida por R es progresiva

Los cambios stepwise fueron aproximadamente:

* R90→R70: **−0.026**
* R70→R50: **−0.028**
* R50→R30: **−0.018**

No aparece un breakpoint universal. La mayoría de las trayectorias R muestran degradación monotónica o casi monotónica.

Por tanto:

> **Removing redundant sequence evidence progressively reduces PU positive recovery rather than producing a single catastrophic threshold transition.**

### El efecto de S es pequeño y principalmente aparece en el régimen más estricto

Stepwise:

* S90→S70: ~0
* S70→S50: pequeño efecto negativo
* S50→S30: **~−0.009**

Además, ~**81%** de las curvas S son nonmonotonic.

Así, la separación homológica tiene un efecto mucho menos uniforme que la reducción de redundancia.

### R domina casi todas las superficies PU

En el análisis de dominance:

* **~90%** de las superficies fueron `R-dominant`
* ~10% balanced
* prácticamente **ninguna S-dominant**
* median dominance index ≈ **0.62**

Esto confirma que:

> **Redundancy is the dominant axis controlling PU generalization.**

### Existe interacción R×S, pero no supera el efecto principal de R

La interacción extrema mediana fue aproximadamente:

$$
I_{R\times S}\approx+0.027
$$

y fue positiva en ~90% de las superficies.

Esto indica cierta **non-additivity / attenuation**: cuando R ya es muy estricto, hacer también S estricto suele añadir menos pérdida de la esperada bajo independencia.

Pero la historia principal sigue siendo `R >> S`.

### Los tres paradigmas PU comparten la misma geometría general

Extreme R effect:

* naive: ~**−0.083**
* BaggingPU: ~**−0.066**
* Elkan–Noto: ~**−0.083**

Extreme S effect:

* naive: ~**−0.008**
* BaggingPU: ~**−0.007**
* Elkan–Noto: ~**−0.007**

Así, cambiar el paradigma PU no altera radicalmente la estructura del R×S landscape.

### BaggingPU sigue siendo casi equivalente al naive en promedio

Median gain U95 vs naive:

$$
\boxed{+0.001}
$$

y ~**76%** de las comparaciones están dentro de ±0.02.

Por tanto, no hay una ventaja global grande de BaggingPU para known-positive recovery.

### Pero BaggingPU se vuelve modestamente más favorable cuando R es estricto

Gain vs naive por R:

* R90: ~0
* R70: ~0
* R50: ~+0.003
* **R30: ~+0.009**

El cambio de gain R30−R90 es ~**+0.008**, positivo en ~75% de los contrastes.

Esto sugiere:

> **BaggingPU becomes modestly more useful as redundant sequence evidence is removed.**

El efecto es real pero pequeño, por lo que no debe presentarse como una gran rescue effect.

### S no aumenta de forma clara el beneficio de BaggingPU

El cambio de gain S90→S30 es solo ~**+0.002**.

Así, no hay evidencia fuerte de que BaggingPU sea especialmente útil para remote-homology generalization.

### Elkan–Noto permanece prácticamente rank-equivalent al naive

Para U95:

* median gain vs naive ≈ **0**
* ~**88%** practical ties
* no cambio relevante del gain con R
* no cambio relevante con S
* prácticamente ninguna interacción del gain

Esto refuerza RQ12: Elkan–Noto actúa principalmente como score correction y no como un cambio sustancial del ranking.

### U99 vuelve a estar dominado por ties de Elkan–Noto

Canonical U99 muestra una aparente ventaja de EN, pero strict U99 cambia completamente la interpretación.

En todas las R/S cells, el U99 threshold de EN alcanza **1.0**, debido al clipping.

Por tanto:

> **Canonical U99 cannot be interpreted as evidence of improved resolved ranking for Elkan–Noto without the tie-aware correction.**

### La reducción de R también reduce fuertemente el soporte positivo

Al pasar R90→R30, el positive-test support cae aproximadamente dos tercios.

Esto se asocia fuertemente con recovery:

$$
\rho(n_P,\mathrm{recovery})\approx0.90-0.93
$$

pero esta relación es descriptiva, porque R cambia simultáneamente soporte, composición, redundancia y topología homológica.

## Take-home message

> **Positive–unlabeled recovery is governed much more strongly by redundancy reduction than by train–test homology separation. Dedicated PU methods do not fundamentally alter this landscape: Elkan–Noto remains essentially equivalent to naive ranking, whereas BaggingPU shows only a modest increase in relative benefit when redundancy control becomes stringent.**

Y la frase conceptual que guardaría:

> **Label uncertainty does not erase the effect of sequence redundancy; even under PU learning, redundancy remains the dominant determinant of the generalization landscape.**

Otra buena para el paper:

> **Uncertainty-aware learning can modestly mitigate, but cannot replace, the information lost when redundant sequence evidence is removed.**

## RQ14 — Effect of the unlabeled-universe definition

RQ14 mostró que ampliar el universo unlabeled desde `U_explicit` hacia un background broad que incorpora `U_other_unresolved` modifica sustancialmente el problema PU, pero principalmente a través de la estructura del espacio de candidatos y no mediante una degradación fuerte del problema original `P_core vs U_explicit`.

En la evaluación native, broad redujo consistentemente el positive recovery a U95, con medianas de aproximadamente **−0.054 para naive**, **−0.056 para BaggingPU** y **−0.047 para Elkan–Noto**. El efecto fue negativo en ~95% de las comparaciones. Esta pérdida persistió al evaluar los mismos positivos bajo las respectivas distribuciones U nativas, descartando que se explicara únicamente por cambios en qué proteínas positivas quedaron asignadas al outer test.

Sin embargo, el análisis common-reference cambió sustancialmente la interpretación. Al evaluar exactamente los mismos `P_core` contra exactamente los mismos `U_explicit`, la penalización broad disminuyó a aproximadamente **−0.009 U95 para naive**, **0 para BaggingPU** y **−0.009 para Elkan–Noto**. Así, aproximadamente 80–100% de la pérdida native desapareció al controlar completamente las entidades utilizadas como referencia. Los percentiles de los positivos sobre el mismo `U_explicit` también permanecieron prácticamente invariantes.

Por tanto, la menor recuperación observada bajo broad no refleja principalmente una pérdida de capacidad para resolver el problema original `P_core vs U_explicit`. En cambio, broad incorpora una población adicional de proteínas unresolved que compite fuertemente en el extremo superior de U.

`U_other_unresolved` representó en mediana aproximadamente **28.8% del broad-U**, pero constituyó ~**53% del top-5% con naive**, ~**58% con BaggingPU** y ~**53% con Elkan–Noto**, correspondientes a enriquecimientos de aproximadamente **1.90×**, **2.07×** y **1.90×**, respectivamente. Además, ~13–16% de estas proteínas superaron el percentil 95 de `U_explicit`. El enriquecimiento relativo aumentó hacia regímenes de redundancy reduction más estrictos, alcanzando aproximadamente **2.08× para naive** y **2.29× para BaggingPU en R30**.

Ampliar U también reorganizó materialmente la priorización incluso entre los mismos `U_explicit`. La correlación de rankings conservative–broad fue ~**0.84 para naive/Elkan–Noto** y ~**0.93 para BaggingPU**, pero el Jaccard del top-1% fue solo ~**0.46**, ~**0.53** y ~**0.48**, respectivamente. Por tanto, una parte importante de las proteínas seleccionadas para follow-up depende de cómo se define el universo de incertidumbre.

BaggingPU mostró una estabilidad relativamente mayor del ranking de `U_explicit` y prácticamente ninguna pérdida common-reference, mientras produjo el mayor enriquecimiento de `U_other_unresolved` entre los candidatos superiores. Sin embargo, ampliar U no aumentó de forma importante la ventaja global de los métodos PU dedicados sobre naive P-vs-U para recuperación de positivos conocidos.

Elkan–Noto mantuvo los problemas de clipping y ties previamente identificados. Su comportamiento canonical U99 no debe interpretarse como una mejora de ranking sin considerar el contraste strict-U99.

### Take-home message

> **Expanding the unlabeled universe does not substantially impair recognition of the original positive-versus-explicit-U problem. Instead, it reveals an additional unresolved protein population that is strongly enriched in the candidate-like tail and materially changes which proteins would be prioritized for discovery.**

La frase conceptual que guardaría para el paper es:

> **In positive–unlabeled protein learning, the unknown class is not a passive background: its biological scope determines which alternatives compete with known positives and, consequently, which unannotated proteins emerge as candidates for experimental investigation.**

## RQ15 — Representation × positive–unlabeled learning paradigm

RQ15 mostró que la interacción entre representación y paradigma PU depende fuertemente del tipo de representación y, de manera crucial, afecta de forma diferente el reconocimiento de positivos conocidos y la priorización de proteínas unlabeled.

Para known-positive recovery, BaggingPU mostró una ganancia sustancial únicamente con one-hot. El gain mediano `BaggingPU−naive` en U95 fue aproximadamente **+0.093** para one-hot, con ~83% de los contrastes positivos. En cambio, la ganancia fue esencialmente nula para Roxy (**~0.000**), ESM2-T6 (**~−0.003**) y ESM2-T30 (**~0.000**), con 88–95% de las comparaciones dentro de la zona de practical tie. Elkan–Noto tampoco produjo una ganancia mediana relevante en ninguna representación.

Por tanto, la dependencia del recovery respecto del paradigma PU no disminuye de forma monotónica con la calidad de representación. En cambio, **one-hot constituye un régimen específicamente sensible al paradigma PU**, mientras Roxy y ambas representaciones ESM2 muestran poca dependencia del learner PU para recuperación de positivos conocidos.

La representación continuó siendo el determinante dominante del recovery. La ganancia `T30−T6` permaneció cercana a **+0.07–0.08 U95** bajo naive, BaggingPU y Elkan–Noto, mientras `T30−Roxy` fue aproximadamente **+0.21** bajo los tres paradigmas. Así, cambiar el paradigma PU no modifica sustancialmente la ventaja proporcionada por una representación pretrained fuerte.

El spread entre paradigmas fue aproximadamente **0.118 U95 para one-hot**, pero solo ~**0.005–0.007** para Roxy, T6 y T30. BaggingPU redujo el representation spread global principalmente mediante su mejora de one-hot, pero no eliminó la gran ventaja de las representaciones ESM2.

La interacción con R y S fue pequeña para las representaciones pretrained. La diferencia de PU gain `T30−T6` aumentó modestamente hacia redundancy regimes estrictos, desde prácticamente cero en R90 hasta ~**+0.009** en R30, pero permaneció pequeña en términos absolutos. No apareció una tendencia clara con S. Por tanto, no existe evidencia de una fuerte interacción adicional `representation × PU paradigm × R/S` para known-positive recovery.

La historia cambió al analizar candidate prioritization. Aunque BaggingPU prácticamente no mejoró U95 recovery con T30, produjo un ranking del universo U sustancialmente diferente del naive. La correlación naive–BaggingPU disminuyó desde ~**0.91 con T6** hasta ~**0.81 con T30**, y el Jaccard tie-aware del top-5% cayó desde ~**0.70 hasta ~0.53**. En comparaciones pareadas, `T30−T6` redujo el Spearman en aproximadamente **−0.083** y el top-5% Jaccard en ~**−0.136**, con reducción del top-5% overlap en todas las 64 comparaciones.

Esta divergencia no parece originarse en inestabilidad de BaggingPU. La variabilidad mediana de scores entre bags disminuyó desde ~**0.015 con T6** hasta ~**0.0068 con T30**, y T30 mostró menor median bag variability en ~**98%** de los contrastes pareados. Por tanto, T30 produce predicciones BaggingPU internamente más estables al mismo tiempo que genera una repriorización más distinta respecto del naive.

El análisis broad-U reforzó esta separación entre recognition y discovery. BaggingPU aumentó el enriquecimiento de `U_other_unresolved` en el top-5% respecto del naive en aproximadamente **+0.206 con T30**, frente a ~**+0.072 con T6** y ~**+0.014 con Roxy**. One-hot mostró el mayor efecto (~**+0.245**), aunque también presentó la mayor variabilidad interna entre bags. Así, la interacción representación–PU en candidate discovery es real pero no monotónica con representation quality.

Elkan–Noto permaneció prácticamente rank-equivalent al naive bajo todas las representaciones. Los Spearman fueron esencialmente 1 y el top-5% tie-aware overlap fue completo. Las diferencias observadas en top-1% para Roxy y one-hot se explicaron por expansión del candidate pool debido a clipping/ties, un efecto que prácticamente desapareció con T30.

### Take-home message

> **Protein representation and PU learning are largely separable for recovery of known positives once informative feature spaces are available, with one-hot forming a distinct PU-sensitive regime. In contrast, representation and learning paradigm remain strongly coupled in the discovery space: BaggingPU can produce stable but substantially different candidate rankings even when it provides essentially no improvement in known-positive recovery.**

La frase conceptual que guardaría para el paper es:

> **A strong protein representation can reduce the need for uncertainty-aware learning to recover what is already known without eliminating its influence on what the model chooses to investigate next.**

De forma compacta, RQ15 identifica un régimen:

> **recovery-substitutive, discovery-coupled.**

## RQ16 — Common-test sensitivity of redundancy and homology effects

RQ16 mostró que los efectos de redundancy reduction (`R`) y homology separation (`S`) presentan sensibilidades muy diferentes a la composición del outer test.

El análisis common-test fue validado recalculando MCC, AP-skill y macro-F1 directamente desde los 2,240 prediction files, reproduciendo los resultados Stage12A con diferencias numéricas únicamente al nivel de precisión de máquina. Las intersecciones pairwise conservaron aproximadamente **28% del test para contrastes R** y **30% para contrastes S**, proporcionando una base razonable para evaluar los mismos proteins bajo dos condiciones experimentales.

### Homology-separation effects are largely preserved on shared proteins

Para `S90→S30`, la pérdida mediana global disminuyó solo moderadamente desde aproximadamente **−0.017 MCC en los native tests** hasta **−0.012 MCC en los common tests**. AP-skill mostró la misma conclusión, pasando de aproximadamente **−0.013 a −0.009**. La dirección del efecto MCC coincidió en ~73% de los contrastes extremos.

Más importante aún, la dependencia del efecto S respecto de R se preservó claramente cuando se evaluaron las mismas proteínas:

```text
             Native ΔMCC      Common-test ΔMCC
R90             −0.047             −0.037
R70             −0.031             −0.020
R50             −0.006             −0.004
R30              ~0                 ~0
```

Por tanto, la desaparición progresiva del penalty de homology separation cuando redundancy control ya es estricto **no puede explicarse simplemente por cambios en qué proteínas entran al test set**.

A nivel de tarea, global easy mantuvo prácticamente intacto su efecto S, mientras global mixed, multiclass y endolysin conservaron penalties negativos algo atenuados. Global hard permaneció prácticamente insensible a S. Depolymerase también mostró un pequeño common-test penalty. VAL no puede interpretarse de forma fuerte debido al reducido número de positivos compartidos.

### Redundancy effects are much more test-composition-sensitive

El contraste `R90→R30` cambió de aproximadamente:

$$
\Delta MCC_{\mathrm{native}}\approx-0.070
$$

a:

$$
\Delta MCC_{\mathrm{common}}\approx-0.012.
$$

A través de todos los contrastes R, la magnitud absoluta mediana disminuyó aproximadamente de **0.043 a 0.018 MCC**, el sign agreement fue ~59%, y en ~72% de las comparaciones el common-test effect fue menor que el native effect.

Así, una fracción importante del efecto R del benchmark no refleja simplemente que exactamente las mismas proteínas se vuelvan más difíciles de predecir. Redundancy reduction modifica el conjunto de evidencia disponible y la población biológica que sobrevive al control de redundancia y alcanza la evaluación.

Este resultado no invalida la dominancia de R observada previamente. En cambio, clarifica su significado:

> **R is the dominant driver of benchmark-wide predictive change because it restructures both the available evidence and the evaluated sequence population, whereas S primarily modifies the train–test generalization relationship within a fixed biological universe.**

Endolysin permaneció como el task más claramente redundancy-sensitive incluso bajo common-test evaluation, aunque su pérdida `R90→R30` disminuyó desde ~**−0.143** hasta ~**−0.032 MCC**. Global easy y global mixed también conservaron pequeñas pérdidas negativas.

En depolymerase, la ligera mejora native observada al reducir R desapareció bajo common-test evaluation y se transformó en un efecto cercano a cero o ligeramente negativo. Por tanto, no existe evidencia robusta de que redundancy reduction mejore intrínsecamente la predicción de depolymerase.

### Exact four-level and four-corner intersections are too restrictive for strong inference

Las intersecciones comunes a los cuatro niveles retuvieron únicamente ~**2–3% del outer test**, y el four-corner `R90S90/R90S30/R30S90/R30S30` retuvo una mediana de apenas ~**74 proteínas**. En algunos tasks, especialmente multiclass, DPO y VAL, el soporte por clase quedó extremadamente bajo o no permitió estimar métricas binarias.

El extreme R×S interaction pasó de una mediana native de ~**+0.051 MCC** a aproximadamente cero en el four-corner common set. Sin embargo, debido a la severa reducción de soporte, este análisis se interpreta como un stress test inconcluso y no como evidencia de que la interacción native sea un artefacto.

La evidencia common-test más informativa para R×S proviene de los contrastes pairwise: el penalty S permanece grande en R90 y desaparece progresivamente hacia R30 sobre las mismas proteínas. Así, **la estructura de atenuación que define la interacción R×S sobrevive al control de test composition**, aunque la interacción extrema no puede ser estimada de forma fiable sobre una única población común a los cuatro corners.

### Take-home message

> **Homology-separated evaluation imposes a generalization penalty that largely persists when the evaluated proteins are held fixed. Redundancy reduction has a different character: much of its benchmark-wide effect arises because it restructures the available evidence and the biological population being evaluated, although task-specific same-protein effects remain, particularly for endolysin.**

La frase conceptual que guardaría para el paper es:

> **Redundancy control and homology separation affect generalization through different mechanisms: the former reshapes the evidence and evaluation population, whereas the latter primarily changes how far models must generalize within that population.**

Y para conectar directamente con RQ4:

> **The attenuation of homology-separation penalties under stringent redundancy control persists on shared test proteins, supporting a genuine dependence between the two evaluation axes beyond differences in test-set composition.**

## RQ17 — Paired finite-test bootstrap uncertainty

RQ17 evaluó la incertidumbre asociada al tamaño finito de los outer test sets mediante 1,000 bootstrap replicates estratificados, manteniendo congeladas todas las predicciones y utilizando paired resampling cuando los modelos fueron evaluados sobre las mismas entidades. Por tanto, los intervalos obtenidos cuantifican finite-test uncertainty condicional a los modelos ajustados y no training-seed variability.

### Representation effects are substantially more resolved than downstream-algorithm effects

La ventaja de ESM2-T30 sobre Roxy fue extremadamente robusta. El gain mediano fue aproximadamente **+0.214 MCC**, y ~**95% de los 140 contrastes** mostraron un 95% CI completamente por encima de cero. La ventaja fue directionally resolved en prácticamente todos los tasks, con cierta reducción de precisión únicamente en VAL.

El contraste T30−T6 fue menor y claramente task-dependent. Su gain mediano global fue aproximadamente **+0.069 MCC**, con ~59% de los contrastes completamente positivos. Sin embargo, global easy, global mixed y endolysin mostraron resolución positiva en **100%** de sus contrastes, mientras global hard presentó un gain mediano de solo ~**+0.006** y ~95% de los intervalos cruzaron cero. Este resultado refuerza la interpretación de representation saturation previamente observada para `global_hard`.

En contraste, las diferencias entre learners competitivos fueron mucho menos resueltas. `MLP−XGB` mostró un gain mediano de solo ~**+0.006 MCC**, con ~**89%** de los intervalos cruzando cero. `MLP−LR` mostró una señal mayor (~**+0.035 MCC**) pero task-dependent. Así, las diferencias entre representaciones fuertes están mucho mejor respaldadas frente a finite-test uncertainty que las diferencias entre clasificadores del upper tier.

### Finite-test uncertainty reproduces the attenuation of the S penalty under stringent R

Sobre las mismas proteínas compartidas, el contraste `S30−S90` presentó una pérdida mediana de aproximadamente **−0.039 MCC en R90**, **−0.020 en R70**, **−0.007 en R50** y aproximadamente **0 en R30**.

La proporción de contrastes cuyo 95% CI quedó completamente por debajo de cero disminuyó de aproximadamente:

$$
48\% \rightarrow 24\% \rightarrow 5\% \rightarrow 2\%
$$

para:

$$
R90 \rightarrow R70 \rightarrow R50 \rightarrow R30.
$$

Por tanto, la atenuación previamente observada no se limita a una reducción de los point estimates. A medida que redundancy control se hace más estricto, también desaparece progresivamente la resolución direccional del efecto de homology separation bajo finite-test resampling.

El common-test R effect mostró una historia diferente. Para `R30−R90` con T30, ~**89%** de los intervalos MCC cruzaron cero y solo ~6% fueron directionally negativos. Esto refuerza la conclusión de RQ16 de que la gran magnitud del native R effect refleja principalmente la reestructuración de evidencia y población inducida por redundancy control, mientras el same-entity component es mucho menor y frecuentemente no resoluble frente a finite-test uncertainty.

### PU uncertainty is strongly representation-dependent

BaggingPU mostró un recovery advantage robusto principalmente bajo one-hot. El gain mediano `BaggingPU−naive` fue ~**+0.093 U95**, con ~**59%** de los contrastes mostrando CIs completamente positivos. Para Roxy, T6 y T30, los gains medianos fueron cercanos a cero y la mayoría de los intervalos incluyeron cero.

En cambio, la ventaja de ESM2-T30 sobre ESM2-T6 fue extraordinariamente robusta bajo todos los paradigmas PU. Los gains medianos fueron ~**+0.072–0.076 U95**, y los **192 contrastes evaluados** mostraron intervalos completamente por encima de cero. Así, la representación T30 proporciona una mejora en known-positive recovery que es robusta frente a finite-test uncertainty y esencialmente independiente del paradigma PU utilizado.

Elkan–Noto mantuvo un comportamiento prácticamente indistinguible del naive para U95 a nivel agregado. Sus resultados continúan interpretándose bajo los problemas de clipping y ties documentados previamente y no proporcionan evidencia de una ventaja general de ranking.

### VAL uncertainty confirms an evidence-limited regime

VAL presentó intervalos mucho más amplios que los tasks con mayor soporte positivo. Con solo **14–26 positivos** en los outer-test corners analizados, los MCC bootstrap CI widths fueron aproximadamente **0.20–0.33**, mientras los intervalos de recall alcanzaron widths de **0.35–0.50**.

El número de positivos mostró una correlación negativa fuerte con la anchura del MCC CI:

$$
\rho \approx -0.71.
$$

Por tanto, la escasez de evidencia en VAL limita no solo el aprendizaje sino también la precisión con que la performance puede ser estimada.

### Take-home message

> **Finite-test uncertainty reinforces the main hierarchy of the benchmark: representation effects, particularly T30 over simpler feature spaces, are substantially more robust than differences among competitive downstream learners. The homology-separation penalty is most clearly resolved when redundancy remains high and loses both magnitude and directional resolution as redundancy control becomes stringent. Dedicated PU learning provides a representation-specific recovery advantage, whereas the T30 representation advantage remains robust across PU paradigms. Low-support tasks such as VAL remain fundamentally uncertainty-limited.**

## RQ18 — Computational burden versus predictive performance

RQ18 mostró que **predictive performance y computational burden no siguen la misma jerarquía**, y que una representación pretrained fuerte puede permitir modelos downstream extremadamente competitivos sin requerir necesariamente los learners más costosos.

El resultado práctico más claro apareció con **ESM2-T30 + KNN**, que alcanzó una MCC mediana de aproximadamente **0.914** con un downstream runtime mediano de solo **~1.8 s**. Esta combinación quedó dentro de **0.02 MCC del mejor modelo de cada `task × R × S` cell en ~79% de los casos** y perteneció al performance/runtime Pareto frontier en ~**87%** de las celdas. En contraste, `T30 + MLP` obtuvo una MCC mediana ligeramente superior (~**0.917**) y fue near-best en ~**90%** de las celdas, pero requirió ~**37 s** medianos, aproximadamente **21× más runtime** que KNN. `T30 + XGB` mostró una MCC mediana de ~**0.911**, con ~76% near-best, pero un runtime mediano de ~**60 s**, aproximadamente **33× el costo downstream de T30+KNN**.

Por tanto, el principal trade-off práctico no fue `T30 + LR`, como se había anticipado inicialmente, sino **T30 + KNN**. LR siguió siendo relativamente liviano frente a MLP/XGB, pero quedó within 0.02 MCC del cell-best en solo ~**26%** de las celdas. Su utilidad fue más clara frente a RF: dentro de T30, LR estuvo within 0.02 MCC de RF en ~**92%** de las comparaciones y fue simultáneamente near-and-faster en ~**80%**. En cambio, frente a MLP y XGB el pequeño costo de LR estuvo acompañado con mayor frecuencia por una pérdida predictiva apreciable.

### KNN provides a strong Pareto solution but shifts cost toward inference

El favorable downstream runtime de KNN debe interpretarse junto con su prediction-time burden. Para T30, el prediction time mediano fue aproximadamente:

```text
KNN   ~0.506 s
RF    ~0.102 s
LR    ~0.093 s
MLP   ~0.060 s
XGB   ~0.036 s
```

Así, KNN prácticamente elimina training complexity pero paga un costo de inferencia mayor. Para los tamaños evaluados aquí sigue siendo un trade-off excelente, pero no debe extrapolarse automáticamente a escenarios de deployment con millones de predicciones.

### One-hot + MLP is computationally dominated

El régimen más costoso del benchmark fue claramente **one-hot + MLP**. Su runtime mediano fue aproximadamente **938 s**, frente a ~37 s para T30+MLP, y algunas ejecuciones superaron los **3,500 s**. A pesar de ello, su MCC mediana fue solo ~**0.545**, muy inferior al ~0.917 obtenido con T30+MLP.

One-hot+MLP presentó:

```text
Pareto rate     = 0%
Near-best rate  = 0%
```

por lo que constituye un régimen claramente dominated: elevado costo computacional sin una compensación predictiva comparable.

### T30 improves prediction without requiring a proportional increase in downstream compute

El aumento de T6 a T30 produjo mejoras importantes de MCC, pero el downstream-runtime premium dependió fuertemente del learner:

```text
Learner    ΔMCC T30−T6    Runtime T30/T6
KNN          +0.073            1.19×
LR           +0.151            0.94×
MLP          +0.066            1.14×
RF           +0.098            1.40×
XGB          +0.066            1.95×
```

Por tanto, una representación pretrained más fuerte **no implica automáticamente un gran incremento del costo downstream**. Para KNN y MLP, T30 produjo ganancias claras con solo ~14–19% de runtime adicional, mientras LR incluso presentó un runtime mediano ligeramente inferior. XGB fue la principal excepción, con casi el doble de downstream runtime.

### The main T30 computational premium lies in representation generation

La generación de embeddings mostró una historia diferente. T6 y T30 fueron extraídos sobre exactamente **434,768 secuencias**, utilizando el mismo workflow, CUDA y precisión FP32:

```text
ESM2-T6
dimensions       320
runtime          ~5,147 s  (~1.43 h)
seconds / 10k    ~118

ESM2-T30
dimensions       640
runtime          ~25,460 s (~7.07 h)
seconds / 10k    ~586
```

Así:

$$
\frac{\mathrm{runtime}_{T30}}
     {\mathrm{runtime}_{T6}}
\approx
\boxed{4.95}
$$

El principal premium computacional de T30 ocurre entonces **durante la generación de la representación**, no durante el ajuste del classifier downstream. Este costo es además amortizable, porque los embeddings pueden reutilizarse entre tasks, condiciones R/S y learners.

### R changes computational scale; S does not

El runtime mediano disminuyó progresivamente al hacer R más estricto:

```text
R90   ~16.4 s
R70   ~13.6 s
R50   ~10.1 s
R30   ~ 8.4 s
```

En comparaciones pareadas:

$$
\mathrm{median}
\left(
\frac{runtime_{R30}}{runtime_{R90}}
\right)
\approx
0.42.
$$

Así, R30 utiliza aproximadamente **42% del runtime de R90**, consistente con la reducción y reestructuración del universo de entrenamiento.

En cambio:

$$
\mathrm{median}
\left(
\frac{runtime_{S30}}{runtime_{S90}}
\right)
\approx
0.995,
$$

por lo que homology separation prácticamente no cambia el workload computacional. Esto refuerza la distinción conceptual de RQ16: **R modifica tanto la estructura estadística como la escala computacional del problema, mientras S modifica principalmente la dificultad de generalización**.

El runtime también mostró fuerte dependencia del training-set size para KNN, RF y MLP, con median Spearman aproximados de **0.99, 0.94 y 0.85**, respectivamente. La asociación fue menor para XGB (~0.46) y LR (~0.39).

### Memory measurements should not support model-level claims

Aunque se registraron `process_peak_rss_kb_start/end`, el `peak_delta` mediano fue cero para prácticamente todas las familias y los mismos high-water RSS values se repitieron entre diferentes learners. Estos campos reflejan un **process-level high-water mark**, no una medición limpia de peak memory incremental por modelo.

Por tanto, RQ18 **no debe utilizarse para establecer una jerarquía de memory consumption entre modelos**.

### PU computational cost is not estimable

Los 768 experimentos PU conservaron `manifest.json` y `pu_diagnostics.tsv`, pero ninguno de estos artifacts contiene un wall-time utilizable y `runtime.json` está ausente en las 768 ejecuciones. En consecuencia:

> **PU compute overhead is not estimable from the retained benchmark artifacts.**

No puede cuantificarse de forma defendible el runtime multiplier de BaggingPU o Elkan–Noto frente a naive, ni el costo broad-versus-conservative. Esto corresponde a una limitación de instrumentation y no afecta los resultados predictivos PU obtenidos en RQ12–RQ17.

### Take-home message

> **Strong protein representations can substantially reduce the need for expensive downstream modeling. ESM2-T30 combined with KNN provides the strongest performance–runtime Pareto trade-off in the executed benchmark, approaching the predictive performance of MLP and XGB at a fraction of their downstream computational burden. The principal computational premium of scaling the pretrained model from T6 to T30 is paid once during representation generation rather than repeatedly during downstream model fitting.**

La frase conceptual que guardaría para el paper es:

> **Representation quality and downstream model complexity are not interchangeable computational burdens: stronger pretrained representations can shift predictive complexity upstream, enabling lightweight downstream models while concentrating most additional computational cost in an amortizable representation-generation step.**

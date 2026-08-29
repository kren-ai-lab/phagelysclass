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

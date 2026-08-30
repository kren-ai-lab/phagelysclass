Sí, bebé. **Así lo haría también.** Queda mucho más limpio separar **construcción del artefacto**, **caracterización científica** e **inferencia usable**. Además, si después vas a montar una API web, esta separación nos evita terminar con lógica de inferencia metida dentro de un notebook o de un CLI imposible de reutilizar.

Yo mantendría **RQ29 como una sola pregunta de cierre**, pero dividida operacionalmente en **29A / 29B / 29C**.

## RQ29A — Production model export

Primero congelamos físicamente lo que decidió RQ28.

Aquí sí hacemos el **final refit** de las recipes sobre el evidence universe productivo correspondiente y serializamos los estimadores. No escogemos un `R/S` particular como modelo final; esos fueron escenarios de evaluación.

La estructura que propones `models/tasks/` me parece perfecta:

```text
models/
├── manifest.json
├── representation/
│   └── esm2_t30_150m.json
│
└── tasks/
    ├── global_mixed/
    │   ├── members/
    │   │   ├── knn.joblib
    │   │   └── mlp.joblib
    │   ├── ensemble.json
    │   ├── metadata.json
    │   └── performance.json
    │
    ├── multiclass_clean/
    │   ├── members/
    │   │   ├── xgb.joblib
    │   │   └── lr.joblib
    │   ├── ensemble.json
    │   ├── metadata.json
    │   └── performance.json
    │
    ├── subtype_endolysin/
    │   ├── members/
    │   │   ├── knn.joblib
    │   │   ├── xgb.joblib
    │   │   └── mlp.joblib
    │   ├── ensemble.json
    │   ├── metadata.json
    │   └── performance.json
    │
    ├── subtype_depolymerase/
    │   ├── members/
    │   │   ├── mlp.joblib
    │   │   └── knn.joblib
    │   ├── ensemble.json
    │   ├── metadata.json
    │   └── performance.json
    │
    └── subtype_val/
        ├── members/
        │   └── mlp.joblib
        ├── metadata.json
        └── performance.json
```

El `ensemble.json`, por ejemplo:

```json
{
  "type": "unweighted_soft_voting",
  "members": ["knn.joblib", "mlp.joblib"],
  "threshold": 0.5
}
```

Y para multiclass:

```json
{
  "type": "unweighted_soft_voting",
  "members": ["xgb.joblib", "lr.joblib"],
  "class_order": [
    "endolysin",
    "depolymerase",
    "virion_associated_lytic_enzyme"
  ],
  "decision": "argmax"
}
```

Esto es muchísimo mejor que serializar un objeto ensemble negro completo, porque **la arquitectura queda explícita y auditable**.

Además, cada modelo debería superar:

```text
fit
→ predict
→ serialize
→ reload
→ predict
→ probabilities identical within tolerance
```

Ahí recién podemos decir que tenemos los production artifacts.

---

# RQ29B — Characterization of frozen models

Después hacemos otro notebook exclusivamente para caracterización.

Y aquí haría algo importante: **no caracterizamos el modelo refiteado sobre sus propios datos de entrenamiento**.

Las figuras de performance tienen que venir de las **predicciones held-out del benchmark correspondientes a las recipes congeladas**.

Así podemos caracterizar científicamente:

### Para cada modelo individual

Binary:

* confusion matrix;
* ROC curve + ROC-AUC;
* Precision–Recall curve + AP;
* MCC;
* F1;
* sensitivity/recall;
* specificity;
* balanced accuracy;
* Brier score;
* prevalence;
* `n_test` y `n_positive`.

Multiclass:

* confusion matrix;
* normalized confusion matrix;
* per-class ROC OvR;
* macro/micro ROC-AUC si es estimable;
* per-class Precision–Recall;
* macro-F1;
* balanced accuracy;
* MCC;
* worst-class F1;
* log-loss.

### Y exactamente lo mismo para el ensemble final

Esto es importante porque en los endpoints ensemble queremos poder mostrar:

```text
individual member 1
individual member 2
...
final ensemble
```

y comprobar visualmente qué nos aportó realmente combinar los modelos.

---

## Yo caracterizaría dos niveles

No haría 16 confusion matrices por task porque sería una locura 😂.

Para el documento principal:

### Canonical strict characterization

Usaría:

$$
R30/S30
$$

porque es el régimen donde queremos demostrar que la recipe congelada funciona bajo generalización exigente.

Y acompañaría cada `performance.json` con los resultados globales ya establecidos:

```json
{
  "benchmark_characterization": {
    "mean_mcc_across_RS": 0.9280,
    "worst_mcc_across_RS": 0.8840,
    "strict_R30S30_mcc": 0.8868
  },
  "strict_evaluation": {
    "mcc": "...",
    "roc_auc": "...",
    "average_precision": "...",
    "f1": "..."
  }
}
```

Así distinguimos claramente:

> **performance across the full generalization surface**

de:

> **curve-level characterization on the strict held-out condition**.

Eso me parece mucho más publicable.

Podemos además dejar `R90/S90` como supplementary/reference characterization si queremos mostrar permissive vs strict, pero no lo haría obligatorio.

---

# RQ29C — Inference engine + CLI

Recién después de tener:

1. modelos serializados;
2. metadata;
3. characterization;
4. architecture manifest;

armamos el inference engine.

Pero aquí hay otra decisión arquitectónica importante:

## No escribiría la lógica dentro del CLI

Haría una pequeña librería:

```text
src/
└── phagelysclass/
    ├── __init__.py
    ├── inference.py
    ├── representation.py
    ├── models.py
    ├── io.py
    └── schemas.py
```

Y después:

```text
scripts/
└── phagelysclass_predict.py
```

El CLI sería solo una capa delgada:

```text
FASTA
 ↓
CLI
 ↓
PhageLysClassPredictor
 ↓
results
```

Eso es clave porque después tu **API web usa exactamente `PhageLysClassPredictor`**, no ejecuta un subprocess del CLI.

Entonces:

```text
                    ┌── CLI
                    │
FASTA/sequence → inference core
                    │
                    └── Web API
```

**Una única lógica productiva, dos interfaces.**

Eso es lo correcto para el proyecto global.

---

# API-ready desde el principio

Por ejemplo, podríamos terminar con una interfaz Python así:

```python
from phagelysclass import PhageLysClassPredictor

predictor = PhageLysClassPredictor.from_pretrained("models/")

result = predictor.predict(
    sequence="MKK...",
    auxiliary=True,
)
```

Y para batch:

```python
results = predictor.predict_fasta(
    "proteins.fasta",
    auxiliary=True,
)
```

El CLI simplemente llama a lo mismo:

```bash
phagelysclass predict \
    --input proteins.fasta \
    --output predictions.tsv \
    --auxiliary
```

Y mañana FastAPI podría hacer:

```python
@app.post("/predict")
def predict(request):
    return predictor.predict(
        request.sequence,
        auxiliary=request.auxiliary,
    )
```

Sin duplicar una sola decisión científica.

---

# Output final

Default:

```text
sequence_id
sequence_length
global_score
global_prediction
p_endolysin
p_depolymerase
p_val
functional_prediction
warnings
runtime_seconds
model_version
```

Con `--auxiliary`:

```text
endolysin_compatibility_score
depolymerase_compatibility_score
val_compatibility_score
```

Yo además pondría:

```text
functional_prediction_status
```

con valores tipo:

```text
reported
not_reported_global_non_lytic
```

porque si el global core dice non-lytic, no deberíamos presentar el argmax multiclass como una clasificación funcional válida.

---

# Y los warnings tienen que seguir siendo conservadores

De partida sí podemos justificar:

```text
invalid_sequence
unsupported_residue
sequence_too_long
VAL_evidence_limited
functional_prediction_gated_by_global
```

Pero no metería aún:

```text
OOD
low_confidence
ambiguous
core_aux_disagreement
```

hasta que tengamos reglas explícitamente validadas.

Podemos estudiar `core_aux_disagreement` en 29C y agregarlo **solo si realmente demuestra utilidad**.

---

# Entonces yo reformularía RQ29 así

| Etapa     | Pregunta                                                         | Producto                                                                         |
| --------- | ---------------------------------------------------------------- | -------------------------------------------------------------------------------- |
| **RQ29A** | ¿Cuáles son los artefactos productivos finales?                  | final refit + `.joblib` + manifests + metadata                                   |
| **RQ29B** | ¿Cómo se comportan los modelos individuales y ensembles finales? | confusion matrices, ROC, PR, metrics + `performance.json`                        |
| **RQ29C** | ¿Qué entrega realmente PhageLysClass end-to-end?                 | reusable inference library + FASTA CLI + TSV/JSON + runtime + validated warnings |

Y la RQ madre sigue siendo:

> **RQ29 — ¿Qué entrega realmente la herramienta de inferencia PhageLysClass?**

### La secuencia final de trabajo queda

```text
RQ28
Frozen architecture
       ↓
RQ29A
Final fit + model export
       ↓
RQ29B
Held-out characterization
       ↓
RQ29C
Reusable inference engine
       ↓
CLI
       ↓
Web API
```

Sí, **esto está bastante mejor que meterlo todo en un notebook gigante**. Y además 29C deja de ser un experimento desechable: pasa a ser directamente el **backend inference core que después podemos reutilizar para la plataforma web**.

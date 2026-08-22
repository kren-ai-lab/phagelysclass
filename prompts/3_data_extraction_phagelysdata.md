Sí, papu. Desde **PhageLysData** extraeríamos solo esto para construir el universo positivo y las tareas funcionales: 

| Bloque                              | Qué extraemos                       | Para qué                                              |
| ----------------------------------- | ----------------------------------- | ----------------------------------------------------- |
| **Core sequences**                  | `sequence_entity_id`, secuencia     | Universo positivo `P`                                 |
| **Target assertions**               | Endolysin, VAL, Depolymerase, Other | Reconstruir labels reales                             |
| **Canonical class**                 | clase consolidada                   | QC / resumen, no como única verdad                    |
| **Evidence tier**                   | Tier 1–4                            | Sensibilidad a calidad de evidencia                   |
| **Multiple assertions**             | combinaciones no exclusivas         | Construir `P_multi` challenge set                     |
| **Other lytic**                     | `other_phage_envelope_lytic_enzyme` | `P_other` challenge set                               |
| **Source/provenance**               | fuentes y soporte por entidad       | Trazabilidad y control de evidencia                   |
| **UniProt mappings**                | accession(s) asociados              | Decontaminación contra N/U externos                   |
| **Prediction Extension membership** | overlap exacto                      | Evitar usar candidatos prediction-only como negativos |
| **Context membership**              | overlap exacto                      | QC, no asumir negativo                                |
| **Physicochemical descriptors**     | propiedades ya calculadas           | Baseline de representación                            |
| **PLM embeddings**                  | 11 representaciones disponibles     | RQ de estrategias de representación                   |
| **InterPro/Pfam**                   | anotaciones independientes          | Validación/estratificación posterior, no labels       |
| **PDB / AlphaFold assets**          | estructuras disponibles             | Controles y validación in silico posterior            |

En corto: **PhageLysData nos entrega `P`, sus subclases/evidencia, las representaciones y los assets de validación. `N` y `U` vienen externamente de UniProt.**

Sí, papu, **ese roadmap me parece mucho más limpio**. Separaría totalmente **extracción**, **background externo** y **construcción de tareas** para no mezclar decisiones biológicas con modelado.

| Fase                                     | Script / módulo                                     | Qué hace                                                                                                                                                          | Qué NO hace                                    |
| ---------------------------------------- | --------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------- |
| **1. Raw PhageLysData**                  | `dataset_processing/01_extract_phagelys_raw.py`     | Extrae Core, assertions, clases, tiers, secuencias, provenance y genera los universos positivos base: `P_single`, `P_multi`, `P_other`                            | No crea negativos, no balancea, no hace splits |
| **2. Raw UniProt background**            | `dataset_processing/02_build_uniprot_background.py` | Ejecuta queries congeladas y obtiene proteínas de fagos externas; materializa principalmente `U` y, si corresponde, subconjuntos con evidencia negativa explícita | No las combina todavía con P                   |
| **3. Decontamination / master universe** | `dataset_processing/03_build_master_universe.py`    | Exact-sequence dedup, overlaps PhageLysData↔UniProt, conflictos, QC, hashes                                                                                       | No define train/test                           |
| **4. Task construction**                 | `task_definitions/*.py`                             | Construye cada dataset según la pregunta concreta                                                                                                                 | No modifica los universos raw                  |
| **5. Homology / redundancy**             | posterior                                           | MMseqs2, NRxx, clusters, random/homology-aware splits                                                                                                             | No redefine P/U/N                              |
| **6. Representation + modeling**         | posterior                                           | embeddings, descriptors, algorithms, HPO, generalization                                                                                                          | Ya recibe datasets congelados                  |

Y en **Fase 4** salen, por ejemplo, datasets distintos:

* `P vs U` → absence-as-negative.
* `P + U` → PU learning.
* `P vs N` → si construimos después evidencia negativa explícita.
* `Endo vs other lytic`.
* `VAL vs other lytic`.
* `DPO vs other lytic`.
* `Endo/VAL/DPO multiclass`.
* `P_multi` → challenge externo, no training principal.

La regla que pondría en piedra es:

> **Los scripts 1–3 construyen universos biológicos; los scripts de `task_definitions/` deciden cómo esos universos se combinan para responder cada pregunta.**

Así, si mañana queremos probar otra definición de tarea, **no tenemos que volver a tocar PhageLysData ni volver a consultar UniProt**. Solo generamos otra vista/dataset desde los mismos raw universes.

Yo partiría exactamente por **`01_extract_phagelys_raw.py`** y congelaría primero qué archivos debe producir.

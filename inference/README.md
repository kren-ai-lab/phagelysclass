# PhageLysClass inference

Minimal single-script inference pipeline for the frozen PhageLysClass production bundle.

```text
phagelysclass_predict.py
README.md
known_positive_depolymerase_demo.fasta
known_positive_depolymerase_demo_metadata.tsv
```

## Representation contract

Inference now generates embeddings through **Sylphy**, matching the backend used by the benchmark:

```text
facebook/esm2_t30_150M_UR50D
final representation
mean sequence pooling
FP32
640 dimensions
```

Install the embedding dependencies with:

```bash
pip install 'sylphy[embeddings]'
```

The remaining fitted-model dependencies must also remain compatible with the versions recorded in:

```text
models/phagelysclass_v1/manifest.json
```

## Run

```bash
python phagelysclass_inference/phagelysclass_predict.py \
    --input known_positive_depolymerase_demo.fasta \
    --models models/phagelysclass_v1 \
    --output known_positive_predictions.tsv \
    --auxiliary
```

The main arguments are:

```text
--input       protein FASTA
--models      exported RQ29A model bundle
--output      .tsv or .csv prediction table
--auxiliary   optional binary Endolysin/DPO/VAL compatibility scores
--device      auto | cpu | cuda
--batch-size  Sylphy embedding batch size
```

## Output semantics

The global core is evaluated first:

```text
global_score
global_prediction = lytic | non_lytic_background
```

The multiclass core answers a different question: **if the sequence is treated as lytic, which of Endolysin / Depolymerase / VAL is most compatible?**

Therefore its columns are explicitly conditional:

```text
p_endolysin_given_lytic
p_depolymerase_given_lytic
p_val_given_lytic
```

These three probabilities are reported only when the global core predicts `lytic`. When the global core predicts `non_lytic_background`, they are written as `NA` and no primary functional label is emitted.

The output also contains:

```text
functional_prediction
functional_prediction_status
warnings
model_version
representation
```

With `--auxiliary`, the three independent binary compatibility scores are also exported:

```text
endolysin_compatibility_score
depolymerase_compatibility_score
val_compatibility_score
```

These are diagnostic outputs and do not override the global or multiclass core.

## Positive-control FASTA

`known_positive_depolymerase_demo.fasta` contains five real sequences selected from the uploaded PhageLysData candidate assets. All five are annotated in the supplied selection tables as:

```text
canonical_target_class = depolymerase
```

and are mapped to UniProt accessions. They are intended as a **positive-control smoke test** for pipeline execution and qualitative sanity checking, not as an independent performance benchmark.

Expected qualitative behavior:

```text
global_prediction: preferably lytic
functional_prediction: preferably depolymerase
```

Because these examples come from the project evidence universe, this test should not be reported as external validation.

## Run metadata

A sidecar JSON is written automatically at:

```text
<output>.metadata.json
```

It records the model bundle, Sylphy backend, ESM2 checkpoint, pooling, precision, device, batch size, valid/invalid sequence counts, and runtimes.

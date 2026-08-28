#!/bin/bash
set -e

TASKS=(
    global_easy
    global_hard
    global_mixed
    subtype_endolysin
    subtype_depolymerase
    subtype_val
    multiclass_clean
)

for TASK in "${TASKS[@]}"; do
    echo "============================================================"
    echo "Running: ${TASK} × esm2_t30_150m"
    echo "============================================================"

    python src/training/11_train_supervised_factorial.py \
        --task "${TASK}" \
        --representation esm2_t30_150m
done

echo "============================================================"
echo "ALL ESM2-t30-150M SUPERVISED TASKS COMPLETED"
echo "============================================================"
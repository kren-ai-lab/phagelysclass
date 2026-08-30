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
    echo "Running: ${TASK} × onehot_1024"
    echo "============================================================"

    python src/training/11_train_supervised_factorial.py \
        --task "${TASK}" \
        --representation onehot_1024
done

echo "============================================================"
echo "ALL ONE HOT SUPERVISED TASKS COMPLETED"
echo "============================================================"
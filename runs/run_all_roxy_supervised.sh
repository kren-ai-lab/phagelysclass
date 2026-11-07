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
    echo "Running: ${TASK} × roxy_48"
    echo "============================================================"

    python src/training/11_train_supervised_factorial.py \
        --task "${TASK}" \
        --representation roxy_48
done

echo "============================================================"
echo "ALL ROXY SUPERVISED TASKS COMPLETED"
echo "============================================================"
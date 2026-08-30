#!/bin/bash
set -e
for TASK in pu_conservative pu_broad; do
    echo "============================================================"
    echo "Running PU: ${TASK} × esm2_t6_8m"
    echo "============================================================"
    python src/training/11_train_pu_factorial.py --task "${TASK}" --representation esm2_t6_8m
done
echo "ALL ESM2-T6 PU TASKS COMPLETED"

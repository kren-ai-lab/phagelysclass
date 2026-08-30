#!/bin/bash
set -e
for TASK in pu_conservative pu_broad; do
    echo "============================================================"
    echo "Running PU: ${TASK} × esm2_t30_150m"
    echo "============================================================"
    python src/training/11_train_pu_factorial.py --task "${TASK}" --representation esm2_t30_150m
done
echo "ALL ESM2-T30 PU TASKS COMPLETED"

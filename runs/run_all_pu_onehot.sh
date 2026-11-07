#!/bin/bash
set -e
for TASK in pu_conservative pu_broad; do
    echo "============================================================"
    echo "Running PU: ${TASK} × onehot_1024"
    echo "============================================================"
    python src/training/11_train_pu_factorial.py --task "${TASK}" --representation onehot_1024
done
echo "ALL ONE-HOT PU TASKS COMPLETED"

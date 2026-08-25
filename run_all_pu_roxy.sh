#!/bin/bash
set -e
for TASK in pu_conservative pu_broad; do
    echo "============================================================"
    echo "Running PU: ${TASK} × roxy_48"
    echo "============================================================"
    python src/training/11_train_pu_factorial.py --task "${TASK}" --representation roxy_48
done
echo "ALL ROXY PU TASKS COMPLETED"

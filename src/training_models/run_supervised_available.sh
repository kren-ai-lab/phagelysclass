#!/bin/bash
# Usage:
#   bash run_supervised_available.sh TASK TASK_TYPE R S SPLIT_FILE [POSITIVE_LABEL]
#
# Example:
#   bash run_supervised_available.sh global_easy binary R90 S90 path/to/split.tsv.gz 1

set -e

TASK="$1"
TASK_TYPE="$2"
R="$3"
S="$4"
SPLIT_FILE="$5"
POSITIVE_LABEL="${6:-}"

REPRESENTATIONS=(
    "roxy_48"
    "onehot_1024"
    "esm2_t6_8m"
)

for REP in "${REPRESENTATIONS[@]}"; do
    echo "============================================================"
    echo "${TASK} | ${R}_${S} | ${REP}"
    echo "============================================================"

    CMD=(
        python src/training_models/11a_train_supervised.py
        --task "$TASK"
        --task-type "$TASK_TYPE"
        --redundancy "$R"
        --split-threshold "$S"
        --split-file "$SPLIT_FILE"
        --representation "$REP"
        --algorithm all
    )

    if [[ -n "$POSITIVE_LABEL" ]]; then
        CMD+=(--positive-label "$POSITIVE_LABEL")
    fi

    "${CMD[@]}"
done

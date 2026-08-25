#!/usr/bin/env python3
"""
10f_freeze_inner_group_splits.py

Freeze one deterministic homology-group-aware inner train/validation split
for every Stage-10E task × R × S cell.

The script is intentionally parameter-free for production execution:
all 9 tasks and all 16 R/S cells are mandatory.

Algorithm
---------
- Input: eligible outer-train entities only.
- Group: Stage-08 S-level cluster_id.
- Candidate partitions: StratifiedGroupKFold(n_splits=5, shuffle=True, seed=42).
- Select the fold whose validation fraction and class/role composition are
  closest to the target 20%, while preserving complete clusters.

Output
------
data/training_splits/<task>/<R>_<S>/
    inner_seed42.tsv.gz
    inner_split_summary.tsv
    inner_split_manifest.json

Global:
data/training_splits/stage10f_summary.tsv
data/training_splits/stage10f_freeze_manifest.json
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedGroupKFold

from _training_common import (
    PU_TASKS,
    R_VALUES,
    S_VALUES,
    SEED,
    TASKS,
    project_root_from_script,
    pu_target,
    read_tsv,
    sha256_file,
    supervised_target,
    write_json,
)


TARGET_VAL = 0.20


def split_one(df: pd.DataFrame, task_id: str, task_type: str):
    outer_train = df.loc[df["split"].eq("train")].copy()
    if len(outer_train) == 0:
        raise RuntimeError("No eligible outer-train entities")

    if task_type == "pu":
        y, target_col = pu_target(outer_train)
        class_names = {0: "U", 1: "P"}
    else:
        y, class_names = supervised_target(outer_train, task_type)
        target_col = "label" if task_type == "binary" else "multiclass_target"

    groups = outer_train["cluster_id"].astype(str).to_numpy()

    if len(np.unique(y)) < 2:
        raise RuntimeError("Inner split requires at least two target states")

    splitter = StratifiedGroupKFold(
        n_splits=5,
        shuffle=True,
        random_state=SEED,
    )

    overall_prop = {
        int(cls): float(np.mean(y == cls)) for cls in np.unique(y)
    }

    candidates = []
    X_dummy = np.zeros((len(outer_train), 1), dtype=np.uint8)

    for fold_id, (train_idx, val_idx) in enumerate(
        splitter.split(X_dummy, y, groups)
    ):
        y_train = y[train_idx]
        y_val = y[val_idx]

        if set(np.unique(y_train)) != set(np.unique(y)):
            continue
        if set(np.unique(y_val)) != set(np.unique(y)):
            continue

        val_frac = len(val_idx) / len(outer_train)
        prevalence_error = 0.0
        for cls, overall in overall_prop.items():
            prevalence_error += abs(float(np.mean(y_val == cls)) - overall)

        score = abs(val_frac - TARGET_VAL) + prevalence_error
        candidates.append(
            (score, fold_id, train_idx, val_idx, val_frac, prevalence_error)
        )

    if not candidates:
        raise RuntimeError(
            f"{task_id}: no 5-fold group-aware candidate retained all classes "
            "in both inner partitions."
        )

    candidates.sort(key=lambda x: (x[0], x[1]))
    _, fold_id, train_idx, val_idx, val_frac, prevalence_error = candidates[0]

    assignment = outer_train[
        ["entity_id", "sequence_uid", "cluster_id"]
    ].copy()
    assignment["inner_split"] = ""
    assignment.iloc[train_idx, assignment.columns.get_loc("inner_split")] = "train"
    assignment.iloc[val_idx, assignment.columns.get_loc("inner_split")] = "validation"

    shared = (
        set(assignment.loc[assignment.inner_split.eq("train"), "cluster_id"])
        & set(
            assignment.loc[
                assignment.inner_split.eq("validation"), "cluster_id"
            ]
        )
    )
    if shared:
        raise RuntimeError(
            f"Inner cluster leakage: {len(shared)} shared clusters"
        )

    return assignment, y, class_names, fold_id, val_frac, prevalence_error


def main():
    root = project_root_from_script(__file__)
    predictive_root = root / "data/predictive_universe"
    output_root = root / "data/training_splits"
    output_root.mkdir(parents=True, exist_ok=True)

    rows = []

    print("=" * 100)
    print("PhageLysClass — Stage 10F group-aware inner split freeze")
    print("=" * 100)

    for task_id, info in TASKS.items():
        task_path = info["path"]
        task_type = info["type"]

        for r in R_VALUES:
            for s in S_VALUES:
                rs = f"R{r}_S{s}"
                input_path = (
                    predictive_root
                    / task_path
                    / rs
                    / "eligible_entities.tsv.gz"
                )
                if not input_path.is_file():
                    raise FileNotFoundError(input_path)

                out_dir = output_root / task_path / rs
                out_dir.mkdir(parents=True, exist_ok=True)

                df = read_tsv(input_path)
                required = {
                    "entity_id",
                    "sequence_uid",
                    "cluster_id",
                    "split",
                }
                missing = required - set(df.columns)
                if missing:
                    raise KeyError(
                        f"{input_path}: missing {sorted(missing)}"
                    )

                assignment, y_outer_train, class_names, fold_id, val_frac, prev_err = split_one(
                    df, task_id, task_type
                )

                train_ids = set(
                    assignment.loc[
                        assignment.inner_split.eq("train"), "entity_id"
                    ]
                )
                val_ids = set(
                    assignment.loc[
                        assignment.inner_split.eq("validation"), "entity_id"
                    ]
                )

                outer_train = df.loc[df["split"].eq("train")].copy()
                target_y = (
                    pu_target(outer_train)[0]
                    if task_type == "pu"
                    else supervised_target(outer_train, task_type)[0]
                )
                target_map = dict(
                    zip(outer_train["entity_id"].astype(str), target_y)
                )

                y_train = np.array([target_map[x] for x in train_ids], dtype=int)
                y_val = np.array([target_map[x] for x in val_ids], dtype=int)

                summary = {
                    "task_id": task_id,
                    "task_type": task_type,
                    "R": r,
                    "S": s,
                    "rs_id": rs,
                    "n_outer_train": len(outer_train),
                    "n_inner_train": len(train_ids),
                    "n_validation": len(val_ids),
                    "validation_fraction": len(val_ids) / len(outer_train),
                    "selected_candidate_fold": fold_id,
                    "composition_error": prev_err,
                    "inner_train_composition": "; ".join(
                        f"{class_names[int(k)]}={int(v)}"
                        for k, v in zip(*np.unique(y_train, return_counts=True))
                    ),
                    "validation_composition": "; ".join(
                        f"{class_names[int(k)]}={int(v)}"
                        for k, v in zip(*np.unique(y_val, return_counts=True))
                    ),
                    "n_inner_train_clusters": int(
                        assignment.loc[
                            assignment.inner_split.eq("train"), "cluster_id"
                        ].nunique()
                    ),
                    "n_validation_clusters": int(
                        assignment.loc[
                            assignment.inner_split.eq("validation"), "cluster_id"
                        ].nunique()
                    ),
                    "shared_inner_train_validation_clusters": 0,
                    "status": "PASS",
                }

                assignment.to_csv(
                    out_dir / "inner_seed42.tsv.gz",
                    sep="\t",
                    index=False,
                    compression="gzip",
                )
                pd.DataFrame([summary]).to_csv(
                    out_dir / "inner_split_summary.tsv",
                    sep="\t",
                    index=False,
                )
                write_json(
                    out_dir / "inner_split_manifest.json",
                    {
                        "stage": "10F",
                        "status": "PASS",
                        "task_id": task_id,
                        "task_type": task_type,
                        "R": r,
                        "S": s,
                        "seed": SEED,
                        "target_validation_fraction": TARGET_VAL,
                        "algorithm": (
                            "StratifiedGroupKFold(n_splits=5, shuffle=True); "
                            "choose candidate minimizing validation-size plus "
                            "class/role-composition deviation."
                        ),
                        "input": str(input_path),
                        "input_sha256": sha256_file(input_path),
                        "summary": summary,
                    },
                )

                rows.append(summary)
                print(
                    f"PASS {task_id} {rs} | "
                    f"outer_train={len(outer_train):,} "
                    f"inner={len(train_ids):,} val={len(val_ids):,} "
                    f"val_frac={len(val_ids)/len(outer_train):.4f}"
                )

    summary_df = pd.DataFrame(rows)
    summary_df.to_csv(
        output_root / "stage10f_summary.tsv",
        sep="\t",
        index=False,
    )

    expected = len(TASKS) * len(R_VALUES) * len(S_VALUES)
    passed = int(summary_df["status"].eq("PASS").sum())
    write_json(
        output_root / "stage10f_freeze_manifest.json",
        {
            "stage": "10F",
            "expected_cells": expected,
            "observed_cells": len(summary_df),
            "passed_cells": passed,
            "seed": SEED,
            "target_validation_fraction": TARGET_VAL,
            "status": "PASS" if passed == expected else "INCOMPLETE",
        },
    )

    print()
    print("=" * 100)
    print(f"Stage 10F cells PASS: {passed}/{expected}")
    print("=" * 100)

    if passed != expected:
        raise SystemExit(2)


if __name__ == "__main__":
    main()

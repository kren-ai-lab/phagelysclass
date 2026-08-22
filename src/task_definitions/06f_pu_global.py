#!/usr/bin/env python3
"""
Global Positive-Unlabeled task.

Creates two PU regimes:
1. conservative: P_core + U_explicit
2. broad:        P_core + U_explicit + U_other_unresolved

U is NOT interpreted as a negative class. `pu_observed_label=0` means
"unlabeled", not "true negative".
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from _task_utils import (
    add_common_args, add_group_counts, clean_output_dir, load_master,
    write_dataset, write_manifest, write_summary,
)


TASK_ID = "pu_global"
REGIMES = {
    "conservative": ["U_explicit"],
    "broad": ["U_explicit", "U_other_unresolved"],
}


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    return add_common_args(
        parser,
        __file__,
        "data/tasks/pu_global",
    ).parse_args()


def main():
    args = parse_args()
    master_path = Path(args.master).expanduser().resolve()
    output_dir = Path(args.output_dir).expanduser().resolve()
    clean_output_dir(output_dir)

    master = load_master(master_path)
    global_rows = []

    for regime, u_states in REGIMES.items():
        subdir = output_dir / regime
        clean_output_dir(subdir)

        mask = (
            master["sequence_state"].eq("P_core")
            | master["sequence_state"].isin(u_states)
        )
        task = master.loc[mask].copy()

        task["pu_role"] = task["sequence_state"].map(
            lambda x: "P" if x == "P_core" else "U"
        )
        task["pu_observed_label"] = task["pu_role"].map({"P": 1, "U": 0})
        task["label_semantics"] = task["pu_role"].map(
            {"P": "known_positive", "U": "unlabeled_not_negative"}
        )

        dataset_path = subdir / "dataset.tsv.gz"
        write_dataset(task, dataset_path)

        rows = []
        add_group_counts(rows, task, "pu_role", "pu_role")
        add_group_counts(rows, task, "sequence_state", "source_state")
        write_summary(rows, subdir / "summary.tsv")

        write_manifest(
            subdir / "manifest.json",
            {
                "task_id": f"{TASK_ID}:{regime}",
                "positive_definition": "All PhageLysData Core exact-sequence entities.",
                "unlabeled_definition": u_states,
                "important": (
                    "pu_observed_label=0 denotes unlabeled observations and must "
                    "not be interpreted as biological negatives."
                ),
                "N_easy_used_for_training": False,
                "N_hard_used_for_training": False,
                "no_balancing": True,
                "no_homology_filtering": True,
                "master": str(master_path),
                "dataset": str(dataset_path),
                "n_examples": len(task),
            },
        )

        counts = task["pu_role"].value_counts()
        global_rows.extend(
            [
                {"regime": regime, "group": "P", "count": int(counts.get("P", 0))},
                {"regime": regime, "group": "U", "count": int(counts.get("U", 0))},
                {"regime": regime, "group": "total", "count": int(len(task))},
            ]
        )

        print(
            f"{regime:16s} "
            f"P={int(counts.get('P', 0)):,} "
            f"U={int(counts.get('U', 0)):,} "
            f"total={len(task):,}"
        )

    pd.DataFrame(global_rows).to_csv(
        output_dir / "summary.tsv",
        sep="\t",
        index=False,
    )


if __name__ == "__main__":
    main()

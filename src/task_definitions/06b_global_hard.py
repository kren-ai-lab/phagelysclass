#!/usr/bin/env python3
"""Global lytic recognition: P_core vs N_hard."""

from __future__ import annotations

import argparse
from pathlib import Path

from _task_utils import (
    add_common_args, add_group_counts, clean_output_dir, load_master,
    print_summary, write_dataset, write_manifest, write_summary,
)


TASK_ID = "global_hard"


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    return add_common_args(
        parser,
        __file__,
        "data/tasks/global_hard",
    ).parse_args()


def main():
    args = parse_args()
    master_path = Path(args.master).expanduser().resolve()
    output_dir = Path(args.output_dir).expanduser().resolve()
    clean_output_dir(output_dir)

    master = load_master(master_path)
    task = master.loc[
        master["sequence_state"].isin(["P_core", "N_hard"])
    ].copy()

    task["label"] = task["sequence_state"].map({"P_core": 1, "N_hard": 0})
    task["class_name"] = task["label"].map({1: "lytic", 0: "non_lytic_hard"})
    task["training_role"] = task["label"].map({1: "P", 0: "N"})

    dataset_path = output_dir / "dataset.tsv.gz"
    write_dataset(task, dataset_path)

    rows = []
    add_group_counts(rows, task, "class_name", "class")
    add_group_counts(rows, task, "sequence_state", "source_state")
    write_summary(rows, output_dir / "summary.tsv")

    write_manifest(
        output_dir / "manifest.json",
        {
            "task_id": TASK_ID,
            "question": "Can phage lytic enzymes be distinguished from lysis-adjacent non-enzymatic phage proteins?",
            "positive_definition": "All PhageLysData Core exact-sequence entities.",
            "negative_definition": "N_hard exact-sequence entities only.",
            "no_balancing": True,
            "no_homology_filtering": True,
            "master": str(master_path),
            "dataset": str(dataset_path),
            "n_examples": len(task),
        },
    )

    print_summary("Task 06b — global hard", task, "class_name")


if __name__ == "__main__":
    main()

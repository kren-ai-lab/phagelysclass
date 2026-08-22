#!/usr/bin/env python3
"""Global lytic recognition: P_core vs N_easy + N_hard."""

from __future__ import annotations

import argparse
from pathlib import Path

from _task_utils import (
    add_common_args, add_group_counts, clean_output_dir, load_master,
    print_summary, write_dataset, write_manifest, write_summary,
)


TASK_ID = "global_mixed"


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    return add_common_args(
        parser,
        __file__,
        "data/tasks/global_mixed",
    ).parse_args()


def main():
    args = parse_args()
    master_path = Path(args.master).expanduser().resolve()
    output_dir = Path(args.output_dir).expanduser().resolve()
    clean_output_dir(output_dir)

    master = load_master(master_path)
    task = master.loc[
        master["sequence_state"].isin(["P_core", "N_easy", "N_hard"])
    ].copy()

    task["label"] = task["sequence_state"].map(
        {"P_core": 1, "N_easy": 0, "N_hard": 0}
    )
    task["class_name"] = task["label"].map({1: "lytic", 0: "non_lytic"})
    task["training_role"] = task["label"].map({1: "P", 0: "N"})
    task["negative_stratum"] = task["sequence_state"].map(
        {"P_core": "", "N_easy": "easy", "N_hard": "hard"}
    )

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
            "question": "Can phage lytic enzymes be recognized against a mixed explicit non-lytic background?",
            "positive_definition": "All PhageLysData Core exact-sequence entities.",
            "negative_definition": "Union of N_easy and N_hard exact-sequence entities.",
            "no_balancing": True,
            "no_homology_filtering": True,
            "master": str(master_path),
            "dataset": str(dataset_path),
            "n_examples": len(task),
        },
    )

    print_summary("Task 06c — global mixed", task, "class_name")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Clean three-class lytic subtype task.

Includes only single-class P_core:
- endolysin
- depolymerase
- virion_associated_lytic_enzyme

Excludes multiple_nonexclusive and other_phage_envelope_lytic_enzyme.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from _task_utils import (
    add_common_args, add_group_counts, clean_output_dir, load_master,
    print_summary, write_dataset, write_manifest, write_summary,
)


TASK_ID = "multiclass_clean"
CLASSES = [
    "endolysin",
    "depolymerase",
    "virion_associated_lytic_enzyme",
]


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    return add_common_args(
        parser,
        __file__,
        "data/tasks/multiclass_clean",
    ).parse_args()


def main():
    args = parse_args()
    master_path = Path(args.master).expanduser().resolve()
    output_dir = Path(args.output_dir).expanduser().resolve()
    clean_output_dir(output_dir)

    master = load_master(master_path)
    task = master.loc[
        master["sequence_state"].eq("P_core")
        & master["canonical_target_class"].isin(CLASSES)
    ].copy()

    label_map = {
        "endolysin": 0,
        "depolymerase": 1,
        "virion_associated_lytic_enzyme": 2,
    }
    task["label"] = task["canonical_target_class"].map(label_map)
    task["class_name"] = task["canonical_target_class"]
    task["training_role"] = "P"

    dataset_path = output_dir / "dataset.tsv.gz"
    write_dataset(task, dataset_path)

    rows = []
    add_group_counts(rows, task, "class_name", "class")
    write_summary(rows, output_dir / "summary.tsv")

    write_manifest(
        output_dir / "manifest.json",
        {
            "task_id": TASK_ID,
            "classes": label_map,
            "multiple_nonexclusive_in_training": False,
            "other_phage_envelope_lytic_enzyme_in_training": False,
            "no_balancing": True,
            "no_homology_filtering": True,
            "master": str(master_path),
            "dataset": str(dataset_path),
            "n_examples": len(task),
        },
    )

    print_summary("Task 06e — multiclass clean", task, "class_name")


if __name__ == "__main__":
    main()

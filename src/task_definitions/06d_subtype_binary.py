#!/usr/bin/env python3
"""
Independent subtype recognition tasks.

Creates one folder per subtype:
- endolysin
- depolymerase
- virion_associated_lytic_enzyme

For each task, positives are P_core sequences of the target subtype and
negatives are N_easy + N_hard. Other lytic Core classes are excluded, not
treated as negatives.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from _task_utils import (
    add_common_args, add_group_counts, clean_output_dir, load_master,
    write_dataset, write_manifest, write_summary,
)


TASK_ID = "subtype_binary"
SUBTYPES = [
    "endolysin",
    "depolymerase",
    "virion_associated_lytic_enzyme",
]


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    return add_common_args(
        parser,
        __file__,
        "data/tasks/subtype_binary",
    ).parse_args()


def main():
    args = parse_args()
    master_path = Path(args.master).expanduser().resolve()
    output_dir = Path(args.output_dir).expanduser().resolve()
    clean_output_dir(output_dir)

    master = load_master(master_path)
    global_rows = []

    for subtype in SUBTYPES:
        subdir = output_dir / subtype
        clean_output_dir(subdir)

        positive_mask = (
            master["sequence_state"].eq("P_core")
            & master["canonical_target_class"].eq(subtype)
        )
        negative_mask = master["sequence_state"].isin(["N_easy", "N_hard"])

        task = master.loc[positive_mask | negative_mask].copy()
        task["label"] = 0
        task.loc[positive_mask.loc[task.index], "label"] = 1
        task["class_name"] = task["label"].map(
            {1: subtype, 0: "non_lytic_background"}
        )
        task["training_role"] = task["label"].map({1: "P", 0: "N"})
        task["negative_stratum"] = task["sequence_state"].map(
            {"N_easy": "easy", "N_hard": "hard"}
        ).fillna("")

        dataset_path = subdir / "dataset.tsv.gz"
        write_dataset(task, dataset_path)

        rows = []
        add_group_counts(rows, task, "class_name", "class")
        add_group_counts(rows, task, "sequence_state", "source_state")
        write_summary(rows, subdir / "summary.tsv")

        write_manifest(
            subdir / "manifest.json",
            {
                "task_id": f"{TASK_ID}:{subtype}",
                "positive_definition": (
                    f"PhageLysData Core canonical_target_class == {subtype}"
                ),
                "negative_definition": "Union of N_easy and N_hard.",
                "other_lytic_classes_are_negatives": False,
                "multiple_nonexclusive_in_training": False,
                "other_phage_envelope_lytic_enzyme_in_training": False,
                "no_balancing": True,
                "no_homology_filtering": True,
                "master": str(master_path),
                "dataset": str(dataset_path),
                "n_examples": len(task),
            },
        )

        n_pos = int((task["label"] == 1).sum())
        n_neg = int((task["label"] == 0).sum())
        global_rows.extend(
            [
                {
                    "subtype": subtype,
                    "group": "positive",
                    "count": n_pos,
                },
                {
                    "subtype": subtype,
                    "group": "negative",
                    "count": n_neg,
                },
                {
                    "subtype": subtype,
                    "group": "total",
                    "count": len(task),
                },
            ]
        )

        print(f"{subtype:34s} P={n_pos:,} N={n_neg:,} total={len(task):,}")

    pd.DataFrame(global_rows).to_csv(
        output_dir / "summary.tsv",
        sep="\t",
        index=False,
    )


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Materialize challenge/evaluation sets excluded from ordinary supervised training.

Sets:
- core_multiple_nonexclusive
- core_other_lytic
- uniprot_lytic_related
- uniprot_ambiguous_structural
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from _task_utils import (
    add_common_args, add_group_counts, clean_output_dir, load_master,
    write_dataset, write_manifest, write_summary,
)


TASK_ID = "challenge_sets"


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    return add_common_args(
        parser,
        __file__,
        "data/tasks/challenge_sets",
    ).parse_args()


def main():
    args = parse_args()
    master_path = Path(args.master).expanduser().resolve()
    output_dir = Path(args.output_dir).expanduser().resolve()
    clean_output_dir(output_dir)

    master = load_master(master_path)

    definitions = {
        "core_multiple_nonexclusive": (
            master["sequence_state"].eq("P_core")
            & master["canonical_target_class"].eq("multiple_nonexclusive")
        ),
        "core_other_lytic": (
            master["sequence_state"].eq("P_core")
            & master["canonical_target_class"].eq(
                "other_phage_envelope_lytic_enzyme"
            )
        ),
        "uniprot_lytic_related": master["sequence_state"].eq(
            "challenge_lytic_related"
        ),
        "uniprot_ambiguous_structural": master["sequence_state"].eq(
            "challenge_ambiguous_structural"
        ),
    }

    global_rows = []

    for challenge_name, mask in definitions.items():
        subdir = output_dir / challenge_name
        clean_output_dir(subdir)

        task = master.loc[mask].copy()
        task["challenge_set"] = challenge_name
        task["training_role"] = "challenge_only"

        dataset_path = subdir / "dataset.tsv.gz"
        write_dataset(task, dataset_path)

        rows = [
            {
                "section": "challenge_set",
                "group": challenge_name,
                "count": len(task),
                "fraction": 1.0 if len(task) else 0.0,
            }
        ]
        if len(task):
            add_group_counts(rows, task, "sequence_state", "source_state")
            if challenge_name.startswith("core_"):
                add_group_counts(
                    rows,
                    task,
                    "canonical_target_class",
                    "canonical_target_class",
                )

        write_summary(rows, subdir / "summary.tsv")

        write_manifest(
            subdir / "manifest.json",
            {
                "task_id": f"{TASK_ID}:{challenge_name}",
                "training_role": "challenge_only",
                "no_training_default": True,
                "no_balancing": True,
                "no_homology_filtering": True,
                "master": str(master_path),
                "dataset": str(dataset_path),
                "n_examples": len(task),
            },
        )

        global_rows.append(
            {
                "challenge_set": challenge_name,
                "count": len(task),
            }
        )
        print(f"{challenge_name:34s} {len(task):>10,}")

    pd.DataFrame(global_rows).to_csv(
        output_dir / "summary.tsv",
        sep="\t",
        index=False,
    )


if __name__ == "__main__":
    main()

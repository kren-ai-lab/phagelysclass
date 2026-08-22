#!/usr/bin/env python3
"""BioSieve/MMseqs2 redundancy reduction for each subtype_binary task."""

from __future__ import annotations
import argparse
from pathlib import Path
from _redundancy_utils import (
    DEFAULT_THRESHOLDS, project_root_from_script, reduce_task
)

SUBTYPES = [
    "endolysin",
    "depolymerase",
    "virion_associated_lytic_enzyme",
]

def main():
    root = project_root_from_script(__file__)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-root",
        default=str(root / "data/tasks/subtype_binary"),
    )
    parser.add_argument(
        "--output-root",
        default=str(root / "data/redundancy/subtype_binary"),
    )
    parser.add_argument(
        "--subtypes",
        nargs="+",
        default=SUBTYPES,
        choices=SUBTYPES,
    )
    parser.add_argument(
        "--thresholds",
        nargs="+",
        type=int,
        default=DEFAULT_THRESHOLDS,
    )
    parser.add_argument("--threshold-param", default="auto")
    parser.add_argument(
        "--threads",
        type=int,
        default=max(1, (__import__("os").cpu_count() or 2) // 2),
    )
    parser.add_argument("--reuse-existing", action="store_true")
    parser.add_argument("--keep-work", action="store_true")
    args = parser.parse_args()

    for subtype in args.subtypes:
        print()
        reduce_task(
            task_name=f"subtype_binary/{subtype}",
            input_path=Path(args.input_root) / subtype / "dataset.tsv.gz",
            output_root=Path(args.output_root) / subtype,
            group_col="sequence_state",
            thresholds=args.thresholds,
            threshold_param=args.threshold_param,
            threads=args.threads,
            reuse_existing=args.reuse_existing,
            keep_work=args.keep_work,
        )

if __name__ == "__main__":
    main()

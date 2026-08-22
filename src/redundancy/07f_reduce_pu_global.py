#!/usr/bin/env python3
"""BioSieve/MMseqs2 redundancy reduction for conservative and broad PU tasks."""

from __future__ import annotations
import argparse
from pathlib import Path
from _redundancy_utils import (
    DEFAULT_THRESHOLDS, project_root_from_script, reduce_task
)

REGIMES = ["conservative", "broad"]

def main():
    root = project_root_from_script(__file__)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-root",
        default=str(root / "data/tasks/pu_global"),
    )
    parser.add_argument(
        "--output-root",
        default=str(root / "data/redundancy/pu_global"),
    )
    parser.add_argument(
        "--regimes",
        nargs="+",
        default=REGIMES,
        choices=REGIMES,
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

    for regime in args.regimes:
        print()
        reduce_task(
            task_name=f"pu_global/{regime}",
            input_path=Path(args.input_root) / regime / "dataset.tsv.gz",
            output_root=Path(args.output_root) / regime,
            group_col="sequence_state",
            thresholds=args.thresholds,
            threshold_param=args.threshold_param,
            threads=args.threads,
            reuse_existing=args.reuse_existing,
            keep_work=args.keep_work,
        )

if __name__ == "__main__":
    main()

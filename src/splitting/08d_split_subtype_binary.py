#!/usr/bin/env python3
"""Stage 08 split for all subtype_binary tasks."""

from __future__ import annotations
import argparse
import os
from pathlib import Path
from _split_utils import DEFAULT_THRESHOLDS, project_root_from_script, split_task_thresholds

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
        default=str(root / "data/redundancy/subtype_binary"),
    )
    parser.add_argument(
        "--output-root",
        default=str(root / "data/splits/subtype_binary"),
    )
    parser.add_argument("--subtypes", nargs="+", default=SUBTYPES, choices=SUBTYPES)
    parser.add_argument("--thresholds", nargs="+", type=int, default=DEFAULT_THRESHOLDS)
    parser.add_argument("--test-size", type=float, default=0.30)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--coverage", type=float, default=0.80)
    parser.add_argument("--cov-mode", type=int, default=0)
    parser.add_argument("--cluster-mode", type=int, default=1)
    parser.add_argument(
        "--threads",
        type=int,
        default=max(1, (os.cpu_count() or 2) // 2),
    )
    parser.add_argument("--reuse-existing", action="store_true")
    parser.add_argument("--keep-work", action="store_true")
    args = parser.parse_args()

    for subtype in args.subtypes:
        print()
        split_task_thresholds(
            task_name=f"subtype_binary/{subtype}",
            input_root=Path(args.input_root) / subtype,
            output_root=Path(args.output_root) / subtype,
            thresholds=args.thresholds,
            summary_cols=["class_name", "sequence_state", "negative_stratum"],
            test_size=args.test_size,
            seed=args.seed,
            coverage=args.coverage,
            cov_mode=args.cov_mode,
            cluster_mode=args.cluster_mode,
            threads=args.threads,
            reuse_existing=args.reuse_existing,
            keep_work=args.keep_work,
        )

if __name__ == "__main__":
    main()

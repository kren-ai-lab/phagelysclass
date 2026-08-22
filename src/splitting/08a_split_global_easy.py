#!/usr/bin/env python3
"""Stage 08 split for global_easy."""

from __future__ import annotations
import argparse
from pathlib import Path
from _split_utils import add_split_args, split_task_thresholds

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    add_split_args(
        parser,
        script_file=__file__,
        default_input_root="data/redundancy/global_easy",
        default_output_root="data/splits/global_easy",
    )
    args = parser.parse_args()
    split_task_thresholds(
        task_name="global_easy",
        input_root=Path(args.input_root),
        output_root=Path(args.output_root),
        thresholds=args.thresholds,
        summary_cols=["class_name", "sequence_state"],
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

#!/usr/bin/env python3
"""BioSieve/MMseqs2 redundancy reduction for the global_mixed task."""

from __future__ import annotations
import argparse
from pathlib import Path
from _redundancy_utils import add_redundancy_args, reduce_task

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    add_redundancy_args(
        parser,
        script_file=__file__,
        default_input="data/tasks/global_mixed/dataset.tsv.gz",
        default_output="data/redundancy/global_mixed",
    )
    args = parser.parse_args()
    reduce_task(
        task_name="global_mixed",
        input_path=Path(args.input),
        output_root=Path(args.output_dir),
        group_col="sequence_state",
        thresholds=args.thresholds,
        threshold_param=args.threshold_param,
        threads=args.threads,
        reuse_existing=args.reuse_existing,
        keep_work=args.keep_work,
    )

if __name__ == "__main__":
    main()

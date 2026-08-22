#!/usr/bin/env python3
"""Shared utilities for PhageLysClass task-definition scripts."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable

import pandas as pd

DEFAULT_MASTER = "data/processed/pnu_universes/pnu_sequence_master.tsv.gz"


def project_root_from_script(script_file: str) -> Path:
    script_dir = Path(script_file).resolve().parent
    return script_dir.parents[1]


def add_common_args(parser: argparse.ArgumentParser, script_file: str, default_output: str):
    root = project_root_from_script(script_file)
    parser.add_argument(
        "--master",
        default=str(root / DEFAULT_MASTER),
        help="Exact-sequence P/N/U master from 05_build_pnu_universes.py",
    )
    parser.add_argument(
        "--output-dir",
        default=str(root / default_output),
        help="Task output directory",
    )
    return parser


def load_master(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"P/N/U master not found: {path}")

    df = pd.read_csv(path, sep="\t", compression="infer", dtype=str).fillna("")

    required = {
        "entity_id",
        "sequence",
        "sequence_sha256",
        "source_universe",
        "pnu_role",
        "sequence_state",
        "stratum",
        "canonical_target_class",
    }
    missing = required - set(df.columns)
    if missing:
        raise KeyError(
            "P/N/U master missing required columns: "
            + ", ".join(sorted(missing))
        )

    if df["sequence"].duplicated().any():
        raise ValueError("P/N/U master contains duplicated exact sequences.")

    return df


def clean_output_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def write_dataset(df: pd.DataFrame, path: Path) -> None:
    if df["sequence"].duplicated().any():
        raise ValueError(f"Task dataset contains duplicated exact sequences: {path}")

    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False, compression="gzip")


def write_summary(
    rows: list[dict],
    path: Path,
) -> pd.DataFrame:
    summary = pd.DataFrame(rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(path, sep="\t", index=False)
    return summary


def add_group_counts(
    rows: list[dict],
    df: pd.DataFrame,
    column: str,
    section: str,
) -> None:
    total = len(df)
    for value, n in df[column].value_counts(dropna=False).items():
        rows.append(
            {
                "section": section,
                "group": str(value),
                "count": int(n),
                "fraction": (n / total) if total else 0.0,
            }
        )


def write_manifest(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def print_summary(title: str, df: pd.DataFrame, label_col: str | None = None) -> None:
    print("=" * 88)
    print(title)
    print("=" * 88)
    print(f"Examples: {len(df):,}")
    if label_col and label_col in df.columns:
        for value, n in df[label_col].value_counts().items():
            print(f"  {str(value):28s} {n:>10,}")
    print()

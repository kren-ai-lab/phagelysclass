#!/usr/bin/env python3
"""
01_extract_phagelys_raw.py

Extract the raw positive PhageLysData universe for PhageLysClass.

This script reads ONLY the released PhageLysData Core master table and
materializes sequence/ID/label datasets. It does not load representations,
perform homology filtering, balance classes, create train/test splits, or
construct external negative/unlabeled sets.

Outputs
-------
1. phagelysdata_core_raw.tsv.gz
   Full evidence-supported Core.

2. phagelysdata_single_class_raw.tsv.gz
   Clean single-class Endolysin / VAL / Depolymerase entities.

3. phagelysdata_multiple_nonexclusive_raw.tsv.gz
   Multiple-nonexclusive challenge universe.

4. phagelysdata_other_lytic_raw.tsv.gz
   Other phage-envelope lytic enzyme challenge universe.

5. phagelysdata_core_raw_summary.tsv
   Counts by canonical class and evidence tier.

The raw extraction preserves PhageLysData labels as released. Task-specific
relabeling is intentionally deferred to src/task_definitions/.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import pandas as pd


PRIMARY_SINGLE_CLASSES = {
    "endolysin",
    "virion_associated_lytic_enzyme",
    "depolymerase",
}

EXPECTED_CANONICAL_CLASSES = PRIMARY_SINGLE_CLASSES | {
    "multiple_nonexclusive",
    "other_phage_envelope_lytic_enzyme",
}

REQUIRED_COLUMNS = [
    "sequence_entity_id",
    "sequence",
    "canonical_target_class",
    "canonical_class_basis",
    "assigned_tier",
]


def resolve_release_root(cli_root: str | None) -> Path:
    """Resolve the PhageLysData release root."""
    candidates: list[Path] = []

    if cli_root:
        candidates.append(Path(cli_root).expanduser())

    env_root = os.environ.get("PHAGELYSDATA_ROOT")
    if env_root:
        candidates.append(Path(env_root).expanduser())

    cwd = Path.cwd()
    candidates.extend(
        [
            cwd / "PhageLysData_v1.0",
            cwd.parent / "PhageLysData_v1.0",
            cwd / "data_release" / "PhageLysData_v1.0",
            cwd.parent / "data_release" / "PhageLysData_v1.0",
            cwd.parent.parent / "data_release" / "PhageLysData_v1.0",
        ]
    )

    seen: set[Path] = set()
    checked: list[Path] = []

    for candidate in candidates:
        resolved = candidate.resolve()
        if resolved in seen:
            continue
        seen.add(resolved)
        checked.append(resolved)

        if (resolved / "core" / "phagelysdata_core.parquet").is_file():
            return resolved

    checked_text = "\n".join(f"  - {path}" for path in checked)
    raise FileNotFoundError(
        "Could not locate PhageLysData.\n\n"
        "Use --phagelysdata-root or set PHAGELYSDATA_ROOT to the directory "
        "containing core/phagelysdata_core.parquet.\n\n"
        f"Checked:\n{checked_text}"
    )


def validate_core(core: pd.DataFrame) -> None:
    """Validate the minimum biological contract required by this extraction."""
    missing_columns = [c for c in REQUIRED_COLUMNS if c not in core.columns]
    if missing_columns:
        raise KeyError(
            "PhageLysData Core is missing required columns: "
            + ", ".join(missing_columns)
        )

    if core.empty:
        raise ValueError("PhageLysData Core is empty.")

    for column in ["sequence_entity_id", "sequence", "canonical_target_class"]:
        if core[column].isna().any():
            n = int(core[column].isna().sum())
            raise ValueError(f"{column} contains {n:,} missing values.")

    if core["sequence_entity_id"].duplicated().any():
        n = int(core["sequence_entity_id"].duplicated().sum())
        raise ValueError(f"Found {n:,} duplicate sequence_entity_id values.")

    # PhageLysData entities are exact-sequence entities; duplicated sequences
    # should therefore not occur in the Core master.
    if core["sequence"].duplicated().any():
        n = int(core["sequence"].duplicated().sum())
        raise ValueError(f"Found {n:,} duplicated exact sequences in the Core.")

    observed_classes = set(core["canonical_target_class"].astype(str))
    unexpected = sorted(observed_classes - EXPECTED_CANONICAL_CLASSES)
    missing_expected = sorted(EXPECTED_CANONICAL_CLASSES - observed_classes)

    if unexpected:
        raise ValueError(
            "Unexpected canonical_target_class values detected: "
            + ", ".join(unexpected)
        )

    if missing_expected:
        raise ValueError(
            "Expected canonical_target_class values are absent: "
            + ", ".join(missing_expected)
        )


def write_tsv_gz(df: pd.DataFrame, path: Path) -> None:
    df.to_csv(path, sep="\t", index=False, compression="gzip")


def build_summary(raw: pd.DataFrame) -> pd.DataFrame:
    """Return a compact class × evidence-tier summary."""
    summary = (
        raw.groupby(
            ["canonical_target_class", "assigned_tier"],
            dropna=False,
            as_index=False,
        )
        .agg(n_sequences=("sequence_entity_id", "nunique"))
        .sort_values(
            ["canonical_target_class", "assigned_tier"],
            kind="stable",
        )
        .reset_index(drop=True)
    )
    return summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Extract sequence IDs, sequences, and released biological labels "
            "from the PhageLysData Core."
        )
    )
    parser.add_argument(
        "--phagelysdata-root",
        default=None,
        help=(
            "Path to the PhageLysData release root. If omitted, the script "
            "uses PHAGELYSDATA_ROOT or common relative locations."
        ),
    )
    parser.add_argument(
        "--output-dir",
        default="data/raw/phagelysdata",
        help="Output directory (default: data/raw/phagelysdata).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    release_root = resolve_release_root(args.phagelysdata_root)
    core_path = release_root / "core" / "phagelysdata_core.parquet"
    output_dir = Path(args.output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 88)
    print("PhageLysClass — raw PhageLysData extraction")
    print("=" * 88)
    print(f"Release root : {release_root}")
    print(f"Core master  : {core_path}")
    print(f"Output dir   : {output_dir}")
    print()

    core = pd.read_parquet(core_path)
    validate_core(core)

    # Keep only sequence identity and biological labels/evidence metadata.
    raw = (
        core[REQUIRED_COLUMNS]
        .copy()
        .sort_values("sequence_entity_id", kind="stable")
        .reset_index(drop=True)
    )

    single = raw.loc[
        raw["canonical_target_class"].isin(PRIMARY_SINGLE_CLASSES)
    ].copy()

    multi = raw.loc[
        raw["canonical_target_class"].eq("multiple_nonexclusive")
    ].copy()

    other = raw.loc[
        raw["canonical_target_class"].eq(
            "other_phage_envelope_lytic_enzyme"
        )
    ].copy()

    # Defensive partition check.
    partition_ids = set(single["sequence_entity_id"])
    partition_ids |= set(multi["sequence_entity_id"])
    partition_ids |= set(other["sequence_entity_id"])

    if partition_ids != set(raw["sequence_entity_id"]):
        missing = len(set(raw["sequence_entity_id"]) - partition_ids)
        raise RuntimeError(
            f"Raw partitioning failed: {missing:,} Core entities were not assigned."
        )

    outputs = {
        "core": output_dir / "phagelysdata_core_raw.tsv.gz",
        "single": output_dir / "phagelysdata_single_class_raw.tsv.gz",
        "multi": output_dir / "phagelysdata_multiple_nonexclusive_raw.tsv.gz",
        "other": output_dir / "phagelysdata_other_lytic_raw.tsv.gz",
        "summary": output_dir / "phagelysdata_core_raw_summary.tsv",
    }

    write_tsv_gz(raw, outputs["core"])
    write_tsv_gz(single, outputs["single"])
    write_tsv_gz(multi, outputs["multi"])
    write_tsv_gz(other, outputs["other"])

    summary = build_summary(raw)
    summary.to_csv(outputs["summary"], sep="\t", index=False)

    print("Extraction complete")
    print("-" * 88)
    print(f"Core total               : {len(raw):,}")
    print(f"Single-class primary     : {len(single):,}")
    print(f"Multiple non-exclusive   : {len(multi):,}")
    print(f"Other lytic / unresolved : {len(other):,}")
    print()

    print("Canonical class counts")
    print("-" * 88)
    counts = raw["canonical_target_class"].value_counts()
    for label, n in counts.items():
        print(f"{label:42s} {n:>8,}")

    print()
    print("Written files")
    print("-" * 88)
    for path in outputs.values():
        print(path)

    print()
    print("PASS — no representations, balancing, homology filtering, or splits applied.")


if __name__ == "__main__":
    main()

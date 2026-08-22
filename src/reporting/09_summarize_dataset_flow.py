#!/usr/bin/env python3
"""
09_summarize_dataset_flow.py

Create a paper-oriented summary of the first major PhageLysClass data-design block:

    Stage 06 — raw task definitions
    Stage 07 — redundancy reduction (H90/H70/H50/H30)
    Stage 08 — homology-aware train/test splitting

The script is intentionally read-only. It never modifies task, redundancy, or
split datasets.

It can be rerun as Stage 08 outputs become available. Missing downstream
outputs are reported but do not stop execution unless --strict is supplied.

Default outputs
---------------
data/reports/stage1_dataset_flow/
├── stage1_task_raw_summary.tsv
├── stage1_dataset_flow_overview.tsv
├── stage1_dataset_flow_long.tsv
├── stage1_challenge_sets.tsv
├── stage1_paper_table.tsv
├── stage1_paper_table.md
├── missing_or_incomplete_inputs.tsv
└── summary_manifest.json
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import pandas as pd


THRESHOLDS = [90, 70, 50, 30]


TASKS = [
    {
        "task_family": "global_binary",
        "task_id": "global_easy",
        "display_name": "Global easy",
        "learning_type": "supervised_binary",
        "task_path": "global_easy",
        "primary_group": "class_name",
        "secondary_group": "sequence_state",
    },
    {
        "task_family": "global_binary",
        "task_id": "global_hard",
        "display_name": "Global hard",
        "learning_type": "supervised_binary",
        "task_path": "global_hard",
        "primary_group": "class_name",
        "secondary_group": "sequence_state",
    },
    {
        "task_family": "global_binary",
        "task_id": "global_mixed",
        "display_name": "Global mixed",
        "learning_type": "supervised_binary",
        "task_path": "global_mixed",
        "primary_group": "class_name",
        "secondary_group": "sequence_state",
    },
    {
        "task_family": "subtype_binary",
        "task_id": "subtype_endolysin",
        "display_name": "Subtype — endolysin",
        "learning_type": "supervised_binary",
        "task_path": "subtype_binary/endolysin",
        "primary_group": "class_name",
        "secondary_group": "sequence_state",
    },
    {
        "task_family": "subtype_binary",
        "task_id": "subtype_depolymerase",
        "display_name": "Subtype — depolymerase",
        "learning_type": "supervised_binary",
        "task_path": "subtype_binary/depolymerase",
        "primary_group": "class_name",
        "secondary_group": "sequence_state",
    },
    {
        "task_family": "subtype_binary",
        "task_id": "subtype_val",
        "display_name": "Subtype — VAL",
        "learning_type": "supervised_binary",
        "task_path": "subtype_binary/virion_associated_lytic_enzyme",
        "primary_group": "class_name",
        "secondary_group": "sequence_state",
    },
    {
        "task_family": "multiclass",
        "task_id": "multiclass_clean",
        "display_name": "Multiclass clean",
        "learning_type": "supervised_multiclass",
        "task_path": "multiclass_clean",
        "primary_group": "class_name",
        "secondary_group": "canonical_target_class",
    },
    {
        "task_family": "pu",
        "task_id": "pu_conservative",
        "display_name": "PU conservative",
        "learning_type": "positive_unlabeled",
        "task_path": "pu_global/conservative",
        "primary_group": "pu_role",
        "secondary_group": "sequence_state",
    },
    {
        "task_family": "pu",
        "task_id": "pu_broad",
        "display_name": "PU broad",
        "learning_type": "positive_unlabeled",
        "task_path": "pu_global/broad",
        "primary_group": "pu_role",
        "secondary_group": "sequence_state",
    },
]


CHALLENGE_SETS = [
    {
        "challenge_id": "core_multiple_nonexclusive",
        "display_name": "Core multiple non-exclusive",
    },
    {
        "challenge_id": "core_other_lytic",
        "display_name": "Core other lytic",
    },
    {
        "challenge_id": "uniprot_lytic_related",
        "display_name": "UniProt lytic-related",
    },
    {
        "challenge_id": "uniprot_ambiguous_structural",
        "display_name": "UniProt ambiguous structural",
    },
]


def parse_args() -> argparse.Namespace:
    script_dir = Path(__file__).resolve().parent
    project_root = script_dir.parents[1]

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--tasks-root",
        default=str(project_root / "data/tasks"),
        help="Stage 06 task root.",
    )
    parser.add_argument(
        "--redundancy-root",
        default=str(project_root / "data/redundancy"),
        help="Stage 07 redundancy root.",
    )
    parser.add_argument(
        "--splits-root",
        default=str(project_root / "data/splits"),
        help="Stage 08 split root.",
    )
    parser.add_argument(
        "--output-dir",
        default=str(project_root / "data/reports/stage1_dataset_flow"),
        help="Output directory.",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help=(
            "Fail if any registered task/redundancy/split input is missing or "
            "if cross-stage counts do not agree."
        ),
    )
    return parser.parse_args()


def read_tsv_gz(path: Path) -> pd.DataFrame:
    return pd.read_csv(
        path,
        sep="\t",
        compression="infer",
        dtype=str,
        low_memory=False,
    ).fillna("")


def count_dict(df: pd.DataFrame, col: str | None) -> dict[str, int]:
    if not col or col not in df.columns:
        return {}
    counts = df[col].astype(str).value_counts(dropna=False)
    return {str(k): int(v) for k, v in counts.items()}


def composition_string(counts: dict[str, int]) -> str:
    if not counts:
        return ""
    items = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))
    return "; ".join(f"{name}={n}" for name, n in items)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def add_long_counts(
    rows: list[dict[str, Any]],
    *,
    task: dict[str, Any],
    stage: str,
    regime: str,
    partition: str,
    df: pd.DataFrame,
    group_type: str,
    group_col: str | None,
) -> None:
    if not group_col or group_col not in df.columns:
        return

    counts = count_dict(df, group_col)
    for group, n in counts.items():
        rows.append(
            {
                "task_family": task["task_family"],
                "task_id": task["task_id"],
                "display_name": task["display_name"],
                "learning_type": task["learning_type"],
                "stage": stage,
                "regime": regime,
                "partition": partition,
                "group_type": group_type,
                "group_column": group_col,
                "group": group,
                "count": n,
                "partition_total": len(df),
                "fraction_within_partition": n / len(df) if len(df) else 0.0,
            }
        )


def get_split_metadata(split_root: Path, threshold: int) -> dict[str, Any]:
    overview_path = split_root / "split_overview.tsv"
    if not overview_path.is_file():
        return {}

    overview = pd.read_csv(overview_path, sep="\t")
    if "threshold" not in overview.columns:
        return {}

    row = overview.loc[
        pd.to_numeric(overview["threshold"], errors="coerce").eq(threshold)
    ]
    if row.empty:
        return {}

    r = row.iloc[0]
    return {
        "n_homology_clusters": _safe_int(r.get("n_homology_clusters")),
        "shared_train_test_clusters": _safe_int(
            r.get("shared_train_test_clusters")
        ),
        "n_warnings": _safe_int(r.get("n_warnings")),
    }


def _safe_int(value: Any) -> int | None:
    if value is None or pd.isna(value) or str(value).strip() == "":
        return None
    try:
        return int(float(value))
    except Exception:
        return None


def markdown_table(df: pd.DataFrame) -> str:
    """
    Small dependency-free Markdown table writer.
    """
    display = df.copy().fillna("")
    columns = list(display.columns)

    def esc(value: Any) -> str:
        return str(value).replace("|", "\\|").replace("\n", " ")

    lines = [
        "| " + " | ".join(esc(c) for c in columns) + " |",
        "| " + " | ".join("---" for _ in columns) + " |",
    ]
    for row in display.itertuples(index=False, name=None):
        lines.append("| " + " | ".join(esc(v) for v in row) + " |")
    return "\n".join(lines)


def main() -> None:
    args = parse_args()

    tasks_root = Path(args.tasks_root).expanduser().resolve()
    redundancy_root = Path(args.redundancy_root).expanduser().resolve()
    splits_root = Path(args.splits_root).expanduser().resolve()
    output_dir = Path(args.output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    overview_rows: list[dict[str, Any]] = []
    raw_rows: list[dict[str, Any]] = []
    long_rows: list[dict[str, Any]] = []
    missing_rows: list[dict[str, Any]] = []

    input_hashes: dict[str, str] = {}

    print("=" * 100)
    print("PhageLysClass — Stage 1 dataset-flow summary")
    print("=" * 100)
    print(f"Tasks root      : {tasks_root}")
    print(f"Redundancy root : {redundancy_root}")
    print(f"Splits root     : {splits_root}")
    print(f"Output dir      : {output_dir}")
    print()

    for task in TASKS:
        task_path = Path(task["task_path"])
        raw_path = tasks_root / task_path / "dataset.tsv.gz"

        if not raw_path.is_file():
            missing_rows.append(
                {
                    "task_id": task["task_id"],
                    "stage": "raw_task",
                    "regime": "RAW",
                    "path": str(raw_path),
                    "issue": "missing",
                }
            )
            continue

        raw = read_tsv_gz(raw_path)
        raw_total = len(raw)
        raw_primary = count_dict(raw, task["primary_group"])
        raw_secondary = count_dict(raw, task["secondary_group"])

        input_hashes[str(raw_path)] = sha256_file(raw_path)

        raw_rows.append(
            {
                "task_family": task["task_family"],
                "task_id": task["task_id"],
                "display_name": task["display_name"],
                "learning_type": task["learning_type"],
                "raw_total": raw_total,
                "primary_group_column": task["primary_group"],
                "primary_composition": composition_string(raw_primary),
                "secondary_group_column": task["secondary_group"],
                "secondary_composition": composition_string(raw_secondary),
            }
        )

        add_long_counts(
            long_rows,
            task=task,
            stage="raw_task",
            regime="RAW",
            partition="all",
            df=raw,
            group_type="primary",
            group_col=task["primary_group"],
        )
        if task["secondary_group"] != task["primary_group"]:
            add_long_counts(
                long_rows,
                task=task,
                stage="raw_task",
                regime="RAW",
                partition="all",
                df=raw,
                group_type="secondary",
                group_col=task["secondary_group"],
            )

        for threshold in THRESHOLDS:
            regime = f"H{threshold}"
            red_path = redundancy_root / task_path / regime / "dataset.tsv.gz"

            if not red_path.is_file():
                missing_rows.append(
                    {
                        "task_id": task["task_id"],
                        "stage": "redundancy",
                        "regime": regime,
                        "path": str(red_path),
                        "issue": "missing",
                    }
                )
                continue

            reduced = read_tsv_gz(red_path)
            reduced_total = len(reduced)
            reduced_primary = count_dict(reduced, task["primary_group"])
            reduced_secondary = count_dict(reduced, task["secondary_group"])

            input_hashes[str(red_path)] = sha256_file(red_path)

            add_long_counts(
                long_rows,
                task=task,
                stage="redundancy",
                regime=regime,
                partition="all",
                df=reduced,
                group_type="primary",
                group_col=task["primary_group"],
            )
            if task["secondary_group"] != task["primary_group"]:
                add_long_counts(
                    long_rows,
                    task=task,
                    stage="redundancy",
                    regime=regime,
                    partition="all",
                    df=reduced,
                    group_type="secondary",
                    group_col=task["secondary_group"],
                )

            split_dir = splits_root / task_path / regime
            train_path = split_dir / "train.tsv.gz"
            test_path = split_dir / "test.tsv.gz"

            train = None
            test = None
            split_complete = train_path.is_file() and test_path.is_file()

            train_total = None
            test_total = None
            actual_test_fraction = None
            train_primary_s = ""
            test_primary_s = ""
            train_secondary_s = ""
            test_secondary_s = ""
            entity_overlap_count = None
            split_count_matches_redundancy = None

            split_meta = get_split_metadata(
                splits_root / task_path,
                threshold,
            )

            if split_complete:
                train = read_tsv_gz(train_path)
                test = read_tsv_gz(test_path)

                train_total = len(train)
                test_total = len(test)
                actual_test_fraction = (
                    test_total / reduced_total if reduced_total else 0.0
                )

                train_primary_s = composition_string(
                    count_dict(train, task["primary_group"])
                )
                test_primary_s = composition_string(
                    count_dict(test, task["primary_group"])
                )
                train_secondary_s = composition_string(
                    count_dict(train, task["secondary_group"])
                )
                test_secondary_s = composition_string(
                    count_dict(test, task["secondary_group"])
                )

                add_long_counts(
                    long_rows,
                    task=task,
                    stage="split",
                    regime=regime,
                    partition="train",
                    df=train,
                    group_type="primary",
                    group_col=task["primary_group"],
                )
                add_long_counts(
                    long_rows,
                    task=task,
                    stage="split",
                    regime=regime,
                    partition="test",
                    df=test,
                    group_type="primary",
                    group_col=task["primary_group"],
                )

                if task["secondary_group"] != task["primary_group"]:
                    add_long_counts(
                        long_rows,
                        task=task,
                        stage="split",
                        regime=regime,
                        partition="train",
                        df=train,
                        group_type="secondary",
                        group_col=task["secondary_group"],
                    )
                    add_long_counts(
                        long_rows,
                        task=task,
                        stage="split",
                        regime=regime,
                        partition="test",
                        df=test,
                        group_type="secondary",
                        group_col=task["secondary_group"],
                    )

                if "entity_id" in train.columns and "entity_id" in test.columns:
                    overlap = set(train["entity_id"]) & set(test["entity_id"])
                    entity_overlap_count = len(overlap)

                split_count_matches_redundancy = (
                    train_total + test_total == reduced_total
                )

                if not split_count_matches_redundancy:
                    missing_rows.append(
                        {
                            "task_id": task["task_id"],
                            "stage": "split",
                            "regime": regime,
                            "path": str(split_dir),
                            "issue": (
                                f"train+test={train_total+test_total} "
                                f"!= reduced={reduced_total}"
                            ),
                        }
                    )

                if entity_overlap_count not in (None, 0):
                    missing_rows.append(
                        {
                            "task_id": task["task_id"],
                            "stage": "split",
                            "regime": regime,
                            "path": str(split_dir),
                            "issue": f"train_test_entity_overlap={entity_overlap_count}",
                        }
                    )
            else:
                missing_rows.append(
                    {
                        "task_id": task["task_id"],
                        "stage": "split",
                        "regime": regime,
                        "path": str(split_dir),
                        "issue": "train/test not yet complete",
                    }
                )

            overview_rows.append(
                {
                    "task_family": task["task_family"],
                    "task_id": task["task_id"],
                    "display_name": task["display_name"],
                    "learning_type": task["learning_type"],
                    "regime": regime,
                    "raw_total": raw_total,
                    "reduced_total": reduced_total,
                    "retained_fraction": (
                        reduced_total / raw_total if raw_total else 0.0
                    ),
                    "removed_by_redundancy": raw_total - reduced_total,
                    "train_total": train_total,
                    "test_total": test_total,
                    "actual_test_fraction": actual_test_fraction,
                    "n_homology_clusters": split_meta.get(
                        "n_homology_clusters"
                    ),
                    "shared_train_test_clusters": split_meta.get(
                        "shared_train_test_clusters"
                    ),
                    "split_warnings": split_meta.get("n_warnings"),
                    "train_test_entity_overlap": entity_overlap_count,
                    "split_count_matches_redundancy": (
                        split_count_matches_redundancy
                    ),
                    "primary_group_column": task["primary_group"],
                    "raw_primary_composition": composition_string(
                        raw_primary
                    ),
                    "reduced_primary_composition": composition_string(
                        reduced_primary
                    ),
                    "train_primary_composition": train_primary_s,
                    "test_primary_composition": test_primary_s,
                    "secondary_group_column": task["secondary_group"],
                    "raw_secondary_composition": composition_string(
                        raw_secondary
                    ),
                    "reduced_secondary_composition": composition_string(
                        reduced_secondary
                    ),
                    "train_secondary_composition": train_secondary_s,
                    "test_secondary_composition": test_secondary_s,
                }
            )

    # Challenge sets: intentionally raw-only.
    challenge_rows = []
    challenge_root = tasks_root / "challenge_sets"

    for challenge in CHALLENGE_SETS:
        path = (
            challenge_root
            / challenge["challenge_id"]
            / "dataset.tsv.gz"
        )

        if not path.is_file():
            missing_rows.append(
                {
                    "task_id": challenge["challenge_id"],
                    "stage": "challenge_raw",
                    "regime": "RAW",
                    "path": str(path),
                    "issue": "missing",
                }
            )
            continue

        df = read_tsv_gz(path)
        challenge_rows.append(
            {
                "challenge_id": challenge["challenge_id"],
                "display_name": challenge["display_name"],
                "count": len(df),
                "sequence_state_composition": composition_string(
                    count_dict(df, "sequence_state")
                ),
                "canonical_target_class_composition": composition_string(
                    count_dict(df, "canonical_target_class")
                ),
                "redundancy_reduced": False,
                "split": False,
            }
        )

    raw_df = pd.DataFrame(raw_rows)
    overview_df = pd.DataFrame(overview_rows)
    long_df = pd.DataFrame(long_rows)
    challenge_df = pd.DataFrame(challenge_rows)
    missing_df = pd.DataFrame(missing_rows)

    # Stable ordering.
    task_order = {task["task_id"]: i for i, task in enumerate(TASKS)}
    threshold_order = {"H90": 0, "H70": 1, "H50": 2, "H30": 3}

    if not raw_df.empty:
        raw_df["_task_order"] = raw_df["task_id"].map(task_order)
        raw_df = raw_df.sort_values("_task_order").drop(columns="_task_order")

    if not overview_df.empty:
        overview_df["_task_order"] = overview_df["task_id"].map(task_order)
        overview_df["_threshold_order"] = overview_df["regime"].map(
            threshold_order
        )
        overview_df = overview_df.sort_values(
            ["_task_order", "_threshold_order"]
        ).drop(columns=["_task_order", "_threshold_order"])

    if not long_df.empty:
        long_df["_task_order"] = long_df["task_id"].map(task_order)
        long_df["_threshold_order"] = long_df["regime"].map(
            {"RAW": -1, **threshold_order}
        )
        long_df = long_df.sort_values(
            [
                "_task_order",
                "_threshold_order",
                "stage",
                "partition",
                "group_type",
                "group",
            ]
        ).drop(columns=["_task_order", "_threshold_order"])

    raw_path = output_dir / "stage1_task_raw_summary.tsv"
    overview_path = output_dir / "stage1_dataset_flow_overview.tsv"
    long_path = output_dir / "stage1_dataset_flow_long.tsv"
    challenge_path = output_dir / "stage1_challenge_sets.tsv"
    missing_path = output_dir / "missing_or_incomplete_inputs.tsv"

    raw_df.to_csv(raw_path, sep="\t", index=False)
    overview_df.to_csv(overview_path, sep="\t", index=False)
    long_df.to_csv(long_path, sep="\t", index=False)
    challenge_df.to_csv(challenge_path, sep="\t", index=False)
    missing_df.to_csv(missing_path, sep="\t", index=False)

    # Concise paper-oriented table.
    paper_columns = [
        "display_name",
        "regime",
        "raw_total",
        "reduced_total",
        "retained_fraction",
        "train_total",
        "test_total",
        "actual_test_fraction",
        "reduced_primary_composition",
        "train_primary_composition",
        "test_primary_composition",
        "n_homology_clusters",
        "shared_train_test_clusters",
    ]

    paper_df = overview_df[paper_columns].copy() if not overview_df.empty else pd.DataFrame()

    if not paper_df.empty:
        paper_df = paper_df.rename(
            columns={
                "display_name": "Task",
                "regime": "Regime",
                "raw_total": "Raw n",
                "reduced_total": "Reduced n",
                "retained_fraction": "Retained fraction",
                "train_total": "Train n",
                "test_total": "Test n",
                "actual_test_fraction": "Test fraction",
                "reduced_primary_composition": "Reduced composition",
                "train_primary_composition": "Train composition",
                "test_primary_composition": "Test composition",
                "n_homology_clusters": "Homology clusters",
                "shared_train_test_clusters": "Shared train/test clusters",
            }
        )

        for col in ["Retained fraction", "Test fraction"]:
            paper_df[col] = pd.to_numeric(
                paper_df[col], errors="coerce"
            ).map(
                lambda x: "" if pd.isna(x) else f"{100*x:.2f}%"
            )

    paper_tsv_path = output_dir / "stage1_paper_table.tsv"
    paper_md_path = output_dir / "stage1_paper_table.md"

    paper_df.to_csv(paper_tsv_path, sep="\t", index=False)

    md_parts = [
        "# PhageLysClass — Stage 1 dataset-design summary",
        "",
        (
            "This table summarizes raw task size, redundancy-controlled "
            "datasets, and available homology-aware outer train/test splits."
        ),
        "",
        markdown_table(paper_df) if not paper_df.empty else "_No rows available._",
        "",
        "## Challenge sets",
        "",
        markdown_table(challenge_df) if not challenge_df.empty else "_No challenge sets available._",
        "",
        "## Notes",
        "",
        "- H90/H70/H50/H30 datasets are generated independently from each raw task.",
        "- Redundancy reduction is performed within frozen biological strata.",
        "- Stage 08 splitting uses whole-task homology components.",
        "- `Shared train/test clusters = 0` is the expected leakage-control result.",
        "- Missing split cells indicate Stage 08 outputs that have not yet been materialized.",
    ]
    paper_md_path.write_text("\n".join(md_parts) + "\n", encoding="utf-8")

    manifest = {
        "report": "PhageLysClass Stage 1 dataset-flow summary",
        "tasks_root": str(tasks_root),
        "redundancy_root": str(redundancy_root),
        "splits_root": str(splits_root),
        "thresholds": THRESHOLDS,
        "n_registered_tasks": len(TASKS),
        "n_raw_tasks_found": len(raw_df),
        "n_regime_rows_found": len(overview_df),
        "n_challenge_sets_found": len(challenge_df),
        "n_missing_or_incomplete": len(missing_df),
        "input_sha256": input_hashes,
        "outputs": {
            "raw_summary": str(raw_path),
            "dataset_flow_overview": str(overview_path),
            "dataset_flow_long": str(long_path),
            "challenge_sets": str(challenge_path),
            "paper_table_tsv": str(paper_tsv_path),
            "paper_table_md": str(paper_md_path),
            "missing_inputs": str(missing_path),
        },
        "design_note": (
            "Read-only aggregation of Stages 06-08. No datasets, labels, "
            "redundancy mappings, or split assignments are modified."
        ),
    }

    manifest_path = output_dir / "summary_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2),
        encoding="utf-8",
    )

    print("Outputs")
    print("-" * 100)
    print(f"Raw tasks          : {raw_path}")
    print(f"Flow overview      : {overview_path}")
    print(f"Long-form counts   : {long_path}")
    print(f"Challenge sets     : {challenge_path}")
    print(f"Paper table TSV    : {paper_tsv_path}")
    print(f"Paper table MD     : {paper_md_path}")
    print(f"Missing/incomplete : {missing_path}")
    print(f"Manifest           : {manifest_path}")
    print()

    if len(missing_df):
        print(
            f"WARNING — {len(missing_df):,} expected inputs are missing or "
            "incomplete. This is allowed for an in-progress Stage 08 run."
        )
        if args.strict:
            raise SystemExit(2)
    else:
        print("PASS — complete Stage 06→07→08 dataset-flow summary.")

    # Strict consistency checks.
    if args.strict and not overview_df.empty:
        bad_split = overview_df.loc[
            overview_df["split_count_matches_redundancy"].eq(False)
        ]
        bad_overlap = overview_df.loc[
            pd.to_numeric(
                overview_df["train_test_entity_overlap"],
                errors="coerce",
            ).fillna(0).gt(0)
        ]
        bad_leak = overview_df.loc[
            pd.to_numeric(
                overview_df["shared_train_test_clusters"],
                errors="coerce",
            ).fillna(0).gt(0)
        ]
        if len(bad_split) or len(bad_overlap) or len(bad_leak):
            raise SystemExit(3)


if __name__ == "__main__":
    main()

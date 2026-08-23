#!/usr/bin/env python3
"""
09_summarize_factorial_splits.py

Read-only audit and paper-oriented summary of the full Stage-08 R × S design.

Expected canonical structure:
data/splits/<task_path>/Rxx_Syy/
    train.tsv.gz
    test.tsv.gz
    dataset_with_split.tsv.gz
    clusters.tsv.gz
    split_summary.tsv
    manifest.json

Outputs:
data/reports/stage08_factorial/
    factorial_split_summary.tsv
    factorial_split_composition_long.tsv
    factorial_split_matrix.tsv
    factorial_split_matrix.md
    missing_or_failed_cells.tsv
    freeze_manifest.json
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import pandas as pd


R_THRESHOLDS = [90, 70, 50, 30]
S_THRESHOLDS = [90, 70, 50, 30]

TASKS = [
    ("global_easy", "global_easy", "class_name"),
    ("global_hard", "global_hard", "class_name"),
    ("global_mixed", "global_mixed", "class_name"),
    ("subtype_endolysin", "subtype_binary/endolysin", "class_name"),
    ("subtype_depolymerase", "subtype_binary/depolymerase", "class_name"),
    ("subtype_val", "subtype_binary/virion_associated_lytic_enzyme", "class_name"),
    ("multiclass_clean", "multiclass_clean", "class_name"),
    ("pu_conservative", "pu_global/conservative", "pu_role"),
    ("pu_broad", "pu_global/broad", "pu_role"),
]


def project_root_from_script(script_file: str) -> Path:
    return Path(script_file).resolve().parent.parents[1]


def parse_args():
    root = project_root_from_script(__file__)
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--splits-root", default=str(root / "data/splits"))
    p.add_argument(
        "--output-dir",
        default=str(root / "data/reports/stage08_factorial"),
    )
    p.add_argument("--strict", action="store_true")
    return p.parse_args()


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def markdown_table(df: pd.DataFrame) -> str:
    x = df.fillna("").copy()
    cols = list(x.columns)
    lines = [
        "| " + " | ".join(cols) + " |",
        "| " + " | ".join("---" for _ in cols) + " |",
    ]
    for row in x.itertuples(index=False, name=None):
        lines.append(
            "| "
            + " | ".join(str(v).replace("|", "\\|") for v in row)
            + " |"
        )
    return "\n".join(lines)


def add_composition_rows(
    rows: list[dict[str, Any]],
    *,
    task_id: str,
    r: int,
    s: int,
    partition: str,
    df: pd.DataFrame,
    group_col: str,
):
    if group_col not in df.columns:
        return
    counts = df[group_col].astype(str).value_counts(dropna=False)
    for group, n in counts.items():
        rows.append({
            "task_id": task_id,
            "redundancy_threshold": r,
            "split_threshold": s,
            "rs_id": f"R{r}_S{s}",
            "partition": partition,
            "group_column": group_col,
            "group": group if str(group).strip() else "NA",
            "count": int(n),
            "partition_total": int(len(df)),
            "fraction": float(n / len(df)) if len(df) else 0.0,
        })


def main():
    args = parse_args()
    splits_root = Path(args.splits_root).expanduser().resolve()
    out = Path(args.output_dir).expanduser().resolve()
    out.mkdir(parents=True, exist_ok=True)

    summary_rows = []
    long_rows = []
    missing = []
    hashes = {}

    for task_id, task_path, group_col in TASKS:
        for r in R_THRESHOLDS:
            for s in S_THRESHOLDS:
                rs_id = f"R{r}_S{s}"
                cell = splits_root / task_path / rs_id
                train_path = cell / "train.tsv.gz"
                test_path = cell / "test.tsv.gz"
                cluster_path = cell / "clusters.tsv.gz"
                manifest_path = cell / "manifest.json"

                required = [train_path, test_path, cluster_path, manifest_path]
                if not all(p.is_file() for p in required):
                    missing.append({
                        "task_id": task_id,
                        "rs_id": rs_id,
                        "issue": "missing required output",
                        "path": str(cell),
                    })
                    continue

                train = pd.read_csv(
                    train_path, sep="\t", compression="infer",
                    dtype=str, low_memory=False,
                ).fillna("")
                test = pd.read_csv(
                    test_path, sep="\t", compression="infer",
                    dtype=str, low_memory=False,
                ).fillna("")
                clusters = pd.read_csv(
                    cluster_path, sep="\t", compression="infer",
                    dtype=str, low_memory=False,
                ).fillna("")
                manifest = json.loads(
                    manifest_path.read_text(encoding="utf-8")
                )

                entity_col = next(
                    (c for c in ["entity_id", "sequence_uid", "sequence_sha256"] if c in train.columns),
                    None,
                )
                if entity_col is None:
                    missing.append({
                        "task_id": task_id,
                        "rs_id": rs_id,
                        "issue": "cannot resolve entity ID column",
                        "path": str(cell),
                    })
                    continue

                entity_overlap = len(
                    set(train[entity_col]) & set(test[entity_col])
                )

                cluster_map = clusters.set_index("entity_id")["cluster_id"].to_dict()
                train_clusters = {
                    cluster_map.get(x) for x in train[entity_col]
                }
                test_clusters = {
                    cluster_map.get(x) for x in test[entity_col]
                }
                train_clusters.discard(None)
                test_clusters.discard(None)
                shared_clusters = len(train_clusters & test_clusters)

                status = (
                    "PASS"
                    if (
                        manifest.get("status") == "PASS"
                        and entity_overlap == 0
                        and shared_clusters == 0
                        and len(train) + len(test) == manifest["audit"]["n_total"]
                    )
                    else "FAIL"
                )

                row = {
                    "task_id": task_id,
                    "redundancy_threshold": r,
                    "split_threshold": s,
                    "rs_id": rs_id,
                    "n_total": len(train) + len(test),
                    "n_train": len(train),
                    "n_test": len(test),
                    "actual_test_fraction": (
                        len(test) / (len(train) + len(test))
                        if len(train) + len(test) else 0.0
                    ),
                    "n_homology_clusters": len(set(cluster_map.values())),
                    "train_test_entity_overlap": entity_overlap,
                    "shared_train_test_clusters": shared_clusters,
                    "status": status,
                }
                summary_rows.append(row)

                add_composition_rows(
                    long_rows,
                    task_id=task_id, r=r, s=s,
                    partition="train", df=train,
                    group_col=group_col,
                )
                add_composition_rows(
                    long_rows,
                    task_id=task_id, r=r, s=s,
                    partition="test", df=test,
                    group_col=group_col,
                )

                for p in required:
                    hashes[str(p)] = sha256_file(p)

                if status != "PASS":
                    missing.append({
                        "task_id": task_id,
                        "rs_id": rs_id,
                        "issue": "audit failed",
                        "path": str(cell),
                    })

    summary = pd.DataFrame(summary_rows)
    long_df = pd.DataFrame(long_rows)
    missing_df = pd.DataFrame(missing)

    summary.to_csv(
        out / "factorial_split_summary.tsv",
        sep="\t", index=False,
    )
    long_df.to_csv(
        out / "factorial_split_composition_long.tsv",
        sep="\t", index=False,
    )
    missing_df.to_csv(
        out / "missing_or_failed_cells.tsv",
        sep="\t", index=False,
    )

    if not summary.empty:
        matrix = summary.pivot_table(
            index=["task_id", "redundancy_threshold"],
            columns="split_threshold",
            values="actual_test_fraction",
            aggfunc="first",
        ).reset_index()
        matrix.columns = [
            "task_id",
            "R",
            *[f"S{int(c)}_test_fraction" for c in matrix.columns[2:]],
        ]
    else:
        matrix = pd.DataFrame()

    matrix.to_csv(
        out / "factorial_split_matrix.tsv",
        sep="\t", index=False,
    )

    md = [
        "# PhageLysClass — Stage 08 factorial split freeze",
        "",
        "Full redundancy × homology-aware split design:",
        "",
        "`R90/R70/R50/R30 × S90/S70/S50/S30`",
        "",
        markdown_table(matrix) if not matrix.empty else "_No complete cells._",
        "",
        "## Validation contract",
        "",
        "- Stage-07 redundancy datasets are the only inputs.",
        "- Clustering is performed across the complete task, not by class.",
        "- No representation or model information is used.",
        "- Train/test exact-entity overlap must be zero.",
        "- Shared MMseqs2 homology clusters across train/test must be zero.",
        "- Challenge sets remain untouched.",
    ]
    (out / "factorial_split_matrix.md").write_text(
        "\n".join(md) + "\n",
        encoding="utf-8",
    )

    expected = len(TASKS) * len(R_THRESHOLDS) * len(S_THRESHOLDS)
    passed = int(summary["status"].eq("PASS").sum()) if not summary.empty else 0

    freeze = {
        "stage": "08-factorial-audit",
        "expected_cells": expected,
        "observed_cells": int(len(summary)),
        "passed_cells": passed,
        "missing_or_failed_cells": int(len(missing_df)),
        "status": "PASS" if passed == expected and len(missing_df) == 0 else "INCOMPLETE",
        "input_output_sha256": hashes,
    }
    (out / "freeze_manifest.json").write_text(
        json.dumps(freeze, indent=2) + "\n",
        encoding="utf-8",
    )

    print(f"PASS cells: {passed}/{expected}")
    print(f"Missing/failed: {len(missing_df)}")
    print(f"Report: {out}")

    if args.strict and freeze["status"] != "PASS":
        raise SystemExit(2)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
10e_freeze_common_predictive_universe.py

Stage 10E — Common predictive-eligibility freeze.

Purpose
-------
Project the frozen Stage-08 factorial splits onto one representation-independent
predictive eligibility universe:

    canonical 20-aa sequence AND length <= 1024 aa

No train/test reassignment is allowed. This stage only removes ineligible
entities while preserving the frozen Stage-08 split and S-level cluster IDs.

The resulting universe is the canonical biological benchmark for ALL predictive
representations. Missing features/embeddings later count as representation
incompleteness; they must never redefine this universe.

Expected inputs
---------------
data/splits/<task_path>/Rxx_Syy/dataset_with_split.tsv.gz
data/representations/sequence_universe/unique_sequences.tsv.gz

Outputs
-------
data/predictive_universe/<task_path>/Rxx_Syy/
    eligible_entities.tsv.gz
    excluded_entities.tsv.gz
    eligibility_summary.tsv
    manifest.json

Global:
data/predictive_universe/
    stage10e_factorial_summary.tsv
    stage10e_exclusion_summary.tsv
    stage10e_freeze_manifest.json
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
from pathlib import Path

import pandas as pd

R_VALUES = [90, 70, 50, 30]
S_VALUES = [90, 70, 50, 30]
CANONICAL_AA = set("ACDEFGHIKLMNPQRSTVWY")
MAX_LENGTH = 1024

TASKS = [
    ("global_easy", "global_easy"),
    ("global_hard", "global_hard"),
    ("global_mixed", "global_mixed"),
    ("subtype_endolysin", "subtype_binary/endolysin"),
    ("subtype_depolymerase", "subtype_binary/depolymerase"),
    ("subtype_val", "subtype_binary/virion_associated_lytic_enzyme"),
    ("multiclass_clean", "multiclass_clean"),
    ("pu_conservative", "pu_global/conservative"),
    ("pu_broad", "pu_global/broad"),
]


def project_root() -> Path:
    # expected location: <repo>/src/representations/10e_...
    return Path(__file__).resolve().parents[2]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def sequence_sha256(seq: str) -> str:
    return hashlib.sha256(seq.encode("utf-8")).hexdigest()


def normalize_sequence(seq: str) -> str:
    return "".join(str(seq).split()).upper()


def detect_sequence_column(df: pd.DataFrame) -> str:
    for c in ("sequence", "sequence_norm", "seq"):
        if c in df.columns:
            return c
    raise KeyError("No sequence column found.")


def detect_uid_column(df: pd.DataFrame) -> str | None:
    for c in ("sequence_uid", "sequence_sha256"):
        if c in df.columns:
            return c
    return None


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(
        path, sep="\t", compression="infer", dtype=str, low_memory=False
    ).fillna("")


def parse_args():
    root = project_root()
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--splits-root", default=str(root / "data/splits"))
    p.add_argument(
        "--sequence-universe",
        default=str(
            root
            / "data/representations/sequence_universe/unique_sequences.tsv.gz"
        ),
    )
    p.add_argument(
        "--output-root",
        default=str(root / "data/predictive_universe"),
    )
    p.add_argument(
        "--only-task",
        action="append",
        choices=[x[0] for x in TASKS],
    )
    p.add_argument("--only-r", action="append", type=int, choices=R_VALUES)
    p.add_argument("--only-s", action="append", type=int, choices=S_VALUES)
    p.add_argument("--overwrite", action="store_true")
    return p.parse_args()


def load_stage10a_universe(path: Path) -> set[str]:
    if not path.is_file():
        raise FileNotFoundError(path)

    df = read_tsv(path)
    seq_col = detect_sequence_column(df)
    uid_col = detect_uid_column(df)

    seqs = df[seq_col].map(normalize_sequence)

    if uid_col is None:
        uids = seqs.map(sequence_sha256)
    else:
        uids = df[uid_col].astype(str)

        # Hard audit that Stage10A IDs actually match the exact sequence.
        expected = seqs.map(sequence_sha256)
        bad = (uids != expected)
        if bad.any():
            raise RuntimeError(
                f"Stage10A UID/sequence mismatch in {int(bad.sum())} rows."
            )

    if uids.duplicated().any():
        raise RuntimeError(
            f"Stage10A contains {int(uids.duplicated().sum())} duplicate UIDs."
        )

    return set(uids)


def choose_role_column(df: pd.DataFrame) -> str:
    for c in ("class_name", "pu_role", "label", "training_role"):
        if c in df.columns:
            return c
    return ""


def partition_counts(df: pd.DataFrame, role_col: str) -> dict[str, str | int | float]:
    out = {
        "n_total": int(len(df)),
        "n_train": int(df["split"].eq("train").sum()),
        "n_test": int(df["split"].eq("test").sum()),
    }

    if role_col:
        for partition in ("train", "test"):
            x = df.loc[df["split"].eq(partition), role_col].astype(str)
            vc = x.value_counts(dropna=False)
            out[f"{partition}_composition"] = "; ".join(
                f"{k if k else 'NA'}={int(v)}" for k, v in vc.items()
            )
    else:
        out["train_composition"] = ""
        out["test_composition"] = ""

    return out


def cell_complete(out_dir: Path, input_hash: str) -> bool:
    manifest = out_dir / "manifest.json"
    req = [
        out_dir / "eligible_entities.tsv.gz",
        out_dir / "excluded_entities.tsv.gz",
        out_dir / "eligibility_summary.tsv",
        manifest,
    ]
    if not all(x.is_file() for x in req):
        return False
    try:
        m = json.loads(manifest.read_text())
        return m.get("status") == "PASS" and m.get("input_sha256") == input_hash
    except Exception:
        return False


def process_cell(
    *,
    task_id: str,
    task_path: str,
    r: int,
    s: int,
    splits_root: Path,
    out_root: Path,
    stage10a_uids: set[str],
    overwrite: bool,
):
    cell_id = f"R{r}_S{s}"
    input_path = splits_root / task_path / cell_id / "dataset_with_split.tsv.gz"
    if not input_path.is_file():
        raise FileNotFoundError(input_path)

    input_hash = sha256_file(input_path)
    out_dir = out_root / task_path / cell_id
    out_dir.mkdir(parents=True, exist_ok=True)

    if not overwrite and cell_complete(out_dir, input_hash):
        row = pd.read_csv(out_dir / "eligibility_summary.tsv", sep="\t").iloc[0]
        print(f"SKIP {task_id} {cell_id} — already PASS")
        return row.to_dict()

    df = read_tsv(input_path)

    required = {"entity_id", "split", "cluster_id"}
    missing = required - set(df.columns)
    if missing:
        raise KeyError(f"{input_path}: missing required columns {sorted(missing)}")

    seq_col = detect_sequence_column(df)
    df["_sequence_norm"] = df[seq_col].map(normalize_sequence)
    df["sequence_length"] = df["_sequence_norm"].str.len().astype(int)
    df["sequence_uid"] = df["_sequence_norm"].map(sequence_sha256)

    # Validate any upstream SHA field instead of assuming equivalence.
    upstream_uid_col = None
    for c in ("sequence_sha256",):
        if c in df.columns:
            upstream_uid_col = c
            bad = df[c].astype(str) != df["sequence_uid"]
            if bad.any():
                raise RuntimeError(
                    f"{task_id} {cell_id}: {int(bad.sum())} upstream "
                    f"{c}/sequence mismatches."
                )

    if df["entity_id"].duplicated().any():
        raise RuntimeError(
            f"{task_id} {cell_id}: duplicate entity_id values detected."
        )
    if df["sequence_uid"].duplicated().any():
        raise RuntimeError(
            f"{task_id} {cell_id}: duplicate exact sequence UIDs detected."
        )
    if set(df["split"]) != {"train", "test"}:
        raise RuntimeError(
            f"{task_id} {cell_id}: expected split values train/test, "
            f"observed {sorted(df['split'].unique())}"
        )

    df["is_stage10a_member"] = df["sequence_uid"].isin(stage10a_uids)
    df["is_canonical20"] = df["_sequence_norm"].map(
        lambda x: bool(x) and set(x).issubset(CANONICAL_AA)
    )
    df["is_length_eligible"] = df["sequence_length"].le(MAX_LENGTH)
    df["is_main_predictive_eligible"] = (
        df["is_stage10a_member"]
        & df["is_canonical20"]
        & df["is_length_eligible"]
    )

    def exclusion_reason(row):
        reasons = []
        if not row["is_stage10a_member"]:
            reasons.append("not_in_stage10a")
        if not row["is_canonical20"]:
            reasons.append("noncanonical")
        if not row["is_length_eligible"]:
            reasons.append("length_gt_1024")
        return ";".join(reasons)

    df["exclusion_reason"] = df.apply(exclusion_reason, axis=1)

    eligible = df.loc[df["is_main_predictive_eligible"]].copy()
    excluded = df.loc[~df["is_main_predictive_eligible"]].copy()

    # Critical preservation audits.
    exact_overlap = len(
        set(eligible.loc[eligible["split"].eq("train"), "entity_id"])
        & set(eligible.loc[eligible["split"].eq("test"), "entity_id"])
    )
    cluster_overlap = len(
        set(eligible.loc[eligible["split"].eq("train"), "cluster_id"])
        & set(eligible.loc[eligible["split"].eq("test"), "cluster_id"])
    )

    if exact_overlap != 0:
        raise RuntimeError(
            f"{task_id} {cell_id}: exact train/test overlap after filtering."
        )
    if cluster_overlap != 0:
        raise RuntimeError(
            f"{task_id} {cell_id}: S-level cluster leakage after filtering."
        )

    # Verify filtering only removed entities and did not alter partition labels.
    original_map = dict(zip(df["entity_id"], df["split"]))
    if any(original_map[e] != sp for e, sp in zip(eligible["entity_id"], eligible["split"])):
        raise RuntimeError(f"{task_id} {cell_id}: split reassignment detected.")

    role_col = choose_role_column(df)
    before = partition_counts(df, role_col)
    after = partition_counts(eligible, role_col)

    summary = {
        "task_id": task_id,
        "task_path": task_path,
        "R": r,
        "S": s,
        "rs_id": cell_id,
        "n_before": before["n_total"],
        "n_after": after["n_total"],
        "n_excluded": int(len(excluded)),
        "n_excluded_noncanonical": int((~df["is_canonical20"]).sum()),
        "n_excluded_length_gt_1024": int((~df["is_length_eligible"]).sum()),
        "n_not_in_stage10a": int((~df["is_stage10a_member"]).sum()),
        "train_before": before["n_train"],
        "train_after": after["n_train"],
        "test_before": before["n_test"],
        "test_after": after["n_test"],
        "role_column": role_col,
        "train_composition_before": before["train_composition"],
        "train_composition_after": after["train_composition"],
        "test_composition_before": before["test_composition"],
        "test_composition_after": after["test_composition"],
        "n_train_clusters_after": int(
            eligible.loc[eligible["split"].eq("train"), "cluster_id"].nunique()
        ),
        "n_test_clusters_after": int(
            eligible.loc[eligible["split"].eq("test"), "cluster_id"].nunique()
        ),
        "shared_train_test_clusters": cluster_overlap,
        "train_test_entity_overlap": exact_overlap,
        "status": "PASS",
    }

    # Canonical output columns: preserve upstream biological metadata plus Stage10E fields.
    drop_private = ["_sequence_norm"]
    eligible.drop(columns=drop_private).to_csv(
        out_dir / "eligible_entities.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    excluded.drop(columns=drop_private).to_csv(
        out_dir / "excluded_entities.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    pd.DataFrame([summary]).to_csv(
        out_dir / "eligibility_summary.tsv", sep="\t", index=False
    )

    manifest = {
        "stage": "10E",
        "status": "PASS",
        "task_id": task_id,
        "task_path": task_path,
        "R": r,
        "S": s,
        "eligibility_contract": {
            "canonical20_only": True,
            "max_length": MAX_LENGTH,
            "must_be_stage10a_member": True,
            "representation_independent": True,
            "representation_missingness_does_not_change_eligibility": True,
        },
        "input": str(input_path),
        "input_sha256": input_hash,
        "stage10a_universe_path": str(args.sequence_universe) if "args" in globals() else "",
        "audit": summary,
        "outputs": {
            "eligible": str(out_dir / "eligible_entities.tsv.gz"),
            "excluded": str(out_dir / "excluded_entities.tsv.gz"),
            "summary": str(out_dir / "eligibility_summary.tsv"),
        },
        "software": {
            "python": platform.python_version(),
            "pandas": pd.__version__,
        },
    }
    (out_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n"
    )

    print(
        f"PASS {task_id} {cell_id} | "
        f"{summary['n_before']:,} -> {summary['n_after']:,} eligible | "
        f"excluded={summary['n_excluded']:,} | "
        f"train={summary['train_after']:,} test={summary['test_after']:,}"
    )
    return summary


def main():
    global args
    args = parse_args()

    splits_root = Path(args.splits_root).expanduser().resolve()
    sequence_universe = Path(args.sequence_universe).expanduser().resolve()
    out_root = Path(args.output_root).expanduser().resolve()
    out_root.mkdir(parents=True, exist_ok=True)

    selected_tasks = [
        x for x in TASKS if not args.only_task or x[0] in args.only_task
    ]
    r_values = args.only_r or R_VALUES
    s_values = args.only_s or S_VALUES

    stage10a_hash = sha256_file(sequence_universe)
    stage10a_uids = load_stage10a_universe(sequence_universe)

    print("=" * 100)
    print("PhageLysClass — Stage 10E common predictive-eligibility freeze")
    print("=" * 100)
    print(f"Stage10A universe : {sequence_universe}")
    print(f"Stage10A UIDs     : {len(stage10a_uids):,}")
    print(f"Splits root       : {splits_root}")
    print(f"Output root       : {out_root}")
    print(f"Eligibility       : canonical20 AND length <= {MAX_LENGTH}")
    print(f"Tasks             : {len(selected_tasks)}")
    print(f"R                 : {r_values}")
    print(f"S                 : {s_values}")
    print()

    rows = []
    for task_id, task_path in selected_tasks:
        for r in r_values:
            for s in s_values:
                rows.append(
                    process_cell(
                        task_id=task_id,
                        task_path=task_path,
                        r=r,
                        s=s,
                        splits_root=splits_root,
                        out_root=out_root,
                        stage10a_uids=stage10a_uids,
                        overwrite=args.overwrite,
                    )
                )

    summary = pd.DataFrame(rows)
    summary.to_csv(
        out_root / "stage10e_factorial_summary.tsv",
        sep="\t",
        index=False,
    )

    exclusion_cols = [
        "task_id", "R", "S", "rs_id", "n_before", "n_after",
        "n_excluded", "n_excluded_noncanonical",
        "n_excluded_length_gt_1024", "n_not_in_stage10a",
    ]
    summary[exclusion_cols].to_csv(
        out_root / "stage10e_exclusion_summary.tsv",
        sep="\t",
        index=False,
    )

    expected = len(selected_tasks) * len(r_values) * len(s_values)
    passed = int(summary["status"].eq("PASS").sum()) if len(summary) else 0

    freeze_manifest = {
        "stage": "10E",
        "status": "PASS" if passed == expected else "INCOMPLETE",
        "expected_cells": expected,
        "observed_cells": int(len(summary)),
        "passed_cells": passed,
        "eligibility": {
            "canonical20": True,
            "max_length": MAX_LENGTH,
            "representation_independent": True,
        },
        "stage10a_universe": str(sequence_universe),
        "stage10a_sha256": stage10a_hash,
        "stage10a_uid_count": len(stage10a_uids),
        "aggregate": {
            "total_rows_before_across_cells": int(summary["n_before"].sum()),
            "total_rows_after_across_cells": int(summary["n_after"].sum()),
            "total_excluded_across_cells": int(summary["n_excluded"].sum()),
        },
    }
    (out_root / "stage10e_freeze_manifest.json").write_text(
        json.dumps(freeze_manifest, indent=2) + "\n"
    )

    print()
    print("=" * 100)
    print(f"Stage 10E cells PASS: {passed}/{expected}")
    print("=" * 100)

    if passed != expected:
        raise SystemExit(2)

    print("PASS — common predictive eligibility is frozen for requested cells.")


if __name__ == "__main__":
    main()

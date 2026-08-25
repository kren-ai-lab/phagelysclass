#!/usr/bin/env python3
"""
12a_aggregate_and_qc_supervised.py

Aggregate and audit all completed supervised PhageLysClass experiments.

The script reads every representation marked READY in:

    config/representation_registry/representation_registry.tsv

For each READY representation it expects:

    7 supervised tasks × 16 R/S cells × 5 algorithms

Outputs
-------
results/analysis_supervised/stage12a/
    supervised_master_results.tsv
    validation_configuration_results.tsv
    per_class_results.tsv
    expected_experiment_registry.tsv
    missing_experiments.tsv
    failed_experiments.tsv
    selection_audit.tsv
    test_prediction_signatures.tsv
    test_set_signature_audit.tsv
    qc_summary.tsv
    stage12a_manifest.json
"""

from __future__ import annotations

import hashlib
import json
from itertools import product
from pathlib import Path

import pandas as pd


TASKS = {
    "global_easy": {
        "path": "global_easy",
        "type": "binary",
    },
    "global_hard": {
        "path": "global_hard",
        "type": "binary",
    },
    "global_mixed": {
        "path": "global_mixed",
        "type": "binary",
    },
    "subtype_endolysin": {
        "path": "subtype_binary/endolysin",
        "type": "binary",
    },
    "subtype_depolymerase": {
        "path": "subtype_binary/depolymerase",
        "type": "binary",
    },
    "subtype_val": {
        "path": "subtype_binary/virion_associated_lytic_enzyme",
        "type": "binary",
    },
    "multiclass_clean": {
        "path": "multiclass_clean",
        "type": "multiclass",
    },
}

R_VALUES = [90, 70, 50, 30]
S_VALUES = [90, 70, 50, 30]
ALGORITHMS = ["LR", "RF", "KNN", "XGB", "MLP"]


def root_from_script() -> Path:
    return Path(__file__).resolve().parents[2]


def read_json(path: Path) -> dict:
    return json.loads(path.read_text())


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", low_memory=False)


def load_ready_representations(root: Path) -> list[str]:
    registry = (
        root
        / "config/representation_registry/representation_registry.tsv"
    )
    if not registry.is_file():
        raise FileNotFoundError(registry)

    df = pd.read_csv(registry, sep="\t", dtype=str).fillna("")

    required = {"representation_id", "status"}
    missing = required - set(df.columns)
    if missing:
        raise KeyError(
            f"{registry}: missing required columns {sorted(missing)}"
        )

    ready = (
        df.loc[
            df["status"].str.upper().eq("READY"),
            "representation_id",
        ]
        .astype(str)
        .tolist()
    )

    if not ready:
        raise RuntimeError("No READY representations in registry.")

    return ready


def experiment_id(
    task: str,
    r: int,
    s: int,
    representation: str,
    algorithm: str,
) -> str:
    return (
        f"SUP__{task}__R{r}__S{s}__"
        f"{representation}__{algorithm}"
    )


def experiment_dir(
    root: Path,
    task: str,
    r: int,
    s: int,
    representation: str,
    algorithm: str,
) -> Path:
    return (
        root
        / "results/training_supervised"
        / TASKS[task]["path"]
        / f"R{r}_S{s}"
        / representation
        / algorithm
    )


def prediction_signature(path: Path) -> tuple[str, int]:
    """
    Hash the exact outer-test entity/label composition.

    All models inside a fixed task × R × S cell should evaluate the same
    outer-test entities and labels.
    """
    if not path.is_file():
        return "", 0

    df = pd.read_parquet(
        path,
        columns=["entity_id", "sequence_uid", "y_true"],
    )

    if df["entity_id"].duplicated().any():
        raise RuntimeError(
            f"{path}: duplicate entity_id in outer-test predictions."
        )

    df = df.sort_values(
        ["entity_id", "sequence_uid"],
        kind="stable",
    )

    h = hashlib.sha256()
    for row in df.itertuples(index=False):
        payload = (
            f"{row.entity_id}\t"
            f"{row.sequence_uid}\t"
            f"{row.y_true}\n"
        )
        h.update(payload.encode("utf-8"))

    return h.hexdigest(), len(df)


def best_config_from_validation(
    df: pd.DataFrame,
    task_type: str,
) -> str:
    """
    Recompute the frozen configuration-selection rule.

    Binary:
        MCC -> AP -> config order

    Multiclass:
        MCC -> macro-F1 -> config order
    """
    if task_type == "binary":
        required = {"mcc", "ap", "config_order", "config_id"}
        missing = required - set(df.columns)
        if missing:
            raise KeyError(
                f"Missing binary selection columns: {sorted(missing)}"
            )

        ranked = df.sort_values(
            ["mcc", "ap", "config_order"],
            ascending=[False, False, True],
            kind="stable",
        )
    else:
        required = {
            "mcc",
            "macro_f1",
            "config_order",
            "config_id",
        }
        missing = required - set(df.columns)
        if missing:
            raise KeyError(
                f"Missing multiclass selection columns: "
                f"{sorted(missing)}"
            )

        ranked = df.sort_values(
            ["mcc", "macro_f1", "config_order"],
            ascending=[False, False, True],
            kind="stable",
        )

    return str(ranked.iloc[0]["config_id"])


def main():
    root = root_from_script()

    out = root / "results/analysis_supervised/stage12a"
    out.mkdir(parents=True, exist_ok=True)

    ready_representations = load_ready_representations(root)

    expected_rows = []
    master_rows = []
    configuration_rows = []
    per_class_rows = []
    missing_rows = []
    failed_rows = []
    selection_rows = []
    signature_rows = []

    for task, r, s, representation, algorithm in product(
        TASKS,
        R_VALUES,
        S_VALUES,
        ready_representations,
        ALGORITHMS,
    ):
        exp_id = experiment_id(
            task,
            r,
            s,
            representation,
            algorithm,
        )

        exp_dir = experiment_dir(
            root,
            task,
            r,
            s,
            representation,
            algorithm,
        )

        expected_rows.append(
            {
                "experiment_id": exp_id,
                "task": task,
                "task_type": TASKS[task]["type"],
                "R": r,
                "S": s,
                "representation": representation,
                "algorithm": algorithm,
                "expected_path": str(exp_dir),
            }
        )

        manifest_path = exp_dir / "manifest.json"

        if not manifest_path.is_file():
            missing_rows.append(
                {
                    "experiment_id": exp_id,
                    "reason": "missing_manifest",
                    "expected_path": str(exp_dir),
                }
            )
            continue

        try:
            manifest = read_json(manifest_path)
        except Exception as exc:
            failed_rows.append(
                {
                    "experiment_id": exp_id,
                    "status": "INVALID_MANIFEST",
                    "error": str(exc),
                    "path": str(exp_dir),
                }
            )
            continue

        status = str(manifest.get("status", "UNKNOWN"))

        if status != "DONE":
            failed_rows.append(
                {
                    "experiment_id": exp_id,
                    "status": status,
                    "error": manifest.get("error", ""),
                    "path": str(exp_dir),
                }
            )
            continue

        required_files = {
            "validation_metrics": (
                exp_dir / "validation_metrics.tsv"
            ),
            "selected_configuration": (
                exp_dir / "selected_configuration.json"
            ),
            "outer_test_metrics": (
                exp_dir / "outer_test_metrics.tsv"
            ),
            "outer_test_predictions": (
                exp_dir / "outer_test_predictions.parquet"
            ),
            "runtime": (
                exp_dir / "runtime.json"
            ),
        }

        absent = [
            name
            for name, path in required_files.items()
            if not path.is_file()
        ]

        if absent:
            failed_rows.append(
                {
                    "experiment_id": exp_id,
                    "status": "INCOMPLETE_DONE",
                    "error": (
                        "missing files: " + ",".join(absent)
                    ),
                    "path": str(exp_dir),
                }
            )
            continue

        validation = read_tsv(
            required_files["validation_metrics"]
        )

        selected_json = read_json(
            required_files["selected_configuration"]
        )

        test_metrics = read_tsv(
            required_files["outer_test_metrics"]
        )

        runtime = read_json(
            required_files["runtime"]
        )

        if len(test_metrics) != 1:
            failed_rows.append(
                {
                    "experiment_id": exp_id,
                    "status": "INVALID_TEST_METRICS",
                    "error": (
                        f"outer_test_metrics rows="
                        f"{len(test_metrics)}"
                    ),
                    "path": str(exp_dir),
                }
            )
            continue

        selected_id = str(
            selected_json.get("config_id", "")
        )

        recomputed_selected_id = best_config_from_validation(
            validation,
            TASKS[task]["type"],
        )

        selection_match = (
            selected_id == recomputed_selected_id
        )

        selection_rows.append(
            {
                "experiment_id": exp_id,
                "task": task,
                "R": r,
                "S": s,
                "representation": representation,
                "algorithm": algorithm,
                "selected_config": selected_id,
                "recomputed_selected_config": (
                    recomputed_selected_id
                ),
                "selection_match": selection_match,
                "n_validation_configs": len(validation),
            }
        )

        for row in validation.to_dict("records"):
            configuration_rows.append(
                {
                    "task": task,
                    "task_type": TASKS[task]["type"],
                    "R": r,
                    "S": s,
                    "representation": representation,
                    "algorithm": algorithm,
                    **row,
                }
            )

        test_row = test_metrics.iloc[0].to_dict()

        master = {
            "experiment_id": exp_id,
            "task": task,
            "task_type": TASKS[task]["type"],
            "R": r,
            "S": s,
            "representation": representation,
            "algorithm": algorithm,
            "selected_config": selected_id,
            **test_row,
        }

        for key, value in runtime.items():
            if (
                isinstance(
                    value,
                    (str, int, float, bool),
                )
                or value is None
            ):
                master[f"runtime__{key}"] = value

        master_rows.append(master)

        per_class_path = (
            exp_dir / "per_class_metrics.tsv"
        )

        if per_class_path.is_file():
            per_class = read_tsv(per_class_path)

            for row in per_class.to_dict("records"):
                per_class_rows.append(
                    {
                        "experiment_id": exp_id,
                        "task": task,
                        "R": r,
                        "S": s,
                        "representation": representation,
                        "algorithm": algorithm,
                        **row,
                    }
                )

        signature, n_predictions = (
            prediction_signature(
                required_files["outer_test_predictions"]
            )
        )

        signature_rows.append(
            {
                "experiment_id": exp_id,
                "task": task,
                "R": r,
                "S": s,
                "representation": representation,
                "algorithm": algorithm,
                "test_signature": signature,
                "n_test_predictions": n_predictions,
            }
        )

    expected = pd.DataFrame(expected_rows)
    master = pd.DataFrame(master_rows)
    configurations = pd.DataFrame(configuration_rows)
    per_class = pd.DataFrame(per_class_rows)
    missing = pd.DataFrame(missing_rows)
    failed = pd.DataFrame(failed_rows)
    selection = pd.DataFrame(selection_rows)
    signatures = pd.DataFrame(signature_rows)

    expected.to_csv(
        out / "expected_experiment_registry.tsv",
        sep="\t",
        index=False,
    )

    master.to_csv(
        out / "supervised_master_results.tsv",
        sep="\t",
        index=False,
    )

    configurations.to_csv(
        out / "validation_configuration_results.tsv",
        sep="\t",
        index=False,
    )

    per_class.to_csv(
        out / "per_class_results.tsv",
        sep="\t",
        index=False,
    )

    missing.to_csv(
        out / "missing_experiments.tsv",
        sep="\t",
        index=False,
    )

    failed.to_csv(
        out / "failed_experiments.tsv",
        sep="\t",
        index=False,
    )

    selection.to_csv(
        out / "selection_audit.tsv",
        sep="\t",
        index=False,
    )

    signatures.to_csv(
        out / "test_prediction_signatures.tsv",
        sep="\t",
        index=False,
    )

    if len(signatures):
        signature_audit = (
            signatures.groupby(
                ["task", "R", "S"],
                as_index=False,
            )
            .agg(
                n_experiments=(
                    "experiment_id",
                    "size",
                ),
                unique_test_signatures=(
                    "test_signature",
                    "nunique",
                ),
                unique_test_sizes=(
                    "n_test_predictions",
                    "nunique",
                ),
                n_test=(
                    "n_test_predictions",
                    "first",
                ),
            )
        )

        signature_audit[
            "same_test_set_across_models"
        ] = (
            signature_audit[
                "unique_test_signatures"
            ].eq(1)
            & signature_audit[
                "unique_test_sizes"
            ].eq(1)
        )
    else:
        signature_audit = pd.DataFrame()

    signature_audit.to_csv(
        out / "test_set_signature_audit.tsv",
        sep="\t",
        index=False,
    )

    expected_n = len(expected)
    done_n = len(master)
    missing_n = len(missing)
    failed_n = len(failed)

    selection_failures = (
        int((~selection["selection_match"]).sum())
        if len(selection)
        else 0
    )

    signature_failures = (
        int(
            (
                ~signature_audit[
                    "same_test_set_across_models"
                ]
            ).sum()
        )
        if len(signature_audit)
        else 0
    )

    qc = pd.DataFrame(
        [
            {
                "check": "ready_representations",
                "value": ",".join(
                    ready_representations
                ),
            },
            {
                "check": "expected_experiment_families",
                "value": expected_n,
            },
            {
                "check": "done_experiment_families",
                "value": done_n,
            },
            {
                "check": "missing_experiment_families",
                "value": missing_n,
            },
            {
                "check": "failed_or_incomplete_families",
                "value": failed_n,
            },
            {
                "check": "selection_audit_failures",
                "value": selection_failures,
            },
            {
                "check": "test_set_signature_failures",
                "value": signature_failures,
            },
        ]
    )

    qc.to_csv(
        out / "qc_summary.tsv",
        sep="\t",
        index=False,
    )

    status = (
        "PASS"
        if (
            done_n == expected_n
            and missing_n == 0
            and failed_n == 0
            and selection_failures == 0
            and signature_failures == 0
        )
        else "INCOMPLETE_OR_FAIL"
    )

    manifest = {
        "stage": "12A",
        "status": status,
        "ready_representations": (
            ready_representations
        ),
        "expected_experiment_families": expected_n,
        "done_experiment_families": done_n,
        "missing_experiment_families": missing_n,
        "failed_or_incomplete_families": failed_n,
        "selection_audit_failures": (
            selection_failures
        ),
        "test_set_signature_failures": (
            signature_failures
        ),
    }

    (
        out / "stage12a_manifest.json"
    ).write_text(
        json.dumps(manifest, indent=2) + "\n"
    )

    print("=" * 100)
    print(
        "PhageLysClass — Stage 12A "
        "supervised aggregation + QC"
    )
    print("=" * 100)
    print(
        "READY representations : "
        + ", ".join(ready_representations)
    )
    print(
        f"Expected families      : {expected_n:,}"
    )
    print(
        f"DONE families          : {done_n:,}"
    )
    print(
        f"Missing families       : {missing_n:,}"
    )
    print(
        f"Failed/incomplete      : {failed_n:,}"
    )
    print(
        f"Selection failures     : "
        f"{selection_failures:,}"
    )
    print(
        f"Test-set mismatches    : "
        f"{signature_failures:,}"
    )
    print(f"STATUS                 : {status}")
    print(f"Outputs                : {out}")
    print("=" * 100)


if __name__ == "__main__":
    main()

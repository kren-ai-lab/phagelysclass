#!/usr/bin/env python3
"""
12a2_audit_pu_training.py

Robust completeness/QC audit for the PhageLysClass PU benchmark.

Unlike the previous version, this auditor DOES NOT assume a specific
task-directory layout under results/training_pu/. It discovers manifests
recursively and validates experiments by their canonical experiment_id.

Expected design
---------------
2 PU tasks × 16 R/S × 4 representations × 3 PU methods × 2 algorithms
= 768 experiment families

Per representation:
2 × 16 × 3 × 2 = 192 families
"""

from __future__ import annotations

import hashlib
import json
import re
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd

TASKS = ["pu_broad", "pu_conservative"]
R_VALUES = [90, 70, 50, 30]
S_VALUES = [90, 70, 50, 30]
REPRESENTATIONS = [
    "roxy_48",
    "onehot_1024",
    "esm2_t6_8m",
    "esm2_t30_150m",
]
METHODS = ["naive", "baggingpu", "elkan_noto"]
ALGORITHMS = ["LR", "XGB"]

EXP_RE = re.compile(
    r"^PU__(pu_(?:broad|conservative))__R(\d+)__S(\d+)__"
    r"(.+)__(naive|baggingpu|elkan_noto)__(LR|XGB)$"
)


def project_root() -> Path:
    return Path(__file__).resolve().parents[2]


def canonical_id(task, r, s, rep, method, alg):
    return f"PU__{task}__R{r}__S{s}__{rep}__{method}__{alg}"


def score_signature(path: Path):
    df = pd.read_parquet(path)
    required = {"entity_id", "sequence_uid", "observed_role"}
    missing = required - set(df.columns)
    if missing:
        raise KeyError(f"{path}: missing columns {sorted(missing)}")

    if df["entity_id"].duplicated().any():
        raise RuntimeError(f"{path}: duplicate entity_id")

    df = df.sort_values(["entity_id", "sequence_uid"], kind="stable")

    h = hashlib.sha256()
    for row in df.itertuples(index=False):
        h.update(
            f"{row.entity_id}\t{row.sequence_uid}\t{row.observed_role}\n".encode()
        )
    return h.hexdigest(), len(df)


def discover_manifests(training_root: Path):
    """
    Discover every PU manifest recursively and index by canonical experiment_id.
    """
    found = {}
    invalid = []

    for manifest_path in training_root.rglob("manifest.json"):
        try:
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        except Exception as exc:
            invalid.append(
                {
                    "manifest_path": str(manifest_path),
                    "reason": f"invalid_json: {exc}",
                }
            )
            continue

        exp_id = str(manifest.get("experiment_id", "")).strip()

        if not exp_id:
            invalid.append(
                {
                    "manifest_path": str(manifest_path),
                    "reason": "missing experiment_id",
                }
            )
            continue

        if not EXP_RE.match(exp_id):
            # Ignore unrelated manifests under the tree.
            continue

        found.setdefault(exp_id, []).append((manifest_path, manifest))

    return found, invalid


def main():
    root = project_root()
    training_root = root / "results/training_pu"
    out = root / "results/analysis_pu/stage12a_qc"
    out.mkdir(parents=True, exist_ok=True)

    if not training_root.is_dir():
        raise FileNotFoundError(training_root)

    discovered, invalid_manifests = discover_manifests(training_root)

    expected_ids = [
        canonical_id(task, r, s, rep, method, alg)
        for task, r, s, rep, method, alg in product(
            TASKS,
            R_VALUES,
            S_VALUES,
            REPRESENTATIONS,
            METHODS,
            ALGORITHMS,
        )
    ]

    completed = []
    problems = []
    method_checks = []
    signatures = []

    for exp_id in expected_ids:
        hits = discovered.get(exp_id, [])

        if len(hits) == 0:
            problems.append(
                {
                    "experiment_id": exp_id,
                    "status": "MISSING_MANIFEST",
                    "reason": "canonical experiment_id not found anywhere under results/training_pu",
                    "path": "",
                }
            )
            continue

        if len(hits) > 1:
            problems.append(
                {
                    "experiment_id": exp_id,
                    "status": "DUPLICATE_MANIFEST",
                    "reason": f"{len(hits)} manifests found",
                    "path": ";".join(str(x[0]) for x in hits),
                }
            )
            continue

        manifest_path, manifest = hits[0]
        d = manifest_path.parent

        m = EXP_RE.match(exp_id)
        task, r, s, rep, method, alg = m.groups()
        r, s = int(r), int(s)

        if str(manifest.get("status", "")) != "DONE":
            problems.append(
                {
                    "experiment_id": exp_id,
                    "status": str(manifest.get("status", "UNKNOWN")),
                    "reason": str(manifest.get("error", "")),
                    "path": str(d),
                }
            )
            continue

        mandatory = [
            d / "pu_diagnostics.tsv",
            d / "outer_test_scores.parquet",
            d / "warnings.log",
        ]
        missing_files = [p.name for p in mandatory if not p.is_file()]

        if missing_files:
            problems.append(
                {
                    "experiment_id": exp_id,
                    "status": "MISSING_OUTPUT",
                    "reason": ",".join(missing_files),
                    "path": str(d),
                }
            )
            continue

        try:
            diagnostics = pd.read_csv(d / "pu_diagnostics.tsv", sep="\t")
            if len(diagnostics) != 1:
                raise RuntimeError(
                    f"pu_diagnostics.tsv expected 1 row, found {len(diagnostics)}"
                )

            signature, n_test = score_signature(d / "outer_test_scores.parquet")

        except Exception as exc:
            problems.append(
                {
                    "experiment_id": exp_id,
                    "status": "INVALID_OUTPUT",
                    "reason": str(exc),
                    "path": str(d),
                }
            )
            continue

        completed.append(
            {
                "experiment_id": exp_id,
                "task": task,
                "R": r,
                "S": s,
                "representation": rep,
                "method": method,
                "algorithm": alg,
                "n_test": n_test,
                "path": str(d),
            }
        )

        signatures.append(
            {
                "experiment_id": exp_id,
                "task": task,
                "R": r,
                "S": s,
                "representation": rep,
                "method": method,
                "algorithm": alg,
                "test_signature": signature,
                "n_test": n_test,
            }
        )

        ok = True
        detail = ""

        if method == "baggingpu":
            p = d / "bagging_diagnostics.tsv"
            if not p.is_file():
                ok = False
                detail = "missing bagging_diagnostics.tsv"
            else:
                n_bags = len(pd.read_csv(p, sep="\t"))
                if n_bags != 30:
                    ok = False
                    detail = f"expected 30 bags, found {n_bags}"

        elif method == "elkan_noto":
            p = d / "propensity_crossfit.tsv"

            if not p.is_file():
                ok = False
                detail = "missing propensity_crossfit.tsv"
            else:
                n_folds = len(pd.read_csv(p, sep="\t"))
                if n_folds != 5:
                    ok = False
                    detail = f"expected 5 folds, found {n_folds}"

            try:
                c = float(manifest.get("c", np.nan))
            except Exception:
                c = np.nan

            if not np.isfinite(c) or c <= 0:
                ok = False
                detail = (detail + f"; invalid c={c}").strip("; ")

        method_checks.append(
            {
                "experiment_id": exp_id,
                "task": task,
                "R": r,
                "S": s,
                "representation": rep,
                "method": method,
                "algorithm": alg,
                "method_specific_ok": ok,
                "detail": detail,
            }
        )

    completed = pd.DataFrame(completed)
    problems = pd.DataFrame(problems)
    method_checks = pd.DataFrame(method_checks)
    signatures = pd.DataFrame(signatures)
    invalid_manifests = pd.DataFrame(invalid_manifests)

    completed.to_csv(out / "pu_completed_registry.tsv", sep="\t", index=False)
    problems.to_csv(out / "pu_missing_or_failed.tsv", sep="\t", index=False)
    method_checks.to_csv(out / "pu_method_specific_audit.tsv", sep="\t", index=False)
    invalid_manifests.to_csv(out / "invalid_discovered_manifests.tsv", sep="\t", index=False)

    if len(signatures):
        signature_audit = (
            signatures.groupby(
                ["task", "R", "S", "representation"],
                as_index=False,
            )
            .agg(
                n_experiments=("experiment_id", "size"),
                unique_test_signatures=("test_signature", "nunique"),
                unique_test_sizes=("n_test", "nunique"),
                n_test=("n_test", "first"),
            )
        )

        signature_audit["same_test_set_across_methods_algorithms"] = (
            signature_audit["n_experiments"].eq(6)
            & signature_audit["unique_test_signatures"].eq(1)
            & signature_audit["unique_test_sizes"].eq(1)
        )
    else:
        signature_audit = pd.DataFrame()

    signature_audit.to_csv(
        out / "pu_test_signature_audit.tsv",
        sep="\t",
        index=False,
    )

    expected_n = len(expected_ids)
    completed_n = len(completed)
    problem_n = len(problems)

    method_failures = (
        int((~method_checks["method_specific_ok"].astype(bool)).sum())
        if len(method_checks)
        else 0
    )

    signature_failures = (
        int(
            (~signature_audit["same_test_set_across_methods_algorithms"].astype(bool)).sum()
        )
        if len(signature_audit)
        else 0
    )

    counts = (
        completed.groupby("representation").size().to_dict()
        if len(completed)
        else {}
    )

    per_rep_expected = 192

    status = (
        "PASS"
        if (
            completed_n == expected_n
            and problem_n == 0
            and method_failures == 0
            and signature_failures == 0
        )
        else "INCOMPLETE_OR_FAIL"
    )

    summary = [
        ("expected_total_families", expected_n),
        ("completed_total_families", completed_n),
        ("missing_or_failed_families", problem_n),
        ("method_specific_failures", method_failures),
        ("test_signature_failures", signature_failures),
    ]

    for rep in REPRESENTATIONS:
        summary.append((f"completed__{rep}", counts.get(rep, 0)))

    pd.DataFrame(summary, columns=["check", "value"]).to_csv(
        out / "pu_qc_summary.tsv",
        sep="\t",
        index=False,
    )

    manifest = {
        "stage": "12A2-PU-QC",
        "status": status,
        "discovery_mode": "recursive_manifest_scan_by_experiment_id",
        "expected_total_families": expected_n,
        "completed_total_families": completed_n,
        "missing_or_failed_families": problem_n,
        "method_specific_failures": method_failures,
        "test_signature_failures": signature_failures,
        "completed_by_representation": counts,
    }

    (out / "stage12a_pu_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n",
        encoding="utf-8",
    )

    print("=" * 90)
    print("PhageLysClass — PU completeness audit v2")
    print("=" * 90)
    print(f"Discovered canonical PU manifests : {len(discovered)}")
    print(f"Expected total families           : {expected_n}")
    print(f"Completed valid families          : {completed_n}")
    print(f"Missing/failed/duplicate          : {problem_n}")
    print(f"Method-specific failures          : {method_failures}")
    print(f"Test-set mismatches               : {signature_failures}")
    for rep in REPRESENTATIONS:
        print(f"{rep:<28}: {counts.get(rep, 0)} / {per_rep_expected}")
    print(f"STATUS                            : {status}")
    print(f"Outputs                           : {out}")


if __name__ == "__main__":
    main()

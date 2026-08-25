#!/usr/bin/env python3
"""
11_train_supervised_factorial.py

Production supervised trainer.

Scientific CLI parameters:
    --task
    --representation

Everything else is mandatory:
- all 16 R × S cells;
- all 5 algorithms;
- all 3 predefined configurations per algorithm;
- validation-only configuration selection;
- full outer-train refit;
- frozen outer-test prediction;
- per-sequence prediction retention.

Operational flags such as --overwrite and --smoke do not change the production
scientific design. --smoke is only for implementation testing.
"""

from __future__ import annotations

import argparse
import json
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd

from _training_common import (
    ALGORITHMS,
    PARAMETER_PANELS,
    R_VALUES,
    S_VALUES,
    SUPERVISED_TASKS,
    TASKS,
    RepresentationStore,
    binary_metrics,
    build_model,
    choose_best_config,
    fit_model,
    load_representation_registry,
    make_binary_prediction_table,
    make_multiclass_prediction_table,
    multiclass_metrics,
    multiclass_scores,
    positive_scores,
    process_peak_rss_kb,
    project_root_from_script,
    read_tsv,
    sha256_file,
    software_versions,
    supervised_target,
    write_json,
)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--task", required=True, choices=SUPERVISED_TASKS)
    p.add_argument("--representation", required=True)
    p.add_argument("--overwrite", action="store_true")
    p.add_argument(
        "--smoke",
        action="store_true",
        help="Implementation test only: R90_S90 + LR. Never use for production freeze.",
    )
    return p.parse_args()


def family_done(out_dir: Path) -> bool:
    manifest = out_dir / "manifest.json"
    if not manifest.is_file():
        return False
    try:
        return json.loads(manifest.read_text()).get("status") == "DONE"
    except Exception:
        return False


def load_cell(root: Path, task_path: str, r: int, s: int):
    rs = f"R{r}_S{s}"
    universe = (
        root / "data/predictive_universe" / task_path / rs
        / "eligible_entities.tsv.gz"
    )
    inner = (
        root / "data/training_splits" / task_path / rs
        / "inner_seed42.tsv.gz"
    )
    if not universe.is_file():
        raise FileNotFoundError(universe)
    if not inner.is_file():
        raise FileNotFoundError(
            f"{inner}\nRun Stage 10F before training."
        )

    df = read_tsv(universe)
    inner_df = read_tsv(inner)
    assignment = dict(
        zip(
            inner_df["sequence_uid"].astype(str),
            inner_df["inner_split"].astype(str),
        )
    )
    df["inner_split"] = df["sequence_uid"].map(assignment).fillna("")
    if not set(df.loc[df["split"].eq("train"), "inner_split"]) <= {
        "train",
        "validation",
    }:
        raise RuntimeError(f"{rs}: invalid/missing inner assignment")
    if (df.loc[df["split"].eq("test"), "inner_split"] != "").any():
        raise RuntimeError(f"{rs}: outer test unexpectedly has inner assignment")
    return df, universe, inner


def train_family(
    *,
    root: Path,
    task_id: str,
    task_type: str,
    task_path: str,
    r: int,
    s: int,
    rep: RepresentationStore,
    algorithm: str,
    overwrite: bool,
):
    rs = f"R{r}_S{s}"
    experiment_id = (
        f"SUP__{task_id}__R{r}__S{s}__"
        f"{rep.spec.representation_id}__{algorithm}"
    )
    out_dir = (
        root / "results/training_supervised" / task_path / rs
        / rep.spec.representation_id / algorithm
    )
    out_dir.mkdir(parents=True, exist_ok=True)

    if family_done(out_dir) and not overwrite:
        print(f"SKIP {experiment_id} — DONE")
        return {"experiment_id": experiment_id, "status": "DONE", "skipped": True}

    start = time.perf_counter()
    start_rss = process_peak_rss_kb()

    try:
        df, universe_path, inner_path = load_cell(
            root, task_path, r, s
        )
        y_all, class_names = supervised_target(df, task_type)

        # Preserve row alignment after target construction.
        df = df.reset_index(drop=True)
        y_all = np.asarray(y_all, dtype=int)

        idx_inner_train = np.where(
            (df["split"].eq("train") & df["inner_split"].eq("train")).to_numpy()
        )[0]
        idx_val = np.where(
            (
                df["split"].eq("train")
                & df["inner_split"].eq("validation")
            ).to_numpy()
        )[0]
        idx_outer_train = np.where(df["split"].eq("train").to_numpy())[0]
        idx_test = np.where(df["split"].eq("test").to_numpy())[0]

        all_uids = df["sequence_uid"].astype(str).tolist()
        rep.require_uids(all_uids)

        X_inner = rep.get(df.iloc[idx_inner_train]["sequence_uid"].astype(str).tolist())
        X_val = rep.get(df.iloc[idx_val]["sequence_uid"].astype(str).tolist())
        X_outer = rep.get(df.iloc[idx_outer_train]["sequence_uid"].astype(str).tolist())
        X_test = rep.get(df.iloc[idx_test]["sequence_uid"].astype(str).tolist())

        y_inner = y_all[idx_inner_train]
        y_val = y_all[idx_val]
        y_outer = y_all[idx_outer_train]
        y_test = y_all[idx_test]

        validation_rows = []
        validation_predictions = []

        for config_order, config in enumerate(PARAMETER_PANELS[algorithm]):
            config_id = config["config_id"]
            model = build_model(
                algorithm,
                config,
                task_type=task_type,
                representation_family=rep.spec.family,
                n_classes=len(class_names),
            )

            t0 = time.perf_counter()
            fit_model(model, X_inner, y_inner, algorithm)
            fit_seconds = time.perf_counter() - t0

            t1 = time.perf_counter()
            if task_type == "binary":
                scores = positive_scores(model, X_val)
                metrics = binary_metrics(y_val, scores)
                pred_df = make_binary_prediction_table(
                    df.iloc[idx_val],
                    y_val,
                    scores,
                    experiment_id=experiment_id,
                    config_id=config_id,
                )
            else:
                prob = multiclass_scores(model, X_val)
                metrics, _, _ = multiclass_metrics(
                    y_val, prob, class_names
                )
                pred_df = make_multiclass_prediction_table(
                    df.iloc[idx_val],
                    y_val,
                    prob,
                    class_names,
                    experiment_id=experiment_id,
                    config_id=config_id,
                )
            pred_seconds = time.perf_counter() - t1

            row = {
                "experiment_id": experiment_id,
                "config_id": config_id,
                "config_order": config_order,
                **{
                    k: json.dumps(v) if isinstance(v, tuple) else v
                    for k, v in config.items()
                    if k != "config_id"
                },
                "n_inner_train": len(idx_inner_train),
                "n_validation": len(idx_val),
                "fit_seconds": fit_seconds,
                "validation_prediction_seconds": pred_seconds,
                **metrics,
            }
            validation_rows.append(row)
            validation_predictions.append(pred_df)

        val_metrics = pd.DataFrame(validation_rows)
        selected_config_id = choose_best_config(val_metrics, task_type)
        val_metrics["selected"] = (
            val_metrics["config_id"].eq(selected_config_id)
        )

        selected_config = next(
            cfg for cfg in PARAMETER_PANELS[algorithm]
            if cfg["config_id"] == selected_config_id
        )

        final_model = build_model(
            algorithm,
            selected_config,
            task_type=task_type,
            representation_family=rep.spec.family,
            n_classes=len(class_names),
        )

        t2 = time.perf_counter()
        fit_model(final_model, X_outer, y_outer, algorithm)
        final_fit_seconds = time.perf_counter() - t2

        t3 = time.perf_counter()
        if task_type == "binary":
            test_scores = positive_scores(final_model, X_test)
            test_metrics = binary_metrics(y_test, test_scores)
            test_predictions = make_binary_prediction_table(
                df.iloc[idx_test],
                y_test,
                test_scores,
                experiment_id=experiment_id,
            )
            per_class = pd.DataFrame(
                [
                    {
                        "class_id": 0,
                        "class_name": "negative",
                        "support": int(np.sum(y_test == 0)),
                    },
                    {
                        "class_id": 1,
                        "class_name": "positive",
                        "support": int(np.sum(y_test == 1)),
                    },
                ]
            )
            cm = np.array(
                [
                    [test_metrics["tn"], test_metrics["fp"]],
                    [test_metrics["fn"], test_metrics["tp"]],
                ]
            )
        else:
            test_prob = multiclass_scores(final_model, X_test)
            test_metrics, per_class, cm = multiclass_metrics(
                y_test, test_prob, class_names
            )
            test_predictions = make_multiclass_prediction_table(
                df.iloc[idx_test],
                y_test,
                test_prob,
                class_names,
                experiment_id=experiment_id,
            )
        test_prediction_seconds = time.perf_counter() - t3

        pd.concat(validation_predictions, ignore_index=True).to_parquet(
            out_dir / "validation_predictions.parquet",
            index=False,
            compression="zstd",
        )
        val_metrics.to_csv(
            out_dir / "validation_metrics.tsv",
            sep="\t",
            index=False,
        )
        write_json(
            out_dir / "selected_configuration.json",
            selected_config,
        )
        pd.DataFrame(
            [
                {
                    "experiment_id": experiment_id,
                    "selected_config": selected_config_id,
                    **test_metrics,
                }
            ]
        ).to_csv(
            out_dir / "outer_test_metrics.tsv",
            sep="\t",
            index=False,
        )
        test_predictions.to_parquet(
            out_dir / "outer_test_predictions.parquet",
            index=False,
            compression="zstd",
        )
        per_class.to_csv(
            out_dir / "per_class_metrics.tsv",
            sep="\t",
            index=False,
        )
        pd.DataFrame(
            cm,
            index=[class_names[i] for i in sorted(class_names)],
            columns=[class_names[i] for i in sorted(class_names)],
        ).to_csv(
            out_dir / "confusion_matrix.tsv",
            sep="\t",
            index=True,
        )

        elapsed = time.perf_counter() - start
        runtime = {
            "experiment_id": experiment_id,
            "validation_and_final_wall_seconds": elapsed,
            "final_fit_seconds": final_fit_seconds,
            "final_test_prediction_seconds": test_prediction_seconds,
            "process_peak_rss_kb_start": start_rss,
            "process_peak_rss_kb_end": process_peak_rss_kb(),
            "representation_dimension": rep.dimension,
            "representation_sparse": rep.sparse,
            "n_outer_train": len(idx_outer_train),
            "n_inner_train": len(idx_inner_train),
            "n_validation": len(idx_val),
            "n_outer_test": len(idx_test),
        }
        write_json(out_dir / "runtime.json", runtime)

        manifest = {
            "stage": "11-supervised",
            "status": "DONE",
            "experiment_id": experiment_id,
            "task": task_id,
            "task_type": task_type,
            "R": r,
            "S": s,
            "representation": rep.spec.representation_id,
            "representation_family": rep.spec.family,
            "representation_dimension": rep.dimension,
            "representation_sparse": rep.sparse,
            "algorithm": algorithm,
            "selected_config": selected_config_id,
            "outer_universe_path": str(universe_path),
            "outer_universe_sha256": sha256_file(universe_path),
            "inner_split_path": str(inner_path),
            "inner_split_sha256": sha256_file(inner_path),
            "representation_path": str(rep.spec.main_path),
            "representation_sha256": sha256_file(rep.spec.main_path),
            "representation_index_path": (
                str(rep.spec.index_path) if rep.spec.index_path else None
            ),
            "software": software_versions(),
            "runtime": runtime,
        }
        write_json(out_dir / "manifest.json", manifest)
        (out_dir / "warnings.log").write_text("", encoding="utf-8")

        print(
            f"DONE {experiment_id} | selected={selected_config_id} | "
            f"primary={test_metrics['mcc']:.4f}"
        )
        return {
            "experiment_id": experiment_id,
            "status": "DONE",
            "selected_config": selected_config_id,
            "primary_metric": test_metrics["mcc"],
        }

    except Exception as exc:
        error = traceback.format_exc()
        (out_dir / "error.txt").write_text(error, encoding="utf-8")
        write_json(
            out_dir / "manifest.json",
            {
                "stage": "11-supervised",
                "status": "FAILED",
                "experiment_id": experiment_id,
                "task": task_id,
                "R": r,
                "S": s,
                "representation": rep.spec.representation_id,
                "algorithm": algorithm,
                "error": str(exc),
                "software": software_versions(),
            },
        )
        print(f"FAILED {experiment_id}: {exc}")
        return {
            "experiment_id": experiment_id,
            "status": "FAILED",
            "error": str(exc),
        }


def main():
    args = parse_args()
    root = project_root_from_script(__file__)

    task_info = TASKS[args.task]
    task_type = task_info["type"]
    task_path = task_info["path"]

    registry_path = root / "config/representation_registry/representation_registry.tsv"
    registry = load_representation_registry(registry_path, root)
    if args.representation not in registry:
        raise KeyError(
            f"Unknown representation {args.representation}. "
            f"Available: {sorted(registry)}"
        )

    print("=" * 100)
    print("PhageLysClass — Stage 11 supervised factorial training")
    print("=" * 100)
    print(f"Task           : {args.task}")
    print(f"Representation : {args.representation}")
    print("Production     : all 16 R/S × all 5 algorithms × all 3 configs")
    if args.smoke:
        print("MODE           : SMOKE ONLY (R90_S90 + LR)")
    print()

    rep = RepresentationStore(registry[args.representation])
    print(
        f"Loaded representation: rows={len(rep.uid_to_row):,} "
        f"dim={rep.dimension:,} sparse={rep.sparse}"
    )

    r_values = [90] if args.smoke else R_VALUES
    s_values = [90] if args.smoke else S_VALUES
    algorithms = ["LR"] if args.smoke else ALGORITHMS

    results = []
    for r in r_values:
        for s in s_values:
            for algorithm in algorithms:
                results.append(
                    train_family(
                        root=root,
                        task_id=args.task,
                        task_type=task_type,
                        task_path=task_path,
                        r=r,
                        s=s,
                        rep=rep,
                        algorithm=algorithm,
                        overwrite=args.overwrite,
                    )
                )

    summary = pd.DataFrame(results)
    out = (
        root / "results/training_supervised" / task_path
        / f"_run_summary__{args.representation}.tsv"
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(out, sep="\t", index=False)

    counts = summary["status"].value_counts().to_dict()
    print()
    print("=" * 100)
    print(f"Run summary: {counts}")
    print("=" * 100)

    if not args.smoke:
        expected = 16 * 5
        if len(summary) != expected:
            raise SystemExit(
                f"Expected {expected} experiment families, got {len(summary)}"
            )


if __name__ == "__main__":
    main()

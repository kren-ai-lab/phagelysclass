#!/usr/bin/env python3
"""
PhageLysClass Stage 11A — supervised predictive benchmark.

One invocation = one task × R × S × representation.
All requested algorithms are trained using the frozen 3-configuration panel.

Example:
    python src/training_models/11a_train_supervised.py \
        --task global_easy \
        --task-type binary \
        --redundancy R90 \
        --split-threshold S90 \
        --split-file data/.../R90_S90/dataset.tsv.gz \
        --representation esm2_t6_8m \
        --algorithm all
"""

from __future__ import annotations

import argparse
import inspect
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import LabelEncoder
from sklearn.utils.class_weight import compute_sample_weight

from _training_utils import (
    REPRESENTATIONS,
    binary_metrics,
    common_eligible_ids,
    create_or_load_inner_split,
    encode_binary_label,
    fit_transform_scaler,
    load_outer_split,
    load_representation,
    multiclass_metrics,
    sha256_file,
    should_standardize,
    transform_with_scaler,
    write_json,
)

ALGORITHMS = ["LR", "RF", "KNN", "XGBoost", "MLP"]

CONFIGS = {
    "LR": [
        ("LR01", {"C": 0.1}),
        ("LR02", {"C": 1.0}),
        ("LR03", {"C": 10.0}),
    ],
    "RF": [
        ("RF01", {"n_estimators": 300, "max_depth": None, "min_samples_leaf": 1}),
        ("RF02", {"n_estimators": 500, "max_depth": None, "min_samples_leaf": 2}),
        ("RF03", {"n_estimators": 500, "max_depth": 30, "min_samples_leaf": 1}),
    ],
    "KNN": [
        ("KNN01", {"n_neighbors": 3}),
        ("KNN02", {"n_neighbors": 5}),
        ("KNN03", {"n_neighbors": 11}),
    ],
    "XGBoost": [
        ("XGB01", {"n_estimators": 300, "max_depth": 4, "learning_rate": 0.10}),
        ("XGB02", {"n_estimators": 500, "max_depth": 6, "learning_rate": 0.05}),
        ("XGB03", {"n_estimators": 800, "max_depth": 4, "learning_rate": 0.03}),
    ],
    "MLP": [
        ("MLP01", {"hidden_layer_sizes": (256,), "alpha": 1e-4, "learning_rate_init": 1e-3}),
        ("MLP02", {"hidden_layer_sizes": (512, 128), "alpha": 1e-4, "learning_rate_init": 5e-4}),
        ("MLP03", {"hidden_layer_sizes": (256, 128, 64), "alpha": 1e-3, "learning_rate_init": 5e-4}),
    ],
}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--project-root", default=".")
    p.add_argument("--task", required=True)
    p.add_argument("--task-type", required=True, choices=["binary", "multiclass"])
    p.add_argument("--redundancy", required=True, choices=["R90", "R70", "R50", "R30"])
    p.add_argument("--split-threshold", required=True, choices=["S90", "S70", "S50", "S30"])
    p.add_argument("--split-file", required=True)
    p.add_argument("--representation", required=True, choices=sorted(REPRESENTATIONS))
    p.add_argument("--algorithm", action="append", choices=ALGORITHMS + ["all"], default=[])
    p.add_argument("--split-col")
    p.add_argument("--label-col")
    p.add_argument("--positive-label")
    p.add_argument(
        "--stage10a",
        default="data/representations/sequence_universe/unique_sequences.tsv.gz",
    )
    p.add_argument("--inner-seed", type=int, default=42)
    p.add_argument("--random-state", type=int, default=42)
    p.add_argument("--n-jobs", type=int, default=-1)
    p.add_argument("--output-root", default="results/training")
    return p.parse_args()


def build_model(algorithm: str, params: dict, *, task_type: str, seed: int, n_jobs: int, y_train):
    if algorithm == "LR":
        return LogisticRegression(
            C=params["C"],
            penalty="l2",
            solver="lbfgs",
            max_iter=2000,
            class_weight="balanced",
            random_state=seed,
        )

    if algorithm == "RF":
        return RandomForestClassifier(
            n_estimators=params["n_estimators"],
            max_depth=params["max_depth"],
            min_samples_leaf=params["min_samples_leaf"],
            max_features="sqrt",
            class_weight="balanced_subsample",
            n_jobs=n_jobs,
            random_state=seed,
        )

    if algorithm == "KNN":
        return KNeighborsClassifier(
            n_neighbors=params["n_neighbors"],
            weights="distance",
            metric="cosine",
            algorithm="brute",
            n_jobs=n_jobs,
        )

    if algorithm == "XGBoost":
        try:
            from xgboost import XGBClassifier
        except ImportError as exc:
            raise RuntimeError("XGBoost requested but xgboost is not installed.") from exc

        common = dict(
            n_estimators=params["n_estimators"],
            max_depth=params["max_depth"],
            learning_rate=params["learning_rate"],
            subsample=0.8,
            colsample_bytree=0.8,
            tree_method="hist",
            n_jobs=n_jobs,
            random_state=seed,
            eval_metric="logloss" if task_type == "binary" else "mlogloss",
        )

        if task_type == "binary":
            positives = max(1, int(np.sum(np.asarray(y_train) == 1)))
            negatives = max(1, int(np.sum(np.asarray(y_train) == 0)))
            common["scale_pos_weight"] = negatives / positives

        return XGBClassifier(**common)

    if algorithm == "MLP":
        return MLPClassifier(
            hidden_layer_sizes=params["hidden_layer_sizes"],
            alpha=params["alpha"],
            learning_rate_init=params["learning_rate_init"],
            max_iter=300,
            batch_size=256,
            random_state=seed,
        )

    raise ValueError(algorithm)


def fit_model(model, X, y, algorithm: str):
    # Balanced sample weights for algorithms that support fit(sample_weight).
    if algorithm in {"XGBoost", "MLP"}:
        try:
            sig = inspect.signature(model.fit)
            if "sample_weight" in sig.parameters:
                weights = compute_sample_weight(class_weight="balanced", y=y)
                return model.fit(X, y, sample_weight=weights)
        except Exception:
            pass
    return model.fit(X, y)


def probability_output(model, X, task_type: str):
    if not hasattr(model, "predict_proba"):
        raise RuntimeError(f"{type(model).__name__} does not expose predict_proba.")

    proba = model.predict_proba(X)
    pred = model.predict(X)

    if task_type == "binary":
        classes = list(model.classes_)
        if 1 not in classes:
            raise RuntimeError(f"Binary classifier classes do not include 1: {classes}")
        pos_idx = classes.index(1)
        return pred.astype(int), proba[:, pos_idx].astype(np.float32)

    return pred, proba.astype(np.float32)


def select_best(validation_df: pd.DataFrame, task_type: str) -> pd.Series:
    if task_type == "binary":
        ranked = validation_df.sort_values(
            ["mcc", "average_precision", "config_id"],
            ascending=[False, False, True],
            kind="stable",
        )
    else:
        ranked = validation_df.sort_values(
            ["mcc", "macro_f1", "config_id"],
            ascending=[False, False, True],
            kind="stable",
        )
    return ranked.iloc[0]


def main():
    args = parse_args()
    root = Path(args.project_root).resolve()
    split_file = (root / args.split_file).resolve() if not Path(args.split_file).is_absolute() else Path(args.split_file)
    stage10a = (root / args.stage10a).resolve() if not Path(args.stage10a).is_absolute() else Path(args.stage10a)

    algorithms = ALGORITHMS if (not args.algorithm or "all" in args.algorithm) else args.algorithm

    outer, resolved_label_col = load_outer_split(
        split_file,
        split_col=args.split_col,
        label_col=args.label_col,
    )

    eligible = common_eligible_ids(stage10a)
    outer = outer.loc[outer["sequence_uid"].isin(eligible)].copy()

    if args.task_type == "binary":
        outer["y"] = encode_binary_label(outer["raw_label"], args.positive_label)
        label_encoder = None
    else:
        label_encoder = LabelEncoder()
        outer["y"] = label_encoder.fit_transform(outer["raw_label"].astype(str))

    outer_train = outer.loc[outer["outer_split"] == "train"].copy()
    outer_test = outer.loc[outer["outer_split"] == "test"].copy()

    if outer_train.empty or outer_test.empty:
        raise RuntimeError("Outer train/test became empty after common representation eligibility.")

    inner_path = (
        root / "data/training_splits" / args.task
        / f"{args.redundancy}_{args.split_threshold}"
        / f"inner_seed{args.inner_seed}.tsv.gz"
    )
    inner = create_or_load_inner_split(
        outer_train=outer_train,
        y=outer_train["y"],
        output_path=inner_path,
        seed=args.inner_seed,
    )

    outer_train = outer_train.merge(inner, on="sequence_uid", how="left", validate="1:1")
    if outer_train["inner_split"].isna().any():
        raise RuntimeError("Missing inner split assignments.")

    # Load representation once for all algorithms/configurations.
    all_ids = outer_train["sequence_uid"].tolist() + outer_test["sequence_uid"].tolist()
    X_all, feature_names, is_sparse, representation_path = load_representation(
        project_root=root,
        representation=args.representation,
        requested_ids=all_ids,
    )

    n_train = len(outer_train)
    X_outer_train = X_all[:n_train]
    X_outer_test = X_all[n_train:]

    inner_train_mask = outer_train["inner_split"].eq("train").to_numpy()
    inner_val_mask = outer_train["inner_split"].eq("validation").to_numpy()

    X_inner_train = X_outer_train[inner_train_mask]
    X_inner_val = X_outer_train[inner_val_mask]
    y_inner_train = outer_train.loc[inner_train_mask, "y"].to_numpy()
    y_inner_val = outer_train.loc[inner_val_mask, "y"].to_numpy()
    y_outer_train = outer_train["y"].to_numpy()
    y_outer_test = outer_test["y"].to_numpy()

    base_out = (
        root / args.output_root / args.task
        / f"{args.redundancy}_{args.split_threshold}"
        / args.representation
    )
    base_out.mkdir(parents=True, exist_ok=True)

    for algorithm in algorithms:
        alg_out = base_out / algorithm
        alg_out.mkdir(parents=True, exist_ok=True)

        validation_rows = []
        fitted_by_config = {}

        for config_id, params in CONFIGS[algorithm]:
            start = time.monotonic()
            standardize = should_standardize(args.representation, algorithm)

            Xi_tr, Xi_val, scaler = fit_transform_scaler(
                X_inner_train,
                X_inner_val,
                standardize=standardize,
            )

            model = build_model(
                algorithm,
                params,
                task_type=args.task_type,
                seed=args.random_state,
                n_jobs=args.n_jobs,
                y_train=y_inner_train,
            )
            fit_model(model, Xi_tr, y_inner_train, algorithm)
            y_pred, y_score = probability_output(model, Xi_val, args.task_type)

            if args.task_type == "binary":
                metrics = binary_metrics(y_inner_val, y_pred, y_score)
            else:
                metrics = multiclass_metrics(y_inner_val, y_pred)

            validation_rows.append({
                "config_id": config_id,
                "algorithm": algorithm,
                "runtime_seconds": time.monotonic() - start,
                **params,
                **metrics,
            })

        validation_df = pd.DataFrame(validation_rows)
        validation_df.to_csv(alg_out / "validation_metrics.tsv", sep="\t", index=False)

        selected = select_best(validation_df, args.task_type)
        selected_id = str(selected["config_id"])
        selected_params = dict(next(params for cid, params in CONFIGS[algorithm] if cid == selected_id))

        # Final refit on complete outer train.
        start = time.monotonic()
        standardize = should_standardize(args.representation, algorithm)
        Xtr, Xte, scaler = fit_transform_scaler(
            X_outer_train,
            X_outer_test,
            standardize=standardize,
        )

        final_model = build_model(
            algorithm,
            selected_params,
            task_type=args.task_type,
            seed=args.random_state,
            n_jobs=args.n_jobs,
            y_train=y_outer_train,
        )
        fit_model(final_model, Xtr, y_outer_train, algorithm)
        y_pred, y_score = probability_output(final_model, Xte, args.task_type)
        final_runtime = time.monotonic() - start

        if args.task_type == "binary":
            metrics = binary_metrics(y_outer_test, y_pred, y_score)
            pred_df = pd.DataFrame({
                "sequence_uid": outer_test["sequence_uid"].to_numpy(),
                "y_true": y_outer_test,
                "y_pred": y_pred,
                "score_positive": y_score,
            })
            cm = confusion_matrix(y_outer_test, y_pred, labels=[0, 1])
            pd.DataFrame(cm, index=["true_0", "true_1"], columns=["pred_0", "pred_1"]).to_csv(
                alg_out / "confusion_matrix.tsv", sep="\t"
            )
        else:
            metrics = multiclass_metrics(y_outer_test, y_pred)
            pred_df = pd.DataFrame({
                "sequence_uid": outer_test["sequence_uid"].to_numpy(),
                "y_true": y_outer_test,
                "y_pred": y_pred,
            })
            for j, cls in enumerate(final_model.classes_):
                name = (
                    str(label_encoder.inverse_transform([int(cls)])[0])
                    if label_encoder is not None else str(cls)
                )
                pred_df[f"p_{name}"] = y_score[:, j]

            labels = np.arange(len(label_encoder.classes_))
            cm = confusion_matrix(y_outer_test, y_pred, labels=labels)
            pd.DataFrame(
                cm,
                index=[f"true_{x}" for x in label_encoder.classes_],
                columns=[f"pred_{x}" for x in label_encoder.classes_],
            ).to_csv(alg_out / "confusion_matrix.tsv", sep="\t")

        pd.DataFrame([{
            "task": args.task,
            "task_type": args.task_type,
            "redundancy": args.redundancy,
            "split_threshold": args.split_threshold,
            "representation": args.representation,
            "algorithm": algorithm,
            "selected_config": selected_id,
            "n_outer_train": len(outer_train),
            "n_outer_test": len(outer_test),
            "runtime_seconds": final_runtime,
            **metrics,
        }]).to_csv(alg_out / "outer_test_metrics.tsv", sep="\t", index=False)

        pred_df.to_parquet(
            alg_out / "outer_test_predictions.parquet",
            index=False,
            compression="zstd",
        )

        write_json(
            {
                "task": args.task,
                "task_type": args.task_type,
                "redundancy": args.redundancy,
                "split_threshold": args.split_threshold,
                "representation": args.representation,
                "representation_file": str(representation_path),
                "representation_sha256": sha256_file(Path(representation_path)),
                "representation_sparse": bool(is_sparse),
                "feature_dimension": int(X_outer_train.shape[1]),
                "algorithm": algorithm,
                "selected_config": selected_id,
                "selected_params": selected_params,
                "inner_seed": args.inner_seed,
                "random_state": args.random_state,
                "split_file": str(split_file),
                "split_file_sha256": sha256_file(split_file),
                "resolved_label_column": resolved_label_col,
                "n_outer_train": int(len(outer_train)),
                "n_outer_test": int(len(outer_test)),
                "status": "completed",
            },
            alg_out / "manifest.json",
        )

        print(
            f"PASS {args.task} {args.redundancy}_{args.split_threshold} "
            f"{args.representation} {algorithm} -> {selected_id}"
        )


if __name__ == "__main__":
    main()

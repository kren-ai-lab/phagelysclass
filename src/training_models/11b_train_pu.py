#!/usr/bin/env python3
"""
PhageLysClass Stage 11B — positive–unlabeled training.

Implemented:
- naive P-vs-U diagnostic baseline
- BaggingPU
- Elkan–Noto

nnPU is intentionally not implemented until the class-prior estimator is frozen.
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.neural_network import MLPClassifier
from sklearn.neighbors import KNeighborsClassifier

from _training_utils import (
    REPRESENTATIONS,
    binary_metrics,
    common_eligible_ids,
    create_or_load_inner_split,
    encode_binary_label,
    fit_transform_scaler,
    load_outer_split,
    load_representation,
    sha256_file,
    should_standardize,
    write_json,
)

BASE_CONFIGS = {
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

METHOD_ALGORITHMS = {
    "naive": ["LR", "RF", "KNN", "XGBoost", "MLP"],
    "bagging": ["LR", "RF", "XGBoost"],
    "elkan_noto": ["LR", "RF", "XGBoost", "MLP"],
    "nnpu": ["MLP"],
}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--project-root", default=".")
    p.add_argument("--task", required=True, choices=["pu_conservative", "pu_broad"])
    p.add_argument("--redundancy", required=True, choices=["R90", "R70", "R50", "R30"])
    p.add_argument("--split-threshold", required=True, choices=["S90", "S70", "S50", "S30"])
    p.add_argument("--split-file", required=True)
    p.add_argument("--representation", required=True, choices=sorted(REPRESENTATIONS))
    p.add_argument("--method", action="append", choices=["naive", "bagging", "elkan_noto", "nnpu", "all"])
    p.add_argument("--algorithm", action="append", choices=["LR", "RF", "KNN", "XGBoost", "MLP"])
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
    p.add_argument("--bags", type=int, default=30)
    p.add_argument("--bag-u-ratio", type=float, default=1.0)
    p.add_argument("--output-root", default="results/training_pu")
    return p.parse_args()


def build_model(algorithm, params, seed, n_jobs, y_train):
    if algorithm == "LR":
        return LogisticRegression(
            C=params["C"], penalty="l2", solver="lbfgs",
            max_iter=2000, class_weight="balanced", random_state=seed,
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
            n_neighbors=params["n_neighbors"], weights="distance",
            metric="cosine", algorithm="brute", n_jobs=n_jobs,
        )
    if algorithm == "XGBoost":
        try:
            from xgboost import XGBClassifier
        except ImportError as exc:
            raise RuntimeError("xgboost is required.") from exc
        positives = max(1, int(np.sum(np.asarray(y_train) == 1)))
        negatives = max(1, int(np.sum(np.asarray(y_train) == 0)))
        return XGBClassifier(
            n_estimators=params["n_estimators"],
            max_depth=params["max_depth"],
            learning_rate=params["learning_rate"],
            subsample=0.8,
            colsample_bytree=0.8,
            tree_method="hist",
            scale_pos_weight=negatives / positives,
            n_jobs=n_jobs,
            random_state=seed,
            eval_metric="logloss",
        )
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


def positive_scores(model, X):
    proba = model.predict_proba(X)
    classes = list(model.classes_)
    return proba[:, classes.index(1)].astype(np.float32)


def pu_diagnostics(y_observed, scores):
    y = np.asarray(y_observed).astype(int)
    s = np.asarray(scores, dtype=float)

    p_scores = s[y == 1]
    u_scores = s[y == 0]

    out = {
        "n_positive": int(len(p_scores)),
        "n_unlabeled": int(len(u_scores)),
        "median_positive_score": float(np.median(p_scores)) if len(p_scores) else np.nan,
        "median_unlabeled_score": float(np.median(u_scores)) if len(u_scores) else np.nan,
        "mean_positive_score": float(np.mean(p_scores)) if len(p_scores) else np.nan,
        "mean_unlabeled_score": float(np.mean(u_scores)) if len(u_scores) else np.nan,
    }

    for q in [0.90, 0.95, 0.99]:
        if len(u_scores):
            threshold = float(np.quantile(u_scores, q))
            out[f"u_q{int(q*100)}_threshold"] = threshold
            out[f"positive_recall_above_u_q{int(q*100)}"] = (
                float(np.mean(p_scores >= threshold)) if len(p_scores) else np.nan
            )
            out[f"unlabeled_fraction_above_u_q{int(q*100)}"] = float(np.mean(u_scores >= threshold))

    # Observed-label diagnostics only.
    observed_pred = (s >= 0.5).astype(int)
    observed = binary_metrics(y, observed_pred, s)
    out.update({f"observed_{k}": v for k, v in observed.items()})
    return out


def fit_naive(algorithm, params, X_train, y_train, seed, n_jobs):
    model = build_model(algorithm, params, seed, n_jobs, y_train)
    model.fit(X_train, y_train)
    return model


def fit_bagging(
    algorithm,
    params,
    X_train,
    y_train,
    *,
    bags,
    u_ratio,
    seed,
    n_jobs,
):
    rng = np.random.default_rng(seed)
    y = np.asarray(y_train)
    p_idx = np.flatnonzero(y == 1)
    u_idx = np.flatnonzero(y == 0)

    if len(p_idx) == 0 or len(u_idx) == 0:
        raise RuntimeError("BaggingPU requires both P and U.")

    n_u = min(len(u_idx), max(1, int(round(len(p_idx) * u_ratio))))
    models = []

    for bag in range(bags):
        sampled_u = rng.choice(u_idx, size=n_u, replace=False)
        idx = np.concatenate([p_idx, sampled_u])
        bag_y = np.concatenate([
            np.ones(len(p_idx), dtype=int),
            np.zeros(len(sampled_u), dtype=int),
        ])
        model = build_model(
            algorithm, params, seed + bag, n_jobs, bag_y
        )
        model.fit(X_train[idx], bag_y)
        models.append(model)

    return models


def bagging_scores(models, X):
    return np.mean([positive_scores(model, X) for model in models], axis=0).astype(np.float32)


def main():
    args = parse_args()
    root = Path(args.project_root).resolve()
    split_file = (root / args.split_file).resolve() if not Path(args.split_file).is_absolute() else Path(args.split_file)
    stage10a = (root / args.stage10a).resolve() if not Path(args.stage10a).is_absolute() else Path(args.stage10a)

    methods = args.method or ["naive", "bagging", "elkan_noto"]
    if "all" in methods:
        methods = ["naive", "bagging", "elkan_noto", "nnpu"]

    if "nnpu" in methods:
        print("nnPU remains PENDING_PRIOR: class-prior estimator has not been frozen.")
        methods = [m for m in methods if m != "nnpu"]

    outer, resolved_label_col = load_outer_split(
        split_file, split_col=args.split_col, label_col=args.label_col
    )
    eligible = common_eligible_ids(stage10a)
    outer = outer.loc[outer["sequence_uid"].isin(eligible)].copy()
    outer["y"] = encode_binary_label(outer["raw_label"], args.positive_label)

    outer_train = outer.loc[outer["outer_split"] == "train"].copy()
    outer_test = outer.loc[outer["outer_split"] == "test"].copy()

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

    all_ids = outer_train["sequence_uid"].tolist() + outer_test["sequence_uid"].tolist()
    X_all, _, is_sparse, representation_path = load_representation(
        project_root=root,
        representation=args.representation,
        requested_ids=all_ids,
    )

    n_train = len(outer_train)
    X_outer_train = X_all[:n_train]
    X_outer_test = X_all[n_train:]
    y_outer_train = outer_train["y"].to_numpy()
    y_outer_test = outer_test["y"].to_numpy()

    tr_mask = outer_train["inner_split"].eq("train").to_numpy()
    val_mask = outer_train["inner_split"].eq("validation").to_numpy()

    X_inner_train = X_outer_train[tr_mask]
    X_inner_val = X_outer_train[val_mask]
    y_inner_train = y_outer_train[tr_mask]
    y_inner_val = y_outer_train[val_mask]

    base_out = (
        root / args.output_root / args.task
        / f"{args.redundancy}_{args.split_threshold}"
        / args.representation
    )

    for method in methods:
        allowed = METHOD_ALGORITHMS[method]
        algorithms = args.algorithm or allowed
        algorithms = [a for a in algorithms if a in allowed]

        for algorithm in algorithms:
            out = base_out / method / algorithm
            out.mkdir(parents=True, exist_ok=True)

            start = time.monotonic()
            standardize = should_standardize(args.representation, algorithm)

            Xi_tr, Xi_val, scaler = fit_transform_scaler(
                X_inner_train, X_inner_val, standardize=standardize
            )
            Xo_tr, Xo_te, final_scaler = fit_transform_scaler(
                X_outer_train, X_outer_test, standardize=standardize
            )

            if method == "naive":
                validation_rows = []

                for config_id, params in BASE_CONFIGS[algorithm]:
                    model = fit_naive(
                        algorithm, params, Xi_tr, y_inner_train,
                        args.random_state, args.n_jobs,
                    )
                    scores = positive_scores(model, Xi_val)
                    diag = pu_diagnostics(y_inner_val, scores)
                    validation_rows.append({
                        "config_id": config_id,
                        **params,
                        **diag,
                    })

                validation_df = pd.DataFrame(validation_rows)
                validation_df.sort_values(
                    ["observed_mcc", "observed_average_precision", "config_id"],
                    ascending=[False, False, True],
                    inplace=True,
                )
                validation_df.to_csv(out / "validation_or_method_diagnostics.tsv", sep="\t", index=False)

                selected_id = str(validation_df.iloc[0]["config_id"])
                params = dict(next(p for cid, p in BASE_CONFIGS[algorithm] if cid == selected_id))
                final_model = fit_naive(
                    algorithm, params, Xo_tr, y_outer_train,
                    args.random_state, args.n_jobs,
                )
                final_scores = positive_scores(final_model, Xo_te)
                method_meta = {"selected_config": selected_id, "selected_params": params}

            elif method == "bagging":
                # Fixed middle base configuration.
                selected_id, params = BASE_CONFIGS[algorithm][1]
                inner_models = fit_bagging(
                    algorithm, params, Xi_tr, y_inner_train,
                    bags=args.bags, u_ratio=args.bag_u_ratio,
                    seed=args.random_state, n_jobs=args.n_jobs,
                )
                val_scores = bagging_scores(inner_models, Xi_val)
                pd.DataFrame([{
                    "config_id": selected_id,
                    "bags": args.bags,
                    "u_ratio": args.bag_u_ratio,
                    **pu_diagnostics(y_inner_val, val_scores),
                }]).to_csv(out / "validation_or_method_diagnostics.tsv", sep="\t", index=False)

                final_models = fit_bagging(
                    algorithm, params, Xo_tr, y_outer_train,
                    bags=args.bags, u_ratio=args.bag_u_ratio,
                    seed=args.random_state, n_jobs=args.n_jobs,
                )
                final_scores = bagging_scores(final_models, Xo_te)
                method_meta = {
                    "selected_config": selected_id,
                    "selected_params": params,
                    "bags": args.bags,
                    "u_ratio": args.bag_u_ratio,
                }

            elif method == "elkan_noto":
                # Fixed middle base configuration. Inner validation positives
                # estimate c; U validation observations are not used to estimate c.
                selected_id, params = BASE_CONFIGS[algorithm][1]
                model = fit_naive(
                    algorithm, params, Xi_tr, y_inner_train,
                    args.random_state, args.n_jobs,
                )
                raw_val_scores = positive_scores(model, Xi_val)
                positive_val = raw_val_scores[y_inner_val == 1]
                if len(positive_val) == 0:
                    raise RuntimeError("Elkan-Noto needs held-out positives to estimate c.")
                c = float(np.mean(positive_val))
                if not (0 < c <= 1):
                    raise RuntimeError(f"Invalid Elkan-Noto c={c}")

                corrected_val = np.clip(raw_val_scores / c, 0.0, 1.0)
                pd.DataFrame([{
                    "config_id": selected_id,
                    "labeling_propensity_c": c,
                    **pu_diagnostics(y_inner_val, corrected_val),
                }]).to_csv(out / "validation_or_method_diagnostics.tsv", sep="\t", index=False)

                final_model = fit_naive(
                    algorithm, params, Xo_tr, y_outer_train,
                    args.random_state, args.n_jobs,
                )
                raw_final = positive_scores(final_model, Xo_te)
                final_scores = np.clip(raw_final / c, 0.0, 1.0).astype(np.float32)
                method_meta = {
                    "selected_config": selected_id,
                    "selected_params": params,
                    "labeling_propensity_c": c,
                }

            else:
                raise ValueError(method)

            diagnostics = pu_diagnostics(y_outer_test, final_scores)
            pd.DataFrame([{
                "task": args.task,
                "redundancy": args.redundancy,
                "split_threshold": args.split_threshold,
                "representation": args.representation,
                "method": method,
                "algorithm": algorithm,
                "runtime_seconds": time.monotonic() - start,
                **diagnostics,
            }]).to_csv(out / "pu_diagnostics.tsv", sep="\t", index=False)

            pd.DataFrame({
                "sequence_uid": outer_test["sequence_uid"].to_numpy(),
                "observed_role": np.where(y_outer_test == 1, "P", "U"),
                "score": final_scores,
            }).to_parquet(
                out / "outer_test_scores.parquet",
                index=False,
                compression="zstd",
            )

            write_json({
                "task": args.task,
                "redundancy": args.redundancy,
                "split_threshold": args.split_threshold,
                "representation": args.representation,
                "representation_file": str(representation_path),
                "representation_sha256": sha256_file(Path(representation_path)),
                "representation_sparse": bool(is_sparse),
                "method": method,
                "algorithm": algorithm,
                "inner_seed": args.inner_seed,
                "random_state": args.random_state,
                "split_file": str(split_file),
                "split_file_sha256": sha256_file(split_file),
                "resolved_label_column": resolved_label_col,
                "n_outer_train": int(len(outer_train)),
                "n_outer_test": int(len(outer_test)),
                "status": "completed",
                **method_meta,
            }, out / "manifest.json")

            print(
                f"PASS {args.task} {args.redundancy}_{args.split_threshold} "
                f"{args.representation} {method}/{algorithm}"
            )


if __name__ == "__main__":
    main()

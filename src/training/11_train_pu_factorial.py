#!/usr/bin/env python3
"""
11_train_pu_factorial.py

Production PU trainer.

Scientific CLI parameters:
    --task            pu_conservative | pu_broad
    --representation

Mandatory production exploration:
- all 16 R × S cells;
- naive P-vs-U with LR/RF/KNN/XGB/MLP and all 3 configs;
- BaggingPU with LR02/RF02/XGB02, 30 bags, sampled U:P=1:1;
- Elkan–Noto with LR02/RF02/XGB02/MLP02 and group-aware 5-fold
  outer-train-only propensity estimation;
- nnPU registered as PENDING_PRIOR until the prior-estimation contract is frozen.
"""

from __future__ import annotations

import argparse
import json
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedGroupKFold

from _training_common import (
    PARAMETER_PANELS,
    PU_FIXED_CONFIG,
    PU_TASKS,
    R_VALUES,
    S_VALUES,
    SEED,
    TASKS,
    RepresentationStore,
    binary_metrics,
    build_model,
    choose_best_config,
    fit_model,
    get_config,
    load_representation_registry,
    positive_scores,
    process_peak_rss_kb,
    project_root_from_script,
    pu_diagnostics,
    pu_target,
    read_tsv,
    sha256_file,
    software_versions,
    write_json,
)


NAIVE_ALGORITHMS = ["LR", "RF", "KNN", "XGB", "MLP"]
BAGGING_ALGORITHMS = ["LR", "RF", "XGB"]
EN_ALGORITHMS = ["LR", "RF", "XGB", "MLP"]
N_BAGS = 30
U_P_RATIO = 1.0
EN_FOLDS = 5


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--task", required=True, choices=PU_TASKS)
    p.add_argument("--representation", required=True)
    p.add_argument("--overwrite", action="store_true")
    p.add_argument(
        "--smoke",
        action="store_true",
        help="Implementation test only: R90_S90 + naive LR.",
    )
    return p.parse_args()


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
        raise FileNotFoundError(f"{inner}\nRun Stage 10F first.")

    df = read_tsv(universe)
    inner_df = read_tsv(inner)
    assignment = dict(
        zip(inner_df["sequence_uid"], inner_df["inner_split"])
    )
    df["inner_split"] = df["sequence_uid"].map(assignment).fillna("")
    return df, universe, inner


def save_pu_scores(
    path: Path,
    metadata: pd.DataFrame,
    y_obs,
    scores,
    *,
    experiment_id: str,
    raw_scores=None,
    score_std=None,
):
    out = metadata[
        ["entity_id", "sequence_uid", "cluster_id"]
    ].reset_index(drop=True).copy()
    out.insert(0, "experiment_id", experiment_id)
    out["observed_role"] = np.where(np.asarray(y_obs) == 1, "P", "U")
    out["score_positive"] = np.asarray(scores, dtype=float)
    if raw_scores is not None:
        out["raw_observed_label_score"] = np.asarray(raw_scores, dtype=float)
    if score_std is not None:
        out["score_std"] = np.asarray(score_std, dtype=float)
    out.to_parquet(path, index=False, compression="zstd")


def run_naive(
    *,
    root,
    task_id,
    task_path,
    r,
    s,
    rep,
    algorithm,
    overwrite,
):
    rs = f"R{r}_S{s}"
    exp_id = (
        f"PU__{task_id}__R{r}__S{s}__"
        f"{rep.spec.representation_id}__naive__{algorithm}"
    )
    out_dir = (
        root / "results/training_pu" / task_path / rs
        / rep.spec.representation_id / "naive" / algorithm
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = out_dir / "manifest.json"
    if manifest_path.is_file() and not overwrite:
        try:
            if json.loads(manifest_path.read_text()).get("status") == "DONE":
                print(f"SKIP {exp_id} — DONE")
                return {"experiment_id": exp_id, "status": "DONE"}
        except Exception:
            pass

    try:
        df, universe, inner = load_cell(root, task_path, r, s)
        y_all, role_col = pu_target(df)
        df = df.reset_index(drop=True)
        y_all = np.asarray(y_all, dtype=int)

        idx_inner = np.where(
            (df["split"].eq("train") & df["inner_split"].eq("train")).to_numpy()
        )[0]
        idx_val = np.where(
            (
                df["split"].eq("train")
                & df["inner_split"].eq("validation")
            ).to_numpy()
        )[0]
        idx_outer = np.where(df["split"].eq("train").to_numpy())[0]
        idx_test = np.where(df["split"].eq("test").to_numpy())[0]

        rep.require_uids(df["sequence_uid"].astype(str).tolist())
        X_inner = rep.get(df.iloc[idx_inner]["sequence_uid"].astype(str).tolist())
        X_val = rep.get(df.iloc[idx_val]["sequence_uid"].astype(str).tolist())
        X_outer = rep.get(df.iloc[idx_outer]["sequence_uid"].astype(str).tolist())
        X_test = rep.get(df.iloc[idx_test]["sequence_uid"].astype(str).tolist())

        rows = []
        pred_parts = []
        for order, cfg in enumerate(PARAMETER_PANELS[algorithm]):
            model = build_model(
                algorithm,
                cfg,
                task_type="binary",
                representation_family=rep.spec.family,
                n_classes=2,
            )
            t0 = time.perf_counter()
            fit_model(model, X_inner, y_all[idx_inner], algorithm)
            fit_sec = time.perf_counter() - t0
            scores = positive_scores(model, X_val)
            obs_metrics = binary_metrics(y_all[idx_val], scores)
            diag = pu_diagnostics(y_all[idx_val], scores)
            rows.append(
                {
                    "experiment_id": exp_id,
                    "config_id": cfg["config_id"],
                    "config_order": order,
                    "fit_seconds": fit_sec,
                    **{f"observed_label_{k}": v for k, v in obs_metrics.items()},
                    **diag,
                }
            )
            pred = df.iloc[idx_val][
                ["entity_id", "sequence_uid", "cluster_id"]
            ].reset_index(drop=True)
            pred["experiment_id"] = exp_id
            pred["config_id"] = cfg["config_id"]
            pred["observed_role"] = np.where(
                y_all[idx_val] == 1, "P", "U"
            )
            pred["score_positive"] = scores
            pred_parts.append(pred)

        metrics_df = pd.DataFrame(rows)
        # Diagnostic baseline: selection uses observed-label validation MCC/AP,
        # explicitly named as such.
        tmp = pd.DataFrame(
            {
                "config_id": metrics_df["config_id"],
                "config_order": metrics_df["config_order"],
                "mcc": metrics_df["observed_label_mcc"],
                "ap": metrics_df["observed_label_ap"],
            }
        )
        selected = choose_best_config(tmp, "binary")
        metrics_df["selected"] = metrics_df["config_id"].eq(selected)
        selected_cfg = get_config(algorithm, selected)

        final_model = build_model(
            algorithm,
            selected_cfg,
            task_type="binary",
            representation_family=rep.spec.family,
            n_classes=2,
        )
        fit_model(final_model, X_outer, y_all[idx_outer], algorithm)
        test_scores = positive_scores(final_model, X_test)
        obs_test = binary_metrics(y_all[idx_test], test_scores)
        diag_test = pu_diagnostics(y_all[idx_test], test_scores)

        metrics_df.to_csv(
            out_dir / "validation_or_method_diagnostics.tsv",
            sep="\t",
            index=False,
        )
        pd.concat(pred_parts, ignore_index=True).to_parquet(
            out_dir / "validation_predictions.parquet",
            index=False,
            compression="zstd",
        )
        write_json(out_dir / "selected_configuration.json", selected_cfg)
        pd.DataFrame(
            [
                {
                    "experiment_id": exp_id,
                    "selected_config": selected,
                    **{f"observed_label_{k}": v for k, v in obs_test.items()},
                    **diag_test,
                }
            ]
        ).to_csv(
            out_dir / "pu_diagnostics.tsv",
            sep="\t",
            index=False,
        )
        save_pu_scores(
            out_dir / "outer_test_scores.parquet",
            df.iloc[idx_test],
            y_all[idx_test],
            test_scores,
            experiment_id=exp_id,
        )
        write_json(
            manifest_path,
            {
                "stage": "11-pu",
                "status": "DONE",
                "method": "naive_p_vs_u",
                "interpretation": "observed-label diagnostic only",
                "experiment_id": exp_id,
                "task": task_id,
                "R": r,
                "S": s,
                "representation": rep.spec.representation_id,
                "algorithm": algorithm,
                "selected_config": selected,
                "outer_universe_sha256": sha256_file(universe),
                "inner_split_sha256": sha256_file(inner),
                "software": software_versions(),
            },
        )
        print(f"DONE {exp_id}")
        return {"experiment_id": exp_id, "status": "DONE"}

    except Exception as exc:
        (out_dir / "error.txt").write_text(
            traceback.format_exc(), encoding="utf-8"
        )
        write_json(
            manifest_path,
            {
                "stage": "11-pu",
                "status": "FAILED",
                "method": "naive_p_vs_u",
                "experiment_id": exp_id,
                "error": str(exc),
            },
        )
        print(f"FAILED {exp_id}: {exc}")
        return {"experiment_id": exp_id, "status": "FAILED", "error": str(exc)}


def run_bagging(
    *,
    root,
    task_id,
    task_path,
    r,
    s,
    rep,
    algorithm,
    overwrite,
):
    rs = f"R{r}_S{s}"
    exp_id = (
        f"PU__{task_id}__R{r}__S{s}__"
        f"{rep.spec.representation_id}__baggingpu__{algorithm}"
    )
    out_dir = (
        root / "results/training_pu" / task_path / rs
        / rep.spec.representation_id / "baggingpu" / algorithm
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = out_dir / "manifest.json"
    if manifest_path.is_file() and not overwrite:
        try:
            if json.loads(manifest_path.read_text()).get("status") == "DONE":
                print(f"SKIP {exp_id} — DONE")
                return {"experiment_id": exp_id, "status": "DONE"}
        except Exception:
            pass

    try:
        df, universe, inner = load_cell(root, task_path, r, s)
        y_all, role_col = pu_target(df)
        df = df.reset_index(drop=True)
        y_all = np.asarray(y_all, dtype=int)

        idx_outer = np.where(df["split"].eq("train").to_numpy())[0]
        idx_test = np.where(df["split"].eq("test").to_numpy())[0]
        y_outer = y_all[idx_outer]
        pos_local = np.where(y_outer == 1)[0]
        unl_local = np.where(y_outer == 0)[0]
        if len(pos_local) == 0 or len(unl_local) == 0:
            raise RuntimeError("BaggingPU requires both P and U")

        rep.require_uids(df["sequence_uid"].astype(str).tolist())
        X_outer = rep.get(df.iloc[idx_outer]["sequence_uid"].astype(str).tolist())
        X_test = rep.get(df.iloc[idx_test]["sequence_uid"].astype(str).tolist())

        cfg = get_config(algorithm, PU_FIXED_CONFIG[algorithm])
        rng = np.random.default_rng(SEED)
        bag_scores = []
        bag_rows = []

        n_u_sample = int(round(len(pos_local) * U_P_RATIO))
        if n_u_sample > len(unl_local):
            raise RuntimeError(
                f"Not enough U for 1:1 BaggingPU without replacement: "
                f"P={len(pos_local)}, U={len(unl_local)}"
            )

        for bag in range(N_BAGS):
            bag_start = time.perf_counter()
            sampled_u = rng.choice(
                unl_local, size=n_u_sample, replace=False
            )
            bag_local = np.concatenate([pos_local, sampled_u])
            bag_y = np.concatenate(
                [
                    np.ones(len(pos_local), dtype=int),
                    np.zeros(len(sampled_u), dtype=int),
                ]
            )
            X_bag = X_outer[bag_local]

            model = build_model(
                algorithm,
                cfg,
                task_type="binary",
                representation_family=rep.spec.family,
                n_classes=2,
            )
            # Sample is already balanced; avoid additional XGB weighting.
            if algorithm == "XGB":
                model.fit(X_bag, bag_y)
            else:
                model.fit(X_bag, bag_y)
            scores = positive_scores(model, X_test)
            bag_scores.append(scores)
            bag_rows.append(
                {
                    "bag": bag,
                    "n_p": len(pos_local),
                    "n_sampled_u": len(sampled_u),
                    "fit_and_predict_seconds": time.perf_counter() - bag_start,
                }
            )

        score_matrix = np.vstack(bag_scores)
        mean_score = score_matrix.mean(axis=0)
        std_score = score_matrix.std(axis=0, ddof=0)
        diag = pu_diagnostics(y_all[idx_test], mean_score)

        pd.DataFrame(bag_rows).to_csv(
            out_dir / "bagging_diagnostics.tsv",
            sep="\t",
            index=False,
        )
        pd.DataFrame([diag]).to_csv(
            out_dir / "pu_diagnostics.tsv",
            sep="\t",
            index=False,
        )
        save_pu_scores(
            out_dir / "outer_test_scores.parquet",
            df.iloc[idx_test],
            y_all[idx_test],
            mean_score,
            experiment_id=exp_id,
            score_std=std_score,
        )
        write_json(
            manifest_path,
            {
                "stage": "11-pu",
                "status": "DONE",
                "method": "BaggingPU",
                "experiment_id": exp_id,
                "task": task_id,
                "R": r,
                "S": s,
                "representation": rep.spec.representation_id,
                "algorithm": algorithm,
                "base_config": cfg["config_id"],
                "n_bags": N_BAGS,
                "u_p_ratio": U_P_RATIO,
                "sampling": "without replacement within each bag",
                "seed": SEED,
                "outer_universe_sha256": sha256_file(universe),
                "software": software_versions(),
            },
        )
        print(f"DONE {exp_id}")
        return {"experiment_id": exp_id, "status": "DONE"}

    except Exception as exc:
        (out_dir / "error.txt").write_text(
            traceback.format_exc(), encoding="utf-8"
        )
        write_json(
            manifest_path,
            {
                "stage": "11-pu",
                "status": "FAILED",
                "method": "BaggingPU",
                "experiment_id": exp_id,
                "error": str(exc),
            },
        )
        print(f"FAILED {exp_id}: {exc}")
        return {"experiment_id": exp_id, "status": "FAILED", "error": str(exc)}


def run_elkan_noto(
    *,
    root,
    task_id,
    task_path,
    r,
    s,
    rep,
    algorithm,
    overwrite,
):
    rs = f"R{r}_S{s}"
    exp_id = (
        f"PU__{task_id}__R{r}__S{s}__"
        f"{rep.spec.representation_id}__elkan_noto__{algorithm}"
    )
    out_dir = (
        root / "results/training_pu" / task_path / rs
        / rep.spec.representation_id / "elkan_noto" / algorithm
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = out_dir / "manifest.json"
    if manifest_path.is_file() and not overwrite:
        try:
            if json.loads(manifest_path.read_text()).get("status") == "DONE":
                print(f"SKIP {exp_id} — DONE")
                return {"experiment_id": exp_id, "status": "DONE"}
        except Exception:
            pass

    try:
        df, universe, inner = load_cell(root, task_path, r, s)
        y_all, role_col = pu_target(df)
        df = df.reset_index(drop=True)
        y_all = np.asarray(y_all, dtype=int)

        idx_outer = np.where(df["split"].eq("train").to_numpy())[0]
        idx_test = np.where(df["split"].eq("test").to_numpy())[0]
        outer_df = df.iloc[idx_outer].reset_index(drop=True)
        y_outer = y_all[idx_outer]

        rep.require_uids(df["sequence_uid"].astype(str).tolist())
        X_outer = rep.get(outer_df["sequence_uid"].astype(str).tolist())
        X_test = rep.get(df.iloc[idx_test]["sequence_uid"].astype(str).tolist())

        cfg = get_config(algorithm, PU_FIXED_CONFIG[algorithm])

        splitter = StratifiedGroupKFold(
            n_splits=EN_FOLDS,
            shuffle=True,
            random_state=SEED,
        )
        groups = outer_df["cluster_id"].astype(str).to_numpy()
        dummy = np.zeros((len(outer_df), 1), dtype=np.uint8)

        positive_oof_scores = []
        fold_rows = []

        for fold, (fit_idx, hold_idx) in enumerate(
            splitter.split(dummy, y_outer, groups)
        ):
            hold_positive = hold_idx[y_outer[hold_idx] == 1]
            if len(hold_positive) == 0:
                raise RuntimeError(
                    f"Elkan-Noto fold {fold} contains no held-out positives"
                )

            model = build_model(
                algorithm,
                cfg,
                task_type="binary",
                representation_family=rep.spec.family,
                n_classes=2,
            )
            fit_model(model, X_outer[fit_idx], y_outer[fit_idx], algorithm)
            scores = positive_scores(model, X_outer[hold_positive])
            positive_oof_scores.extend(scores.tolist())
            fold_rows.append(
                {
                    "fold": fold,
                    "n_fit": len(fit_idx),
                    "n_holdout": len(hold_idx),
                    "n_holdout_p": len(hold_positive),
                    "c_fold_mean_positive_score": float(np.mean(scores)),
                }
            )

        c = float(np.mean(positive_oof_scores))
        if not np.isfinite(c) or c <= 0:
            raise RuntimeError(f"Invalid Elkan-Noto c={c}")

        final_model = build_model(
            algorithm,
            cfg,
            task_type="binary",
            representation_family=rep.spec.family,
            n_classes=2,
        )
        fit_model(final_model, X_outer, y_outer, algorithm)
        raw_score = positive_scores(final_model, X_test)
        corrected = np.clip(raw_score / c, 0.0, 1.0)
        clipping_fraction = float(np.mean((raw_score / c) > 1.0))

        diag = pu_diagnostics(y_all[idx_test], corrected)
        diag["elkan_noto_c"] = c
        diag["corrected_score_clipping_fraction"] = clipping_fraction

        pd.DataFrame(fold_rows).to_csv(
            out_dir / "propensity_crossfit.tsv",
            sep="\t",
            index=False,
        )
        pd.DataFrame([diag]).to_csv(
            out_dir / "pu_diagnostics.tsv",
            sep="\t",
            index=False,
        )
        save_pu_scores(
            out_dir / "outer_test_scores.parquet",
            df.iloc[idx_test],
            y_all[idx_test],
            corrected,
            experiment_id=exp_id,
            raw_scores=raw_score,
        )
        write_json(
            manifest_path,
            {
                "stage": "11-pu",
                "status": "DONE",
                "method": "Elkan-Noto",
                "experiment_id": exp_id,
                "task": task_id,
                "R": r,
                "S": s,
                "representation": rep.spec.representation_id,
                "algorithm": algorithm,
                "base_config": cfg["config_id"],
                "crossfit_folds": EN_FOLDS,
                "group_aware_crossfit": True,
                "seed": SEED,
                "c": c,
                "outer_universe_sha256": sha256_file(universe),
                "software": software_versions(),
            },
        )
        print(f"DONE {exp_id} | c={c:.4f}")
        return {"experiment_id": exp_id, "status": "DONE"}

    except Exception as exc:
        (out_dir / "error.txt").write_text(
            traceback.format_exc(), encoding="utf-8"
        )
        write_json(
            manifest_path,
            {
                "stage": "11-pu",
                "status": "FAILED",
                "method": "Elkan-Noto",
                "experiment_id": exp_id,
                "error": str(exc),
            },
        )
        print(f"FAILED {exp_id}: {exc}")
        return {"experiment_id": exp_id, "status": "FAILED", "error": str(exc)}


def register_nnpu_pending(root, task_id, task_path, r, s, rep):
    rs = f"R{r}_S{s}"
    exp_id = (
        f"PU__{task_id}__R{r}__S{s}__"
        f"{rep.spec.representation_id}__nnpu__MLP"
    )
    out_dir = (
        root / "results/training_pu" / task_path / rs
        / rep.spec.representation_id / "nnpu" / "MLP"
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    write_json(
        out_dir / "manifest.json",
        {
            "stage": "11-pu",
            "status": "PENDING_PRIOR",
            "method": "nnPU",
            "experiment_id": exp_id,
            "task": task_id,
            "R": r,
            "S": s,
            "representation": rep.spec.representation_id,
            "algorithm": "MLP",
            "reason": "Positive-class prior estimator not yet frozen.",
        },
    )
    return {"experiment_id": exp_id, "status": "PENDING_PRIOR"}


def main():
    args = parse_args()
    root = project_root_from_script(__file__)
    task_path = TASKS[args.task]["path"]

    registry = load_representation_registry(
        root / "config/representation_registry/representation_registry.tsv",
        root,
    )
    if args.representation not in registry:
        raise KeyError(args.representation)

    rep = RepresentationStore(registry[args.representation])

    print("=" * 100)
    print("PhageLysClass — Stage 11 PU factorial training")
    print("=" * 100)
    print(f"Task           : {args.task}")
    print(f"Representation : {args.representation}")
    if args.smoke:
        print("MODE           : SMOKE ONLY (R90_S90 + naive LR)")
    else:
        print("Production     : all 16 R/S + naive + BaggingPU + Elkan-Noto + nnPU pending")
    print()

    results = []
    r_values = [90] if args.smoke else R_VALUES
    s_values = [90] if args.smoke else S_VALUES

    for r in r_values:
        for s in s_values:
            if args.smoke:
                results.append(
                    run_naive(
                        root=root,
                        task_id=args.task,
                        task_path=task_path,
                        r=r,
                        s=s,
                        rep=rep,
                        algorithm="LR",
                        overwrite=args.overwrite,
                    )
                )
                continue

            for algorithm in NAIVE_ALGORITHMS:
                results.append(
                    run_naive(
                        root=root,
                        task_id=args.task,
                        task_path=task_path,
                        r=r,
                        s=s,
                        rep=rep,
                        algorithm=algorithm,
                        overwrite=args.overwrite,
                    )
                )
            for algorithm in BAGGING_ALGORITHMS:
                results.append(
                    run_bagging(
                        root=root,
                        task_id=args.task,
                        task_path=task_path,
                        r=r,
                        s=s,
                        rep=rep,
                        algorithm=algorithm,
                        overwrite=args.overwrite,
                    )
                )
            for algorithm in EN_ALGORITHMS:
                results.append(
                    run_elkan_noto(
                        root=root,
                        task_id=args.task,
                        task_path=task_path,
                        r=r,
                        s=s,
                        rep=rep,
                        algorithm=algorithm,
                        overwrite=args.overwrite,
                    )
                )
            results.append(
                register_nnpu_pending(
                    root, args.task, task_path, r, s, rep
                )
            )

    summary = pd.DataFrame(results)
    out = (
        root / "results/training_pu" / task_path
        / f"_run_summary__{args.representation}.tsv"
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(out, sep="\t", index=False)

    print()
    print("=" * 100)
    print(f"Run summary: {summary['status'].value_counts().to_dict()}")
    print("=" * 100)


if __name__ == "__main__":
    main()

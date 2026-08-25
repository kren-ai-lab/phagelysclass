#!/usr/bin/env python3
"""
11_train_pu_factorial.py

Compact Protein Science PU benchmark for PhageLysClass.

Scientific CLI parameters
-------------------------
--task
    pu_conservative
    pu_broad

--representation
    any READY representation in:
    config/representation_registry/representation_registry.tsv

Mandatory exploration
---------------------
For the selected PU task and representation, the script runs:

    16 R/S cells
    × 2 base learners: LR02, XGB02
    × 3 PU strategies:
        1. naive P-vs-U diagnostic
        2. BaggingPU
        3. Elkan–Noto

No R, S, algorithm, method, or hyperparameter option is exposed.

Why only these PU tasks?
------------------------
global_easy/global_hard/global_mixed remain supervised reference tasks because
their non-positive classes are evidence-defined negatives. They are combined
with pu_conservative/pu_broad later at the analysis stage, yielding five
global analytical conditions without scientifically mislabeling explicit
negatives as unlabeled.
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
    R_VALUES,
    S_VALUES,
    SEED,
    TASKS,
    RepresentationStore,
    build_model,
    fit_model,
    get_config,
    load_representation_registry,
    positive_scores,
    project_root_from_script,
    pu_diagnostics,
    pu_target,
    read_tsv,
    sha256_file,
    software_versions,
    write_json,
)

PU_TASKS = ["pu_conservative", "pu_broad"]
BASE_ALGORITHMS = ["LR", "XGB"]
BASE_CONFIGS = {"LR": "LR02", "XGB": "XGB02"}
METHODS = ["naive", "baggingpu", "elkan_noto"]

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
        help="Operational smoke test only: R90_S90, LR02, naive P-vs-U.",
    )
    return p.parse_args()


def load_cell(root: Path, task_path: str, r: int, s: int):
    rs = f"R{r}_S{s}"
    universe = root / "data/predictive_universe" / task_path / rs / "eligible_entities.tsv.gz"
    inner = root / "data/training_splits" / task_path / rs / "inner_seed42.tsv.gz"
    if not universe.is_file():
        raise FileNotFoundError(universe)
    if not inner.is_file():
        raise FileNotFoundError(f"{inner}\nStage 10F must be frozen before PU training.")
    df = read_tsv(universe).reset_index(drop=True)
    required = {"entity_id", "sequence_uid", "cluster_id", "split"}
    missing = required - set(df.columns)
    if missing:
        raise KeyError(f"{universe}: missing required columns {sorted(missing)}")
    return df, universe, inner


def output_dir(root, task_path, r, s, representation, method, algorithm):
    return root / "results/training_pu" / task_path / f"R{r}_S{s}" / representation / method / algorithm


def done(out_dir: Path) -> bool:
    manifest = out_dir / "manifest.json"
    if not manifest.is_file():
        return False
    try:
        return json.loads(manifest.read_text()).get("status") == "DONE"
    except Exception:
        return False


def save_scores(path, metadata, y_obs, score, *, experiment_id, raw_score=None, score_std=None):
    out = metadata[["entity_id", "sequence_uid", "cluster_id"]].reset_index(drop=True).copy()
    out.insert(0, "experiment_id", experiment_id)
    out["observed_role"] = np.where(np.asarray(y_obs) == 1, "P", "U")
    out["score_positive"] = np.asarray(score, dtype=float)
    if raw_score is not None:
        out["raw_observed_label_score"] = np.asarray(raw_score, dtype=float)
    if score_std is not None:
        out["score_std"] = np.asarray(score_std, dtype=float)
    out.to_parquet(path, index=False, compression="zstd")


def common_manifest(*, method, experiment_id, task_id, r, s, representation, algorithm, config_id, universe, inner):
    return {
        "stage": "11-pu",
        "status": "DONE",
        "method": method,
        "experiment_id": experiment_id,
        "task": task_id,
        "R": r,
        "S": s,
        "representation": representation.spec.representation_id,
        "representation_family": representation.spec.family,
        "representation_dimension": representation.dimension,
        "representation_sparse": representation.sparse,
        "algorithm": algorithm,
        "base_config": config_id,
        "seed": SEED,
        "outer_universe_path": str(universe),
        "outer_universe_sha256": sha256_file(universe),
        "inner_split_path": str(inner),
        "inner_split_sha256": sha256_file(inner),
        "representation_path": str(representation.spec.main_path),
        "representation_sha256": sha256_file(representation.spec.main_path),
        "software": software_versions(),
    }


def run_naive(*, root, task_id, task_path, r, s, rep, algorithm, overwrite):
    method = "naive"
    config_id = BASE_CONFIGS[algorithm]
    exp_id = f"PU__{task_id}__R{r}__S{s}__{rep.spec.representation_id}__{method}__{algorithm}"
    out = output_dir(root, task_path, r, s, rep.spec.representation_id, method, algorithm)
    out.mkdir(parents=True, exist_ok=True)
    if done(out) and not overwrite:
        print(f"SKIP {exp_id} — DONE")
        return {"experiment_id": exp_id, "status": "DONE"}
    try:
        df, universe, inner = load_cell(root, task_path, r, s)
        y_all, _ = pu_target(df)
        idx_train = np.where(df["split"].eq("train").to_numpy())[0]
        idx_test = np.where(df["split"].eq("test").to_numpy())[0]
        rep.require_uids(df["sequence_uid"].astype(str).tolist())
        X_train = rep.get(df.iloc[idx_train]["sequence_uid"].astype(str).tolist())
        X_test = rep.get(df.iloc[idx_test]["sequence_uid"].astype(str).tolist())
        y_train, y_test = y_all[idx_train], y_all[idx_test]
        cfg = get_config(algorithm, config_id)
        model = build_model(algorithm, cfg, task_type="binary", representation_family=rep.spec.family, n_classes=2)
        t0 = time.perf_counter()
        fit_model(model, X_train, y_train, algorithm)
        fit_seconds = time.perf_counter() - t0
        t1 = time.perf_counter()
        test_score = positive_scores(model, X_test)
        pred_seconds = time.perf_counter() - t1
        diag = pu_diagnostics(y_test, test_score)
        diag.update({"method": method, "algorithm": algorithm, "base_config": config_id, "fit_seconds": fit_seconds, "prediction_seconds": pred_seconds})
        pd.DataFrame([diag]).to_csv(out / "pu_diagnostics.tsv", sep="\t", index=False)
        save_scores(out / "outer_test_scores.parquet", df.iloc[idx_test], y_test, test_score, experiment_id=exp_id)
        manifest = common_manifest(method=method, experiment_id=exp_id, task_id=task_id, r=r, s=s, representation=rep, algorithm=algorithm, config_id=config_id, universe=universe, inner=inner)
        manifest["interpretation"] = "Observed-label P-vs-U diagnostic only. U is not a confirmed negative class."
        write_json(out / "manifest.json", manifest)
        (out / "warnings.log").write_text("", encoding="utf-8")
        print(f"DONE {exp_id}")
        return {"experiment_id": exp_id, "status": "DONE"}
    except Exception as exc:
        (out / "error.txt").write_text(traceback.format_exc(), encoding="utf-8")
        write_json(out / "manifest.json", {"stage": "11-pu", "status": "FAILED", "experiment_id": exp_id, "method": method, "error": str(exc)})
        print(f"FAILED {exp_id}: {exc}")
        return {"experiment_id": exp_id, "status": "FAILED", "error": str(exc)}


def run_baggingpu(*, root, task_id, task_path, r, s, rep, algorithm, overwrite):
    method = "baggingpu"
    config_id = BASE_CONFIGS[algorithm]
    exp_id = f"PU__{task_id}__R{r}__S{s}__{rep.spec.representation_id}__{method}__{algorithm}"
    out = output_dir(root, task_path, r, s, rep.spec.representation_id, method, algorithm)
    out.mkdir(parents=True, exist_ok=True)
    if done(out) and not overwrite:
        print(f"SKIP {exp_id} — DONE")
        return {"experiment_id": exp_id, "status": "DONE"}
    try:
        df, universe, inner = load_cell(root, task_path, r, s)
        y_all, _ = pu_target(df)
        idx_train = np.where(df["split"].eq("train").to_numpy())[0]
        idx_test = np.where(df["split"].eq("test").to_numpy())[0]
        rep.require_uids(df["sequence_uid"].astype(str).tolist())
        X_train = rep.get(df.iloc[idx_train]["sequence_uid"].astype(str).tolist())
        X_test = rep.get(df.iloc[idx_test]["sequence_uid"].astype(str).tolist())
        y_train, y_test = y_all[idx_train], y_all[idx_test]
        p_local = np.where(y_train == 1)[0]
        u_local = np.where(y_train == 0)[0]
        if len(p_local) == 0 or len(u_local) == 0:
            raise RuntimeError("BaggingPU requires both P and U.")
        n_sample_u = int(round(len(p_local) * U_P_RATIO))
        if n_sample_u > len(u_local):
            raise RuntimeError(f"Not enough U for frozen 1:1 BaggingPU: P={len(p_local):,}, U={len(u_local):,}")
        cfg = get_config(algorithm, config_id)
        rng = np.random.default_rng(SEED)
        score_parts, bag_rows = [], []
        for bag_id in range(N_BAGS):
            sampled_u = rng.choice(u_local, size=n_sample_u, replace=False)
            fit_idx = np.concatenate([p_local, sampled_u])
            fit_y = np.concatenate([np.ones(len(p_local), dtype=int), np.zeros(len(sampled_u), dtype=int)])
            model = build_model(algorithm, cfg, task_type="binary", representation_family=rep.spec.family, n_classes=2)
            t0 = time.perf_counter()
            model.fit(X_train[fit_idx], fit_y)
            score = positive_scores(model, X_test)
            elapsed = time.perf_counter() - t0
            score_parts.append(score)
            bag_rows.append({"bag_id": bag_id, "n_p": len(p_local), "n_sampled_u": len(sampled_u), "fit_plus_predict_seconds": elapsed})
        score_matrix = np.vstack(score_parts)
        mean_score = score_matrix.mean(axis=0)
        std_score = score_matrix.std(axis=0, ddof=0)
        diag = pu_diagnostics(y_test, mean_score)
        diag.update({"method": method, "algorithm": algorithm, "base_config": config_id, "n_bags": N_BAGS, "u_p_ratio": U_P_RATIO, "median_test_score_std": float(np.median(std_score)), "mean_test_score_std": float(np.mean(std_score))})
        pd.DataFrame(bag_rows).to_csv(out / "bagging_diagnostics.tsv", sep="\t", index=False)
        pd.DataFrame([diag]).to_csv(out / "pu_diagnostics.tsv", sep="\t", index=False)
        save_scores(out / "outer_test_scores.parquet", df.iloc[idx_test], y_test, mean_score, experiment_id=exp_id, score_std=std_score)
        manifest = common_manifest(method="BaggingPU", experiment_id=exp_id, task_id=task_id, r=r, s=s, representation=rep, algorithm=algorithm, config_id=config_id, universe=universe, inner=inner)
        manifest.update({"n_bags": N_BAGS, "u_p_ratio": U_P_RATIO, "u_sampling": "without replacement within each bag", "bag_score_aggregation": "mean"})
        write_json(out / "manifest.json", manifest)
        (out / "warnings.log").write_text("", encoding="utf-8")
        print(f"DONE {exp_id}")
        return {"experiment_id": exp_id, "status": "DONE"}
    except Exception as exc:
        (out / "error.txt").write_text(traceback.format_exc(), encoding="utf-8")
        write_json(out / "manifest.json", {"stage": "11-pu", "status": "FAILED", "experiment_id": exp_id, "method": method, "error": str(exc)})
        print(f"FAILED {exp_id}: {exc}")
        return {"experiment_id": exp_id, "status": "FAILED", "error": str(exc)}


def run_elkan_noto(*, root, task_id, task_path, r, s, rep, algorithm, overwrite):
    method = "elkan_noto"
    config_id = BASE_CONFIGS[algorithm]
    exp_id = f"PU__{task_id}__R{r}__S{s}__{rep.spec.representation_id}__{method}__{algorithm}"
    out = output_dir(root, task_path, r, s, rep.spec.representation_id, method, algorithm)
    out.mkdir(parents=True, exist_ok=True)
    if done(out) and not overwrite:
        print(f"SKIP {exp_id} — DONE")
        return {"experiment_id": exp_id, "status": "DONE"}
    try:
        df, universe, inner = load_cell(root, task_path, r, s)
        y_all, _ = pu_target(df)
        idx_train = np.where(df["split"].eq("train").to_numpy())[0]
        idx_test = np.where(df["split"].eq("test").to_numpy())[0]
        outer_df = df.iloc[idx_train].reset_index(drop=True)
        y_train, y_test = y_all[idx_train], y_all[idx_test]
        rep.require_uids(df["sequence_uid"].astype(str).tolist())
        X_train = rep.get(outer_df["sequence_uid"].astype(str).tolist())
        X_test = rep.get(df.iloc[idx_test]["sequence_uid"].astype(str).tolist())
        cfg = get_config(algorithm, config_id)
        splitter = StratifiedGroupKFold(n_splits=EN_FOLDS, shuffle=True, random_state=SEED)
        groups = outer_df["cluster_id"].astype(str).to_numpy()
        dummy = np.zeros((len(outer_df), 1), dtype=np.uint8)
        held_positive_scores, fold_rows = [], []
        for fold_id, (fit_idx, hold_idx) in enumerate(splitter.split(dummy, y_train, groups)):
            hold_positive = hold_idx[y_train[hold_idx] == 1]
            if len(hold_positive) == 0:
                raise RuntimeError(f"Elkan–Noto fold {fold_id} contains no held-out P.")
            model = build_model(algorithm, cfg, task_type="binary", representation_family=rep.spec.family, n_classes=2)
            fit_model(model, X_train[fit_idx], y_train[fit_idx], algorithm)
            p_score = positive_scores(model, X_train[hold_positive])
            held_positive_scores.extend(p_score.tolist())
            fold_rows.append({"fold_id": fold_id, "n_fit": len(fit_idx), "n_holdout": len(hold_idx), "n_holdout_p": len(hold_positive), "mean_holdout_p_score": float(np.mean(p_score))})
        c = float(np.mean(held_positive_scores))
        if not np.isfinite(c) or c <= 0:
            raise RuntimeError(f"Invalid Elkan–Noto propensity c={c}")
        final_model = build_model(algorithm, cfg, task_type="binary", representation_family=rep.spec.family, n_classes=2)
        fit_model(final_model, X_train, y_train, algorithm)
        raw_score = positive_scores(final_model, X_test)
        unbounded = raw_score / c
        corrected = np.clip(unbounded, 0.0, 1.0)
        clipping_fraction = float(np.mean(unbounded > 1.0))
        diag = pu_diagnostics(y_test, corrected)
        diag.update({"method": method, "algorithm": algorithm, "base_config": config_id, "c": c, "crossfit_folds": EN_FOLDS, "corrected_score_clipping_fraction": clipping_fraction})
        pd.DataFrame(fold_rows).to_csv(out / "propensity_crossfit.tsv", sep="\t", index=False)
        pd.DataFrame([diag]).to_csv(out / "pu_diagnostics.tsv", sep="\t", index=False)
        save_scores(out / "outer_test_scores.parquet", df.iloc[idx_test], y_test, corrected, experiment_id=exp_id, raw_score=raw_score)
        manifest = common_manifest(method="Elkan-Noto", experiment_id=exp_id, task_id=task_id, r=r, s=s, representation=rep, algorithm=algorithm, config_id=config_id, universe=universe, inner=inner)
        manifest.update({"crossfit_folds": EN_FOLDS, "group_aware_crossfit": True, "c": c, "score_correction": "clip(raw_score / c, 0, 1)", "clipping_fraction": clipping_fraction})
        write_json(out / "manifest.json", manifest)
        (out / "warnings.log").write_text("", encoding="utf-8")
        print(f"DONE {exp_id} | c={c:.4f}")
        return {"experiment_id": exp_id, "status": "DONE"}
    except Exception as exc:
        (out / "error.txt").write_text(traceback.format_exc(), encoding="utf-8")
        write_json(out / "manifest.json", {"stage": "11-pu", "status": "FAILED", "experiment_id": exp_id, "method": method, "error": str(exc)})
        print(f"FAILED {exp_id}: {exc}")
        return {"experiment_id": exp_id, "status": "FAILED", "error": str(exc)}


def main():
    args = parse_args()
    root = project_root_from_script(__file__)
    task_path = TASKS[args.task]["path"]
    registry_path = root / "config/representation_registry/representation_registry.tsv"
    registry = load_representation_registry(registry_path, root)
    if args.representation not in registry:
        raise KeyError(f"Unknown representation {args.representation}. Available: {sorted(registry)}")
    rep = RepresentationStore(registry[args.representation])

    print("=" * 100)
    print("PhageLysClass — compact PU factorial benchmark")
    print("=" * 100)
    print(f"Task           : {args.task}")
    print(f"Representation : {args.representation}")
    print("Production     : all 16 R/S × LR02/XGB02 × naive/BaggingPU/Elkan-Noto")
    if args.smoke:
        print("MODE           : SMOKE ONLY (R90_S90 × LR02 × naive)")
    print()

    r_values = [90] if args.smoke else R_VALUES
    s_values = [90] if args.smoke else S_VALUES
    algorithms = ["LR"] if args.smoke else BASE_ALGORITHMS
    methods = ["naive"] if args.smoke else METHODS

    results = []
    for r in r_values:
        for s in s_values:
            for algorithm in algorithms:
                for method in methods:
                    kwargs = dict(root=root, task_id=args.task, task_path=task_path, r=r, s=s, rep=rep, algorithm=algorithm, overwrite=args.overwrite)
                    if method == "naive":
                        result = run_naive(**kwargs)
                    elif method == "baggingpu":
                        result = run_baggingpu(**kwargs)
                    else:
                        result = run_elkan_noto(**kwargs)
                    results.append(result)

    summary = pd.DataFrame(results)
    out_path = root / "results/training_pu" / task_path / f"_run_summary__{args.representation}.tsv"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(out_path, sep="\t", index=False)

    print()
    print("=" * 100)
    print(f"Run summary: {summary['status'].value_counts().to_dict()}")
    print("=" * 100)

    if not args.smoke:
        expected = 16 * len(BASE_ALGORITHMS) * len(METHODS)
        if len(summary) != expected:
            raise SystemExit(f"Expected {expected} PU experiment families, observed {len(summary)}")


if __name__ == "__main__":
    main()

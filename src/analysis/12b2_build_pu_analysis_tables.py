#!/usr/bin/env python3
"""
12b2_build_pu_analysis_tables.py

Comprehensive analysis-table builder for the final PhageLysClass PU benchmark.

Prerequisite
------------
Run and PASS:

    python src/analysis/12a2_audit_pu_training.py

This script intentionally uses the QC-produced completed registry, including
the discovered experiment paths, so it does not assume any particular
results/training_pu directory layout.

Final PU design
---------------
2 PU tasks
× 16 R/S cells
× 4 representations
× 3 PU methods
× 2 algorithms
= 768 experiment families

Important interpretation
------------------------
U is unlabeled, not confirmed negative.

Therefore MCC/AP/AUROC or related metrics calculated against observed P/U
labels are treated here as OBSERVED-LABEL DIAGNOSTICS. They are useful for
controlled comparisons among methods but must not be described as unbiased
biological classification performance.

Outputs
-------
results/analysis_pu/stage12b2/

Core aggregation:
    pu_master_results.tsv
    pu_score_distribution_summary.tsv
    pu_bagging_stability_summary.tsv
    pu_elkan_noto_summary.tsv

Model comparisons on identical native test sets:
    pu_method_pairwise_deltas.tsv
    pu_representation_pairwise_deltas.tsv
    pu_algorithm_pairwise_deltas.tsv

Complete R×S analysis:
    pu_rs_native_surfaces.tsv
    pu_rs_all_pairwise_cell_deltas.tsv
    pu_rs_redundancy_deltas_fixed_S.tsv
    pu_rs_split_deltas_fixed_R.tsv
    pu_rs_all_interaction_rectangles.tsv
    pu_rs_extreme_summary_contrasts.tsv

Score-level agreement on shared U entities:
    pu_method_score_concordance.tsv
    pu_representation_score_concordance.tsv
    pu_algorithm_score_concordance.tsv

Robustness / ranking:
    pu_robustness_summary.tsv
    pu_diagnostic_winner_per_cell.tsv
    pu_diagnostic_winner_counts.tsv

QC:
    pu_cross_representation_test_audit.tsv
    stage12b2_summary.tsv
    stage12b2_manifest.json
"""

from __future__ import annotations

import hashlib
import json
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd


R_LEVELS = [90, 70, 50, 30]
S_LEVELS = [90, 70, 50, 30]

EXPECTED_TOTAL = 768
EXPECTED_REPRESENTATIONS = [
    "roxy_48",
    "onehot_1024",
    "esm2_t6_8m",
    "esm2_t30_150m",
]
EXPECTED_METHODS = [
    "naive",
    "baggingpu",
    "elkan_noto",
]
EXPECTED_ALGORITHMS = [
    "LR",
    "XGB",
]


def project_root() -> Path:
    return Path(__file__).resolve().parents[2]


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def safe_numeric_metrics(df: pd.DataFrame) -> list[str]:
    """
    Detect PU observed-label diagnostic metrics robustly.

    The PU trainer may prefix diagnostics as, e.g.:
        observed_label_mcc
        observed_label_ap
        observed_label_roc_auc

    Earlier versions of this script only recognized names such as
    ``observed_mcc`` or ``mcc``, which could incorrectly return an empty
    metric list even when valid diagnostics were present.

    We therefore identify numeric columns whose names correspond to standard
    classification diagnostics, while excluding runtime/count/PU nuisance
    fields.
    """
    metric_tokens = (
        "mcc",
        "average_precision",
        "roc_auc",
        "auroc",
        "balanced_accuracy",
        "macro_f1",
        "weighted_f1",
        "precision",
        "recall",
        "specificity",
        "f1",
        "brier",
        "log_loss",
    )

    exact_ap_names = {
        "ap",
        "observed_ap",
        "observed_label_ap",
    }

    metrics = []

    for column in df.columns:
        if not pd.api.types.is_numeric_dtype(df[column]):
            continue

        name = str(column).lower()

        # Explicitly exclude non-performance numeric fields.
        if (
            name in {"r", "s", "c", "n_bags", "u_p_ratio", "crossfit_folds"}
            or name.startswith("manifest__")
            or "seconds" in name
            or "runtime" in name
            or name.startswith("n_")
            or name.endswith("_count")
            or "clipping_fraction" in name
            or "score_std" in name
        ):
            continue

        is_metric = (
            name in exact_ap_names
            or any(token in name for token in metric_tokens)
        )

        if is_metric:
            metrics.append(column)

    # Stable, scientifically useful ordering.
    priority = [
        "observed_label_mcc",
        "observed_mcc",
        "mcc",
        "observed_label_ap",
        "observed_ap",
        "ap",
        "average_precision",
        "observed_label_roc_auc",
        "observed_roc_auc",
        "roc_auc",
        "auroc",
        "observed_label_balanced_accuracy",
        "balanced_accuracy",
        "observed_label_f1",
        "observed_f1",
        "f1",
        "precision",
        "recall",
        "specificity",
        "brier",
        "log_loss",
    ]

    rank = {name: i for i, name in enumerate(priority)}

    return sorted(
        dict.fromkeys(metrics),
        key=lambda c: (rank.get(str(c).lower(), len(priority)), str(c)),
    )


def primary_metric(metrics: list[str]) -> str:
    """
    Prefer MCC as the primary observed-label diagnostic, regardless of
    whether the trainer used an ``observed_label_`` prefix.
    """
    preferred = [
        "observed_label_mcc",
        "observed_mcc",
        "mcc",
        "observed_label_ap",
        "observed_ap",
        "ap",
        "average_precision",
        "observed_label_roc_auc",
        "observed_roc_auc",
        "roc_auc",
        "auroc",
    ]

    lower_to_original = {
        str(metric).lower(): metric
        for metric in metrics
    }

    for candidate in preferred:
        if candidate in lower_to_original:
            return lower_to_original[candidate]

    # Generic fallback: any detected MCC-like metric, then AP-like, then AUC.
    for token in ["mcc", "average_precision", "_ap", "roc_auc", "auroc"]:
        for metric in metrics:
            if token in str(metric).lower():
                return metric

    if metrics:
        return metrics[0]

    raise RuntimeError(
        "No PU diagnostic metric was detected. "
        "Inspect results/analysis_pu/stage12b2/pu_master_results.tsv "
        "or one pu_diagnostics.tsv to confirm the output schema."
    )


def score_signature(path: Path) -> tuple[str, int]:
    df = pd.read_parquet(
        path,
        columns=["entity_id", "sequence_uid", "observed_role"],
    )
    if df["entity_id"].duplicated().any():
        raise RuntimeError(f"{path}: duplicate entity_id")

    df = df.sort_values(
        ["entity_id", "sequence_uid"],
        kind="stable",
    )

    h = hashlib.sha256()
    for row in df.itertuples(index=False):
        h.update(
            (
                f"{row.entity_id}\t"
                f"{row.sequence_uid}\t"
                f"{row.observed_role}\n"
            ).encode("utf-8")
        )
    return h.hexdigest(), len(df)


def lookup_group(
    group: pd.DataFrame,
    metric: str,
) -> dict[tuple[int, int], float]:
    out = {}
    for row in group.itertuples(index=False):
        value = getattr(row, metric)
        if pd.notna(value):
            out[(int(row.R), int(row.S))] = float(value)
    return out


def pairwise_model_deltas(
    df: pd.DataFrame,
    *,
    compare_col: str,
    matched_cols: list[str],
    metrics: list[str],
) -> pd.DataFrame:
    """
    Compare methods / representations / algorithms while all matched columns
    are identical. These are metric-level comparisons on the same native test.
    """
    rows = []
    values = sorted(
        df[compare_col].dropna().astype(str).unique()
    )

    for a, b in combinations(values, 2):
        left = df.loc[
            df[compare_col].astype(str).eq(a),
            matched_cols + metrics,
        ].copy()

        right = df.loc[
            df[compare_col].astype(str).eq(b),
            matched_cols + metrics,
        ].copy()

        left = left.rename(
            columns={m: f"{m}__A" for m in metrics}
        )
        right = right.rename(
            columns={m: f"{m}__B" for m in metrics}
        )

        merged = left.merge(
            right,
            on=matched_cols,
            how="inner",
            validate="one_to_one",
        )

        for rec in merged.to_dict("records"):
            out = {
                **{c: rec[c] for c in matched_cols},
                f"{compare_col}_A": a,
                f"{compare_col}_B": b,
            }

            for metric in metrics:
                va = rec.get(f"{metric}__A")
                vb = rec.get(f"{metric}__B")

                out[f"{metric}_A"] = va
                out[f"{metric}_B"] = vb
                out[f"delta_{metric}_B_minus_A"] = (
                    float(vb) - float(va)
                    if pd.notna(va) and pd.notna(vb)
                    else np.nan
                )

            rows.append(out)

    return pd.DataFrame(rows)


def all_cell_pair_deltas(
    df: pd.DataFrame,
    metric: str,
) -> pd.DataFrame:
    """
    All 120 unordered comparisons among the 16 R×S cells for each:
        task × representation × method × algorithm
    """
    rows = []

    group_cols = [
        "task",
        "representation",
        "method",
        "algorithm",
    ]

    for keys, group in df.groupby(group_cols):
        task, rep, method, alg = keys
        lookup = lookup_group(group, metric)

        cells = [
            (r, s)
            for r in R_LEVELS
            for s in S_LEVELS
            if (r, s) in lookup
        ]

        for c1, c2 in combinations(cells, 2):
            r1, s1 = c1
            r2, s2 = c2

            if s1 == s2:
                a, b = (
                    (c1, c2)
                    if r1 > r2
                    else (c2, c1)
                )
                relation = "same_S_redundancy"

            elif r1 == r2:
                a, b = (
                    (c1, c2)
                    if s1 > s2
                    else (c2, c1)
                )
                relation = "same_R_split"

            elif r2 < r1 and s2 < s1:
                a, b = c1, c2
                relation = "both_factors_stricter"

            elif r1 < r2 and s1 < s2:
                a, b = c2, c1
                relation = "both_factors_stricter"

            else:
                a, b = c1, c2
                relation = "tradeoff_incomparable"

            ra, sa = a
            rb, sb = b
            va = lookup[a]
            vb = lookup[b]

            rows.append(
                {
                    "task": task,
                    "representation": rep,
                    "method": method,
                    "algorithm": alg,
                    "metric": metric,
                    "relation": relation,
                    "R_A": ra,
                    "S_A": sa,
                    "R_B": rb,
                    "S_B": sb,
                    "value_A": va,
                    "value_B": vb,
                    "delta_B_minus_A": vb - va,
                    "native_test_note": (
                        "PU observed-label native-test contrast; "
                        "not necessarily entity-level paired across R/S."
                    ),
                }
            )

    return pd.DataFrame(rows)


def redundancy_deltas_fixed_s(
    df: pd.DataFrame,
    metric: str,
) -> pd.DataFrame:
    rows = []

    for keys, group in df.groupby(
        ["task", "representation", "method", "algorithm"]
    ):
        task, rep, method, alg = keys
        lookup = lookup_group(group, metric)

        for s in S_LEVELS:
            for r_loose, r_strict in combinations(R_LEVELS, 2):
                if (
                    (r_loose, s) not in lookup
                    or (r_strict, s) not in lookup
                ):
                    continue

                va = lookup[(r_loose, s)]
                vb = lookup[(r_strict, s)]

                rows.append(
                    {
                        "task": task,
                        "representation": rep,
                        "method": method,
                        "algorithm": alg,
                        "metric": metric,
                        "S_fixed": s,
                        "R_from": r_loose,
                        "R_to": r_strict,
                        "threshold_step": (
                            f"R{r_loose}->R{r_strict}"
                        ),
                        "value_from": va,
                        "value_to": vb,
                        "delta": vb - va,
                    }
                )

    return pd.DataFrame(rows)


def split_deltas_fixed_r(
    df: pd.DataFrame,
    metric: str,
) -> pd.DataFrame:
    rows = []

    for keys, group in df.groupby(
        ["task", "representation", "method", "algorithm"]
    ):
        task, rep, method, alg = keys
        lookup = lookup_group(group, metric)

        for r in R_LEVELS:
            for s_loose, s_strict in combinations(S_LEVELS, 2):
                if (
                    (r, s_loose) not in lookup
                    or (r, s_strict) not in lookup
                ):
                    continue

                va = lookup[(r, s_loose)]
                vb = lookup[(r, s_strict)]

                rows.append(
                    {
                        "task": task,
                        "representation": rep,
                        "method": method,
                        "algorithm": alg,
                        "metric": metric,
                        "R_fixed": r,
                        "S_from": s_loose,
                        "S_to": s_strict,
                        "threshold_step": (
                            f"S{s_loose}->S{s_strict}"
                        ),
                        "value_from": va,
                        "value_to": vb,
                        "delta": vb - va,
                    }
                )

    return pd.DataFrame(rows)


def all_interaction_rectangles(
    df: pd.DataFrame,
    metric: str,
) -> pd.DataFrame:
    rows = []

    for keys, group in df.groupby(
        ["task", "representation", "method", "algorithm"]
    ):
        task, rep, method, alg = keys
        lookup = lookup_group(group, metric)

        for r_loose, r_strict in combinations(R_LEVELS, 2):
            for s_loose, s_strict in combinations(S_LEVELS, 2):
                corners = [
                    (r_loose, s_loose),
                    (r_strict, s_loose),
                    (r_loose, s_strict),
                    (r_strict, s_strict),
                ]

                if not all(c in lookup for c in corners):
                    continue

                m00 = lookup[(r_loose, s_loose)]
                m10 = lookup[(r_strict, s_loose)]
                m01 = lookup[(r_loose, s_strict)]
                m11 = lookup[(r_strict, s_strict)]

                delta_r_loose_s = m10 - m00
                delta_r_strict_s = m11 - m01

                delta_s_loose_r = m01 - m00
                delta_s_strict_r = m11 - m10

                interaction = (
                    delta_r_strict_s
                    - delta_r_loose_s
                )

                rows.append(
                    {
                        "task": task,
                        "representation": rep,
                        "method": method,
                        "algorithm": alg,
                        "metric": metric,
                        "R_loose": r_loose,
                        "R_strict": r_strict,
                        "S_loose": s_loose,
                        "S_strict": s_strict,
                        "M_Rloose_Sloose": m00,
                        "M_Rstrict_Sloose": m10,
                        "M_Rloose_Sstrict": m01,
                        "M_Rstrict_Sstrict": m11,
                        "delta_R_at_Sloose": delta_r_loose_s,
                        "delta_R_at_Sstrict": delta_r_strict_s,
                        "delta_S_at_Rloose": delta_s_loose_r,
                        "delta_S_at_Rstrict": delta_s_strict_r,
                        "interaction": interaction,
                        "interaction_equivalence_check": (
                            delta_s_strict_r
                            - delta_s_loose_r
                        ),
                    }
                )

    return pd.DataFrame(rows)


def extreme_summary(
    df: pd.DataFrame,
    metric: str,
) -> pd.DataFrame:
    rows = []

    for keys, group in df.groupby(
        ["task", "representation", "method", "algorithm"]
    ):
        task, rep, method, alg = keys
        lookup = lookup_group(group, metric)

        for s in S_LEVELS:
            if (90, s) in lookup and (30, s) in lookup:
                rows.append(
                    {
                        "task": task,
                        "representation": rep,
                        "method": method,
                        "algorithm": alg,
                        "metric": metric,
                        "contrast": "R30_minus_R90",
                        "fixed_factor": "S",
                        "fixed_value": s,
                        "value": (
                            lookup[(30, s)]
                            - lookup[(90, s)]
                        ),
                    }
                )

        for r in R_LEVELS:
            if (r, 90) in lookup and (r, 30) in lookup:
                rows.append(
                    {
                        "task": task,
                        "representation": rep,
                        "method": method,
                        "algorithm": alg,
                        "metric": metric,
                        "contrast": "S30_minus_S90",
                        "fixed_factor": "R",
                        "fixed_value": r,
                        "value": (
                            lookup[(r, 30)]
                            - lookup[(r, 90)]
                        ),
                    }
                )

        corners = [(90, 90), (30, 90), (90, 30), (30, 30)]
        if all(c in lookup for c in corners):
            interaction = (
                (
                    lookup[(30, 30)]
                    - lookup[(90, 30)]
                )
                -
                (
                    lookup[(30, 90)]
                    - lookup[(90, 90)]
                )
            )

            rows.append(
                {
                    "task": task,
                    "representation": rep,
                    "method": method,
                    "algorithm": alg,
                    "metric": metric,
                    "contrast": "extreme_RxS_interaction",
                    "fixed_factor": "",
                    "fixed_value": np.nan,
                    "value": interaction,
                }
            )

    return pd.DataFrame(rows)


def top_fraction_ids(
    df: pd.DataFrame,
    fraction: float,
) -> set:
    if not len(df):
        return set()

    n = max(
        1,
        int(np.ceil(len(df) * fraction)),
    )

    return set(
        df.nlargest(
            n,
            "score_positive",
        )["entity_id"].astype(str)
    )


def jaccard(a: set, b: set) -> float:
    union = a | b
    if not union:
        return np.nan
    return len(a & b) / len(union)


def score_concordance_pair(
    path_a: Path,
    path_b: Path,
) -> dict:
    """
    Compare two score files on shared U entities only.
    """
    a = pd.read_parquet(
        path_a,
        columns=[
            "entity_id",
            "observed_role",
            "score_positive",
        ],
    )

    b = pd.read_parquet(
        path_b,
        columns=[
            "entity_id",
            "observed_role",
            "score_positive",
        ],
    )

    a = a[a["observed_role"].eq("U")].copy()
    b = b[b["observed_role"].eq("U")].copy()

    a = a[
        ["entity_id", "score_positive"]
    ].rename(
        columns={
            "score_positive": "score_A"
        }
    )

    b = b[
        ["entity_id", "score_positive"]
    ].rename(
        columns={
            "score_positive": "score_B"
        }
    )

    merged = a.merge(
        b,
        on="entity_id",
        how="inner",
        validate="one_to_one",
    )

    if not len(merged):
        return {
            "n_shared_U": 0,
            "pearson": np.nan,
            "spearman": np.nan,
            "mae": np.nan,
            "top_1pct_jaccard": np.nan,
            "top_5pct_jaccard": np.nan,
        }

    pearson = merged[
        ["score_A", "score_B"]
    ].corr(method="pearson").iloc[0, 1]

    spearman = merged[
        ["score_A", "score_B"]
    ].corr(method="spearman").iloc[0, 1]

    mae = float(
        np.mean(
            np.abs(
                merged["score_A"]
                - merged["score_B"]
            )
        )
    )

    ma = merged[
        ["entity_id", "score_A"]
    ].rename(
        columns={"score_A": "score_positive"}
    )

    mb = merged[
        ["entity_id", "score_B"]
    ].rename(
        columns={"score_B": "score_positive"}
    )

    return {
        "n_shared_U": len(merged),
        "pearson": pearson,
        "spearman": spearman,
        "mae": mae,
        "top_1pct_jaccard": jaccard(
            top_fraction_ids(ma, 0.01),
            top_fraction_ids(mb, 0.01),
        ),
        "top_5pct_jaccard": jaccard(
            top_fraction_ids(ma, 0.05),
            top_fraction_ids(mb, 0.05),
        ),
    }


def build_score_concordance(
    registry: pd.DataFrame,
    *,
    compare_col: str,
    matched_cols: list[str],
) -> pd.DataFrame:
    rows = []

    values = sorted(
        registry[compare_col]
        .dropna()
        .astype(str)
        .unique()
    )

    for a, b in combinations(values, 2):
        left = registry.loc[
            registry[compare_col].astype(str).eq(a),
            matched_cols + ["path"],
        ].rename(
            columns={"path": "path_A"}
        )

        right = registry.loc[
            registry[compare_col].astype(str).eq(b),
            matched_cols + ["path"],
        ].rename(
            columns={"path": "path_B"}
        )

        merged = left.merge(
            right,
            on=matched_cols,
            how="inner",
            validate="one_to_one",
        )

        for rec in merged.to_dict("records"):
            stats = score_concordance_pair(
                Path(rec["path_A"])
                / "outer_test_scores.parquet",
                Path(rec["path_B"])
                / "outer_test_scores.parquet",
            )

            rows.append(
                {
                    **{
                        c: rec[c]
                        for c in matched_cols
                    },
                    f"{compare_col}_A": a,
                    f"{compare_col}_B": b,
                    **stats,
                }
            )

    return pd.DataFrame(rows)


def main():
    root = project_root()

    qc_dir = (
        root
        / "results/analysis_pu/stage12a_qc"
    )

    qc_manifest_path = (
        qc_dir / "stage12a_pu_manifest.json"
    )

    registry_path = (
        qc_dir / "pu_completed_registry.tsv"
    )

    if not qc_manifest_path.is_file():
        raise FileNotFoundError(
            f"{qc_manifest_path}\nRun Stage 12A2 PU QC first."
        )

    qc_manifest = read_json(qc_manifest_path)

    if qc_manifest.get("status") != "PASS":
        raise RuntimeError(
            "Stage 12A2 PU QC did not PASS. "
            "Do not build PU analysis tables yet."
        )

    if not registry_path.is_file():
        raise FileNotFoundError(registry_path)

    registry = pd.read_csv(
        registry_path,
        sep="\t",
        low_memory=False,
    )

    if len(registry) != EXPECTED_TOTAL:
        raise RuntimeError(
            f"Expected {EXPECTED_TOTAL} completed PU families, "
            f"found {len(registry)}."
        )

    required_registry = {
        "experiment_id",
        "task",
        "R",
        "S",
        "representation",
        "method",
        "algorithm",
        "path",
    }

    missing = (
        required_registry
        - set(registry.columns)
    )

    if missing:
        raise KeyError(
            f"PU completed registry missing columns: {sorted(missing)}"
        )

    out = (
        root
        / "results/analysis_pu/stage12b2"
    )

    out.mkdir(
        parents=True,
        exist_ok=True,
    )

    # --------------------------------------------------
    # A. Aggregate all PU diagnostics + manifests
    # --------------------------------------------------

    master_rows = []
    score_summary_rows = []
    bagging_rows = []
    en_rows = []
    signature_rows = []

    for row in registry.itertuples(index=False):
        d = Path(row.path)

        diagnostics = pd.read_csv(
            d / "pu_diagnostics.tsv",
            sep="\t",
            low_memory=False,
        )

        if len(diagnostics) != 1:
            raise RuntimeError(
                f"{d}/pu_diagnostics.tsv expected 1 row."
            )

        manifest = read_json(
            d / "manifest.json"
        )

        diag = diagnostics.iloc[0].to_dict()

        master_row = {
            "experiment_id": row.experiment_id,
            "task": row.task,
            "R": int(row.R),
            "S": int(row.S),
            "representation": row.representation,
            "method": row.method,
            "algorithm": row.algorithm,
            **diag,
        }

        # Useful method-specific manifest fields.
        for key in [
            "c",
            "clipping_fraction",
            "n_bags",
            "u_p_ratio",
            "crossfit_folds",
        ]:
            if key in manifest:
                master_row[
                    f"manifest__{key}"
                ] = manifest[key]

        master_rows.append(master_row)

        scores = pd.read_parquet(
            d / "outer_test_scores.parquet"
        )

        for role, g in scores.groupby(
            "observed_role"
        ):
            s = g["score_positive"].astype(float)

            score_summary_rows.append(
                {
                    "experiment_id": row.experiment_id,
                    "task": row.task,
                    "R": int(row.R),
                    "S": int(row.S),
                    "representation": row.representation,
                    "method": row.method,
                    "algorithm": row.algorithm,
                    "observed_role": role,
                    "n": len(g),
                    "score_mean": s.mean(),
                    "score_median": s.median(),
                    "score_std": s.std(ddof=0),
                    "score_q05": s.quantile(0.05),
                    "score_q25": s.quantile(0.25),
                    "score_q75": s.quantile(0.75),
                    "score_q95": s.quantile(0.95),
                    "score_min": s.min(),
                    "score_max": s.max(),
                }
            )

        sig, n_test = score_signature(
            d / "outer_test_scores.parquet"
        )

        signature_rows.append(
            {
                "experiment_id": row.experiment_id,
                "task": row.task,
                "R": int(row.R),
                "S": int(row.S),
                "representation": row.representation,
                "method": row.method,
                "algorithm": row.algorithm,
                "test_signature": sig,
                "n_test": n_test,
            }
        )

        if row.method == "baggingpu":
            bags = pd.read_csv(
                d / "bagging_diagnostics.tsv",
                sep="\t",
                low_memory=False,
            )

            bagging_rows.append(
                {
                    "experiment_id": row.experiment_id,
                    "task": row.task,
                    "R": int(row.R),
                    "S": int(row.S),
                    "representation": row.representation,
                    "algorithm": row.algorithm,
                    "n_bags": len(bags),
                    "median_bag_seconds": (
                        bags[
                            "fit_plus_predict_seconds"
                        ].median()
                        if "fit_plus_predict_seconds"
                        in bags.columns
                        else np.nan
                    ),
                    "total_bag_seconds": (
                        bags[
                            "fit_plus_predict_seconds"
                        ].sum()
                        if "fit_plus_predict_seconds"
                        in bags.columns
                        else np.nan
                    ),
                    "median_test_score_std": (
                        diag.get(
                            "median_test_score_std",
                            np.nan,
                        )
                    ),
                    "mean_test_score_std": (
                        diag.get(
                            "mean_test_score_std",
                            np.nan,
                        )
                    ),
                }
            )

        elif row.method == "elkan_noto":
            folds = pd.read_csv(
                d / "propensity_crossfit.tsv",
                sep="\t",
                low_memory=False,
            )

            en_rows.append(
                {
                    "experiment_id": row.experiment_id,
                    "task": row.task,
                    "R": int(row.R),
                    "S": int(row.S),
                    "representation": row.representation,
                    "algorithm": row.algorithm,
                    "c": manifest.get(
                        "c",
                        diag.get("c", np.nan),
                    ),
                    "clipping_fraction": manifest.get(
                        "clipping_fraction",
                        diag.get(
                            "corrected_score_clipping_fraction",
                            np.nan,
                        ),
                    ),
                    "n_crossfit_folds": len(folds),
                    "mean_fold_positive_score": (
                        folds[
                            "mean_holdout_p_score"
                        ].mean()
                        if "mean_holdout_p_score"
                        in folds.columns
                        else np.nan
                    ),
                    "min_fold_positive_score": (
                        folds[
                            "mean_holdout_p_score"
                        ].min()
                        if "mean_holdout_p_score"
                        in folds.columns
                        else np.nan
                    ),
                    "max_fold_positive_score": (
                        folds[
                            "mean_holdout_p_score"
                        ].max()
                        if "mean_holdout_p_score"
                        in folds.columns
                        else np.nan
                    ),
                }
            )

    master = pd.DataFrame(master_rows)
    score_summary = pd.DataFrame(
        score_summary_rows
    )
    bagging = pd.DataFrame(bagging_rows)
    en = pd.DataFrame(en_rows)
    signatures = pd.DataFrame(
        signature_rows
    )

    master.to_csv(
        out / "pu_master_results.tsv",
        sep="\t",
        index=False,
    )

    score_summary.to_csv(
        out / "pu_score_distribution_summary.tsv",
        sep="\t",
        index=False,
    )

    bagging.to_csv(
        out / "pu_bagging_stability_summary.tsv",
        sep="\t",
        index=False,
    )

    en.to_csv(
        out / "pu_elkan_noto_summary.tsv",
        sep="\t",
        index=False,
    )

    # --------------------------------------------------
    # B. Cross-representation test-set QC
    # --------------------------------------------------

    cross_rep_audit = (
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
                "n_test",
                "nunique",
            ),
        )
    )

    cross_rep_audit[
        "same_test_set_across_all_PU_models"
    ] = (
        cross_rep_audit[
            "n_experiments"
        ].eq(
            len(EXPECTED_REPRESENTATIONS)
            * len(EXPECTED_METHODS)
            * len(EXPECTED_ALGORITHMS)
        )
        & cross_rep_audit[
            "unique_test_signatures"
        ].eq(1)
        & cross_rep_audit[
            "unique_test_sizes"
        ].eq(1)
    )

    cross_rep_audit.to_csv(
        out
        / "pu_cross_representation_test_audit.tsv",
        sep="\t",
        index=False,
    )

    cross_rep_failures = int(
        (
            ~cross_rep_audit[
                "same_test_set_across_all_PU_models"
            ]
        ).sum()
    )

    # --------------------------------------------------
    # C. Identify diagnostic metrics
    # --------------------------------------------------

    metrics = safe_numeric_metrics(master)
    primary = primary_metric(metrics)

    # --------------------------------------------------
    # D. Pairwise method / representation / algorithm
    # --------------------------------------------------

    method_deltas = pairwise_model_deltas(
        master,
        compare_col="method",
        matched_cols=[
            "task",
            "R",
            "S",
            "representation",
            "algorithm",
        ],
        metrics=metrics,
    )

    representation_deltas = (
        pairwise_model_deltas(
            master,
            compare_col="representation",
            matched_cols=[
                "task",
                "R",
                "S",
                "method",
                "algorithm",
            ],
            metrics=metrics,
        )
    )

    algorithm_deltas = (
        pairwise_model_deltas(
            master,
            compare_col="algorithm",
            matched_cols=[
                "task",
                "R",
                "S",
                "representation",
                "method",
            ],
            metrics=metrics,
        )
    )

    method_deltas.to_csv(
        out / "pu_method_pairwise_deltas.tsv",
        sep="\t",
        index=False,
    )

    representation_deltas.to_csv(
        out
        / "pu_representation_pairwise_deltas.tsv",
        sep="\t",
        index=False,
    )

    algorithm_deltas.to_csv(
        out
        / "pu_algorithm_pairwise_deltas.tsv",
        sep="\t",
        index=False,
    )

    # --------------------------------------------------
    # E. Complete R×S factorial contrasts
    # --------------------------------------------------

    surface_cols = [
        "experiment_id",
        "task",
        "R",
        "S",
        "representation",
        "method",
        "algorithm",
    ] + metrics

    master[
        surface_cols
    ].to_csv(
        out / "pu_rs_native_surfaces.tsv",
        sep="\t",
        index=False,
    )

    all_pair_parts = []
    r_parts = []
    s_parts = []
    interaction_parts = []
    extreme_parts = []

    for metric in metrics:
        all_pair_parts.append(
            all_cell_pair_deltas(
                master,
                metric,
            )
        )

        r_parts.append(
            redundancy_deltas_fixed_s(
                master,
                metric,
            )
        )

        s_parts.append(
            split_deltas_fixed_r(
                master,
                metric,
            )
        )

        interaction_parts.append(
            all_interaction_rectangles(
                master,
                metric,
            )
        )

        extreme_parts.append(
            extreme_summary(
                master,
                metric,
            )
        )

    all_pairs = pd.concat(
        all_pair_parts,
        ignore_index=True,
    )

    r_deltas = pd.concat(
        r_parts,
        ignore_index=True,
    )

    s_deltas = pd.concat(
        s_parts,
        ignore_index=True,
    )

    interactions = pd.concat(
        interaction_parts,
        ignore_index=True,
    )

    extremes = pd.concat(
        extreme_parts,
        ignore_index=True,
    )

    all_pairs.to_csv(
        out
        / "pu_rs_all_pairwise_cell_deltas.tsv",
        sep="\t",
        index=False,
    )

    r_deltas.to_csv(
        out
        / "pu_rs_redundancy_deltas_fixed_S.tsv",
        sep="\t",
        index=False,
    )

    s_deltas.to_csv(
        out
        / "pu_rs_split_deltas_fixed_R.tsv",
        sep="\t",
        index=False,
    )

    interactions.to_csv(
        out
        / "pu_rs_all_interaction_rectangles.tsv",
        sep="\t",
        index=False,
    )

    extremes.to_csv(
        out
        / "pu_rs_extreme_summary_contrasts.tsv",
        sep="\t",
        index=False,
    )

    # --------------------------------------------------
    # F. Score-level concordance on shared U entities
    # --------------------------------------------------

    method_concordance = build_score_concordance(
        registry,
        compare_col="method",
        matched_cols=[
            "task",
            "R",
            "S",
            "representation",
            "algorithm",
        ],
    )

    representation_concordance = (
        build_score_concordance(
            registry,
            compare_col="representation",
            matched_cols=[
                "task",
                "R",
                "S",
                "method",
                "algorithm",
            ],
        )
    )

    algorithm_concordance = (
        build_score_concordance(
            registry,
            compare_col="algorithm",
            matched_cols=[
                "task",
                "R",
                "S",
                "representation",
                "method",
            ],
        )
    )

    method_concordance.to_csv(
        out
        / "pu_method_score_concordance.tsv",
        sep="\t",
        index=False,
    )

    representation_concordance.to_csv(
        out
        / "pu_representation_score_concordance.tsv",
        sep="\t",
        index=False,
    )

    algorithm_concordance.to_csv(
        out
        / "pu_algorithm_score_concordance.tsv",
        sep="\t",
        index=False,
    )

    # --------------------------------------------------
    # G. Robustness over the complete 16-cell R×S space
    # --------------------------------------------------

    robustness_rows = []

    for keys, group in master.groupby(
        [
            "task",
            "representation",
            "method",
            "algorithm",
        ]
    ):
        task, rep, method, alg = keys

        row = {
            "task": task,
            "representation": rep,
            "method": method,
            "algorithm": alg,
            "n_rs_cells": len(group),
        }

        for metric in metrics:
            if group[metric].notna().any():
                row[
                    f"{metric}_mean"
                ] = group[metric].mean()

                row[
                    f"{metric}_median"
                ] = group[metric].median()

                row[
                    f"{metric}_min"
                ] = group[metric].min()

                row[
                    f"{metric}_max"
                ] = group[metric].max()

                row[
                    f"{metric}_range"
                ] = (
                    group[metric].max()
                    - group[metric].min()
                )

        robustness_rows.append(row)

    robustness = pd.DataFrame(
        robustness_rows
    )

    robustness.to_csv(
        out / "pu_robustness_summary.tsv",
        sep="\t",
        index=False,
    )

    # --------------------------------------------------
    # H. Diagnostic winners per native cell
    # --------------------------------------------------

    winners = (
        master.sort_values(
            [
                "task",
                "R",
                "S",
                primary,
            ],
            ascending=[
                True,
                True,
                True,
                False,
            ],
            kind="stable",
        )
        .groupby(
            ["task", "R", "S"],
            as_index=False,
        )
        .first()
    )

    winner_cols = [
        "task",
        "R",
        "S",
        "representation",
        "method",
        "algorithm",
        primary,
    ]

    winners[
        winner_cols
    ].to_csv(
        out
        / "pu_diagnostic_winner_per_cell.tsv",
        sep="\t",
        index=False,
    )

    winner_counts = (
        winners.groupby(
            [
                "task",
                "representation",
                "method",
                "algorithm",
            ],
            as_index=False,
        )
        .size()
        .rename(
            columns={"size": "n_cells_won"}
        )
        .sort_values(
            ["task", "n_cells_won"],
            ascending=[True, False],
        )
    )

    winner_counts.to_csv(
        out
        / "pu_diagnostic_winner_counts.tsv",
        sep="\t",
        index=False,
    )

    # --------------------------------------------------
    # I. Summary + manifest
    # --------------------------------------------------

    summary = pd.DataFrame(
        [
            {
                "item": "experiment_families",
                "value": len(master),
            },
            {
                "item": "diagnostic_metrics",
                "value": ",".join(metrics),
            },
            {
                "item": "primary_diagnostic_metric",
                "value": primary,
            },
            {
                "item": "method_pairwise_rows",
                "value": len(method_deltas),
            },
            {
                "item": "representation_pairwise_rows",
                "value": len(
                    representation_deltas
                ),
            },
            {
                "item": "algorithm_pairwise_rows",
                "value": len(
                    algorithm_deltas
                ),
            },
            {
                "item": "all_RxS_pairwise_rows",
                "value": len(all_pairs),
            },
            {
                "item": "redundancy_fixed_S_rows",
                "value": len(r_deltas),
            },
            {
                "item": "split_fixed_R_rows",
                "value": len(s_deltas),
            },
            {
                "item": "interaction_rectangle_rows",
                "value": len(interactions),
            },
            {
                "item": "method_score_concordance_rows",
                "value": len(
                    method_concordance
                ),
            },
            {
                "item": (
                    "representation_score_concordance_rows"
                ),
                "value": len(
                    representation_concordance
                ),
            },
            {
                "item": "cross_representation_test_failures",
                "value": cross_rep_failures,
            },
        ]
    )

    summary.to_csv(
        out / "stage12b2_summary.tsv",
        sep="\t",
        index=False,
    )

    status = (
        "PASS"
        if (
            len(master) == EXPECTED_TOTAL
            and cross_rep_failures == 0
        )
        else "FAIL"
    )

    manifest = {
        "stage": "12B2-PU",
        "status": status,
        "experiment_families": len(master),
        "diagnostic_metrics": metrics,
        "primary_diagnostic_metric": primary,
        "cross_representation_test_failures": (
            cross_rep_failures
        ),
        "interpretation": (
            "PU performance metrics are observed-label diagnostics. "
            "U is not a confirmed negative class."
        ),
    }

    (
        out / "stage12b2_manifest.json"
    ).write_text(
        json.dumps(
            manifest,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    print("=" * 100)
    print(
        "PhageLysClass — Stage 12B2 comprehensive PU analysis"
    )
    print("=" * 100)
    print(
        f"Experiment families              : {len(master):,}"
    )
    print(
        f"Diagnostic metrics               : {', '.join(metrics)}"
    )
    print(
        f"Primary diagnostic metric        : {primary}"
    )
    print(
        f"Method pairwise rows             : {len(method_deltas):,}"
    )
    print(
        f"Representation pairwise rows     : {len(representation_deltas):,}"
    )
    print(
        f"Algorithm pairwise rows          : {len(algorithm_deltas):,}"
    )
    print(
        f"All R×S pairwise rows            : {len(all_pairs):,}"
    )
    print(
        f"All R deltas at fixed S          : {len(r_deltas):,}"
    )
    print(
        f"All S deltas at fixed R          : {len(s_deltas):,}"
    )
    print(
        f"All 2×2 interaction rectangles   : {len(interactions):,}"
    )
    print(
        f"Method score-concordance rows    : {len(method_concordance):,}"
    )
    print(
        f"Representation concordance rows  : {len(representation_concordance):,}"
    )
    print(
        f"Cross-representation test failures: {cross_rep_failures:,}"
    )
    print(
        f"STATUS                           : {status}"
    )
    print(
        f"Outputs                          : {out}"
    )
    print("=" * 100)


if __name__ == "__main__":
    main()

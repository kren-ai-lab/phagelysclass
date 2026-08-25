#!/usr/bin/env python3
"""
12b_build_supervised_analysis_tables.py

Comprehensive factorial analysis tables for PhageLysClass supervised results.

This implementation evaluates the COMPLETE 4×4 redundancy × split surface.

For every:

    task × representation × algorithm × metric

the script generates:

1. ALL pairwise comparisons among the 16 R×S cells:
       C(16,2) = 120 cell-pair deltas

2. ALL redundancy contrasts at fixed S:
       C(4,2) × 4 = 24 deltas

3. ALL split-stringency contrasts at fixed R:
       C(4,2) × 4 = 24 deltas

4. ALL 2×2 R×S interaction rectangles:
       C(4,2) × C(4,2) = 36 interaction contrasts

5. Extreme contrasts remain only as compact summaries.

Important
---------
Native test membership can differ across R and/or S. Therefore R/S deltas
are factorial native-test contrasts, not automatically entity-level paired
comparisons. Common-test sensitivity analysis is performed in a later stage.

Inputs
------
results/analysis_supervised/stage12a/
    supervised_master_results.tsv
    validation_configuration_results.tsv

Outputs
-------
results/analysis_supervised/stage12b/
    representation_pairwise_deltas.tsv
    algorithm_pairwise_deltas.tsv
    rs_native_surfaces.tsv
    rs_all_pairwise_cell_deltas.tsv
    rs_redundancy_deltas_fixed_S.tsv
    rs_split_deltas_fixed_R.tsv
    rs_all_interaction_rectangles.tsv
    rs_extreme_summary_contrasts.tsv
    marginal_by_R.tsv
    marginal_by_S.tsv
    negative_definition_contrasts.tsv
    configuration_selection_frequency.tsv
    best_model_per_cell.tsv
    robustness_summary.tsv
    stage12b_summary.tsv
"""

from __future__ import annotations

from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd


R_LEVELS = [90, 70, 50, 30]
S_LEVELS = [90, 70, 50, 30]


def root_from_script() -> Path:
    return Path(__file__).resolve().parents[2]


def numeric_metric_columns(
    df: pd.DataFrame,
) -> list[str]:
    preferred = [
        "mcc",
        "ap",
        "ap_gain",
        "roc_auc",
        "f1",
        "precision",
        "recall",
        "specificity",
        "balanced_accuracy",
        "macro_f1",
        "macro_precision",
        "macro_recall",
        "weighted_f1",
        "log_loss",
        "brier",
    ]

    return [
        column
        for column in preferred
        if column in df.columns
    ]


def pairwise_model_deltas(
    df: pd.DataFrame,
    *,
    compare_col: str,
    matched_cols: list[str],
    metrics: list[str],
) -> pd.DataFrame:
    """
    Pairwise model comparisons within identical native test cells.

    Representation comparison:
        fixed task × R × S × algorithm

    Algorithm comparison:
        fixed task × R × S × representation
    """
    rows = []

    values = sorted(
        df[compare_col]
        .dropna()
        .astype(str)
        .unique()
    )

    for value_a, value_b in combinations(
        values,
        2,
    ):
        left = df.loc[
            df[compare_col].astype(str).eq(
                value_a
            ),
            matched_cols + metrics,
        ].copy()

        right = df.loc[
            df[compare_col].astype(str).eq(
                value_b
            ),
            matched_cols + metrics,
        ].copy()

        left = left.rename(
            columns={
                metric: f"{metric}__A"
                for metric in metrics
            }
        )

        right = right.rename(
            columns={
                metric: f"{metric}__B"
                for metric in metrics
            }
        )

        merged = left.merge(
            right,
            on=matched_cols,
            how="inner",
        )

        for record in merged.to_dict("records"):
            out = {
                **{
                    column: record[column]
                    for column in matched_cols
                },
                f"{compare_col}_A": value_a,
                f"{compare_col}_B": value_b,
            }

            for metric in metrics:
                value_metric_a = record.get(
                    f"{metric}__A"
                )
                value_metric_b = record.get(
                    f"{metric}__B"
                )

                out[f"{metric}_A"] = (
                    value_metric_a
                )
                out[f"{metric}_B"] = (
                    value_metric_b
                )

                out[
                    f"delta_{metric}_B_minus_A"
                ] = (
                    float(value_metric_b)
                    - float(value_metric_a)
                    if (
                        pd.notna(value_metric_a)
                        and pd.notna(value_metric_b)
                    )
                    else np.nan
                )

            rows.append(out)

    return pd.DataFrame(rows)


def lookup_group(
    group: pd.DataFrame,
    metric: str,
) -> dict[tuple[int, int], float]:
    return {
        (
            int(row.R),
            int(row.S),
        ): float(
            getattr(row, metric)
        )
        for row in group.itertuples()
        if pd.notna(
            getattr(row, metric)
        )
    }


def all_cell_pair_deltas(
    df: pd.DataFrame,
    metric: str,
) -> pd.DataFrame:
    """
    All unordered pairwise comparisons among the 16 R×S cells.

    Pair orientation:
    - same S:
        B = stricter redundancy condition
    - same R:
        B = stricter split condition
    - both factors stricter:
        B = stricter in both
    - tradeoff pairs:
        no natural ordering exists; retain canonical pair order
    """
    rows = []

    group_cols = [
        "task",
        "representation",
        "algorithm",
    ]

    for (
        task,
        representation,
        algorithm,
    ), group in df.groupby(group_cols):

        lookup = lookup_group(
            group,
            metric,
        )

        cells = [
            (r, s)
            for r in R_LEVELS
            for s in S_LEVELS
            if (r, s) in lookup
        ]

        for cell_1, cell_2 in combinations(
            cells,
            2,
        ):
            r1, s1 = cell_1
            r2, s2 = cell_2

            if s1 == s2:
                if r1 > r2:
                    cell_a, cell_b = (
                        cell_1,
                        cell_2,
                    )
                else:
                    cell_a, cell_b = (
                        cell_2,
                        cell_1,
                    )

                relation = (
                    "same_S_redundancy"
                )

            elif r1 == r2:
                if s1 > s2:
                    cell_a, cell_b = (
                        cell_1,
                        cell_2,
                    )
                else:
                    cell_a, cell_b = (
                        cell_2,
                        cell_1,
                    )

                relation = "same_R_split"

            else:
                cell_2_stricter_both = (
                    r2 < r1
                    and s2 < s1
                )

                cell_1_stricter_both = (
                    r1 < r2
                    and s1 < s2
                )

                if cell_2_stricter_both:
                    cell_a, cell_b = (
                        cell_1,
                        cell_2,
                    )
                    relation = (
                        "both_factors_stricter"
                    )

                elif cell_1_stricter_both:
                    cell_a, cell_b = (
                        cell_2,
                        cell_1,
                    )
                    relation = (
                        "both_factors_stricter"
                    )

                else:
                    cell_a, cell_b = (
                        cell_1,
                        cell_2,
                    )
                    relation = (
                        "tradeoff_incomparable"
                    )

            r_a, s_a = cell_a
            r_b, s_b = cell_b

            value_a = lookup[cell_a]
            value_b = lookup[cell_b]

            rows.append(
                {
                    "task": task,
                    "representation": (
                        representation
                    ),
                    "algorithm": algorithm,
                    "metric": metric,
                    "relation": relation,
                    "R_A": r_a,
                    "S_A": s_a,
                    "R_B": r_b,
                    "S_B": s_b,
                    "value_A": value_a,
                    "value_B": value_b,
                    "delta_B_minus_A": (
                        value_b - value_a
                    ),
                    "delta_R_threshold": (
                        r_b - r_a
                    ),
                    "delta_S_threshold": (
                        s_b - s_a
                    ),
                    "native_test_note": (
                        "Factorial native-test "
                        "contrast; not necessarily "
                        "entity-level paired."
                    ),
                }
            )

    return pd.DataFrame(rows)


def redundancy_deltas_fixed_s(
    df: pd.DataFrame,
    metric: str,
) -> pd.DataFrame:
    """
    All six R-pair deltas at each fixed S.

    6 R pairs × 4 S levels = 24 contrasts
    per task × representation × algorithm × metric.
    """
    rows = []

    for (
        task,
        representation,
        algorithm,
    ), group in df.groupby(
        [
            "task",
            "representation",
            "algorithm",
        ]
    ):
        lookup = lookup_group(
            group,
            metric,
        )

        for s in S_LEVELS:
            for (
                r_loose,
                r_strict,
            ) in combinations(
                R_LEVELS,
                2,
            ):
                if (
                    (r_loose, s)
                    not in lookup
                    or (r_strict, s)
                    not in lookup
                ):
                    continue

                value_from = lookup[
                    (r_loose, s)
                ]

                value_to = lookup[
                    (r_strict, s)
                ]

                rows.append(
                    {
                        "task": task,
                        "representation": (
                            representation
                        ),
                        "algorithm": algorithm,
                        "metric": metric,
                        "S_fixed": s,
                        "R_from": r_loose,
                        "R_to": r_strict,
                        "threshold_step": (
                            f"R{r_loose}->"
                            f"R{r_strict}"
                        ),
                        "value_from": (
                            value_from
                        ),
                        "value_to": value_to,
                        "delta": (
                            value_to
                            - value_from
                        ),
                        "absolute_threshold_change": (
                            r_loose
                            - r_strict
                        ),
                        "native_test_note": (
                            "Native-test "
                            "redundancy contrast "
                            "at fixed S; not "
                            "entity-level paired."
                        ),
                    }
                )

    return pd.DataFrame(rows)


def split_deltas_fixed_r(
    df: pd.DataFrame,
    metric: str,
) -> pd.DataFrame:
    """
    All six S-pair deltas at each fixed R.

    6 S pairs × 4 R levels = 24 contrasts
    per task × representation × algorithm × metric.
    """
    rows = []

    for (
        task,
        representation,
        algorithm,
    ), group in df.groupby(
        [
            "task",
            "representation",
            "algorithm",
        ]
    ):
        lookup = lookup_group(
            group,
            metric,
        )

        for r in R_LEVELS:
            for (
                s_loose,
                s_strict,
            ) in combinations(
                S_LEVELS,
                2,
            ):
                if (
                    (r, s_loose)
                    not in lookup
                    or (r, s_strict)
                    not in lookup
                ):
                    continue

                value_from = lookup[
                    (r, s_loose)
                ]

                value_to = lookup[
                    (r, s_strict)
                ]

                rows.append(
                    {
                        "task": task,
                        "representation": (
                            representation
                        ),
                        "algorithm": algorithm,
                        "metric": metric,
                        "R_fixed": r,
                        "S_from": s_loose,
                        "S_to": s_strict,
                        "threshold_step": (
                            f"S{s_loose}->"
                            f"S{s_strict}"
                        ),
                        "value_from": (
                            value_from
                        ),
                        "value_to": value_to,
                        "delta": (
                            value_to
                            - value_from
                        ),
                        "absolute_threshold_change": (
                            s_loose
                            - s_strict
                        ),
                        "native_test_note": (
                            "Native-test "
                            "split-stringency "
                            "contrast at fixed R; "
                            "not entity-level paired."
                        ),
                    }
                )

    return pd.DataFrame(rows)


def all_interaction_rectangles(
    df: pd.DataFrame,
    metric: str,
) -> pd.DataFrame:
    """
    Every 2×2 rectangle in the 4×4 factorial.

    C(4,2) × C(4,2) = 36 interactions
    per task × representation × algorithm × metric.

    For:

        R_loose > R_strict
        S_loose > S_strict

    interaction =

        [M(R_strict,S_strict)
         - M(R_loose,S_strict)]

        -

        [M(R_strict,S_loose)
         - M(R_loose,S_loose)]
    """
    rows = []

    for (
        task,
        representation,
        algorithm,
    ), group in df.groupby(
        [
            "task",
            "representation",
            "algorithm",
        ]
    ):
        lookup = lookup_group(
            group,
            metric,
        )

        for (
            r_loose,
            r_strict,
        ) in combinations(
            R_LEVELS,
            2,
        ):
            for (
                s_loose,
                s_strict,
            ) in combinations(
                S_LEVELS,
                2,
            ):
                corners = [
                    (
                        r_loose,
                        s_loose,
                    ),
                    (
                        r_strict,
                        s_loose,
                    ),
                    (
                        r_loose,
                        s_strict,
                    ),
                    (
                        r_strict,
                        s_strict,
                    ),
                ]

                if not all(
                    corner in lookup
                    for corner in corners
                ):
                    continue

                m_rloose_sloose = lookup[
                    (
                        r_loose,
                        s_loose,
                    )
                ]

                m_rstrict_sloose = lookup[
                    (
                        r_strict,
                        s_loose,
                    )
                ]

                m_rloose_sstrict = lookup[
                    (
                        r_loose,
                        s_strict,
                    )
                ]

                m_rstrict_sstrict = lookup[
                    (
                        r_strict,
                        s_strict,
                    )
                ]

                delta_r_at_s_loose = (
                    m_rstrict_sloose
                    - m_rloose_sloose
                )

                delta_r_at_s_strict = (
                    m_rstrict_sstrict
                    - m_rloose_sstrict
                )

                interaction = (
                    delta_r_at_s_strict
                    - delta_r_at_s_loose
                )

                delta_s_at_r_loose = (
                    m_rloose_sstrict
                    - m_rloose_sloose
                )

                delta_s_at_r_strict = (
                    m_rstrict_sstrict
                    - m_rstrict_sloose
                )

                rows.append(
                    {
                        "task": task,
                        "representation": (
                            representation
                        ),
                        "algorithm": algorithm,
                        "metric": metric,
                        "R_loose": r_loose,
                        "R_strict": r_strict,
                        "S_loose": s_loose,
                        "S_strict": s_strict,
                        "M_Rloose_Sloose": (
                            m_rloose_sloose
                        ),
                        "M_Rstrict_Sloose": (
                            m_rstrict_sloose
                        ),
                        "M_Rloose_Sstrict": (
                            m_rloose_sstrict
                        ),
                        "M_Rstrict_Sstrict": (
                            m_rstrict_sstrict
                        ),
                        "delta_R_at_Sloose": (
                            delta_r_at_s_loose
                        ),
                        "delta_R_at_Sstrict": (
                            delta_r_at_s_strict
                        ),
                        "delta_S_at_Rloose": (
                            delta_s_at_r_loose
                        ),
                        "delta_S_at_Rstrict": (
                            delta_s_at_r_strict
                        ),
                        "interaction": (
                            interaction
                        ),
                        "interaction_equivalence_check": (
                            delta_s_at_r_strict
                            - delta_s_at_r_loose
                        ),
                        "native_test_note": (
                            "2x2 factorial "
                            "interaction over "
                            "native test cells; "
                            "not an entity-level "
                            "paired interaction."
                        ),
                    }
                )

    return pd.DataFrame(rows)


def extreme_summary(
    df: pd.DataFrame,
    metric: str,
) -> pd.DataFrame:
    """
    Extreme contrasts retained only as compact summaries.
    """
    rows = []

    for (
        task,
        representation,
        algorithm,
    ), group in df.groupby(
        [
            "task",
            "representation",
            "algorithm",
        ]
    ):
        lookup = lookup_group(
            group,
            metric,
        )

        for s in S_LEVELS:
            if (
                (90, s) in lookup
                and (30, s) in lookup
            ):
                rows.append(
                    {
                        "task": task,
                        "representation": (
                            representation
                        ),
                        "algorithm": algorithm,
                        "metric": metric,
                        "contrast": (
                            "R30_minus_R90"
                        ),
                        "fixed_factor": "S",
                        "fixed_value": s,
                        "value": (
                            lookup[(30, s)]
                            - lookup[(90, s)]
                        ),
                    }
                )

        for r in R_LEVELS:
            if (
                (r, 90) in lookup
                and (r, 30) in lookup
            ):
                rows.append(
                    {
                        "task": task,
                        "representation": (
                            representation
                        ),
                        "algorithm": algorithm,
                        "metric": metric,
                        "contrast": (
                            "S30_minus_S90"
                        ),
                        "fixed_factor": "R",
                        "fixed_value": r,
                        "value": (
                            lookup[(r, 30)]
                            - lookup[(r, 90)]
                        ),
                    }
                )

        corners = [
            (90, 90),
            (30, 90),
            (90, 30),
            (30, 30),
        ]

        if all(
            corner in lookup
            for corner in corners
        ):
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
                    "representation": (
                        representation
                    ),
                    "algorithm": algorithm,
                    "metric": metric,
                    "contrast": (
                        "extreme_RxS_interaction"
                    ),
                    "fixed_factor": "",
                    "fixed_value": np.nan,
                    "value": interaction,
                }
            )

    return pd.DataFrame(rows)


def main():
    root = root_from_script()

    stage12a = (
        root
        / "results/analysis_supervised/stage12a"
    )

    master_path = (
        stage12a
        / "supervised_master_results.tsv"
    )

    configuration_path = (
        stage12a
        / "validation_configuration_results.tsv"
    )

    if not master_path.is_file():
        raise FileNotFoundError(
            f"{master_path}\nRun Stage 12A first."
        )

    if not configuration_path.is_file():
        raise FileNotFoundError(
            configuration_path
        )

    df = pd.read_csv(
        master_path,
        sep="\t",
        low_memory=False,
    )

    configurations = pd.read_csv(
        configuration_path,
        sep="\t",
        low_memory=False,
    )

    out = (
        root
        / "results/analysis_supervised/stage12b"
    )
    out.mkdir(
        parents=True,
        exist_ok=True,
    )

    metrics = numeric_metric_columns(df)

    factorial_metrics = [
        metric
        for metric in [
            "mcc",
            "ap",
            "ap_gain",
            "macro_f1",
        ]
        if metric in df.columns
    ]

    # --------------------------------------------------
    # A. Representation comparisons
    # --------------------------------------------------

    representation_deltas = (
        pairwise_model_deltas(
            df,
            compare_col="representation",
            matched_cols=[
                "task",
                "R",
                "S",
                "algorithm",
            ],
            metrics=factorial_metrics,
        )
    )

    representation_deltas.to_csv(
        out
        / "representation_pairwise_deltas.tsv",
        sep="\t",
        index=False,
    )

    # --------------------------------------------------
    # B. Algorithm comparisons
    # --------------------------------------------------

    algorithm_deltas = (
        pairwise_model_deltas(
            df,
            compare_col="algorithm",
            matched_cols=[
                "task",
                "R",
                "S",
                "representation",
            ],
            metrics=factorial_metrics,
        )
    )

    algorithm_deltas.to_csv(
        out
        / "algorithm_pairwise_deltas.tsv",
        sep="\t",
        index=False,
    )

    # --------------------------------------------------
    # C. Native R×S surfaces
    # --------------------------------------------------

    surface_columns = [
        "experiment_id",
        "task",
        "task_type",
        "R",
        "S",
        "representation",
        "algorithm",
        "selected_config",
    ] + metrics

    surfaces = df[
        [
            column
            for column in surface_columns
            if column in df.columns
        ]
    ].copy()

    surfaces.to_csv(
        out / "rs_native_surfaces.tsv",
        sep="\t",
        index=False,
    )

    # --------------------------------------------------
    # D. Complete R×S contrast system
    # --------------------------------------------------

    all_pair_parts = []
    redundancy_parts = []
    split_parts = []
    interaction_parts = []
    extreme_parts = []

    for metric in factorial_metrics:
        all_pair_parts.append(
            all_cell_pair_deltas(
                df,
                metric,
            )
        )

        redundancy_parts.append(
            redundancy_deltas_fixed_s(
                df,
                metric,
            )
        )

        split_parts.append(
            split_deltas_fixed_r(
                df,
                metric,
            )
        )

        interaction_parts.append(
            all_interaction_rectangles(
                df,
                metric,
            )
        )

        extreme_parts.append(
            extreme_summary(
                df,
                metric,
            )
        )

    all_pairs = pd.concat(
        all_pair_parts,
        ignore_index=True,
    )

    redundancy_deltas = pd.concat(
        redundancy_parts,
        ignore_index=True,
    )

    split_deltas = pd.concat(
        split_parts,
        ignore_index=True,
    )

    interactions = pd.concat(
        interaction_parts,
        ignore_index=True,
    )

    extreme_summaries = pd.concat(
        extreme_parts,
        ignore_index=True,
    )

    all_pairs.to_csv(
        out
        / "rs_all_pairwise_cell_deltas.tsv",
        sep="\t",
        index=False,
    )

    redundancy_deltas.to_csv(
        out
        / "rs_redundancy_deltas_fixed_S.tsv",
        sep="\t",
        index=False,
    )

    split_deltas.to_csv(
        out
        / "rs_split_deltas_fixed_R.tsv",
        sep="\t",
        index=False,
    )

    interactions.to_csv(
        out
        / "rs_all_interaction_rectangles.tsv",
        sep="\t",
        index=False,
    )

    extreme_summaries.to_csv(
        out
        / "rs_extreme_summary_contrasts.tsv",
        sep="\t",
        index=False,
    )

    # --------------------------------------------------
    # E. Marginal summaries
    # --------------------------------------------------

    aggregation_spec = {}

    for metric in factorial_metrics:
        aggregation_spec[
            f"{metric}_mean"
        ] = (
            metric,
            "mean",
        )

        aggregation_spec[
            f"{metric}_median"
        ] = (
            metric,
            "median",
        )

        aggregation_spec[
            f"{metric}_min"
        ] = (
            metric,
            "min",
        )

        aggregation_spec[
            f"{metric}_max"
        ] = (
            metric,
            "max",
        )

    marginal_by_r = (
        df.groupby(
            [
                "task",
                "representation",
                "algorithm",
                "R",
            ],
            as_index=False,
        )
        .agg(**aggregation_spec)
    )

    marginal_by_s = (
        df.groupby(
            [
                "task",
                "representation",
                "algorithm",
                "S",
            ],
            as_index=False,
        )
        .agg(**aggregation_spec)
    )

    marginal_by_r.to_csv(
        out / "marginal_by_R.tsv",
        sep="\t",
        index=False,
    )

    marginal_by_s.to_csv(
        out / "marginal_by_S.tsv",
        sep="\t",
        index=False,
    )

    # --------------------------------------------------
    # F. Negative-definition contrasts
    # --------------------------------------------------

    global_tasks = df.loc[
        df["task"].isin(
            [
                "global_easy",
                "global_hard",
                "global_mixed",
            ]
        )
    ].copy()

    negative_rows = []

    for metric in [
        metric
        for metric in [
            "mcc",
            "ap",
            "ap_gain",
        ]
        if metric in df.columns
    ]:
        pivot = global_tasks.pivot_table(
            index=[
                "R",
                "S",
                "representation",
                "algorithm",
            ],
            columns="task",
            values=metric,
            aggfunc="first",
        ).reset_index()

        required = {
            "global_easy",
            "global_hard",
            "global_mixed",
        }

        if not required <= set(pivot.columns):
            continue

        for row in pivot.itertuples(
            index=False
        ):
            negative_rows.append(
                {
                    "R": row.R,
                    "S": row.S,
                    "representation": (
                        row.representation
                    ),
                    "algorithm": row.algorithm,
                    "metric": metric,
                    "global_easy": (
                        row.global_easy
                    ),
                    "global_hard": (
                        row.global_hard
                    ),
                    "global_mixed": (
                        row.global_mixed
                    ),
                    "hard_minus_easy": (
                        row.global_hard
                        - row.global_easy
                    ),
                    "mixed_minus_easy": (
                        row.global_mixed
                        - row.global_easy
                    ),
                    "hard_minus_mixed": (
                        row.global_hard
                        - row.global_mixed
                    ),
                    "comparison_note": (
                        "Matched-factor "
                        "native-test contrast; "
                        "not entity-level paired."
                    ),
                }
            )

    pd.DataFrame(
        negative_rows
    ).to_csv(
        out
        / "negative_definition_contrasts.tsv",
        sep="\t",
        index=False,
    )

    # --------------------------------------------------
    # G. Configuration-selection frequency
    # --------------------------------------------------

    if "selected" in configurations.columns:
        selected_mask = (
            configurations["selected"]
            .astype(str)
            .str.lower()
            .isin(
                {
                    "true",
                    "1",
                    "yes",
                }
            )
        )

        selected_configurations = (
            configurations.loc[
                selected_mask
            ].copy()
        )
    else:
        selected_configurations = (
            pd.DataFrame()
        )

    if len(selected_configurations):
        configuration_frequency = (
            selected_configurations
            .groupby(
                [
                    "task",
                    "representation",
                    "algorithm",
                    "config_id",
                ],
                as_index=False,
            )
            .size()
            .rename(
                columns={
                    "size": "n_selected"
                }
            )
        )
    else:
        configuration_frequency = (
            pd.DataFrame()
        )

    configuration_frequency.to_csv(
        out
        / "configuration_selection_frequency.tsv",
        sep="\t",
        index=False,
    )

    # --------------------------------------------------
    # H. Best model per native cell
    # --------------------------------------------------

    if "mcc" in df.columns:
        best_model = (
            df.sort_values(
                [
                    "task",
                    "R",
                    "S",
                    "mcc",
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
                [
                    "task",
                    "R",
                    "S",
                ],
                as_index=False,
            )
            .first()
        )
    else:
        best_model = pd.DataFrame()

    best_model.to_csv(
        out / "best_model_per_cell.tsv",
        sep="\t",
        index=False,
    )

    # --------------------------------------------------
    # I. Robustness summaries
    # --------------------------------------------------

    robustness_rows = []

    for (
        task,
        representation,
        algorithm,
    ), group in df.groupby(
        [
            "task",
            "representation",
            "algorithm",
        ]
    ):
        row = {
            "task": task,
            "representation": (
                representation
            ),
            "algorithm": algorithm,
            "n_rs_cells": len(group),
        }

        for metric in factorial_metrics:
            if group[metric].notna().any():
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
        out / "robustness_summary.tsv",
        sep="\t",
        index=False,
    )

    # --------------------------------------------------
    # J. Stage summary
    # --------------------------------------------------

    n_factorial_groups = (
        df[
            [
                "task",
                "representation",
                "algorithm",
            ]
        ]
        .drop_duplicates()
        .shape[0]
    )

    summary = pd.DataFrame(
        [
            {
                "item": "experiment_families",
                "value": len(df),
            },
            {
                "item": "factorial_model_groups",
                "value": n_factorial_groups,
            },
            {
                "item": "factorial_metrics",
                "value": ",".join(
                    factorial_metrics
                ),
            },
            {
                "item": "all_RxS_pairwise_rows",
                "value": len(all_pairs),
            },
            {
                "item": "redundancy_fixed_S_rows",
                "value": len(
                    redundancy_deltas
                ),
            },
            {
                "item": "split_fixed_R_rows",
                "value": len(
                    split_deltas
                ),
            },
            {
                "item": (
                    "interaction_rectangle_rows"
                ),
                "value": len(interactions),
            },
            {
                "item": "extreme_summary_rows",
                "value": len(
                    extreme_summaries
                ),
            },
            {
                "item": (
                    "representation_pairwise_rows"
                ),
                "value": len(
                    representation_deltas
                ),
            },
            {
                "item": (
                    "algorithm_pairwise_rows"
                ),
                "value": len(
                    algorithm_deltas
                ),
            },
        ]
    )

    summary.to_csv(
        out / "stage12b_summary.tsv",
        sep="\t",
        index=False,
    )

    print("=" * 100)
    print(
        "PhageLysClass — Stage 12B "
        "comprehensive R×S analysis"
    )
    print("=" * 100)
    print(
        f"Experiment families           : "
        f"{len(df):,}"
    )
    print(
        f"Task×rep×algorithm groups     : "
        f"{n_factorial_groups:,}"
    )
    print(
        "Factorial metrics             : "
        + ", ".join(factorial_metrics)
    )
    print(
        f"All 16-cell pairwise deltas   : "
        f"{len(all_pairs):,}"
    )
    print(
        f"All R-pair deltas at fixed S  : "
        f"{len(redundancy_deltas):,}"
    )
    print(
        f"All S-pair deltas at fixed R  : "
        f"{len(split_deltas):,}"
    )
    print(
        f"All 2×2 interaction rectangles: "
        f"{len(interactions):,}"
    )
    print(
        f"Extreme summaries only        : "
        f"{len(extreme_summaries):,}"
    )
    print(f"Outputs                       : {out}")
    print("=" * 100)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
A4 — Prediction × property × geometry integration.

Purpose
-------
Integrate the frozen PhageLysClass predictions (A0/A1), the independent
26-descriptor property characterization (A2), and the full-space/reference
geometry analysis (A3).

A4 does NOT re-train models, re-select descriptors, re-fit decision thresholds,
or use PCA/UMAP/t-SNE coordinates for quantitative inference.

Core integration logic
----------------------
For every globally lytic sequence, three independent subtype-support axes are
retained:

1. Property support
   Fraction of the 26 frozen descriptors inside the predicted subtype's
   annotated training-reference 95% interval.

2. Geometric coverage
   Whether the sequence falls inside the predicted subtype's empirical 95%
   full-space reference coverage, based on leave-one-out nearest-neighbor
   distances in the standardized 26D descriptor space.

3. Local geometric agreement
   Majority agreement of three interpretable full-space rules:
       - kNN majority subtype
       - raw nearest-neighbor subtype
       - nearest centroid subtype
   Density-normalized NN remains in the master as a diagnostic but is not used
   in the primary concordance call because its behavior is strongly affected
   by subtype-specific reference density.

Operational subtype concordance
-------------------------------
property_supported:
    fraction_within_reference_95 >= --property-support-threshold (default 0.80)

geometry_supported:
    predicted_subtype_coverage95 == True

local_geometry_supported:
    >= --geometry-votes-required of the 3 primary geometry rules agree
    with the PhageLysClass functional prediction (default 2 of 3)

functional_support_count:
    0..3 support axes satisfied

functional_concordance:
    3/3  -> concordant
    2/3  -> mixed
    0-1  -> discordant

This is an explicit operational evidence-integration rule, not a probability
calibration or a new classifier.

Visualization
-------------
- Publication candidate multipanel integrating A1–A4.
- Geometry matrix with PCA, reference-fitted UMAP and joint t-SNE for all
  three external sources.
- t-SNE is visualization only and is never used quantitatively.
- PCA and UMAP are reference-fitted transforms generated in A3.

Default execution
-----------------
python analysis_explored_seqs/A4_prediction_geometry_integration.py
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import spearmanr


SOURCE_ORDER = [
    "phagelys_prediction_extension",
    "zymctrl",
    "ambiguous_structural",
]

SOURCE_LABELS = {
    "phagelys_prediction_extension": "Prediction Extension",
    "zymctrl": "ZymCTRL",
    "ambiguous_structural": "Ambiguous structural",
}

FUNCTION_ORDER = [
    "endolysin",
    "depolymerase",
    "val",
]

FUNCTION_LABELS = {
    "endolysin": "Endolysin",
    "depolymerase": "Depolymerase",
    "val": "VAL",
}

CONCORDANCE_ORDER = [
    "concordant",
    "mixed",
    "discordant",
]

sns.set_theme(style="whitegrid", context="paper")


# ---------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="A4 prediction/property/geometry integration."
    )

    p.add_argument(
        "--master",
        type=Path,
        default=Path(
            "analysis_explored_seqs/"
            "A0_integrate_external_sources/"
            "external_sources_master.tsv.gz"
        ),
    )
    p.add_argument(
        "--a1-dir",
        type=Path,
        default=Path(
            "analysis_explored_seqs/"
            "A1_prediction_landscape/"
            "tables"
        ),
    )
    p.add_argument(
        "--a2-dir",
        type=Path,
        default=Path(
            "analysis_explored_seqs/"
            "A2_property_characterization/"
            "tables"
        ),
    )
    p.add_argument(
        "--a3-dir",
        type=Path,
        default=Path(
            "analysis_explored_seqs/"
            "A3_reference_geometry/"
            "tables"
        ),
    )
    p.add_argument(
        "--output-dir",
        type=Path,
        default=Path(
            "analysis_explored_seqs/"
            "A4_prediction_geometry_integration"
        ),
    )
    p.add_argument(
        "--property-support-threshold",
        type=float,
        default=0.80,
        help=(
            "Minimum fraction of 26 descriptors within the predicted "
            "subtype reference 95% interval to call property support."
        ),
    )
    p.add_argument(
        "--geometry-votes-required",
        type=int,
        default=2,
        choices=[1, 2, 3],
        help=(
            "Number of primary geometry rules that must agree with the "
            "predicted subtype to call local geometric support."
        ),
    )
    p.add_argument(
        "--plot-max-reference",
        type=int,
        default=2500,
    )
    p.add_argument(
        "--plot-max-external-per-source",
        type=int,
        default=2500,
    )
    p.add_argument(
        "--seed",
        type=int,
        default=42,
    )

    return p.parse_args()


# ---------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------

def require_file(path: Path) -> Path:
    if not path.is_file():
        raise FileNotFoundError(path)
    return path


def optional_file(path: Path) -> Path | None:
    return path if path.is_file() else None


def canonical_function(value: object) -> str | None:
    if pd.isna(value):
        return None

    x = str(value).strip().lower()

    mapping = {
        "endolysin": "endolysin",
        "depolymerase": "depolymerase",
        "val": "val",
        "virion_associated_lytic_enzyme": "val",
        "virion-associated lytic enzyme": "val",
    }
    return mapping.get(x)


def save_figure(fig: plt.Figure, stem: Path) -> None:
    fig.savefig(
        stem.with_suffix(".png"),
        dpi=300,
        bbox_inches="tight",
    )
    fig.savefig(
        stem.with_suffix(".pdf"),
        bbox_inches="tight",
    )
    plt.close(fig)


def panel_label(
    ax: plt.Axes,
    label: str,
) -> None:
    ax.text(
        -0.10,
        1.06,
        label,
        transform=ax.transAxes,
        fontsize=14,
        fontweight="bold",
        va="top",
        ha="left",
    )


def stratified_sample(
    df: pd.DataFrame,
    group_col: str,
    max_total: int,
    seed: int,
) -> pd.DataFrame:
    if len(df) <= max_total:
        return df.copy()

    counts = df[group_col].value_counts(dropna=False)

    if counts.empty:
        return df.sample(
            n=max_total,
            random_state=seed,
        )

    allocations = {}
    for group, n in counts.items():
        allocations[group] = max(
            1,
            int(round(max_total * n / counts.sum())),
        )

    while sum(allocations.values()) > max_total:
        largest = max(
            allocations,
            key=allocations.get,
        )
        if allocations[largest] > 1:
            allocations[largest] -= 1
        else:
            break

    pieces = []

    for i, (group, n_take) in enumerate(
        allocations.items()
    ):
        if pd.isna(group):
            g = df[df[group_col].isna()]
        else:
            g = df[df[group_col] == group]

        n_take = min(n_take, len(g))

        pieces.append(
            g.sample(
                n=n_take,
                random_state=seed + i,
            )
        )

    return pd.concat(
        pieces,
        ignore_index=True,
    )


def prefix_columns(
    df: pd.DataFrame,
    prefix: str,
    key: str = "analysis_id",
) -> pd.DataFrame:
    rename = {
        c: f"{prefix}{c}"
        for c in df.columns
        if c != key
    }
    return df.rename(columns=rename)


def safe_spearman(
    x: pd.Series,
    y: pd.Series,
) -> tuple[float, float, int]:
    pair = pd.DataFrame(
        {
            "x": pd.to_numeric(
                x,
                errors="coerce",
            ),
            "y": pd.to_numeric(
                y,
                errors="coerce",
            ),
        }
    ).dropna()

    if len(pair) < 3:
        return np.nan, np.nan, int(len(pair))

    if (
        pair["x"].nunique() < 2
        or pair["y"].nunique() < 2
    ):
        return np.nan, np.nan, int(len(pair))

    rho, p = spearmanr(
        pair["x"],
        pair["y"],
    )

    return float(rho), float(p), int(len(pair))


# ---------------------------------------------------------------------
# Load + integrate
# ---------------------------------------------------------------------

def load_inputs(
    args: argparse.Namespace,
) -> dict[str, pd.DataFrame]:
    master_path = require_file(
        args.master.resolve()
    )

    a2 = args.a2_dir.resolve()
    a3 = args.a3_dir.resolve()

    paths = {
        "master": master_path,
        "property_global": require_file(
            a2 / "individual_global_reference_fit.tsv.gz"
        ),
        "property_functional": require_file(
            a2 / "individual_functional_reference_fit.tsv.gz"
        ),
        "property_group_stats": require_file(
            a2 / "group_property_comparisons.tsv"
        ),
        "geometry_global": require_file(
            a3 / "global_reference_coverage.tsv.gz"
        ),
        "geometry_functional": require_file(
            a3 / "functional_reference_geometry.tsv.gz"
        ),
        "pca_reference": require_file(
            a3 / "pca_reference_coordinates.tsv.gz"
        ),
        "pca_external": require_file(
            a3 / "pca_external_coordinates.tsv.gz"
        ),
        "umap_reference": require_file(
            a3 / "umap_reference_coordinates.tsv.gz"
        ),
        "umap_external": require_file(
            a3 / "umap_external_coordinates.tsv.gz"
        ),
    }

    frames = {}

    for name, path in paths.items():
        frames[name] = pd.read_csv(
            path,
            sep="\t",
            compression="infer",
            low_memory=False,
        )

    tsne_parts = []

    for source in SOURCE_ORDER:
        path = require_file(
            a3 / f"tsne_{source}_coordinates.tsv.gz"
        )
        x = pd.read_csv(
            path,
            sep="\t",
            compression="infer",
            low_memory=False,
        )
        x["tsne_source"] = source
        tsne_parts.append(x)

    frames["tsne"] = pd.concat(
        tsne_parts,
        ignore_index=True,
    )

    return frames


def build_integrated_master(
    frames: dict[str, pd.DataFrame],
    property_threshold: float,
    geometry_votes_required: int,
) -> pd.DataFrame:
    master = frames["master"].copy()

    required = {
        "analysis_id",
        "sequence_id",
        "source",
        "global_score",
        "global_prediction",
        "functional_prediction",
        "p_endolysin_given_lytic",
        "p_depolymerase_given_lytic",
        "p_val_given_lytic",
    }

    missing = sorted(
        required - set(master.columns)
    )
    if missing:
        raise ValueError(
            f"A0 master missing columns: {missing}"
        )

    master["functional_class_a4"] = (
        master["functional_prediction"].map(
            canonical_function
        )
    )

    # --------------------------------------------------------------
    # A2 property evidence
    # --------------------------------------------------------------
    prop_global = frames[
        "property_global"
    ][
        [
            "analysis_id",
            "fraction_within_reference_95",
            "fraction_within_reference_iqr",
            "median_abs_robust_z",
            "max_abs_robust_z",
            "n_abs_robust_z_gt2",
            "n_abs_robust_z_gt3",
        ]
    ].copy()

    prop_global = prefix_columns(
        prop_global,
        "prop_global_",
    )

    prop_func = frames[
        "property_functional"
    ][
        [
            "analysis_id",
            "reference_cohort",
            "fraction_within_reference_95",
            "fraction_within_reference_iqr",
            "median_abs_robust_z",
            "max_abs_robust_z",
            "n_abs_robust_z_gt2",
            "n_abs_robust_z_gt3",
        ]
    ].copy()

    prop_func = prefix_columns(
        prop_func,
        "prop_func_",
    )

    # --------------------------------------------------------------
    # A3 geometry evidence
    # --------------------------------------------------------------
    global_geom_cols = [
        "analysis_id",
        "exact_global_reference_overlap",
        "global_nn_distance",
        "global_mean_k_distance",
        "global_nn_distance_percentile",
        "global_mean_k_distance_percentile",
        "global_coverage95",
        "global_coverage99",
        "global_ood95",
        "global_ood99",
    ]

    geom_global = frames[
        "geometry_global"
    ][global_geom_cols].copy()

    geom_global = prefix_columns(
        geom_global,
        "geom_global_",
    )

    functional_geom_cols = [
        "analysis_id",
        "exact_multiclass_reference_overlap",
        "knn_majority_subtype",
        "knn_neighborhood_entropy",
        "nearest_raw_nn_subtype",
        "nearest_density_normalized_subtype",
        "nearest_centroid_subtype",
        "agreement_knn_majority",
        "agreement_raw_nn",
        "agreement_density_normalized_nn",
        "agreement_centroid",
        "predicted_subtype_nn_distance",
        "predicted_subtype_density_normalized_nn",
        "predicted_subtype_centroid_distance",
        "predicted_subtype_nn_percentile",
        "predicted_subtype_coverage95",
        "predicted_subtype_coverage99",
        "predicted_subtype_ood95",
        "predicted_subtype_ood99",
        "predicted_subtype_knn_fraction",
    ]

    geom_func = frames[
        "geometry_functional"
    ][functional_geom_cols].copy()

    geom_func = prefix_columns(
        geom_func,
        "geom_func_",
    )

    # --------------------------------------------------------------
    # PCA/UMAP external coordinates
    # --------------------------------------------------------------
    pca = frames["pca_external"][
        [
            "analysis_id",
            "PC1",
            "PC2",
        ]
    ].copy()

    umap = frames["umap_external"][
        [
            "analysis_id",
            "UMAP1",
            "UMAP2",
        ]
    ].copy()

    # External t-SNE coordinates exist only for the visualization sample.
    tsne_ext = frames["tsne"][
        frames["tsne"]["origin"] == "External"
    ][
        [
            "analysis_id",
            "TSNE1",
            "TSNE2",
        ]
    ].copy()

    tsne_ext = tsne_ext.drop_duplicates(
        "analysis_id"
    )

    # --------------------------------------------------------------
    # Merge
    # --------------------------------------------------------------
    out = master.merge(
        prop_global,
        on="analysis_id",
        how="left",
        validate="one_to_one",
    )

    out = out.merge(
        prop_func,
        on="analysis_id",
        how="left",
        validate="one_to_one",
    )

    out = out.merge(
        geom_global,
        on="analysis_id",
        how="left",
        validate="one_to_one",
    )

    out = out.merge(
        geom_func,
        on="analysis_id",
        how="left",
        validate="one_to_one",
    )

    out = out.merge(
        pca,
        on="analysis_id",
        how="left",
        validate="one_to_one",
    )

    out = out.merge(
        umap,
        on="analysis_id",
        how="left",
        validate="one_to_one",
    )

    out = out.merge(
        tsne_ext,
        on="analysis_id",
        how="left",
        validate="one_to_one",
    )

    # --------------------------------------------------------------
    # Prediction-strength descriptors
    # --------------------------------------------------------------
    probs = out[
        [
            "p_endolysin_given_lytic",
            "p_depolymerase_given_lytic",
            "p_val_given_lytic",
        ]
    ].apply(
        pd.to_numeric,
        errors="coerce",
    )

    prob_array = probs.to_numpy(dtype=float)

    with np.errstate(
        invalid="ignore",
    ):
        sorted_probs = np.sort(
            prob_array,
            axis=1,
        )

    # Use pandas' skip-NA semantics to avoid warnings on non-lytic rows,
    # where all conditional subtype probabilities are intentionally NA.
    out["functional_max_probability"] = probs.max(
        axis=1,
        skipna=True,
    )

    out["functional_probability_margin"] = (
        sorted_probs[:, -1]
        - sorted_probs[:, -2]
    )

    non_lytic = (
        out["global_prediction"]
        != "lytic"
    )
    out.loc[
        non_lytic,
        [
            "functional_max_probability",
            "functional_probability_margin",
        ],
    ] = np.nan

    # --------------------------------------------------------------
    # Exact training overlap audit
    # --------------------------------------------------------------
    global_overlap = out[
        "geom_global_exact_global_reference_overlap"
    ].fillna(False).astype(bool)

    functional_overlap = out[
        "geom_func_exact_multiclass_reference_overlap"
    ].fillna(False).astype(bool)

    out["exact_training_reference_overlap"] = (
        global_overlap
        | functional_overlap
    )
    out["novel_vs_training_reference"] = (
        ~out["exact_training_reference_overlap"]
    )

    # --------------------------------------------------------------
    # Global lytic support: 2 transparent independent axes
    # --------------------------------------------------------------
    is_lytic = (
        out["global_prediction"] == "lytic"
    )

    out["global_property_supported"] = (
        out[
            "prop_global_fraction_within_reference_95"
        ]
        >= property_threshold
    )

    out["global_geometry_supported"] = (
        out[
            "geom_global_global_coverage95"
        ]
        .fillna(False)
        .astype(bool)
    )

    out["global_support_count"] = np.where(
        is_lytic,
        out[
            [
                "global_property_supported",
                "global_geometry_supported",
            ]
        ]
        .astype(int)
        .sum(axis=1),
        np.nan,
    )

    out["global_support_class"] = pd.NA

    out.loc[
        is_lytic
        & (out["global_support_count"] == 2),
        "global_support_class",
    ] = "concordant"

    out.loc[
        is_lytic
        & (out["global_support_count"] == 1),
        "global_support_class",
    ] = "mixed"

    out.loc[
        is_lytic
        & (out["global_support_count"] == 0),
        "global_support_class",
    ] = "discordant"

    # --------------------------------------------------------------
    # Functional subtype support: 3 independent axes
    # --------------------------------------------------------------
    out["functional_property_supported"] = (
        out[
            "prop_func_fraction_within_reference_95"
        ]
        >= property_threshold
    )

    out["functional_geometry_coverage_supported"] = (
        out[
            "geom_func_predicted_subtype_coverage95"
        ]
        .fillna(False)
        .astype(bool)
    )

    primary_agreement_cols = [
        "geom_func_agreement_knn_majority",
        "geom_func_agreement_raw_nn",
        "geom_func_agreement_centroid",
    ]

    out["functional_local_geometry_votes"] = (
        out[
            primary_agreement_cols
        ]
        .fillna(False)
        .astype(bool)
        .astype(int)
        .sum(axis=1)
    )

    out["functional_local_geometry_supported"] = (
        out["functional_local_geometry_votes"]
        >= geometry_votes_required
    )

    out["functional_support_count"] = np.where(
        is_lytic,
        out[
            [
                "functional_property_supported",
                "functional_geometry_coverage_supported",
                "functional_local_geometry_supported",
            ]
        ]
        .astype(int)
        .sum(axis=1),
        np.nan,
    )

    out["functional_concordance"] = pd.NA

    out.loc[
        is_lytic
        & (out["functional_support_count"] == 3),
        "functional_concordance",
    ] = "concordant"

    out.loc[
        is_lytic
        & (out["functional_support_count"] == 2),
        "functional_concordance",
    ] = "mixed"

    out.loc[
        is_lytic
        & (out["functional_support_count"] <= 1),
        "functional_concordance",
    ] = "discordant"

    # Useful explicit discordance flags.
    out["prediction_property_discordant"] = (
        is_lytic
        & (~out["functional_property_supported"])
    )

    out["prediction_geometry_discordant"] = (
        is_lytic
        & (
            ~out[
                "functional_geometry_coverage_supported"
            ]
        )
    )

    out["prediction_local_geometry_discordant"] = (
        is_lytic
        & (
            ~out[
                "functional_local_geometry_supported"
            ]
        )
    )

    return out


# ---------------------------------------------------------------------
# Summaries
# ---------------------------------------------------------------------

def build_source_summary(
    master: pd.DataFrame,
) -> pd.DataFrame:
    rows = []

    for source in SOURCE_ORDER:
        g = master[
            (master["source"] == source)
            & (
                master["global_prediction"]
                == "lytic"
            )
            & (
                master[
                    "novel_vs_training_reference"
                ]
            )
        ]

        if g.empty:
            continue

        rows.append(
            {
                "source": source,
                "n_lytic_novel": int(len(g)),
                "median_global_score": float(
                    g["global_score"].median()
                ),
                "global_property_support_fraction": float(
                    g[
                        "global_property_supported"
                    ].mean()
                ),
                "global_geometry_support_fraction": float(
                    g[
                        "global_geometry_supported"
                    ].mean()
                ),
                "global_concordant_fraction": float(
                    (
                        g[
                            "global_support_class"
                        ]
                        == "concordant"
                    ).mean()
                ),
                "functional_property_support_fraction": float(
                    g[
                        "functional_property_supported"
                    ].mean()
                ),
                "functional_geometry_coverage_fraction": float(
                    g[
                        "functional_geometry_coverage_supported"
                    ].mean()
                ),
                "functional_local_geometry_support_fraction": float(
                    g[
                        "functional_local_geometry_supported"
                    ].mean()
                ),
                "functional_concordant_fraction": float(
                    (
                        g[
                            "functional_concordance"
                        ]
                        == "concordant"
                    ).mean()
                ),
                "functional_mixed_fraction": float(
                    (
                        g[
                            "functional_concordance"
                        ]
                        == "mixed"
                    ).mean()
                ),
                "functional_discordant_fraction": float(
                    (
                        g[
                            "functional_concordance"
                        ]
                        == "discordant"
                    ).mean()
                ),
            }
        )

    return pd.DataFrame(rows)


def build_source_class_summary(
    master: pd.DataFrame,
) -> pd.DataFrame:
    rows = []

    for source in SOURCE_ORDER:
        for cls in FUNCTION_ORDER:
            g = master[
                (master["source"] == source)
                & (
                    master[
                        "functional_class_a4"
                    ]
                    == cls
                )
                & (
                    master["global_prediction"]
                    == "lytic"
                )
                & (
                    master[
                        "novel_vs_training_reference"
                    ]
                )
            ]

            if g.empty:
                continue

            rows.append(
                {
                    "source": source,
                    "predicted_class": cls,
                    "n": int(len(g)),
                    "median_global_score": float(
                        g["global_score"].median()
                    ),
                    "median_functional_max_probability": float(
                        g[
                            "functional_max_probability"
                        ].median()
                    ),
                    "median_functional_probability_margin": float(
                        g[
                            "functional_probability_margin"
                        ].median()
                    ),
                    "median_property_fit_95": float(
                        g[
                            "prop_func_fraction_within_reference_95"
                        ].median()
                    ),
                    "property_supported_fraction": float(
                        g[
                            "functional_property_supported"
                        ].mean()
                    ),
                    "geometry_coverage95_fraction": float(
                        g[
                            "functional_geometry_coverage_supported"
                        ].mean()
                    ),
                    "local_geometry_supported_fraction": float(
                        g[
                            "functional_local_geometry_supported"
                        ].mean()
                    ),
                    "median_predicted_subtype_nn_percentile": float(
                        g[
                            "geom_func_predicted_subtype_nn_percentile"
                        ].median()
                    ),
                    "median_predicted_subtype_knn_fraction": float(
                        g[
                            "geom_func_predicted_subtype_knn_fraction"
                        ].median()
                    ),
                    "concordant_fraction": float(
                        (
                            g[
                                "functional_concordance"
                            ]
                            == "concordant"
                        ).mean()
                    ),
                    "mixed_fraction": float(
                        (
                            g[
                                "functional_concordance"
                            ]
                            == "mixed"
                        ).mean()
                    ),
                    "discordant_fraction": float(
                        (
                            g[
                                "functional_concordance"
                            ]
                            == "discordant"
                        ).mean()
                    ),
                }
            )

    return pd.DataFrame(rows)


def build_correlations(
    master: pd.DataFrame,
) -> pd.DataFrame:
    rows = []

    specs = [
        (
            "global_score_vs_global_property_fit",
            "global_score",
            "prop_global_fraction_within_reference_95",
        ),
        (
            "global_score_vs_global_nn_percentile",
            "global_score",
            "geom_global_global_nn_distance_percentile",
        ),
        (
            "functional_probability_vs_property_fit",
            "functional_max_probability",
            "prop_func_fraction_within_reference_95",
        ),
        (
            "functional_probability_vs_subtype_nn_percentile",
            "functional_max_probability",
            "geom_func_predicted_subtype_nn_percentile",
        ),
        (
            "functional_probability_vs_subtype_knn_fraction",
            "functional_max_probability",
            "geom_func_predicted_subtype_knn_fraction",
        ),
    ]

    for source in SOURCE_ORDER:
        for cls in [None, *FUNCTION_ORDER]:
            g = master[
                (master["source"] == source)
                & (
                    master["global_prediction"]
                    == "lytic"
                )
                & (
                    master[
                        "novel_vs_training_reference"
                    ]
                )
            ]

            if cls is not None:
                g = g[
                    g[
                        "functional_class_a4"
                    ]
                    == cls
                ]

            if g.empty:
                continue

            for name, xcol, ycol in specs:
                rho, p, n = safe_spearman(
                    g[xcol],
                    g[ycol],
                )

                rows.append(
                    {
                        "source": source,
                        "predicted_class": (
                            "ALL"
                            if cls is None
                            else cls
                        ),
                        "relationship": name,
                        "n": n,
                        "spearman_rho": rho,
                        "p_value": p,
                    }
                )

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------
# Plot data
# ---------------------------------------------------------------------

def support_heatmap_matrix(
    source_class_summary: pd.DataFrame,
) -> pd.DataFrame:
    df = source_class_summary.copy()

    df["row"] = (
        df["source"].map(SOURCE_LABELS)
        + " — "
        + df["predicted_class"].map(
            FUNCTION_LABELS
        )
    )

    metrics = {
        "Property\nsupport": "property_supported_fraction",
        "Subtype\ncoverage": "geometry_coverage95_fraction",
        "Local geometry\nagreement": "local_geometry_supported_fraction",
        "Full\nconcordance": "concordant_fraction",
    }

    rows = []

    for _, r in df.iterrows():
        for label, col in metrics.items():
            rows.append(
                {
                    "row": r["row"],
                    "metric": label,
                    "value": 100.0 * r[col],
                }
            )

    long = pd.DataFrame(rows)

    matrix = long.pivot(
        index="row",
        columns="metric",
        values="value",
    )

    ordered_rows = []

    for source in SOURCE_ORDER:
        for cls in FUNCTION_ORDER:
            label = (
                f"{SOURCE_LABELS[source]} — "
                f"{FUNCTION_LABELS[cls]}"
            )
            if label in matrix.index:
                ordered_rows.append(label)

    ordered_cols = list(metrics.keys())

    return matrix.reindex(
        index=ordered_rows,
        columns=ordered_cols,
    )


# ---------------------------------------------------------------------
# Standalone figures
# ---------------------------------------------------------------------

def plot_concordance_heatmap(
    source_class_summary: pd.DataFrame,
    figures: Path,
) -> None:
    matrix = support_heatmap_matrix(
        source_class_summary
    )

    fig, ax = plt.subplots(
        figsize=(8.0, max(5.5, 0.55 * len(matrix)))
    )

    sns.heatmap(
        matrix,
        vmin=0,
        vmax=100,
        annot=True,
        fmt=".1f",
        linewidths=0.25,
        cbar_kws={
            "label": "Sequences (%)"
        },
        ax=ax,
    )

    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title(
        "Independent support for PhageLysClass functional predictions"
    )
    ax.tick_params(
        axis="x",
        rotation=0,
    )
    ax.tick_params(
        axis="y",
        rotation=0,
    )

    fig.tight_layout()

    save_figure(
        fig,
        figures / "G1_prediction_evidence_concordance",
    )


def plot_concordance_by_source(
    master: pd.DataFrame,
    figures: Path,
) -> None:
    lytic = master[
        (master["global_prediction"] == "lytic")
        & (
            master[
                "novel_vs_training_reference"
            ]
        )
    ].copy()

    tab = (
        lytic.groupby(
            [
                "source",
                "functional_concordance",
            ]
        )
        .size()
        .rename("n")
        .reset_index()
    )

    tab["fraction"] = (
        tab["n"]
        / tab.groupby("source")[
            "n"
        ].transform("sum")
    )
    tab["percent"] = (
        100.0 * tab["fraction"]
    )
    tab["source_label"] = tab[
        "source"
    ].map(SOURCE_LABELS)

    fig, ax = plt.subplots(
        figsize=(7.5, 4.8)
    )

    sns.barplot(
        data=tab,
        x="source_label",
        y="percent",
        hue="functional_concordance",
        order=[
            SOURCE_LABELS[x]
            for x in SOURCE_ORDER
        ],
        hue_order=CONCORDANCE_ORDER,
        errorbar=None,
        ax=ax,
    )

    ax.set_xlabel("")
    ax.set_ylabel(
        "Predicted lytic sequences (%)"
    )
    ax.set_ylim(0, 100)
    ax.set_title(
        "Prediction–evidence concordance"
    )
    ax.legend(
        title="",
        frameon=False,
    )
    ax.tick_params(
        axis="x",
        rotation=15,
    )
    sns.despine(ax=ax)

    fig.tight_layout()

    save_figure(
        fig,
        figures / "G2_concordance_by_source",
    )


def plot_geometry_matrix(
    frames: dict[str, pd.DataFrame],
    figures: Path,
    max_reference: int,
    max_external: int,
    seed: int,
) -> None:
    pca_ref = frames[
        "pca_reference"
    ].copy()
    pca_ext = frames[
        "pca_external"
    ].copy()

    umap_ref = frames[
        "umap_reference"
    ].copy()
    umap_ext = frames[
        "umap_external"
    ].copy()

    tsne = frames["tsne"].copy()

    fig, axes = plt.subplots(
        3,
        3,
        figsize=(18, 16),
    )

    methods = [
        (
            "PCA",
            pca_ref,
            pca_ext,
            "PC1",
            "PC2",
        ),
        (
            "UMAP",
            umap_ref,
            umap_ext,
            "UMAP1",
            "UMAP2",
        ),
    ]

    legend_handles = None
    legend_labels = None

    for row_idx, source in enumerate(
        SOURCE_ORDER
    ):
        for col_idx, (
            method,
            ref,
            ext,
            xcol,
            ycol,
        ) in enumerate(methods):
            ax = axes[
                row_idx,
                col_idx,
            ]

            ref_plot = stratified_sample(
                ref,
                group_col="reference_class",
                max_total=max_reference,
                seed=seed + row_idx,
            ).copy()

            ext_plot = ext[
                (ext["source"] == source)
                & (
                    ext["global_prediction"]
                    == "lytic"
                )
            ].copy()

            if len(ext_plot) > max_external:
                ext_plot = stratified_sample(
                    ext_plot,
                    group_col="functional_class_a2",
                    max_total=max_external,
                    seed=seed + 100 + row_idx,
                )

            ref_plot["origin"] = "Reference"
            ref_plot["functional_class"] = (
                ref_plot["reference_class"]
            )

            ext_plot["origin"] = "External"
            ext_plot["functional_class"] = (
                ext_plot[
                    "functional_class_a2"
                ]
            )

            combo = pd.concat(
                [
                    ref_plot[
                        [
                            xcol,
                            ycol,
                            "origin",
                            "functional_class",
                        ]
                    ],
                    ext_plot[
                        [
                            xcol,
                            ycol,
                            "origin",
                            "functional_class",
                        ]
                    ],
                ],
                ignore_index=True,
            )

            combo["functional_label"] = (
                combo[
                    "functional_class"
                ].map(FUNCTION_LABELS)
            )

            sns.scatterplot(
                data=combo,
                x=xcol,
                y=ycol,
                hue="functional_label",
                style="origin",
                alpha=0.45,
                s=16,
                linewidth=0,
                legend=(
                    "full"
                    if row_idx == 0
                    and col_idx == 0
                    else False
                ),
                ax=ax,
            )

            ax.set_title(
                f"{SOURCE_LABELS[source]} — {method}"
            )
            ax.set_xlabel(xcol)
            ax.set_ylabel(ycol)
            sns.despine(ax=ax)

            if (
                row_idx == 0
                and col_idx == 0
            ):
                legend_handles, legend_labels = (
                    ax.get_legend_handles_labels()
                )
                legend = ax.get_legend()
                if legend is not None:
                    legend.remove()

        # t-SNE column
        ax = axes[row_idx, 2]

        t = tsne[
            tsne["tsne_source"]
            == source
        ].copy()

        t["functional_label"] = (
            t["functional_class"].map(
                FUNCTION_LABELS
            )
        )

        sns.scatterplot(
            data=t,
            x="TSNE1",
            y="TSNE2",
            hue="functional_label",
            style="origin",
            alpha=0.50,
            s=16,
            linewidth=0,
            legend=False,
            ax=ax,
        )

        ax.set_title(
            f"{SOURCE_LABELS[source]} — t-SNE"
        )
        ax.set_xlabel("t-SNE1")
        ax.set_ylabel("t-SNE2")
        sns.despine(ax=ax)

    if legend_handles is not None:
        fig.legend(
            legend_handles,
            legend_labels,
            loc="upper center",
            bbox_to_anchor=(0.5, 1.01),
            ncol=5,
            frameon=False,
        )

    fig.suptitle(
        "Reference geometry and external PhageLysClass predictions",
        y=1.025,
        fontsize=16,
    )

    fig.tight_layout()

    save_figure(
        fig,
        figures
        / "G3_geometry_PCA_UMAP_tSNE_matrix",
    )


# ---------------------------------------------------------------------
# Full publication-candidate figure
# ---------------------------------------------------------------------

def draw_source_overview(
    ax: plt.Axes,
    master: pd.DataFrame,
) -> None:
    tab = (
        master.groupby(
            [
                "source",
                "global_prediction",
            ]
        )
        .size()
        .rename("n")
        .reset_index()
    )

    totals = (
        tab.groupby("source")["n"]
        .sum()
        .rename("total")
        .reset_index()
    )

    lytic = tab[
        tab["global_prediction"] == "lytic"
    ].merge(
        totals,
        on="source",
        how="left",
    )

    lytic["percent"] = (
        100.0
        * lytic["n"]
        / lytic["total"]
    )
    lytic["source_label"] = (
        lytic["source"].map(
            SOURCE_LABELS
        )
    )

    sns.barplot(
        data=lytic,
        x="source_label",
        y="percent",
        order=[
            SOURCE_LABELS[x]
            for x in SOURCE_ORDER
        ],
        errorbar=None,
        ax=ax,
    )

    for i, source in enumerate(
        SOURCE_ORDER
    ):
        row = lytic[
            lytic["source"] == source
        ].iloc[0]

        ax.text(
            i,
            row["percent"] + 2,
            f"n={int(row['total']):,}",
            ha="center",
            va="bottom",
            fontsize=8,
        )

    ax.set_ylim(0, 100)
    ax.set_xlabel("")
    ax.set_ylabel("Predicted lytic (%)")
    ax.set_title(
        "Global prediction landscape"
    )
    ax.tick_params(
        axis="x",
        rotation=15,
    )
    sns.despine(ax=ax)


def draw_global_scores(
    ax: plt.Axes,
    master: pd.DataFrame,
) -> None:
    x = master.copy()

    x["source_label"] = (
        x["source"].map(
            SOURCE_LABELS
        )
    )

    sns.boxplot(
        data=x,
        x="source_label",
        y="global_score",
        order=[
            SOURCE_LABELS[s]
            for s in SOURCE_ORDER
        ],
        showfliers=False,
        ax=ax,
    )

    ax.set_ylim(0, 1)
    ax.set_xlabel("")
    ax.set_ylabel("Global lytic score")
    ax.set_title(
        "Global-score distributions"
    )
    ax.tick_params(
        axis="x",
        rotation=15,
    )
    sns.despine(ax=ax)


def draw_functional_composition(
    ax: plt.Axes,
    master: pd.DataFrame,
) -> None:
    lytic = master[
        master["global_prediction"] == "lytic"
    ].copy()

    tab = (
        lytic.groupby(
            [
                "source",
                "functional_class_a4",
            ]
        )
        .size()
        .rename("n")
        .reset_index()
    )

    tab["percent"] = (
        100.0
        * tab["n"]
        / tab.groupby("source")[
            "n"
        ].transform("sum")
    )

    tab["source_label"] = (
        tab["source"].map(
            SOURCE_LABELS
        )
    )
    tab["functional_label"] = (
        tab[
            "functional_class_a4"
        ].map(FUNCTION_LABELS)
    )

    sns.barplot(
        data=tab,
        x="source_label",
        y="percent",
        hue="functional_label",
        order=[
            SOURCE_LABELS[s]
            for s in SOURCE_ORDER
        ],
        hue_order=[
            FUNCTION_LABELS[c]
            for c in FUNCTION_ORDER
        ],
        errorbar=None,
        ax=ax,
    )

    ax.set_ylim(0, 100)
    ax.set_xlabel("")
    ax.set_ylabel(
        "Predicted lytic sequences (%)"
    )
    ax.set_title(
        "Functional predictions"
    )
    ax.legend(
        title="",
        frameon=False,
        fontsize=8,
    )
    ax.tick_params(
        axis="x",
        rotation=15,
    )
    sns.despine(ax=ax)


def draw_zymctrl_ec(
    ax: plt.Axes,
    a1_dir: Path,
) -> None:
    path = optional_file(
        a1_dir
        / "D1_zymctrl_by_ec.tsv"
    )

    if path is None:
        ax.text(
            0.5,
            0.5,
            "ZymCTRL EC table unavailable",
            ha="center",
            va="center",
        )
        ax.axis("off")
        return

    df = pd.read_csv(
        path,
        sep="\t",
    ).copy()

    df["percent"] = (
        100.0 * df["fraction_lytic"]
    )

    df = df.sort_values(
        "percent"
    )

    sns.barplot(
        data=df,
        y="ec_number",
        x="percent",
        order=df["ec_number"].tolist(),
        errorbar=None,
        ax=ax,
    )

    ax.set_xlim(0, 100)
    ax.set_xlabel("Predicted lytic (%)")
    ax.set_ylabel("EC number")
    ax.set_title(
        "ZymCTRL by EC"
    )
    sns.despine(ax=ax)


def draw_ambiguous_categories(
    ax: plt.Axes,
    a1_dir: Path,
) -> None:
    path = optional_file(
        a1_dir
        / "D4_ambiguous_structural_by_category.tsv"
    )

    if path is None:
        ax.text(
            0.5,
            0.5,
            "Structural-category table unavailable",
            ha="center",
            va="center",
        )
        ax.axis("off")
        return

    df = pd.read_csv(
        path,
        sep="\t",
    ).copy()

    df["percent"] = (
        100.0 * df["fraction_lytic"]
    )

    df = df.sort_values(
        "percent"
    )

    sns.barplot(
        data=df,
        y="structural_category",
        x="percent",
        order=df[
            "structural_category"
        ].tolist(),
        errorbar=None,
        ax=ax,
    )

    ax.set_xlim(0, 100)
    ax.set_xlabel("Predicted lytic (%)")
    ax.set_ylabel("")
    ax.set_title(
        "Ambiguous structural challenge"
    )
    sns.despine(ax=ax)


def draw_property_fit_heatmap(
    ax: plt.Axes,
    source_class_summary: pd.DataFrame,
) -> None:
    df = source_class_summary.copy()

    df["row"] = (
        df["source"].map(
            SOURCE_LABELS
        )
        + " — "
        + df["predicted_class"].map(
            FUNCTION_LABELS
        )
    )

    matrix = (
        df.set_index("row")[
            [
                "median_property_fit_95",
            ]
        ]
        * 100.0
    )

    ordered = []

    for source in SOURCE_ORDER:
        for cls in FUNCTION_ORDER:
            label = (
                f"{SOURCE_LABELS[source]} — "
                f"{FUNCTION_LABELS[cls]}"
            )
            if label in matrix.index:
                ordered.append(label)

    matrix = matrix.reindex(
        ordered
    )
    matrix.columns = [
        "Descriptors within\nreference 95% range"
    ]

    sns.heatmap(
        matrix,
        vmin=0,
        vmax=100,
        annot=True,
        fmt=".1f",
        linewidths=0.25,
        cbar=False,
        ax=ax,
    )

    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title(
        "Property compatibility"
    )
    ax.tick_params(
        axis="y",
        labelsize=7,
        rotation=0,
    )


def draw_tsne_panel(
    ax: plt.Axes,
    tsne: pd.DataFrame,
    source: str,
) -> None:
    t = tsne[
        tsne["tsne_source"] == source
    ].copy()

    t["functional_label"] = (
        t["functional_class"].map(
            FUNCTION_LABELS
        )
    )

    sns.scatterplot(
        data=t,
        x="TSNE1",
        y="TSNE2",
        hue="functional_label",
        style="origin",
        alpha=0.50,
        s=13,
        linewidth=0,
        legend=False,
        ax=ax,
    )

    ax.set_title(
        SOURCE_LABELS[source]
    )
    ax.set_xlabel("t-SNE1")
    ax.set_ylabel("t-SNE2")
    sns.despine(ax=ax)


def draw_concordance_heatmap(
    ax: plt.Axes,
    source_class_summary: pd.DataFrame,
) -> None:
    matrix = support_heatmap_matrix(
        source_class_summary
    )

    sns.heatmap(
        matrix,
        vmin=0,
        vmax=100,
        annot=True,
        fmt=".1f",
        linewidths=0.25,
        cbar_kws={
            "label": "Sequences (%)"
        },
        ax=ax,
    )

    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title(
        "Prediction × property × geometry concordance"
    )
    ax.tick_params(
        axis="x",
        rotation=0,
    )
    ax.tick_params(
        axis="y",
        rotation=0,
        labelsize=8,
    )


def build_full_figure(
    master: pd.DataFrame,
    source_class_summary: pd.DataFrame,
    frames: dict[str, pd.DataFrame],
    a1_dir: Path,
    figures: Path,
) -> None:
    fig = plt.figure(
        figsize=(21, 25)
    )

    gs = fig.add_gridspec(
        4,
        3,
        height_ratios=[
            1.0,
            1.15,
            1.35,
            1.45,
        ],
        hspace=0.48,
        wspace=0.38,
    )

    ax_a = fig.add_subplot(
        gs[0, 0]
    )
    ax_b = fig.add_subplot(
        gs[0, 1]
    )
    ax_c = fig.add_subplot(
        gs[0, 2]
    )

    draw_source_overview(
        ax_a,
        master,
    )
    draw_global_scores(
        ax_b,
        master,
    )
    draw_functional_composition(
        ax_c,
        master,
    )

    panel_label(
        ax_a,
        "A",
    )
    panel_label(
        ax_b,
        "B",
    )
    panel_label(
        ax_c,
        "C",
    )

    ax_d1 = fig.add_subplot(
        gs[1, 0]
    )
    ax_d2 = fig.add_subplot(
        gs[1, 1]
    )
    ax_e = fig.add_subplot(
        gs[1, 2]
    )

    draw_zymctrl_ec(
        ax_d1,
        a1_dir,
    )
    draw_ambiguous_categories(
        ax_d2,
        a1_dir,
    )
    draw_property_fit_heatmap(
        ax_e,
        source_class_summary,
    )

    panel_label(
        ax_d1,
        "D1",
    )
    panel_label(
        ax_d2,
        "D2",
    )
    panel_label(
        ax_e,
        "E",
    )

    tsne = frames["tsne"]

    for j, source in enumerate(
        SOURCE_ORDER
    ):
        ax = fig.add_subplot(
            gs[2, j]
        )

        draw_tsne_panel(
            ax,
            tsne,
            source,
        )

        panel_label(
            ax,
            f"F{j + 1}",
        )

    ax_g = fig.add_subplot(
        gs[3, :]
    )

    draw_concordance_heatmap(
        ax_g,
        source_class_summary,
    )

    panel_label(
        ax_g,
        "G",
    )

    fig.suptitle(
        "Application of PhageLysClass to heterogeneous external sequence spaces",
        fontsize=18,
        y=0.995,
    )

    save_figure(
        fig,
        figures
        / "Figure_final_external_application",
    )


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def main() -> None:
    args = parse_args()

    if not (
        0.0
        <= args.property_support_threshold
        <= 1.0
    ):
        raise ValueError(
            "--property-support-threshold must be between 0 and 1."
        )

    outdir = args.output_dir.resolve()
    tables = outdir / "tables"
    figures = outdir / "figures"

    tables.mkdir(
        parents=True,
        exist_ok=True,
    )
    figures.mkdir(
        parents=True,
        exist_ok=True,
    )

    frames = load_inputs(
        args
    )

    print(
        "A4: integrating predictions, properties and geometry...",
        flush=True,
    )

    master = build_integrated_master(
        frames,
        property_threshold=(
            args.property_support_threshold
        ),
        geometry_votes_required=(
            args.geometry_votes_required
        ),
    )

    final_master_path = (
        tables
        / "external_prediction_property_geometry_master.tsv.gz"
    )

    master.to_csv(
        final_master_path,
        sep="\t",
        index=False,
        compression="gzip",
    )

    source_summary = (
        build_source_summary(
            master
        )
    )
    source_class_summary = (
        build_source_class_summary(
            master
        )
    )
    correlations = (
        build_correlations(
            master
        )
    )

    source_summary.to_csv(
        tables
        / "integration_summary_by_source.tsv",
        sep="\t",
        index=False,
    )

    source_class_summary.to_csv(
        tables
        / "integration_summary_by_source_class.tsv",
        sep="\t",
        index=False,
    )

    correlations.to_csv(
        tables
        / "prediction_support_correlations.tsv",
        sep="\t",
        index=False,
    )

    lytic_novel = master[
        (master["global_prediction"] == "lytic")
        & (
            master[
                "novel_vs_training_reference"
            ]
        )
    ].copy()

    concordant = lytic_novel[
        lytic_novel[
            "functional_concordance"
        ]
        == "concordant"
    ].copy()

    discordant = lytic_novel[
        lytic_novel[
            "functional_concordance"
        ]
        == "discordant"
    ].copy()

    concordant = concordant.sort_values(
        [
            "global_score",
            "functional_max_probability",
        ],
        ascending=False,
    )

    discordant = discordant.sort_values(
        [
            "global_score",
            "functional_max_probability",
        ],
        ascending=False,
    )

    concordant.to_csv(
        tables
        / "concordant_predictions.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    discordant.to_csv(
        tables
        / "discordant_predictions.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    overlap = master[
        master[
            "exact_training_reference_overlap"
        ]
    ].copy()

    overlap.to_csv(
        tables
        / "exact_training_reference_overlap_audit.tsv",
        sep="\t",
        index=False,
    )

    print(
        "A4: concordance figures...",
        flush=True,
    )

    plot_concordance_heatmap(
        source_class_summary,
        figures,
    )

    plot_concordance_by_source(
        master,
        figures,
    )

    print(
        "A4: PCA / UMAP / t-SNE geometry matrix...",
        flush=True,
    )

    plot_geometry_matrix(
        frames,
        figures,
        max_reference=args.plot_max_reference,
        max_external=(
            args.plot_max_external_per_source
        ),
        seed=args.seed,
    )

    print(
        "A4: final multipanel figure...",
        flush=True,
    )

    build_full_figure(
        master,
        source_class_summary,
        frames,
        args.a1_dir.resolve(),
        figures,
    )

    manifest = {
        "stage": "A4",
        "stage_name": (
            "prediction_geometry_integration"
        ),
        "generated_utc": datetime.now(
            timezone.utc
        ).isoformat(),
        "n_sequences": int(len(master)),
        "n_lytic": int(
            (
                master[
                    "global_prediction"
                ]
                == "lytic"
            ).sum()
        ),
        "n_exact_training_reference_overlap": int(
            master[
                "exact_training_reference_overlap"
            ].sum()
        ),
        "concordance_contract": {
            "property_support_threshold": float(
                args.property_support_threshold
            ),
            "geometry_votes_required": int(
                args.geometry_votes_required
            ),
            "primary_geometry_rules": [
                "kNN majority subtype",
                "raw nearest-neighbor subtype",
                "nearest centroid subtype",
            ],
            "density_normalized_nn": (
                "Retained as a diagnostic only; excluded from the "
                "primary local-concordance vote."
            ),
            "functional_concordance": {
                "3": "concordant",
                "2": "mixed",
                "0-1": "discordant",
            },
        },
        "quantitative_geometry": (
            "Full 26D standardized descriptor space from A3."
        ),
        "visualization_geometry": {
            "PCA": (
                "Reference-fitted; external transform."
            ),
            "UMAP": (
                "Reference-fitted; external transform."
            ),
            "t-SNE": (
                "Joint stratified visualization sample only; "
                "not used quantitatively."
            ),
        },
        "outputs": {
            "master": str(
                final_master_path
            ),
            "source_summary": str(
                tables
                / "integration_summary_by_source.tsv"
            ),
            "source_class_summary": str(
                tables
                / "integration_summary_by_source_class.tsv"
            ),
            "correlations": str(
                tables
                / "prediction_support_correlations.tsv"
            ),
            "concordant_predictions": str(
                tables
                / "concordant_predictions.tsv.gz"
            ),
            "discordant_predictions": str(
                tables
                / "discordant_predictions.tsv.gz"
            ),
            "full_figure": str(
                figures
                / "Figure_final_external_application"
            ),
            "geometry_matrix": str(
                figures
                / "G3_geometry_PCA_UMAP_tSNE_matrix"
            ),
        },
        "status": "PASS",
    }

    (
        outdir / "manifest.json"
    ).write_text(
        json.dumps(
            manifest,
            indent=2,
        ),
        encoding="utf-8",
    )

    print("\nA4 — PASS")
    print("=" * 118)
    print(
        f"Integrated sequences:                  "
        f"{len(master):,}"
    )
    print(
        f"Predicted lytic:                      "
        f"{(master['global_prediction'] == 'lytic').sum():,}"
    )
    print(
        f"Exact training-reference overlaps:    "
        f"{master['exact_training_reference_overlap'].sum():,}"
    )
    print(
        f"Concordant functional predictions:    "
        f"{len(concordant):,}"
    )
    print(
        f"Discordant functional predictions:    "
        f"{len(discordant):,}"
    )
    print("=" * 118)
    print("\nSOURCE-LEVEL INTEGRATION")
    print(
        source_summary.to_string(
            index=False
        )
    )
    print("=" * 118)
    print("\nSOURCE × FUNCTIONAL CLASS")
    print(
        source_class_summary.to_string(
            index=False
        )
    )
    print("=" * 118)
    print(f"Master:  {final_master_path}")
    print(f"Tables:  {tables}")
    print(f"Figures: {figures}")
    print(f"Output:  {outdir}")


if __name__ == "__main__":
    main()

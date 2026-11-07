#!/usr/bin/env python3
"""
A3 — Reference geometry for explored PhageLysClass sequences.

Scientific objective
--------------------
Determine whether external PhageLysClass predictions occupy the multivariate
property geometry of annotated training references.

This stage uses the exact 26 frozen Roxy descriptors calculated in A2.
No feature selection, model fitting, threshold optimization, or joint
reference+external manifold fitting is performed.

Geometry contract
-----------------
1. GLOBAL POSITIVE GEOMETRY
   Reference:
       eligible positive class from all-data global_mixed training universe
       already materialized by A2 as reference_cohort == global_positive

   Scaling:
       StandardScaler fitted on the global positive reference only

   Quantities:
       - full 26D standardized nearest-reference distance
       - mean kNN distance (k=25 by default)
       - empirical distance percentile relative to reference leave-one-out
         nearest-neighbor distances
       - 95% and 99% reference-coverage calls

2. FUNCTIONAL GEOMETRY
   Reference:
       multiclass_clean Endolysin / Depolymerase / VAL training classes

   Scaling:
       one StandardScaler fitted on the complete multiclass reference only

   Quantities:
       - nearest neighbor distance to each subtype
       - distance to each subtype centroid
       - subtype-specific distance normalized by that subtype's own reference
         nearest-neighbor density
       - nearest subtype under raw NN, density-normalized NN and centroid rules
       - kNN composition over the complete multiclass reference
       - neighborhood entropy
       - prediction-vs-neighborhood agreement
       - class-conditional 95% / 99% coverage and OOD status

3. MANIFOLD VISUALIZATION
   PCA:
       fitted only on the standardized multiclass reference;
       external sequences transformed out-of-sample.

   UMAP:
       fitted only on the standardized multiclass reference;
       external sequences transformed out-of-sample with UMAP.transform().

   t-SNE:
       no out-of-sample projection is claimed. A stratified reference sample
       and an external-source sample are embedded jointly for visualization
       only. t-SNE coordinates are never used for quantitative geometry.

   All quantitative distances remain calculated in the full 26-dimensional
   standardized descriptor space.

Default execution
-----------------
python analysis_explored_seqs/A3_reference_geometry.py

Inputs
------
analysis_explored_seqs/A2_property_characterization/tables/
    external_frozen26_descriptors.tsv.gz
    reference_frozen26_descriptors.tsv.gz
    frozen26_descriptor_registry.tsv

Outputs
-------
analysis_explored_seqs/A3_reference_geometry/
├── tables/
│   ├── geometry_feature_contract.tsv
│   ├── global_reference_self_distance_profile.tsv
│   ├── subtype_reference_self_distance_profiles.tsv
│   ├── global_reference_coverage.tsv.gz
│   ├── functional_reference_geometry.tsv.gz
│   ├── pca_reference_coordinates.tsv.gz
│   ├── pca_external_coordinates.tsv.gz
│   ├── global_geometry_summary.tsv
│   ├── functional_geometry_summary.tsv
│   └── neighborhood_agreement_summary.tsv
├── figures/
│   ├── F1a_pca_prediction_extension.{png,pdf}
│   ├── F1b_pca_zymctrl.{png,pdf}
│   ├── F1c_pca_ambiguous_structural.{png,pdf}
│   ├── F1d_umap_prediction_extension.{png,pdf}
│   ├── F1e_umap_zymctrl.{png,pdf}
│   ├── F1f_umap_ambiguous_structural.{png,pdf}
│   ├── F1g_tsne_prediction_extension.{png,pdf}
│   ├── F1h_tsne_zymctrl.{png,pdf}
│   ├── F1i_tsne_ambiguous_structural.{png,pdf}
│   ├── F2_reference_coverage_by_source.{png,pdf}
│   ├── F3_subtype_geometry_agreement.{png,pdf}
│   └── F4_predicted_subtype_distance_percentile.{png,pdf}
└── manifest.json
"""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler

try:
    import umap
except ImportError:
    umap = None


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

sns.set_theme(style="whitegrid", context="paper")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="A3 reference geometry.")

    p.add_argument(
        "--external-descriptors",
        type=Path,
        default=Path(
            "analysis_explored_seqs/"
            "A2_property_characterization/"
            "tables/"
            "external_frozen26_descriptors.tsv.gz"
        ),
    )
    p.add_argument(
        "--reference-descriptors",
        type=Path,
        default=Path(
            "analysis_explored_seqs/"
            "A2_property_characterization/"
            "tables/"
            "reference_frozen26_descriptors.tsv.gz"
        ),
    )
    p.add_argument(
        "--descriptor-registry",
        type=Path,
        default=Path(
            "analysis_explored_seqs/"
            "A2_property_characterization/"
            "tables/"
            "frozen26_descriptor_registry.tsv"
        ),
    )
    p.add_argument(
        "--output-dir",
        type=Path,
        default=Path("analysis_explored_seqs/A3_reference_geometry"),
    )
    p.add_argument(
        "--k",
        type=int,
        default=25,
        help="Number of neighbors for local-reference composition.",
    )
    p.add_argument(
        "--pca-components",
        type=int,
        default=10,
        help="Number of reference-fitted PCA components to retain.",
    )
    p.add_argument(
        "--plot-max-reference",
        type=int,
        default=5000,
    )
    p.add_argument(
        "--plot-max-external-per-source",
        type=int,
        default=5000,
    )

    p.add_argument(
        "--umap-n-neighbors",
        type=int,
        default=30,
    )
    p.add_argument(
        "--umap-min-dist",
        type=float,
        default=0.2,
    )
    p.add_argument(
        "--tsne-perplexity",
        type=float,
        default=30.0,
    )
    p.add_argument(
        "--tsne-max-reference",
        type=int,
        default=3000,
    )
    p.add_argument(
        "--tsne-max-external-per-source",
        type=int,
        default=3000,
    )
    p.add_argument(
        "--seed",
        type=int,
        default=42,
    )

    return p.parse_args()


def require_file(path: Path) -> Path:
    if not path.is_file():
        raise FileNotFoundError(f"Required file not found: {path}")
    return path


def save_figure(fig: plt.Figure, stem: Path) -> None:
    fig.tight_layout()
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def empirical_percentile(
    values: np.ndarray,
    reference_distribution: np.ndarray,
) -> np.ndarray:
    ref = np.asarray(reference_distribution, dtype=float)
    ref = np.sort(ref[np.isfinite(ref)])

    x = np.asarray(values, dtype=float)

    if len(ref) == 0:
        return np.full(len(x), np.nan)

    return (
        np.searchsorted(ref, x, side="right")
        / len(ref)
    )


def entropy3(probabilities: np.ndarray) -> np.ndarray:
    p = np.asarray(probabilities, dtype=float)
    p = np.clip(p, 0.0, 1.0)

    with np.errstate(divide="ignore", invalid="ignore"):
        terms = np.where(p > 0, p * np.log(p), 0.0)

    h = -np.sum(terms, axis=1)
    return h / np.log(3.0)


def second_smallest_margin(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Returns smallest distance and the gap between the second-smallest and
    smallest value. Larger gap = more geometrically unambiguous.
    """
    ordered = np.sort(matrix, axis=1)
    return ordered[:, 0], ordered[:, 1] - ordered[:, 0]


def stratified_sample(
    df: pd.DataFrame,
    group_col: str,
    max_total: int,
    seed: int,
) -> pd.DataFrame:
    if len(df) <= max_total:
        return df.copy()

    groups = list(df[group_col].dropna().unique())
    if not groups:
        return df.sample(n=max_total, random_state=seed)

    counts = df[group_col].value_counts()
    allocations = {}

    for group in groups:
        frac = counts[group] / counts.sum()
        allocations[group] = max(1, int(round(max_total * frac)))

    # Keep total <= max_total.
    while sum(allocations.values()) > max_total:
        largest = max(allocations, key=allocations.get)
        if allocations[largest] > 1:
            allocations[largest] -= 1
        else:
            break

    pieces = []
    for i, group in enumerate(groups):
        g = df[df[group_col] == group]
        n = min(len(g), allocations[group])
        pieces.append(
            g.sample(
                n=n,
                random_state=seed + i,
            )
        )

    return pd.concat(pieces, ignore_index=True)


def load_feature_contract(
    registry_path: Path,
    reference: pd.DataFrame,
    external: pd.DataFrame,
) -> list[str]:
    registry = pd.read_csv(registry_path, sep="\t")

    if "feature" not in registry.columns:
        raise ValueError("Descriptor registry lacks 'feature' column.")

    features = registry["feature"].astype(str).tolist()

    if len(features) != 26:
        raise ValueError(
            f"Expected 26 frozen descriptors, found {len(features)}."
        )

    for name, frame in [
        ("reference", reference),
        ("external", external),
    ]:
        missing = [f for f in features if f not in frame.columns]
        if missing:
            raise ValueError(
                f"{name} table missing frozen descriptors: {missing}"
            )

        x = (
            frame[features]
            .apply(pd.to_numeric, errors="coerce")
            .replace([np.inf, -np.inf], np.nan)
        )

        nonfinite = x.isna().sum()
        bad = nonfinite[nonfinite > 0]

        if len(bad):
            raise ValueError(
                f"{name} geometry contains non-finite values. "
                f"No imputation is allowed. Counts: {bad.to_dict()}"
            )

    return features


def prepare_inputs(
    external_path: Path,
    reference_path: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    external = pd.read_csv(
        external_path,
        sep="\t",
        compression="infer",
        low_memory=False,
    )
    reference = pd.read_csv(
        reference_path,
        sep="\t",
        compression="infer",
        low_memory=False,
    )

    ext_required = {
        "analysis_id",
        "sequence_id",
        "source",
        "sequence_sha256_a2",
        "global_score",
        "global_prediction",
        "functional_prediction",
        "functional_class_a2",
    }
    ref_required = {
        "entity_id",
        "sequence_sha256_a2",
        "reference_cohort",
        "reference_class",
    }

    missing_ext = sorted(ext_required - set(external.columns))
    missing_ref = sorted(ref_required - set(reference.columns))

    if missing_ext:
        raise ValueError(
            f"External descriptor table missing columns: {missing_ext}"
        )
    if missing_ref:
        raise ValueError(
            f"Reference descriptor table missing columns: {missing_ref}"
        )

    unexpected_sources = sorted(
        set(external["source"].dropna()) - set(SOURCE_ORDER)
    )
    if unexpected_sources:
        raise ValueError(
            f"Unexpected external source labels: {unexpected_sources}"
        )

    return external, reference


def reference_self_geometry(
    X: np.ndarray,
    k: int,
) -> dict:
    n = len(X)
    if n < 2:
        raise ValueError("Reference geometry requires at least 2 sequences.")

    k_eff = min(k, n - 1)

    nn = NearestNeighbors(
        n_neighbors=k_eff + 1,
        metric="euclidean",
        algorithm="auto",
        n_jobs=-1,
    )
    nn.fit(X)

    distances, _ = nn.kneighbors(X)

    # First neighbor is self at distance 0.
    loo = distances[:, 1:]
    nearest = loo[:, 0]
    mean_k = loo.mean(axis=1)

    return {
        "n_reference": int(n),
        "k_effective": int(k_eff),
        "nearest_distances": nearest,
        "mean_k_distances": mean_k,
        "nearest_median": float(np.median(nearest)),
        "nearest_q95": float(np.quantile(nearest, 0.95)),
        "nearest_q99": float(np.quantile(nearest, 0.99)),
        "mean_k_median": float(np.median(mean_k)),
        "mean_k_q95": float(np.quantile(mean_k, 0.95)),
        "mean_k_q99": float(np.quantile(mean_k, 0.99)),
    }


def query_geometry(
    X_query: np.ndarray,
    X_reference: np.ndarray,
    k: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    k_eff = min(k, len(X_reference))

    nn = NearestNeighbors(
        n_neighbors=k_eff,
        metric="euclidean",
        algorithm="auto",
        n_jobs=-1,
    )
    nn.fit(X_reference)

    distances, indices = nn.kneighbors(X_query)

    nearest = distances[:, 0]
    mean_k = distances.mean(axis=1)

    return nearest, mean_k, indices


def build_global_geometry(
    external: pd.DataFrame,
    reference: pd.DataFrame,
    features: list[str],
    k: int,
) -> tuple[pd.DataFrame, pd.DataFrame, StandardScaler]:
    ref = (
        reference[
            reference["reference_cohort"] == "global_positive"
        ]
        .drop_duplicates("sequence_sha256_a2")
        .copy()
    )

    if ref.empty:
        raise RuntimeError("global_positive reference cohort is empty.")

    scaler = StandardScaler()
    X_ref = scaler.fit_transform(
        ref[features].to_numpy(dtype=float)
    )
    X_query = scaler.transform(
        external[features].to_numpy(dtype=float)
    )

    self_geom = reference_self_geometry(X_ref, k=k)

    nearest, mean_k, _ = query_geometry(
        X_query,
        X_ref,
        k=k,
    )

    ref_hashes = set(ref["sequence_sha256_a2"].astype(str))

    result = external[
        [
            "analysis_id",
            "sequence_id",
            "source",
            "sequence_sha256_a2",
            "global_score",
            "global_prediction",
            "functional_prediction",
            "functional_class_a2",
        ]
    ].copy()

    result["exact_global_reference_overlap"] = (
        result["sequence_sha256_a2"].astype(str).isin(ref_hashes)
    )

    result["global_nn_distance"] = nearest
    result["global_mean_k_distance"] = mean_k

    result["global_nn_distance_percentile"] = empirical_percentile(
        nearest,
        self_geom["nearest_distances"],
    )
    result["global_mean_k_distance_percentile"] = empirical_percentile(
        mean_k,
        self_geom["mean_k_distances"],
    )

    result["global_coverage95"] = (
        nearest <= self_geom["nearest_q95"]
    )
    result["global_coverage99"] = (
        nearest <= self_geom["nearest_q99"]
    )

    result["global_ood95"] = ~result["global_coverage95"]
    result["global_ood99"] = ~result["global_coverage99"]

    profile = pd.DataFrame(
        [
            {
                "reference_cohort": "global_positive",
                "n_reference": self_geom["n_reference"],
                "k_effective": self_geom["k_effective"],
                "nearest_median": self_geom["nearest_median"],
                "nearest_q95": self_geom["nearest_q95"],
                "nearest_q99": self_geom["nearest_q99"],
                "mean_k_median": self_geom["mean_k_median"],
                "mean_k_q95": self_geom["mean_k_q95"],
                "mean_k_q99": self_geom["mean_k_q99"],
            }
        ]
    )

    return result, profile, scaler


def build_functional_geometry(
    external: pd.DataFrame,
    reference: pd.DataFrame,
    features: list[str],
    k: int,
) -> tuple[
    pd.DataFrame,
    pd.DataFrame,
    StandardScaler,
    np.ndarray,
    pd.DataFrame,
]:
    ref = (
        reference[
            reference["reference_cohort"].isin(
                [f"multiclass_{x}" for x in FUNCTION_ORDER]
            )
        ]
        .copy()
    )

    ref["reference_class"] = ref["reference_class"].astype(str)

    unexpected = sorted(
        set(ref["reference_class"]) - set(FUNCTION_ORDER)
    )
    if unexpected:
        raise ValueError(
            f"Unexpected functional reference classes: {unexpected}"
        )

    ref = ref.drop_duplicates(
        ["reference_class", "sequence_sha256_a2"]
    ).copy()

    lytic = external[
        external["global_prediction"] == "lytic"
    ].copy()

    scaler = StandardScaler()
    X_ref = scaler.fit_transform(
        ref[features].to_numpy(dtype=float)
    )
    X_query = scaler.transform(
        lytic[features].to_numpy(dtype=float)
    )

    ref_labels = ref["reference_class"].to_numpy(dtype=str)

    # --------------------------------------------------------------
    # Complete-reference kNN composition
    # --------------------------------------------------------------
    k_eff = min(k, len(ref))

    knn = NearestNeighbors(
        n_neighbors=k_eff,
        metric="euclidean",
        algorithm="auto",
        n_jobs=-1,
    )
    knn.fit(X_ref)

    distances_all, indices_all = knn.kneighbors(X_query)
    neighbor_labels = ref_labels[indices_all]

    knn_fracs = {}
    for cls in FUNCTION_ORDER:
        knn_fracs[cls] = (
            neighbor_labels == cls
        ).mean(axis=1)

    fraction_matrix = np.column_stack(
        [knn_fracs[c] for c in FUNCTION_ORDER]
    )

    majority_idx = np.argmax(fraction_matrix, axis=1)
    majority_class = np.array(FUNCTION_ORDER, dtype=object)[
        majority_idx
    ]

    neighborhood_entropy = entropy3(fraction_matrix)

    # --------------------------------------------------------------
    # Class-specific nearest distances and density calibration
    # --------------------------------------------------------------
    class_profiles = []
    raw_nn_matrix = []
    normalized_nn_matrix = []
    centroid_matrix = []
    class_percentile_matrix = []

    exact_multiclass_hashes = set(
        ref["sequence_sha256_a2"].astype(str)
    )

    class_results = {}

    for cls in FUNCTION_ORDER:
        mask = ref_labels == cls
        X_cls = X_ref[mask]

        if len(X_cls) < 2:
            raise RuntimeError(
                f"Functional reference class {cls} has <2 sequences."
            )

        self_geom = reference_self_geometry(
            X_cls,
            k=min(k, len(X_cls) - 1),
        )

        nearest, mean_k_cls, _ = query_geometry(
            X_query,
            X_cls,
            k=min(k, len(X_cls)),
        )

        centroid = X_cls.mean(axis=0)
        centroid_distance = np.linalg.norm(
            X_query - centroid,
            axis=1,
        )

        median_self_nn = self_geom["nearest_median"]
        if median_self_nn <= 0:
            density_normalized = np.full(
                len(nearest),
                np.nan,
            )
        else:
            density_normalized = nearest / median_self_nn

        percentiles = empirical_percentile(
            nearest,
            self_geom["nearest_distances"],
        )

        raw_nn_matrix.append(nearest)
        normalized_nn_matrix.append(density_normalized)
        centroid_matrix.append(centroid_distance)
        class_percentile_matrix.append(percentiles)

        class_results[cls] = {
            "nearest": nearest,
            "mean_k": mean_k_cls,
            "centroid": centroid_distance,
            "density_normalized": density_normalized,
            "percentile": percentiles,
            "coverage95": nearest <= self_geom["nearest_q95"],
            "coverage99": nearest <= self_geom["nearest_q99"],
        }

        class_profiles.append(
            {
                "reference_class": cls,
                "n_reference": self_geom["n_reference"],
                "k_effective": self_geom["k_effective"],
                "nearest_median": self_geom["nearest_median"],
                "nearest_q95": self_geom["nearest_q95"],
                "nearest_q99": self_geom["nearest_q99"],
                "mean_k_median": self_geom["mean_k_median"],
                "mean_k_q95": self_geom["mean_k_q95"],
                "mean_k_q99": self_geom["mean_k_q99"],
            }
        )

    raw_nn_matrix = np.column_stack(raw_nn_matrix)
    normalized_nn_matrix = np.column_stack(
        normalized_nn_matrix
    )
    centroid_matrix = np.column_stack(centroid_matrix)
    class_percentile_matrix = np.column_stack(
        class_percentile_matrix
    )

    raw_idx = np.argmin(raw_nn_matrix, axis=1)
    norm_idx = np.argmin(normalized_nn_matrix, axis=1)
    centroid_idx = np.argmin(centroid_matrix, axis=1)

    nearest_raw_class = np.array(
        FUNCTION_ORDER,
        dtype=object,
    )[raw_idx]
    nearest_norm_class = np.array(
        FUNCTION_ORDER,
        dtype=object,
    )[norm_idx]
    nearest_centroid_class = np.array(
        FUNCTION_ORDER,
        dtype=object,
    )[centroid_idx]

    _, raw_margin = second_smallest_margin(
        raw_nn_matrix
    )
    _, normalized_margin = second_smallest_margin(
        normalized_nn_matrix
    )
    _, centroid_margin = second_smallest_margin(
        centroid_matrix
    )

    # --------------------------------------------------------------
    # Assemble per-sequence geometry table
    # --------------------------------------------------------------
    result = lytic[
        [
            "analysis_id",
            "sequence_id",
            "source",
            "sequence_sha256_a2",
            "global_score",
            "functional_prediction",
            "functional_class_a2",
            "p_endolysin_given_lytic",
            "p_depolymerase_given_lytic",
            "p_val_given_lytic",
            "endolysin_compatibility_score",
            "depolymerase_compatibility_score",
            "val_compatibility_score",
        ]
    ].copy()

    result["exact_multiclass_reference_overlap"] = (
        result["sequence_sha256_a2"]
        .astype(str)
        .isin(exact_multiclass_hashes)
    )

    for cls in FUNCTION_ORDER:
        result[f"nn_distance_{cls}"] = class_results[cls][
            "nearest"
        ]
        result[f"mean_k_distance_{cls}"] = class_results[cls][
            "mean_k"
        ]
        result[f"centroid_distance_{cls}"] = class_results[cls][
            "centroid"
        ]
        result[f"density_normalized_nn_{cls}"] = class_results[
            cls
        ]["density_normalized"]
        result[f"nn_percentile_{cls}"] = class_results[cls][
            "percentile"
        ]
        result[f"coverage95_{cls}"] = class_results[cls][
            "coverage95"
        ]
        result[f"coverage99_{cls}"] = class_results[cls][
            "coverage99"
        ]

        result[f"knn_fraction_{cls}"] = knn_fracs[cls]

    result["knn_k_effective"] = int(k_eff)
    result["knn_majority_subtype"] = majority_class
    result["knn_neighborhood_entropy"] = neighborhood_entropy

    result["nearest_raw_nn_subtype"] = nearest_raw_class
    result["nearest_density_normalized_subtype"] = nearest_norm_class
    result["nearest_centroid_subtype"] = nearest_centroid_class

    result["raw_nn_subtype_margin"] = raw_margin
    result["density_normalized_subtype_margin"] = normalized_margin
    result["centroid_subtype_margin"] = centroid_margin

    predicted = result["functional_class_a2"].astype(str).to_numpy()

    result["agreement_knn_majority"] = (
        predicted == majority_class
    )
    result["agreement_raw_nn"] = (
        predicted == nearest_raw_class
    )
    result["agreement_density_normalized_nn"] = (
        predicted == nearest_norm_class
    )
    result["agreement_centroid"] = (
        predicted == nearest_centroid_class
    )

    # Extract the predicted class-specific geometry for each row.
    predicted_nn = np.full(len(result), np.nan)
    predicted_norm = np.full(len(result), np.nan)
    predicted_centroid = np.full(len(result), np.nan)
    predicted_percentile = np.full(len(result), np.nan)
    predicted_cov95 = np.zeros(len(result), dtype=bool)
    predicted_cov99 = np.zeros(len(result), dtype=bool)
    predicted_knn_frac = np.full(len(result), np.nan)

    for j, cls in enumerate(FUNCTION_ORDER):
        mask = predicted == cls

        predicted_nn[mask] = raw_nn_matrix[mask, j]
        predicted_norm[mask] = normalized_nn_matrix[mask, j]
        predicted_centroid[mask] = centroid_matrix[mask, j]
        predicted_percentile[mask] = class_percentile_matrix[mask, j]
        predicted_cov95[mask] = class_results[cls][
            "coverage95"
        ][mask]
        predicted_cov99[mask] = class_results[cls][
            "coverage99"
        ][mask]
        predicted_knn_frac[mask] = knn_fracs[cls][mask]

    result["predicted_subtype_nn_distance"] = predicted_nn
    result[
        "predicted_subtype_density_normalized_nn"
    ] = predicted_norm
    result["predicted_subtype_centroid_distance"] = (
        predicted_centroid
    )
    result["predicted_subtype_nn_percentile"] = (
        predicted_percentile
    )
    result["predicted_subtype_coverage95"] = predicted_cov95
    result["predicted_subtype_coverage99"] = predicted_cov99
    result["predicted_subtype_ood95"] = ~predicted_cov95
    result["predicted_subtype_ood99"] = ~predicted_cov99
    result["predicted_subtype_knn_fraction"] = predicted_knn_frac

    profile = pd.DataFrame(class_profiles)

    return result, profile, scaler, X_ref, ref


def build_reference_pca(
    external: pd.DataFrame,
    reference: pd.DataFrame,
    features: list[str],
    pca_components: int,
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    ref = (
        reference[
            reference["reference_cohort"].isin(
                [f"multiclass_{x}" for x in FUNCTION_ORDER]
            )
        ]
        .drop_duplicates(
            ["reference_class", "sequence_sha256_a2"]
        )
        .copy()
    )

    scaler = StandardScaler()
    X_ref = scaler.fit_transform(
        ref[features].to_numpy(dtype=float)
    )
    X_ext = scaler.transform(
        external[features].to_numpy(dtype=float)
    )

    n_components = min(
        pca_components,
        X_ref.shape[1],
        X_ref.shape[0] - 1,
    )

    pca = PCA(
        n_components=n_components,
        svd_solver="auto",
        random_state=42,
    )
    Z_ref = pca.fit_transform(X_ref)
    Z_ext = pca.transform(X_ext)

    ref_out = ref[
        [
            "entity_id",
            "sequence_sha256_a2",
            "reference_class",
            "reference_cohort",
        ]
    ].copy()

    ext_out = external[
        [
            "analysis_id",
            "sequence_id",
            "source",
            "sequence_sha256_a2",
            "global_score",
            "global_prediction",
            "functional_prediction",
            "functional_class_a2",
        ]
    ].copy()

    for i in range(n_components):
        ref_out[f"PC{i + 1}"] = Z_ref[:, i]
        ext_out[f"PC{i + 1}"] = Z_ext[:, i]

    metadata = {
        "n_components": int(n_components),
        "explained_variance_ratio": [
            float(x)
            for x in pca.explained_variance_ratio_
        ],
        "cumulative_explained_variance_ratio": [
            float(x)
            for x in np.cumsum(
                pca.explained_variance_ratio_
            )
        ],
    }

    return ref_out, ext_out, metadata


def build_reference_umap(
    external: pd.DataFrame,
    reference: pd.DataFrame,
    features: list[str],
    n_neighbors: int,
    min_dist: float,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    if umap is None:
        raise RuntimeError(
            "UMAP visualization requested but 'umap-learn' is not installed. "
            "Install it with: pip install umap-learn"
        )

    ref = (
        reference[
            reference["reference_cohort"].isin(
                [f"multiclass_{x}" for x in FUNCTION_ORDER]
            )
        ]
        .drop_duplicates(
            ["reference_class", "sequence_sha256_a2"]
        )
        .copy()
    )

    scaler = StandardScaler()
    X_ref = scaler.fit_transform(
        ref[features].to_numpy(dtype=float)
    )
    X_ext = scaler.transform(
        external[features].to_numpy(dtype=float)
    )

    reducer = umap.UMAP(
        n_components=2,
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        metric="euclidean",
        random_state=seed,
        transform_seed=seed,
    )

    Z_ref = reducer.fit_transform(X_ref)
    Z_ext = reducer.transform(X_ext)

    ref_out = ref[
        [
            "entity_id",
            "sequence_sha256_a2",
            "reference_class",
            "reference_cohort",
        ]
    ].copy()
    ref_out["UMAP1"] = Z_ref[:, 0]
    ref_out["UMAP2"] = Z_ref[:, 1]

    ext_out = external[
        [
            "analysis_id",
            "sequence_id",
            "source",
            "sequence_sha256_a2",
            "global_score",
            "global_prediction",
            "functional_prediction",
            "functional_class_a2",
        ]
    ].copy()
    ext_out["UMAP1"] = Z_ext[:, 0]
    ext_out["UMAP2"] = Z_ext[:, 1]

    metadata = {
        "n_neighbors": int(n_neighbors),
        "min_dist": float(min_dist),
        "metric": "euclidean",
        "fit_policy": (
            "StandardScaler and UMAP fitted on multiclass reference only; "
            "external sequences transformed out-of-sample."
        ),
    }

    return ref_out, ext_out, metadata


def build_joint_tsne_source(
    external: pd.DataFrame,
    reference: pd.DataFrame,
    features: list[str],
    source: str,
    max_reference: int,
    max_external: int,
    perplexity: float,
    seed: int,
) -> tuple[pd.DataFrame, dict]:
    ref = (
        reference[
            reference["reference_cohort"].isin(
                [f"multiclass_{x}" for x in FUNCTION_ORDER]
            )
        ]
        .drop_duplicates(
            ["reference_class", "sequence_sha256_a2"]
        )
        .copy()
    )

    ext = external[
        (external["source"] == source)
        & (external["global_prediction"] == "lytic")
    ].copy()

    ref_sample = stratified_sample(
        ref,
        group_col="reference_class",
        max_total=max_reference,
        seed=seed,
    ).copy()

    if len(ext) > max_external:
        # Preserve predicted subtype composition when sampling.
        ext = stratified_sample(
            ext,
            group_col="functional_class_a2",
            max_total=max_external,
            seed=seed + 100,
        )

    scaler = StandardScaler()
    X_ref_all = scaler.fit(
        ref[features].to_numpy(dtype=float)
    )

    X_ref = scaler.transform(
        ref_sample[features].to_numpy(dtype=float)
    )
    X_ext = scaler.transform(
        ext[features].to_numpy(dtype=float)
    )

    meta_ref = pd.DataFrame(
        {
            "origin": "Reference",
            "source": "reference",
            "entity_id": ref_sample["entity_id"].astype(str).to_numpy(),
            "analysis_id": pd.NA,
            "functional_class": ref_sample[
                "reference_class"
            ].astype(str).to_numpy(),
        }
    )

    meta_ext = pd.DataFrame(
        {
            "origin": "External",
            "source": ext["source"].astype(str).to_numpy(),
            "entity_id": pd.NA,
            "analysis_id": ext["analysis_id"].astype(str).to_numpy(),
            "functional_class": ext[
                "functional_class_a2"
            ].astype(str).to_numpy(),
        }
    )

    X = np.vstack([X_ref, X_ext])
    meta = pd.concat(
        [meta_ref, meta_ext],
        ignore_index=True,
    )

    effective_perplexity = min(
        float(perplexity),
        max(5.0, (len(X) - 1) / 3.0),
    )

    tsne = TSNE(
        n_components=2,
        perplexity=effective_perplexity,
        learning_rate="auto",
        init="pca",
        metric="euclidean",
        random_state=seed,
    )

    Z = tsne.fit_transform(X)
    meta["TSNE1"] = Z[:, 0]
    meta["TSNE2"] = Z[:, 1]

    metadata = {
        "source": source,
        "n_reference_sample": int(len(ref_sample)),
        "n_external_sample": int(len(ext)),
        "perplexity": float(effective_perplexity),
        "fit_policy": (
            "Joint t-SNE on a stratified reference sample plus a sampled "
            "external source. Visualization only; no out-of-sample projection "
            "or quantitative distance interpretation is claimed."
        ),
    }

    return meta, metadata


def plot_reference_external_embedding(
    reference_coords: pd.DataFrame,
    external_coords: pd.DataFrame,
    source: str,
    x_col: str,
    y_col: str,
    method_name: str,
    figures: Path,
    max_reference: int,
    max_external: int,
    seed: int,
    stem_name: str,
) -> None:
    ref_plot = stratified_sample(
        reference_coords,
        group_col="reference_class",
        max_total=max_reference,
        seed=seed,
    ).copy()

    ext = external_coords[
        (external_coords["source"] == source)
        & (external_coords["global_prediction"] == "lytic")
    ].copy()

    if len(ext) > max_external:
        ext = stratified_sample(
            ext,
            group_col="functional_class_a2",
            max_total=max_external,
            seed=seed + 1,
        )

    ref_plot["origin"] = "Reference"
    ref_plot["functional_class"] = ref_plot[
        "reference_class"
    ]
    ext["origin"] = "External"
    ext["functional_class"] = ext[
        "functional_class_a2"
    ]

    combined = pd.concat(
        [
            ref_plot[
                [x_col, y_col, "origin", "functional_class"]
            ],
            ext[
                [x_col, y_col, "origin", "functional_class"]
            ],
        ],
        ignore_index=True,
    )

    combined["functional_label"] = combined[
        "functional_class"
    ].map(FUNCTION_LABELS)

    fig, ax = plt.subplots(figsize=(7.2, 5.5))

    sns.scatterplot(
        data=combined,
        x=x_col,
        y=y_col,
        hue="functional_label",
        style="origin",
        alpha=0.48,
        s=20,
        linewidth=0,
        ax=ax,
    )

    ax.set_title(
        f"{method_name} — {SOURCE_LABELS[source]}"
    )
    ax.set_xlabel(x_col)
    ax.set_ylabel(y_col)
    ax.legend(
        title="",
        frameon=False,
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
    )
    sns.despine(ax=ax)

    save_figure(
        fig,
        figures / stem_name,
    )


def plot_joint_tsne(
    tsne_coords: pd.DataFrame,
    source: str,
    figures: Path,
    stem_name: str,
) -> None:
    plot_df = tsne_coords.copy()
    plot_df["functional_label"] = plot_df[
        "functional_class"
    ].map(FUNCTION_LABELS)

    fig, ax = plt.subplots(figsize=(7.2, 5.5))

    sns.scatterplot(
        data=plot_df,
        x="TSNE1",
        y="TSNE2",
        hue="functional_label",
        style="origin",
        alpha=0.50,
        s=20,
        linewidth=0,
        ax=ax,
    )

    ax.set_title(
        f"Joint t-SNE — {SOURCE_LABELS[source]}"
    )
    ax.set_xlabel("t-SNE1")
    ax.set_ylabel("t-SNE2")
    ax.legend(
        title="",
        frameon=False,
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
    )
    sns.despine(ax=ax)

    save_figure(
        fig,
        figures / stem_name,
    )


def build_global_summary(
    global_geometry: pd.DataFrame,
) -> pd.DataFrame:
    rows = []

    for source in SOURCE_ORDER:
        for prediction in [
            "lytic",
            "non_lytic_background",
        ]:
            g = global_geometry[
                (global_geometry["source"] == source)
                & (
                    global_geometry["global_prediction"]
                    == prediction
                )
            ]

            if g.empty:
                continue

            rows.append(
                {
                    "source": source,
                    "global_prediction": prediction,
                    "n": int(len(g)),
                    "exact_reference_overlap_n": int(
                        g[
                            "exact_global_reference_overlap"
                        ].sum()
                    ),
                    "coverage95_fraction": float(
                        g["global_coverage95"].mean()
                    ),
                    "coverage99_fraction": float(
                        g["global_coverage99"].mean()
                    ),
                    "median_nn_distance": float(
                        g["global_nn_distance"].median()
                    ),
                    "median_nn_percentile": float(
                        g[
                            "global_nn_distance_percentile"
                        ].median()
                    ),
                    "median_mean_k_distance": float(
                        g[
                            "global_mean_k_distance"
                        ].median()
                    ),
                }
            )

    return pd.DataFrame(rows)


def build_functional_summary(
    geometry: pd.DataFrame,
) -> pd.DataFrame:
    rows = []

    for source in SOURCE_ORDER:
        for cls in FUNCTION_ORDER:
            g = geometry[
                (geometry["source"] == source)
                & (
                    geometry["functional_class_a2"]
                    == cls
                )
            ]

            if g.empty:
                continue

            rows.append(
                {
                    "source": source,
                    "predicted_class": cls,
                    "n": int(len(g)),
                    "exact_reference_overlap_n": int(
                        g[
                            "exact_multiclass_reference_overlap"
                        ].sum()
                    ),
                    "predicted_subtype_coverage95_fraction": float(
                        g[
                            "predicted_subtype_coverage95"
                        ].mean()
                    ),
                    "predicted_subtype_coverage99_fraction": float(
                        g[
                            "predicted_subtype_coverage99"
                        ].mean()
                    ),
                    "median_predicted_subtype_nn_percentile": float(
                        g[
                            "predicted_subtype_nn_percentile"
                        ].median()
                    ),
                    "median_density_normalized_nn": float(
                        g[
                            "predicted_subtype_density_normalized_nn"
                        ].median()
                    ),
                    "median_predicted_subtype_knn_fraction": float(
                        g[
                            "predicted_subtype_knn_fraction"
                        ].median()
                    ),
                    "agreement_knn_majority_fraction": float(
                        g["agreement_knn_majority"].mean()
                    ),
                    "agreement_raw_nn_fraction": float(
                        g["agreement_raw_nn"].mean()
                    ),
                    "agreement_density_normalized_nn_fraction": float(
                        g[
                            "agreement_density_normalized_nn"
                        ].mean()
                    ),
                    "agreement_centroid_fraction": float(
                        g["agreement_centroid"].mean()
                    ),
                    "median_knn_entropy": float(
                        g[
                            "knn_neighborhood_entropy"
                        ].median()
                    ),
                }
            )

    return pd.DataFrame(rows)


def build_agreement_summary(
    geometry: pd.DataFrame,
) -> pd.DataFrame:
    metrics = [
        "agreement_knn_majority",
        "agreement_raw_nn",
        "agreement_density_normalized_nn",
        "agreement_centroid",
    ]

    rows = []

    for source in SOURCE_ORDER:
        for cls in FUNCTION_ORDER:
            g = geometry[
                (geometry["source"] == source)
                & (
                    geometry["functional_class_a2"]
                    == cls
                )
            ]

            if g.empty:
                continue

            for metric in metrics:
                rows.append(
                    {
                        "source": source,
                        "predicted_class": cls,
                        "geometry_rule": metric,
                        "n": int(len(g)),
                        "agreement_fraction": float(
                            g[metric].mean()
                        ),
                    }
                )

    return pd.DataFrame(rows)


def plot_pca_source(
    reference_pca: pd.DataFrame,
    external_pca: pd.DataFrame,
    source: str,
    figures: Path,
    max_reference: int,
    max_external: int,
    seed: int,
    stem_name: str,
) -> None:
    ref_plot = stratified_sample(
        reference_pca,
        group_col="reference_class",
        max_total=max_reference,
        seed=seed,
    ).copy()

    ext = external_pca[
        (external_pca["source"] == source)
        & (
            external_pca["global_prediction"]
            == "lytic"
        )
    ].copy()

    if len(ext) > max_external:
        ext = ext.sample(
            n=max_external,
            random_state=seed,
        )

    ref_plot["plot_group"] = (
        "Reference "
        + ref_plot["reference_class"].map(
            FUNCTION_LABELS
        )
    )

    ext["plot_group"] = (
        "External "
        + ext["functional_class_a2"].map(
            FUNCTION_LABELS
        )
    )

    combined = pd.concat(
        [
            ref_plot[
                ["PC1", "PC2", "plot_group"]
            ],
            ext[
                ["PC1", "PC2", "plot_group"]
            ],
        ],
        ignore_index=True,
    )

    fig, ax = plt.subplots(figsize=(7.2, 5.5))

    sns.scatterplot(
        data=combined,
        x="PC1",
        y="PC2",
        hue="plot_group",
        alpha=0.45,
        s=18,
        linewidth=0,
        ax=ax,
    )

    ax.set_title(
        f"Reference-fitted PCA — {SOURCE_LABELS[source]}"
    )
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.legend(
        title="",
        frameon=False,
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
    )
    sns.despine(ax=ax)

    save_figure(
        fig,
        figures / stem_name,
    )


def plot_coverage(
    global_summary: pd.DataFrame,
    functional_summary: pd.DataFrame,
    figures: Path,
) -> None:
    global_lytic = global_summary[
        global_summary["global_prediction"] == "lytic"
    ][
        [
            "source",
            "coverage95_fraction",
        ]
    ].copy()
    global_lytic["geometry_level"] = "Global positive"

    functional = (
        functional_summary
        .groupby("source", as_index=False)
        .apply(
            lambda g: pd.Series(
                {
                    "coverage95_fraction": np.average(
                        g[
                            "predicted_subtype_coverage95_fraction"
                        ],
                        weights=g["n"],
                    )
                }
            ),
            include_groups=False,
        )
        .reset_index()
    )

    if "level_0" in functional.columns:
        functional = functional.drop(
            columns=["level_0"]
        )

    functional["geometry_level"] = "Predicted subtype"

    plot_df = pd.concat(
        [
            global_lytic,
            functional[
                [
                    "source",
                    "coverage95_fraction",
                    "geometry_level",
                ]
            ],
        ],
        ignore_index=True,
    )

    plot_df["source_label"] = plot_df[
        "source"
    ].map(SOURCE_LABELS)
    plot_df["coverage_percent"] = (
        100.0 * plot_df["coverage95_fraction"]
    )

    fig, ax = plt.subplots(figsize=(7.5, 4.8))

    sns.barplot(
        data=plot_df,
        x="source_label",
        y="coverage_percent",
        hue="geometry_level",
        order=[SOURCE_LABELS[x] for x in SOURCE_ORDER],
        errorbar=None,
        ax=ax,
    )

    ax.set_xlabel("")
    ax.set_ylabel("Within reference 95% coverage (%)")
    ax.set_ylim(0, 100)
    ax.set_title("Reference-space coverage")
    ax.legend(title="", frameon=False)
    ax.tick_params(axis="x", rotation=15)
    sns.despine(ax=ax)

    save_figure(
        fig,
        figures / "F2_reference_coverage_by_source",
    )


def plot_agreement_heatmap(
    agreement: pd.DataFrame,
    figures: Path,
) -> None:
    if agreement.empty:
        return

    df = agreement.copy()

    df["row"] = (
        df["source"].map(SOURCE_LABELS)
        + " — "
        + df["predicted_class"].map(FUNCTION_LABELS)
    )

    rule_labels = {
        "agreement_knn_majority": "kNN majority",
        "agreement_raw_nn": "Raw NN",
        "agreement_density_normalized_nn": "Density-normalized NN",
        "agreement_centroid": "Centroid",
    }
    df["rule"] = df["geometry_rule"].map(
        rule_labels
    )

    matrix = df.pivot(
        index="row",
        columns="rule",
        values="agreement_fraction",
    )

    desired_rows = []
    for source in SOURCE_ORDER:
        for cls in FUNCTION_ORDER:
            label = (
                f"{SOURCE_LABELS[source]} — "
                f"{FUNCTION_LABELS[cls]}"
            )
            if label in matrix.index:
                desired_rows.append(label)

    desired_cols = [
        "kNN majority",
        "Raw NN",
        "Density-normalized NN",
        "Centroid",
    ]

    matrix = matrix.reindex(
        index=desired_rows,
        columns=desired_cols,
    )

    fig, ax = plt.subplots(
        figsize=(7.0, max(4.8, 0.48 * len(matrix)))
    )

    sns.heatmap(
        100.0 * matrix,
        annot=True,
        fmt=".1f",
        linewidths=0.25,
        cbar_kws={"label": "Agreement (%)"},
        ax=ax,
    )

    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title(
        "Prediction agreement with reference geometry"
    )
    ax.tick_params(axis="x", rotation=30)
    ax.tick_params(axis="y", rotation=0)

    save_figure(
        fig,
        figures / "F3_subtype_geometry_agreement",
    )


def plot_distance_percentile(
    geometry: pd.DataFrame,
    figures: Path,
) -> None:
    plot_df = geometry[
        [
            "source",
            "functional_class_a2",
            "predicted_subtype_nn_percentile",
        ]
    ].copy()

    plot_df["source_label"] = plot_df[
        "source"
    ].map(SOURCE_LABELS)
    plot_df["functional_label"] = plot_df[
        "functional_class_a2"
    ].map(FUNCTION_LABELS)
    plot_df["distance_percentile"] = (
        100.0
        * plot_df[
            "predicted_subtype_nn_percentile"
        ]
    )

    fig, ax = plt.subplots(figsize=(9.0, 5.2))

    sns.boxplot(
        data=plot_df,
        x="functional_label",
        y="distance_percentile",
        hue="source_label",
        order=[
            FUNCTION_LABELS[x]
            for x in FUNCTION_ORDER
        ],
        hue_order=[
            SOURCE_LABELS[x]
            for x in SOURCE_ORDER
        ],
        showfliers=False,
        ax=ax,
    )

    ax.axhline(
        95,
        linestyle="--",
        linewidth=1,
    )
    ax.set_xlabel("")
    ax.set_ylabel(
        "Predicted-subtype NN distance percentile"
    )
    ax.set_ylim(0, 100)
    ax.set_title(
        "Class-conditional geometric remoteness"
    )
    ax.legend(title="", frameon=False)
    sns.despine(ax=ax)

    save_figure(
        fig,
        figures
        / "F4_predicted_subtype_distance_percentile",
    )


def main() -> None:
    args = parse_args()

    external_path = require_file(
        args.external_descriptors.resolve()
    )
    reference_path = require_file(
        args.reference_descriptors.resolve()
    )
    registry_path = require_file(
        args.descriptor_registry.resolve()
    )

    outdir = args.output_dir.resolve()
    tables = outdir / "tables"
    figures = outdir / "figures"

    tables.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)

    external, reference = prepare_inputs(
        external_path,
        reference_path,
    )

    features = load_feature_contract(
        registry_path,
        reference,
        external,
    )

    registry = pd.read_csv(
        registry_path,
        sep="\t",
    )
    registry.to_csv(
        tables / "geometry_feature_contract.tsv",
        sep="\t",
        index=False,
    )

    print(
        "A3: global positive geometry...",
        flush=True,
    )

    (
        global_geometry,
        global_profile,
        _,
    ) = build_global_geometry(
        external,
        reference,
        features,
        k=args.k,
    )

    global_profile.to_csv(
        tables
        / "global_reference_self_distance_profile.tsv",
        sep="\t",
        index=False,
    )

    global_geometry.to_csv(
        tables / "global_reference_coverage.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    print(
        "A3: functional subtype geometry...",
        flush=True,
    )

    (
        functional_geometry,
        subtype_profiles,
        _,
        _,
        _,
    ) = build_functional_geometry(
        external,
        reference,
        features,
        k=args.k,
    )

    subtype_profiles.to_csv(
        tables
        / "subtype_reference_self_distance_profiles.tsv",
        sep="\t",
        index=False,
    )

    functional_geometry.to_csv(
        tables / "functional_reference_geometry.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    print(
        "A3: reference-fitted PCA projection...",
        flush=True,
    )

    (
        reference_pca,
        external_pca,
        pca_metadata,
    ) = build_reference_pca(
        external,
        reference,
        features,
        pca_components=args.pca_components,
    )

    reference_pca.to_csv(
        tables / "pca_reference_coordinates.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    external_pca.to_csv(
        tables / "pca_external_coordinates.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    print(
        "A3: reference-fitted UMAP projection...",
        flush=True,
    )

    (
        reference_umap,
        external_umap,
        umap_metadata,
    ) = build_reference_umap(
        external,
        reference,
        features,
        n_neighbors=args.umap_n_neighbors,
        min_dist=args.umap_min_dist,
        seed=args.seed,
    )

    reference_umap.to_csv(
        tables / "umap_reference_coordinates.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    external_umap.to_csv(
        tables / "umap_external_coordinates.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    tsne_metadata = {}
    tsne_tables = {}

    for i, source in enumerate(SOURCE_ORDER):
        print(
            f"A3: joint t-SNE visualization — {source}...",
            flush=True,
        )

        tsne_coords, tsne_meta = build_joint_tsne_source(
            external,
            reference,
            features,
            source=source,
            max_reference=args.tsne_max_reference,
            max_external=args.tsne_max_external_per_source,
            perplexity=args.tsne_perplexity,
            seed=args.seed + i,
        )

        tsne_path = (
            tables / f"tsne_{source}_coordinates.tsv.gz"
        )
        tsne_coords.to_csv(
            tsne_path,
            sep="\t",
            index=False,
            compression="gzip",
        )

        tsne_metadata[source] = tsne_meta
        tsne_tables[source] = tsne_path

    global_summary = build_global_summary(
        global_geometry
    )
    functional_summary = build_functional_summary(
        functional_geometry
    )
    agreement_summary = build_agreement_summary(
        functional_geometry
    )

    global_summary.to_csv(
        tables / "global_geometry_summary.tsv",
        sep="\t",
        index=False,
    )

    functional_summary.to_csv(
        tables / "functional_geometry_summary.tsv",
        sep="\t",
        index=False,
    )

    agreement_summary.to_csv(
        tables / "neighborhood_agreement_summary.tsv",
        sep="\t",
        index=False,
    )

    print(
        "A3: figures...",
        flush=True,
    )

    # PCA remains as a linear control visualization.
    for source, stem in [
        (
            "phagelys_prediction_extension",
            "F1a_pca_prediction_extension",
        ),
        (
            "zymctrl",
            "F1b_pca_zymctrl",
        ),
        (
            "ambiguous_structural",
            "F1c_pca_ambiguous_structural",
        ),
    ]:
        plot_reference_external_embedding(
            reference_pca,
            external_pca,
            source=source,
            x_col="PC1",
            y_col="PC2",
            method_name="Reference-fitted PCA",
            figures=figures,
            max_reference=args.plot_max_reference,
            max_external=args.plot_max_external_per_source,
            seed=args.seed,
            stem_name=stem,
        )

    # UMAP is the principal nonlinear out-of-sample projection.
    for source, stem in [
        (
            "phagelys_prediction_extension",
            "F1d_umap_prediction_extension",
        ),
        (
            "zymctrl",
            "F1e_umap_zymctrl",
        ),
        (
            "ambiguous_structural",
            "F1f_umap_ambiguous_structural",
        ),
    ]:
        plot_reference_external_embedding(
            reference_umap,
            external_umap,
            source=source,
            x_col="UMAP1",
            y_col="UMAP2",
            method_name="Reference-fitted UMAP",
            figures=figures,
            max_reference=args.plot_max_reference,
            max_external=args.plot_max_external_per_source,
            seed=args.seed,
            stem_name=stem,
        )

    # t-SNE is visualization only and is explicitly not called a projection.
    tsne_stems = {
        "phagelys_prediction_extension":
            "F1g_tsne_prediction_extension",
        "zymctrl":
            "F1h_tsne_zymctrl",
        "ambiguous_structural":
            "F1i_tsne_ambiguous_structural",
    }

    for source in SOURCE_ORDER:
        coords = pd.read_csv(
            tsne_tables[source],
            sep="\t",
            compression="infer",
            low_memory=False,
        )
        plot_joint_tsne(
            coords,
            source=source,
            figures=figures,
            stem_name=tsne_stems[source],
        )

    plot_coverage(
        global_summary,
        functional_summary,
        figures,
    )

    plot_agreement_heatmap(
        agreement_summary,
        figures,
    )

    plot_distance_percentile(
        functional_geometry,
        figures,
    )

    exact_global_overlap = int(
        global_geometry[
            "exact_global_reference_overlap"
        ].sum()
    )
    exact_multiclass_overlap = int(
        functional_geometry[
            "exact_multiclass_reference_overlap"
        ].sum()
    )

    manifest = {
        "stage": "A3",
        "stage_name": "reference_geometry",
        "generated_utc": datetime.now(
            timezone.utc
        ).isoformat(),
        "input_external_descriptors": str(
            external_path
        ),
        "input_reference_descriptors": str(
            reference_path
        ),
        "descriptor_registry": str(
            registry_path
        ),
        "n_frozen_features": int(len(features)),
        "frozen_features": features,
        "distance_space": (
            "Euclidean distance in full 26D descriptor space after "
            "reference-only StandardScaler transformation."
        ),
        "global_reference": (
            "A2 global_positive cohort from production global_mixed "
            "positive training class."
        ),
        "functional_reference": (
            "A2 multiclass_clean Endolysin / Depolymerase / VAL "
            "training classes."
        ),
        "reference_coverage_definition": (
            "A query is within 95% reference coverage when its nearest "
            "reference distance is <= the 95th percentile of the "
            "leave-one-out nearest-neighbor distance distribution of "
            "the corresponding annotated reference cohort."
        ),
        "knn_k_requested": int(args.k),
        "pca": {
            **pca_metadata,
            "fit_policy": (
                "StandardScaler and PCA fitted on multiclass reference "
                "only; external sequences transformed out-of-sample."
            ),
            "interpretation": (
                "Linear visualization only; quantitative geometry uses the "
                "full 26-dimensional standardized descriptor space."
            ),
        },
        "umap": {
            **umap_metadata,
            "interpretation": (
                "Nonlinear reference-fitted visualization only; quantitative "
                "geometry uses the full 26-dimensional standardized space."
            ),
        },
        "tsne": tsne_metadata,
        "exact_overlap_audit": {
            "external_vs_global_positive": exact_global_overlap,
            "external_lytic_vs_multiclass_reference": (
                exact_multiclass_overlap
            ),
        },
        "notes": {
            "manifold_visualization": (
                "UMAP is fitted on reference only and external sequences are "
                "transformed out-of-sample. t-SNE is fitted jointly only on "
                "stratified visualization samples and is never used as a "
                "quantitative geometry or described as a projection."
            ),
            "class_imbalance": (
                "Raw kNN composition is retained, but subtype assignment "
                "is also evaluated with class-specific nearest distances "
                "normalized by each subtype's own leave-one-out reference "
                "density to reduce dependence on unequal class sizes."
            ),
            "no_retraining": True,
            "no_feature_reselection": True,
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

    print("\nA3 — PASS")
    print("=" * 112)
    print(
        f"External sequences:                    "
        f"{len(external):,}"
    )
    print(
        f"External lytic sequences:              "
        f"{len(functional_geometry):,}"
    )
    print(
        f"Frozen descriptor dimensions:          "
        f"{len(features)}"
    )
    print(
        f"Exact overlap vs global reference:      "
        f"{exact_global_overlap:,}"
    )
    print(
        f"Exact overlap vs multiclass reference:  "
        f"{exact_multiclass_overlap:,}"
    )
    print(
        f"PCA explained variance PC1+PC2:         "
        f"{100 * sum(pca_metadata['explained_variance_ratio'][:2]):.2f}%"
    )
    print("=" * 112)
    print("\nGLOBAL REFERENCE COVERAGE")
    print(
        global_summary.to_string(
            index=False
        )
    )
    print("=" * 112)
    print("\nFUNCTIONAL REFERENCE GEOMETRY")
    print(
        functional_summary.to_string(
            index=False
        )
    )
    print("=" * 112)
    print(f"Tables:  {tables}")
    print(f"Figures: {figures}")
    print(f"Output:  {outdir}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
A2 — Property characterization of explored PhageLysClass predictions.

Scientific objective
--------------------
Test whether externally predicted proteins are physicochemically consistent
with the annotated sequence populations used to train the production models.

Two complementary levels are evaluated:

1. Population level
   - Predicted lytic sequences from each external source
     vs the positive class of the all-data global_mixed training universe.
   - Predicted Endolysin / Depolymerase / VAL sequences
     vs the corresponding class of the all-data multiclass_clean training universe.
   - 26 descriptors are frozen from the previous E2/E3 characterization.
   - Median/IQR, Cliff's delta, Mann–Whitney U, BH-FDR and bootstrap median-delta CI.

2. Individual level
   - Every predicted lytic sequence is placed against the relevant empirical
     training distribution descriptor by descriptor.
   - Empirical percentile, robust z-score, reference 95% coverage,
     fraction of descriptors inside the training range, and strongest
     descriptor deviations are retained.

No descriptor selection is performed in A2.
No model is re-trained or re-calibrated.

Default execution from repository root
--------------------------------------
python analysis_explored_seqs/A2_property_characterization.py

Expected inputs
---------------
analysis_explored_seqs/A0_integrate_external_sources/external_sources_master.tsv.gz
data/tasks/global_mixed/dataset.tsv.gz
data/tasks/multiclass_clean/dataset.tsv.gz

Outputs
-------
analysis_explored_seqs/A2_property_characterization/
├── tables/
│   ├── frozen26_descriptor_registry.tsv
│   ├── frozen26_descriptor_cache.tsv.gz
│   ├── external_frozen26_descriptors.tsv.gz
│   ├── reference_frozen26_descriptors.tsv.gz
│   ├── reference_distribution_profiles.tsv
│   ├── group_property_comparisons.tsv
│   ├── individual_global_reference_fit.tsv.gz
│   ├── individual_global_top_deviations.tsv.gz
│   ├── individual_functional_reference_fit.tsv.gz
│   ├── individual_functional_top_deviations.tsv.gz
│   └── source_property_fit_summary.tsv
├── figures/
│   ├── E1_cliffs_delta_heatmap.{png,pdf}
│   ├── E2_median_iqr_shift_heatmap.{png,pdf}
│   ├── E3_global_individual_reference_fit.{png,pdf}
│   └── E4_functional_individual_reference_fit.{png,pdf}
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
from scipy.stats import mannwhitneyu


# ---------------------------------------------------------------------
# Frozen descriptor contract
# ---------------------------------------------------------------------

FROZEN_FEATURES = [
    "global_basic_polarity_std",
    "global_basic_fcr",
    "global_basic_charged_fraction",
    "grouped_charged_frac",
    "grouped_ratio_charged_uncharged",
    "global_basic_polarity_mean",
    "grouped_positive_frac",
    "aac_freq_G",
    "aac_freq_S",
    "global_basic_acceptors_per_residue",
    "aaindex_FAUJ880104_mean",
    "structural_propensity_helix_mean",
    "aac_freq_E",
    "compositional_bias_positive_negative_skew",
    "aac_freq_T",
    "grouped_hydroxyl_frac",
    "entropy_complexity_unique_kmer_fraction_k2",
    "grouped_negative_frac",
    "grouped_small_frac",
    "aac_freq_N",
    "grouped_disorder_promoting_frac",
    "aac_freq_I",
    "aac_freq_W",
    "entropy_complexity_repeated_kmer_fraction_k2",
    "entropy_complexity_linguistic_complexity_k3",
    "entropy_complexity_unique_kmer_fraction_k3",
]

FEATURE_TO_DESCRIPTOR = {
    # global_basic
    "global_basic_polarity_std": "global_basic",
    "global_basic_fcr": "global_basic",
    "global_basic_charged_fraction": "global_basic",
    "global_basic_polarity_mean": "global_basic",
    "global_basic_acceptors_per_residue": "global_basic",

    # grouped
    "grouped_charged_frac": "grouped",
    "grouped_ratio_charged_uncharged": "grouped",
    "grouped_positive_frac": "grouped",
    "grouped_hydroxyl_frac": "grouped",
    "grouped_negative_frac": "grouped",
    "grouped_small_frac": "grouped",
    "grouped_disorder_promoting_frac": "grouped",

    # AAC
    "aac_freq_G": "aac",
    "aac_freq_S": "aac",
    "aac_freq_E": "aac",
    "aac_freq_T": "aac",
    "aac_freq_N": "aac",
    "aac_freq_I": "aac",
    "aac_freq_W": "aac",

    # other frozen families
    "aaindex_FAUJ880104_mean": "aaindex",
    "structural_propensity_helix_mean": "structural_propensity",
    "compositional_bias_positive_negative_skew": "compositional_bias",
    "entropy_complexity_unique_kmer_fraction_k2": "entropy_complexity",
    "entropy_complexity_repeated_kmer_fraction_k2": "entropy_complexity",
    "entropy_complexity_linguistic_complexity_k3": "entropy_complexity",
    "entropy_complexity_unique_kmer_fraction_k3": "entropy_complexity",
}

DESCRIPTOR_KEYS = [
    "global_basic",
    "grouped",
    "aac",
    "aaindex",
    "structural_propensity",
    "compositional_bias",
    "entropy_complexity",
]

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

FUNCTION_ORDER = ["endolysin", "depolymerase", "val"]

FUNCTION_LABELS = {
    "endolysin": "Endolysin",
    "depolymerase": "Depolymerase",
    "val": "VAL",
}

AA20 = set("ACDEFGHIKLMNPQRSTVWY")

sns.set_theme(style="whitegrid", context="paper")


# ---------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="A2 property characterization.")

    p.add_argument(
        "--external-master",
        type=Path,
        default=Path(
            "analysis_explored_seqs/"
            "A0_integrate_external_sources/"
            "external_sources_master.tsv.gz"
        ),
    )
    p.add_argument(
        "--global-training",
        type=Path,
        default=Path("data/tasks/global_mixed/dataset.tsv.gz"),
    )
    p.add_argument(
        "--multiclass-training",
        type=Path,
        default=Path("data/tasks/multiclass_clean/dataset.tsv.gz"),
    )
    p.add_argument(
        "--output-dir",
        type=Path,
        default=Path("analysis_explored_seqs/A2_property_characterization"),
    )
    p.add_argument(
        "--batch-size",
        type=int,
        default=2000,
        help="Roxy descriptor batch size.",
    )
    p.add_argument(
        "--bootstrap-reps",
        type=int,
        default=500,
        help="Bootstrap replicates for median-difference confidence intervals.",
    )
    p.add_argument(
        "--bootstrap-max-n",
        type=int,
        default=1000,
        help=(
            "Maximum observations sampled from each group per bootstrap "
            "replicate. Full observed medians/effect sizes still use all rows."
        ),
    )
    p.add_argument(
        "--seed",
        type=int,
        default=42,
    )

    return p.parse_args()


# ---------------------------------------------------------------------
# General utilities
# ---------------------------------------------------------------------

def require_file(path: Path) -> Path:
    if not path.is_file():
        raise FileNotFoundError(f"Required file not found: {path}")
    return path


def sha256_sequence(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("utf-8")).hexdigest()


def normalize_sequence(value: object) -> str:
    return "".join(str(value).split()).upper()


def eligible(sequence: str) -> bool:
    return (
        bool(sequence)
        and len(sequence) <= 1024
        and set(sequence).issubset(AA20)
    )


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


def bh_fdr(pvalues: pd.Series) -> pd.Series:
    p = pd.to_numeric(pvalues, errors="coerce").to_numpy(dtype=float)
    out = np.full(len(p), np.nan, dtype=float)

    mask = np.isfinite(p)
    if not mask.any():
        return pd.Series(out, index=pvalues.index)

    vals = p[mask]
    order = np.argsort(vals)
    ranked = vals[order]

    m = len(ranked)
    adjusted = ranked * m / np.arange(1, m + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)

    restored = np.empty(m, dtype=float)
    restored[order] = adjusted
    out[np.where(mask)[0]] = restored

    return pd.Series(out, index=pvalues.index)


def deterministic_seed(text: str, base_seed: int) -> int:
    digest = hashlib.sha256(
        f"{base_seed}|{text}".encode("utf-8")
    ).digest()
    return int.from_bytes(digest[:4], "little")


def save_figure(fig: plt.Figure, stem: Path) -> None:
    fig.tight_layout()
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------
# Training reference construction
# ---------------------------------------------------------------------

def prepare_external(path: Path) -> pd.DataFrame:
    df = pd.read_csv(
        path,
        sep="\t",
        compression="infer",
        low_memory=False,
    )

    required = {
        "analysis_id",
        "sequence_id",
        "source",
        "sequence",
        "input_status",
        "global_score",
        "global_prediction",
        "functional_prediction",
    }
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(f"External master missing columns: {missing}")

    df = df[df["input_status"] == "valid"].copy()
    df["sequence"] = df["sequence"].map(normalize_sequence)

    bad = ~df["sequence"].map(eligible)
    if bad.any():
        raise ValueError(
            f"External master contains {int(bad.sum())} sequences outside "
            "canonical20/<=1024 eligibility."
        )

    df["sequence_sha256_a2"] = df["sequence"].map(sha256_sequence)
    df["functional_class_a2"] = df["functional_prediction"].map(
        canonical_function
    )

    return df


def prepare_global_reference(path: Path) -> tuple[pd.DataFrame, dict]:
    df = pd.read_csv(
        path,
        sep="\t",
        compression="infer",
        low_memory=False,
    )

    required = {"entity_id", "sequence", "label"}
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(f"global_mixed missing columns: {missing}")

    n_input = len(df)

    # Positive class actually seen by the production global model.
    df = df[pd.to_numeric(df["label"], errors="coerce") == 1].copy()
    n_positive_before_eligibility = len(df)

    df["sequence"] = df["sequence"].map(normalize_sequence)
    df = df[df["sequence"].map(eligible)].copy()
    df["sequence_sha256_a2"] = df["sequence"].map(sha256_sequence)

    df = df.drop_duplicates("sequence_sha256_a2").copy()

    df["reference_cohort"] = "global_positive"
    df["reference_class"] = "lytic"

    report = {
        "n_input_global_mixed": int(n_input),
        "n_positive_before_eligibility": int(n_positive_before_eligibility),
        "n_positive_eligible_unique": int(len(df)),
    }

    return df, report


def prepare_multiclass_reference(path: Path) -> tuple[pd.DataFrame, dict]:
    df = pd.read_csv(
        path,
        sep="\t",
        compression="infer",
        low_memory=False,
    )

    required = {"entity_id", "sequence", "class_name"}
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(f"multiclass_clean missing columns: {missing}")

    n_input = len(df)

    df["sequence"] = df["sequence"].map(normalize_sequence)
    df = df[df["sequence"].map(eligible)].copy()

    df["reference_class"] = df["class_name"].map(canonical_function)
    df = df[df["reference_class"].isin(FUNCTION_ORDER)].copy()

    df["sequence_sha256_a2"] = df["sequence"].map(sha256_sequence)
    df = df.drop_duplicates(
        ["reference_class", "sequence_sha256_a2"]
    ).copy()

    df["reference_cohort"] = (
        "multiclass_" + df["reference_class"].astype(str)
    )

    report = {
        "n_input_multiclass_clean": int(n_input),
        "n_eligible_recognized": int(len(df)),
        "class_counts": {
            k: int(v)
            for k, v in df["reference_class"].value_counts().to_dict().items()
        },
    }

    return df, report


# ---------------------------------------------------------------------
# Roxy frozen descriptor calculation
# ---------------------------------------------------------------------

def descriptor_registry_table() -> pd.DataFrame:
    rows = []
    for i, feature in enumerate(FROZEN_FEATURES, start=1):
        rows.append(
            {
                "frozen_order": i,
                "feature": feature,
                "roxy_descriptor": FEATURE_TO_DESCRIPTOR[feature],
            }
        )
    return pd.DataFrame(rows)


def compute_roxy_descriptors(
    sequence_table: pd.DataFrame,
    batch_size: int,
) -> tuple[pd.DataFrame, str]:
    """
    sequence_table requires:
      sequence_sha256_a2
      sequence

    Returns one row per unique sequence hash with the exact 26 frozen features.
    """
    try:
        import roxy
        from roxy import DESCRIPTOR_REGISTRY
    except ImportError as exc:
        raise RuntimeError(
            "Roxy is required for A2. Install/activate the same Roxy "
            "environment used for the frozen descriptor analysis."
        ) from exc

    roxy_version = getattr(roxy, "__version__", "unknown")

    missing_keys = [k for k in DESCRIPTOR_KEYS if k not in DESCRIPTOR_REGISTRY]
    if missing_keys:
        raise RuntimeError(
            "Required Roxy descriptors are not registered: "
            + ", ".join(missing_keys)
        )

    seqs = (
        sequence_table[
            ["sequence_sha256_a2", "sequence"]
        ]
        .drop_duplicates("sequence_sha256_a2")
        .sort_values("sequence_sha256_a2")
        .reset_index(drop=True)
    )

    descriptors = {
        key: DESCRIPTOR_REGISTRY[key]()
        for key in DESCRIPTOR_KEYS
    }

    pieces = []

    for start in range(0, len(seqs), batch_size):
        stop = min(start + batch_size, len(seqs))
        batch = seqs.iloc[start:stop]

        ids = batch["sequence_sha256_a2"].astype(str).tolist()
        sequences = batch["sequence"].astype(str).tolist()

        combined = pd.DataFrame(
            {"sequence_sha256_a2": ids}
        )

        for descriptor_key, descriptor in descriptors.items():
            out = descriptor.compute(sequences, ids=ids)

            # Roxy returns Polars DataFrames in the frozen release.
            if hasattr(out, "to_dicts"):
                out_pd = pd.DataFrame(out.to_dicts())
            else:
                out_pd = pd.DataFrame(out)

            if "id" not in out_pd.columns:
                raise RuntimeError(
                    f"Roxy descriptor {descriptor_key} returned no id column."
                )

            expected = [
                f
                for f in FROZEN_FEATURES
                if FEATURE_TO_DESCRIPTOR[f] == descriptor_key
            ]

            missing_features = [
                f for f in expected if f not in out_pd.columns
            ]
            if missing_features:
                available = ", ".join(map(str, out_pd.columns[:30]))
                raise RuntimeError(
                    f"Roxy descriptor {descriptor_key} does not expose the "
                    f"frozen features {missing_features}. "
                    f"First available columns: {available}"
                )

            selected = out_pd[["id", *expected]].copy()
            selected = selected.rename(
                columns={"id": "sequence_sha256_a2"}
            )
            selected["sequence_sha256_a2"] = selected[
                "sequence_sha256_a2"
            ].astype(str)

            combined = combined.merge(
                selected,
                on="sequence_sha256_a2",
                how="left",
                validate="one_to_one",
            )

        pieces.append(combined)

        print(
            f"Roxy frozen26: {stop:,}/{len(seqs):,}",
            flush=True,
        )

    result = pd.concat(pieces, ignore_index=True)

    missing_final = [f for f in FROZEN_FEATURES if f not in result.columns]
    if missing_final:
        raise RuntimeError(
            f"Frozen descriptor calculation incomplete: {missing_final}"
        )

    for feature in FROZEN_FEATURES:
        result[feature] = pd.to_numeric(
            result[feature],
            errors="coerce",
        )

    return result, roxy_version


def build_descriptor_cache(
    all_sequences: pd.DataFrame,
    cache_path: Path,
    batch_size: int,
) -> tuple[pd.DataFrame, str]:
    needed = (
        all_sequences[
            ["sequence_sha256_a2", "sequence"]
        ]
        .drop_duplicates("sequence_sha256_a2")
        .copy()
    )

    if cache_path.is_file():
        cache = pd.read_csv(
            cache_path,
            sep="\t",
            compression="infer",
            low_memory=False,
        )

        required_cache = {"sequence_sha256_a2", *FROZEN_FEATURES}
        if not required_cache.issubset(cache.columns):
            print(
                "Existing descriptor cache is incompatible; rebuilding.",
                flush=True,
            )
            cache = pd.DataFrame()
    else:
        cache = pd.DataFrame()

    cached_ids = (
        set(cache["sequence_sha256_a2"].astype(str))
        if not cache.empty
        else set()
    )

    missing = needed[
        ~needed["sequence_sha256_a2"].astype(str).isin(cached_ids)
    ].copy()

    if len(missing):
        print(
            f"Descriptor cache: computing {len(missing):,} missing unique "
            f"sequences.",
            flush=True,
        )
        new, roxy_version = compute_roxy_descriptors(
            missing,
            batch_size=batch_size,
        )

        if cache.empty:
            cache = new
        else:
            cache = pd.concat([cache, new], ignore_index=True)
            cache = cache.drop_duplicates(
                "sequence_sha256_a2",
                keep="last",
            )

        cache.to_csv(
            cache_path,
            sep="\t",
            index=False,
            compression="gzip",
        )
    else:
        print("Descriptor cache: all required sequences already present.")
        try:
            import roxy
            roxy_version = getattr(roxy, "__version__", "unknown")
        except ImportError:
            roxy_version = "not_imported_cache_only"

    needed_ids = set(needed["sequence_sha256_a2"].astype(str))
    final = cache[
        cache["sequence_sha256_a2"].astype(str).isin(needed_ids)
    ].copy()

    if len(set(final["sequence_sha256_a2"])) != len(needed_ids):
        raise RuntimeError("Descriptor cache does not cover all requested sequences.")

    return final, roxy_version


# ---------------------------------------------------------------------
# Reference distributions
# ---------------------------------------------------------------------

def describe_reference_distribution(
    values: np.ndarray,
) -> dict[str, float]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]

    if len(values) == 0:
        return {
            "n": 0,
            "q025": np.nan,
            "q25": np.nan,
            "median": np.nan,
            "q75": np.nan,
            "q975": np.nan,
            "iqr": np.nan,
            "mad": np.nan,
            "robust_scale": np.nan,
        }

    q025, q25, median, q75, q975 = np.quantile(
        values,
        [0.025, 0.25, 0.50, 0.75, 0.975],
    )
    iqr = q75 - q25
    mad = np.median(np.abs(values - median))

    # IQR/1.349 approximates sigma under a Gaussian distribution.
    robust_scale = iqr / 1.349 if iqr > 0 else mad * 1.4826
    if robust_scale <= 0:
        robust_scale = np.nan

    return {
        "n": int(len(values)),
        "q025": float(q025),
        "q25": float(q25),
        "median": float(median),
        "q75": float(q75),
        "q975": float(q975),
        "iqr": float(iqr),
        "mad": float(mad),
        "robust_scale": float(robust_scale),
    }


def build_reference_profiles(
    reference: pd.DataFrame,
) -> pd.DataFrame:
    rows = []

    cohort_specs = [
        ("global_positive", "lytic"),
        ("multiclass_endolysin", "endolysin"),
        ("multiclass_depolymerase", "depolymerase"),
        ("multiclass_val", "val"),
    ]

    for cohort, cls in cohort_specs:
        g = reference[
            reference["reference_cohort"] == cohort
        ]

        if g.empty:
            raise RuntimeError(f"Reference cohort is empty: {cohort}")

        for feature in FROZEN_FEATURES:
            stats = describe_reference_distribution(
                pd.to_numeric(g[feature], errors="coerce").to_numpy()
            )
            rows.append(
                {
                    "reference_cohort": cohort,
                    "reference_class": cls,
                    "feature": feature,
                    **stats,
                }
            )

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------
# Group-level statistics
# ---------------------------------------------------------------------

def bootstrap_median_difference(
    query: np.ndarray,
    reference: np.ndarray,
    reps: int,
    max_n: int,
    seed: int,
) -> tuple[float, float, int, int]:
    q = np.asarray(query, dtype=float)
    r = np.asarray(reference, dtype=float)

    q = q[np.isfinite(q)]
    r = r[np.isfinite(r)]

    if reps <= 0 or len(q) == 0 or len(r) == 0:
        return np.nan, np.nan, 0, 0

    qn = min(len(q), max_n)
    rn = min(len(r), max_n)

    rng = np.random.default_rng(seed)

    # Deterministic cap for tractable bootstrap on very large cohorts.
    q_base = (
        rng.choice(q, size=qn, replace=False)
        if len(q) > qn
        else q
    )
    r_base = (
        rng.choice(r, size=rn, replace=False)
        if len(r) > rn
        else r
    )

    deltas = np.empty(reps, dtype=float)

    for i in range(reps):
        qb = rng.choice(q_base, size=len(q_base), replace=True)
        rb = rng.choice(r_base, size=len(r_base), replace=True)
        deltas[i] = np.median(qb) - np.median(rb)

    low, high = np.quantile(deltas, [0.025, 0.975])
    return float(low), float(high), int(len(q_base)), int(len(r_base))


def compare_groups(
    query: pd.DataFrame,
    reference: pd.DataFrame,
    comparison_id: str,
    comparison_level: str,
    source: str,
    predicted_class: str,
    reference_cohort: str,
    bootstrap_reps: int,
    bootstrap_max_n: int,
    base_seed: int,
) -> pd.DataFrame:
    rows = []

    for feature in FROZEN_FEATURES:
        q = pd.to_numeric(query[feature], errors="coerce").to_numpy(dtype=float)
        r = pd.to_numeric(reference[feature], errors="coerce").to_numpy(dtype=float)

        q = q[np.isfinite(q)]
        r = r[np.isfinite(r)]

        if len(q) == 0 or len(r) == 0:
            continue

        q25_q, med_q, q75_q = np.quantile(q, [0.25, 0.50, 0.75])
        q25_r, med_r, q75_r = np.quantile(r, [0.25, 0.50, 0.75])

        test = mannwhitneyu(
            q,
            r,
            alternative="two-sided",
            method="asymptotic",
        )

        # For Mann-Whitney U defined on query first:
        # Cliff's delta = 2U/(n_query*n_ref) - 1.
        cliffs_delta = (
            2.0 * float(test.statistic) / (len(q) * len(r)) - 1.0
        )

        ref_iqr = q75_r - q25_r
        median_delta = med_q - med_r
        median_iqr_shift = (
            median_delta / ref_iqr
            if ref_iqr > 0
            else np.nan
        )

        boot_seed = deterministic_seed(
            f"{comparison_id}|{feature}",
            base_seed,
        )

        ci_low, ci_high, boot_qn, boot_rn = bootstrap_median_difference(
            q,
            r,
            reps=bootstrap_reps,
            max_n=bootstrap_max_n,
            seed=boot_seed,
        )

        rows.append(
            {
                "comparison_id": comparison_id,
                "comparison_level": comparison_level,
                "source": source,
                "predicted_class": predicted_class,
                "reference_cohort": reference_cohort,
                "feature": feature,
                "n_query": int(len(q)),
                "n_reference": int(len(r)),
                "query_q25": float(q25_q),
                "query_median": float(med_q),
                "query_q75": float(q75_q),
                "reference_q25": float(q25_r),
                "reference_median": float(med_r),
                "reference_q75": float(q75_r),
                "median_delta": float(median_delta),
                "median_iqr_shift": float(median_iqr_shift),
                "cliffs_delta": float(cliffs_delta),
                "mannwhitney_u": float(test.statistic),
                "p_value": float(test.pvalue),
                "bootstrap_median_delta_ci_low": ci_low,
                "bootstrap_median_delta_ci_high": ci_high,
                "bootstrap_query_n": boot_qn,
                "bootstrap_reference_n": boot_rn,
            }
        )

    out = pd.DataFrame(rows)
    if not out.empty:
        out["q_value_bh"] = bh_fdr(out["p_value"])
    return out


def build_group_comparisons(
    external: pd.DataFrame,
    reference: pd.DataFrame,
    bootstrap_reps: int,
    bootstrap_max_n: int,
    seed: int,
) -> pd.DataFrame:
    comparisons = []

    global_ref = reference[
        reference["reference_cohort"] == "global_positive"
    ]

    # Global predicted lytic vs global positive training class.
    for source in SOURCE_ORDER:
        q = external[
            (external["source"] == source)
            & (external["global_prediction"] == "lytic")
        ]

        if q.empty:
            continue

        comparison_id = f"{source}__lytic_vs_global_positive"

        comparisons.append(
            compare_groups(
                q,
                global_ref,
                comparison_id=comparison_id,
                comparison_level="global",
                source=source,
                predicted_class="lytic",
                reference_cohort="global_positive",
                bootstrap_reps=bootstrap_reps,
                bootstrap_max_n=bootstrap_max_n,
                base_seed=seed,
            )
        )

    # Predicted functional subtype vs multiclass training subtype.
    for source in SOURCE_ORDER:
        for cls in FUNCTION_ORDER:
            q = external[
                (external["source"] == source)
                & (external["global_prediction"] == "lytic")
                & (external["functional_class_a2"] == cls)
            ]

            ref_cohort = f"multiclass_{cls}"
            r = reference[
                reference["reference_cohort"] == ref_cohort
            ]

            if q.empty or r.empty:
                continue

            comparison_id = f"{source}__{cls}_vs_{ref_cohort}"

            comparisons.append(
                compare_groups(
                    q,
                    r,
                    comparison_id=comparison_id,
                    comparison_level="functional",
                    source=source,
                    predicted_class=cls,
                    reference_cohort=ref_cohort,
                    bootstrap_reps=bootstrap_reps,
                    bootstrap_max_n=bootstrap_max_n,
                    base_seed=seed,
                )
            )

    if not comparisons:
        return pd.DataFrame()

    return pd.concat(comparisons, ignore_index=True)


# ---------------------------------------------------------------------
# Individual sequence vs reference distribution
# ---------------------------------------------------------------------

def individual_reference_fit(
    query: pd.DataFrame,
    reference_profiles: pd.DataFrame,
    reference_cohort: str,
    comparison_level: str,
    top_k: int = 3,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    profile = (
        reference_profiles[
            reference_profiles["reference_cohort"] == reference_cohort
        ]
        .set_index("feature")
    )

    missing = [f for f in FROZEN_FEATURES if f not in profile.index]
    if missing:
        raise RuntimeError(
            f"Reference profile {reference_cohort} missing features: {missing}"
        )

    rows = []
    top_rows = []

    for row in query.itertuples(index=False):
        robust_z = []
        within95 = []
        within_iqr = []
        feature_details = []

        for feature in FROZEN_FEATURES:
            value = getattr(row, feature)
            ref = profile.loc[feature]

            if pd.isna(value):
                continue

            value = float(value)
            scale = float(ref["robust_scale"])
            median = float(ref["median"])

            z = (
                (value - median) / scale
                if np.isfinite(scale) and scale > 0
                else np.nan
            )

            in95 = bool(
                value >= float(ref["q025"])
                and value <= float(ref["q975"])
            )
            iniqr = bool(
                value >= float(ref["q25"])
                and value <= float(ref["q75"])
            )

            robust_z.append(z)
            within95.append(in95)
            within_iqr.append(iniqr)

            feature_details.append(
                {
                    "feature": feature,
                    "value": value,
                    "reference_median": median,
                    "reference_q025": float(ref["q025"]),
                    "reference_q975": float(ref["q975"]),
                    "robust_z": z,
                    "abs_robust_z": abs(z) if np.isfinite(z) else np.nan,
                }
            )

        z_arr = np.asarray(robust_z, dtype=float)
        finite_z = z_arr[np.isfinite(z_arr)]

        summary = {
            "analysis_id": row.analysis_id,
            "sequence_id": row.sequence_id,
            "source": row.source,
            "global_score": row.global_score,
            "global_prediction": row.global_prediction,
            "functional_prediction": row.functional_prediction,
            "functional_class": row.functional_class_a2,
            "comparison_level": comparison_level,
            "reference_cohort": reference_cohort,
            "n_features_available": int(len(feature_details)),
            "fraction_within_reference_95": (
                float(np.mean(within95)) if within95 else np.nan
            ),
            "fraction_within_reference_iqr": (
                float(np.mean(within_iqr)) if within_iqr else np.nan
            ),
            "median_abs_robust_z": (
                float(np.median(np.abs(finite_z)))
                if len(finite_z)
                else np.nan
            ),
            "max_abs_robust_z": (
                float(np.max(np.abs(finite_z)))
                if len(finite_z)
                else np.nan
            ),
            "n_abs_robust_z_gt2": (
                int(np.sum(np.abs(finite_z) > 2))
                if len(finite_z)
                else 0
            ),
            "n_abs_robust_z_gt3": (
                int(np.sum(np.abs(finite_z) > 3))
                if len(finite_z)
                else 0
            ),
        }
        rows.append(summary)

        ranked = sorted(
            feature_details,
            key=lambda x: (
                -x["abs_robust_z"]
                if np.isfinite(x["abs_robust_z"])
                else np.inf
            ),
        )[:top_k]

        for rank, detail in enumerate(ranked, start=1):
            top_rows.append(
                {
                    "analysis_id": row.analysis_id,
                    "sequence_id": row.sequence_id,
                    "source": row.source,
                    "global_score": row.global_score,
                    "functional_prediction": row.functional_prediction,
                    "functional_class": row.functional_class_a2,
                    "comparison_level": comparison_level,
                    "reference_cohort": reference_cohort,
                    "deviation_rank": rank,
                    **detail,
                }
            )

    return pd.DataFrame(rows), pd.DataFrame(top_rows)


def build_individual_fits(
    external: pd.DataFrame,
    profiles: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    lytic = external[
        external["global_prediction"] == "lytic"
    ].copy()

    global_fit, global_top = individual_reference_fit(
        lytic,
        profiles,
        reference_cohort="global_positive",
        comparison_level="global",
    )

    functional_parts = []
    functional_top_parts = []

    for cls in FUNCTION_ORDER:
        q = lytic[
            lytic["functional_class_a2"] == cls
        ].copy()

        if q.empty:
            continue

        fit, top = individual_reference_fit(
            q,
            profiles,
            reference_cohort=f"multiclass_{cls}",
            comparison_level="functional",
        )

        functional_parts.append(fit)
        functional_top_parts.append(top)

    functional_fit = (
        pd.concat(functional_parts, ignore_index=True)
        if functional_parts
        else pd.DataFrame()
    )
    functional_top = (
        pd.concat(functional_top_parts, ignore_index=True)
        if functional_top_parts
        else pd.DataFrame()
    )

    return global_fit, global_top, functional_fit, functional_top


# ---------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------

def readable_feature(feature: str) -> str:
    replacements = {
        "global_basic_": "",
        "grouped_": "",
        "aac_freq_": "AAC ",
        "aaindex_": "AAIndex ",
        "structural_propensity_": "",
        "compositional_bias_": "",
        "entropy_complexity_": "",
        "_fraction": " fraction",
        "_frac": " fraction",
        "_mean": " mean",
        "_std": " std",
        "_": " ",
    }

    label = feature
    for old, new in replacements.items():
        label = label.replace(old, new)

    return label


def comparison_label(row: pd.Series) -> str:
    source = SOURCE_LABELS.get(row["source"], row["source"])
    cls = row["predicted_class"]

    if cls == "lytic":
        return f"{source} — Lytic"

    return f"{source} — {FUNCTION_LABELS.get(cls, cls)}"


def plot_effect_heatmap(
    stats: pd.DataFrame,
    value_col: str,
    title: str,
    output_stem: Path,
    center: float = 0.0,
) -> None:
    meta = (
        stats[
            [
                "comparison_id",
                "source",
                "predicted_class",
            ]
        ]
        .drop_duplicates()
        .copy()
    )
    meta["comparison_label"] = meta.apply(comparison_label, axis=1)

    label_map = dict(
        zip(meta["comparison_id"], meta["comparison_label"])
    )

    matrix = stats.pivot(
        index="comparison_id",
        columns="feature",
        values=value_col,
    )

    desired_rows = []

    # global first
    for source in SOURCE_ORDER:
        cid = f"{source}__lytic_vs_global_positive"
        if cid in matrix.index:
            desired_rows.append(cid)

    # then subtypes
    for cls in FUNCTION_ORDER:
        for source in SOURCE_ORDER:
            cid = f"{source}__{cls}_vs_multiclass_{cls}"
            if cid in matrix.index:
                desired_rows.append(cid)

    matrix = matrix.reindex(
        index=desired_rows,
        columns=FROZEN_FEATURES,
    )

    matrix.index = [label_map.get(x, x) for x in matrix.index]
    matrix.columns = [readable_feature(x) for x in matrix.columns]

    fig_width = max(13.0, 0.55 * len(FROZEN_FEATURES))
    fig_height = max(5.5, 0.48 * len(matrix))

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    sns.heatmap(
        matrix,
        center=center,
        cmap="vlag",
        linewidths=0.25,
        linecolor="white",
        cbar_kws={"shrink": 0.8},
        ax=ax,
    )

    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title(title)
    ax.tick_params(axis="x", rotation=70)
    ax.tick_params(axis="y", rotation=0)

    save_figure(fig, output_stem)


def plot_global_fit(
    global_fit: pd.DataFrame,
    output_stem: Path,
) -> None:
    plot_df = global_fit.copy()
    plot_df["source_label"] = plot_df["source"].map(SOURCE_LABELS)

    fig, ax = plt.subplots(figsize=(7.2, 4.8))

    sns.boxplot(
        data=plot_df,
        x="source_label",
        y="fraction_within_reference_95",
        order=[SOURCE_LABELS[s] for s in SOURCE_ORDER],
        showfliers=False,
        ax=ax,
    )

    ax.set_xlabel("")
    ax.set_ylabel("Fraction of 26 descriptors within reference 95% range")
    ax.set_ylim(0, 1.02)
    ax.set_title("Individual compatibility with global positive training space")
    ax.tick_params(axis="x", rotation=15)
    sns.despine(ax=ax)

    save_figure(fig, output_stem)


def plot_functional_fit(
    functional_fit: pd.DataFrame,
    output_stem: Path,
) -> None:
    plot_df = functional_fit.copy()
    plot_df["source_label"] = plot_df["source"].map(SOURCE_LABELS)
    plot_df["function_label"] = plot_df["functional_class"].map(
        FUNCTION_LABELS
    )

    fig, ax = plt.subplots(figsize=(9.0, 5.2))

    sns.boxplot(
        data=plot_df,
        x="function_label",
        y="fraction_within_reference_95",
        hue="source_label",
        order=[FUNCTION_LABELS[x] for x in FUNCTION_ORDER],
        hue_order=[SOURCE_LABELS[s] for s in SOURCE_ORDER],
        showfliers=False,
        ax=ax,
    )

    ax.set_xlabel("")
    ax.set_ylabel("Fraction of 26 descriptors within subtype reference 95% range")
    ax.set_ylim(0, 1.02)
    ax.set_title("Individual compatibility with functional training classes")
    ax.legend(title="", frameon=False)
    sns.despine(ax=ax)

    save_figure(fig, output_stem)


# ---------------------------------------------------------------------
# Compact summaries
# ---------------------------------------------------------------------

def source_fit_summary(
    global_fit: pd.DataFrame,
    functional_fit: pd.DataFrame,
) -> pd.DataFrame:
    rows = []

    for source in SOURCE_ORDER:
        g = global_fit[global_fit["source"] == source]

        if len(g):
            rows.append(
                {
                    "comparison_level": "global",
                    "source": source,
                    "predicted_class": "lytic",
                    "n_sequences": int(len(g)),
                    "median_fraction_within_reference_95": float(
                        g["fraction_within_reference_95"].median()
                    ),
                    "median_fraction_within_reference_iqr": float(
                        g["fraction_within_reference_iqr"].median()
                    ),
                    "median_of_median_abs_robust_z": float(
                        g["median_abs_robust_z"].median()
                    ),
                    "fraction_with_zero_extreme_features_z_gt3": float(
                        (g["n_abs_robust_z_gt3"] == 0).mean()
                    ),
                }
            )

    for source in SOURCE_ORDER:
        for cls in FUNCTION_ORDER:
            g = functional_fit[
                (functional_fit["source"] == source)
                & (functional_fit["functional_class"] == cls)
            ]

            if len(g):
                rows.append(
                    {
                        "comparison_level": "functional",
                        "source": source,
                        "predicted_class": cls,
                        "n_sequences": int(len(g)),
                        "median_fraction_within_reference_95": float(
                            g["fraction_within_reference_95"].median()
                        ),
                        "median_fraction_within_reference_iqr": float(
                            g["fraction_within_reference_iqr"].median()
                        ),
                        "median_of_median_abs_robust_z": float(
                            g["median_abs_robust_z"].median()
                        ),
                        "fraction_with_zero_extreme_features_z_gt3": float(
                            (g["n_abs_robust_z_gt3"] == 0).mean()
                        ),
                    }
                )

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def main() -> None:
    args = parse_args()

    external_path = require_file(args.external_master.resolve())
    global_path = require_file(args.global_training.resolve())
    multiclass_path = require_file(args.multiclass_training.resolve())

    outdir = args.output_dir.resolve()
    tables = outdir / "tables"
    figures = outdir / "figures"

    tables.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)

    # Frozen descriptor registry.
    registry = descriptor_registry_table()
    registry.to_csv(
        tables / "frozen26_descriptor_registry.tsv",
        sep="\t",
        index=False,
    )

    # Inputs and training-reference cohorts.
    external = prepare_external(external_path)
    global_ref, global_report = prepare_global_reference(global_path)
    multiclass_ref, multiclass_report = prepare_multiclass_reference(
        multiclass_path
    )

    reference_meta = pd.concat(
        [
            global_ref[
                [
                    "entity_id",
                    "sequence",
                    "sequence_sha256_a2",
                    "reference_cohort",
                    "reference_class",
                ]
            ],
            multiclass_ref[
                [
                    "entity_id",
                    "sequence",
                    "sequence_sha256_a2",
                    "reference_cohort",
                    "reference_class",
                ]
            ],
        ],
        ignore_index=True,
    )

    # Build one unique-sequence descriptor cache across external + reference.
    all_sequences = pd.concat(
        [
            external[
                ["sequence_sha256_a2", "sequence"]
            ],
            reference_meta[
                ["sequence_sha256_a2", "sequence"]
            ],
        ],
        ignore_index=True,
    )

    cache_path = tables / "frozen26_descriptor_cache.tsv.gz"

    descriptor_cache, roxy_version = build_descriptor_cache(
        all_sequences,
        cache_path=cache_path,
        batch_size=args.batch_size,
    )

    # Join external metadata -> descriptors.
    external_desc = external.merge(
        descriptor_cache,
        on="sequence_sha256_a2",
        how="left",
        validate="many_to_one",
    )

    # Join reference cohorts -> descriptors.
    reference_desc = reference_meta.merge(
        descriptor_cache,
        on="sequence_sha256_a2",
        how="left",
        validate="many_to_one",
    )

    # Hard closure: every sequence must have all frozen features columns.
    for name, frame in [
        ("external", external_desc),
        ("reference", reference_desc),
    ]:
        absent = [f for f in FROZEN_FEATURES if f not in frame.columns]
        if absent:
            raise RuntimeError(f"{name}: missing frozen features: {absent}")

    external_desc.to_csv(
        tables / "external_frozen26_descriptors.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    reference_desc.to_csv(
        tables / "reference_frozen26_descriptors.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    # Reference distribution profiles.
    profiles = build_reference_profiles(reference_desc)
    profiles.to_csv(
        tables / "reference_distribution_profiles.tsv",
        sep="\t",
        index=False,
    )

    # Population-level statistical comparisons.
    print("Running population-level property comparisons...", flush=True)

    comparisons = build_group_comparisons(
        external_desc,
        reference_desc,
        bootstrap_reps=args.bootstrap_reps,
        bootstrap_max_n=args.bootstrap_max_n,
        seed=args.seed,
    )

    comparisons.to_csv(
        tables / "group_property_comparisons.tsv",
        sep="\t",
        index=False,
    )

    # Individual sequence ↔ training-distribution compatibility.
    print("Running individual sequence/reference comparisons...", flush=True)

    (
        global_fit,
        global_top,
        functional_fit,
        functional_top,
    ) = build_individual_fits(
        external_desc,
        profiles,
    )

    global_fit.to_csv(
        tables / "individual_global_reference_fit.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    global_top.to_csv(
        tables / "individual_global_top_deviations.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    functional_fit.to_csv(
        tables / "individual_functional_reference_fit.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    functional_top.to_csv(
        tables / "individual_functional_top_deviations.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    summary = source_fit_summary(
        global_fit,
        functional_fit,
    )
    summary.to_csv(
        tables / "source_property_fit_summary.tsv",
        sep="\t",
        index=False,
    )

    # Candidate panel-E figures.
    if not comparisons.empty:
        plot_effect_heatmap(
            comparisons,
            value_col="cliffs_delta",
            title="Descriptor effect sizes relative to annotated training references",
            output_stem=figures / "E1_cliffs_delta_heatmap",
        )

        plot_effect_heatmap(
            comparisons,
            value_col="median_iqr_shift",
            title="Median descriptor shifts relative to training-reference IQR",
            output_stem=figures / "E2_median_iqr_shift_heatmap",
        )

    if not global_fit.empty:
        plot_global_fit(
            global_fit,
            figures / "E3_global_individual_reference_fit",
        )

    if not functional_fit.empty:
        plot_functional_fit(
            functional_fit,
            figures / "E4_functional_individual_reference_fit",
        )

    # Manifest.
    manifest = {
        "stage": "A2",
        "stage_name": "property_characterization",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "external_master": str(external_path),
            "global_training": str(global_path),
            "multiclass_training": str(multiclass_path),
        },
        "frozen_descriptor_count": len(FROZEN_FEATURES),
        "frozen_descriptors": FROZEN_FEATURES,
        "roxy_descriptor_keys": DESCRIPTOR_KEYS,
        "roxy_version": roxy_version,
        "reference_contract": {
            "global": (
                "Eligible label==1 sequences from the all-data global_mixed "
                "training universe."
            ),
            "functional": (
                "Eligible Endolysin, Depolymerase and VAL classes from the "
                "all-data multiclass_clean training universe."
            ),
        },
        "reference_counts": {
            "global": global_report,
            "multiclass": multiclass_report,
        },
        "n_external_sequences": int(len(external_desc)),
        "n_external_lytic": int(
            external_desc["global_prediction"].eq("lytic").sum()
        ),
        "bootstrap": {
            "replicates": int(args.bootstrap_reps),
            "maximum_rows_per_group_per_replicate": int(
                args.bootstrap_max_n
            ),
            "seed": int(args.seed),
            "note": (
                "Observed medians, Mann-Whitney tests and Cliff's delta use "
                "all finite observations. Bootstrap CIs use a deterministic "
                "cap for computational tractability when cohorts are large."
            ),
        },
        "individual_fit_metrics": [
            "fraction_within_reference_95",
            "fraction_within_reference_iqr",
            "median_abs_robust_z",
            "max_abs_robust_z",
            "n_abs_robust_z_gt2",
            "n_abs_robust_z_gt3",
            "top_3_descriptor_deviations",
        ],
        "statistics": [
            "median",
            "IQR",
            "Cliff's delta",
            "Mann-Whitney U",
            "Benjamini-Hochberg FDR",
            "bootstrap median-difference 95% CI",
        ],
        "status": "PASS",
    }

    (outdir / "manifest.json").write_text(
        json.dumps(manifest, indent=2),
        encoding="utf-8",
    )

    print("\nA2 — PASS")
    print("=" * 100)
    print(f"External valid sequences:   {len(external_desc):,}")
    print(
        f"External predicted lytic:   "
        f"{external_desc['global_prediction'].eq('lytic').sum():,}"
    )
    print(
        f"Global positive reference:  "
        f"{global_report['n_positive_eligible_unique']:,}"
    )
    print("Multiclass references:")
    for cls in FUNCTION_ORDER:
        n = multiclass_report["class_counts"].get(cls, 0)
        print(f"  {FUNCTION_LABELS[cls]:<15} {n:,}")
    print(f"Roxy version:               {roxy_version}")
    print(f"Group comparisons:          {comparisons['comparison_id'].nunique():,}")
    print(f"Individual global fits:     {len(global_fit):,}")
    print(f"Individual subtype fits:    {len(functional_fit):,}")
    print("=" * 100)
    print(summary.to_string(index=False))
    print("=" * 100)
    print(f"Tables:  {tables}")
    print(f"Figures: {figures}")
    print(f"Output:  {outdir}")


if __name__ == "__main__":
    main()

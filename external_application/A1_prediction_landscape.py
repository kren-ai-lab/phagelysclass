#!/usr/bin/env python3
"""
A1 — Prediction landscape across explored sequence sources.

Input
-----
analysis_explored_seqs/A0_integrate_external_sources/external_sources_master.tsv.gz

Outputs
-------
analysis_explored_seqs/A1_prediction_landscape/
├── tables/
├── figures/
└── manifest.json

Scope
-----
B. Global lytic/non-lytic landscape and global scores
C. Functional predictions and auxiliary compatibility scores
D. Source-specific analyses
   - Prediction Extension: BioSieve/MMseqs2 cluster-size strata, when available
   - ZymCTRL: EC number and catalytic block
   - Ambiguous Structural: protein-name-derived structural categories

This stage is descriptive. It does not re-fit models, alter thresholds, or
perform property/geometry analyses.
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

FUNCTION_ORDER = ["endolysin", "depolymerase", "VAL"]

AUX_COLUMNS = [
    "endolysin_compatibility_score",
    "depolymerase_compatibility_score",
    "val_compatibility_score",
]

AUX_LABELS = {
    "endolysin_compatibility_score": "Endolysin compatibility",
    "depolymerase_compatibility_score": "Depolymerase compatibility",
    "val_compatibility_score": "VAL compatibility",
}


sns.set_theme(style="whitegrid", context="paper")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="A1 prediction landscape.")
    p.add_argument(
        "--input",
        type=Path,
        default=Path(
            "analysis_explored_seqs/"
            "A0_integrate_external_sources/"
            "external_sources_master.tsv.gz"
        ),
    )
    p.add_argument(
        "--output-dir",
        type=Path,
        default=Path("analysis_explored_seqs/A1_prediction_landscape"),
    )
    return p.parse_args()


def require_columns(df: pd.DataFrame, columns: list[str]) -> None:
    missing = [c for c in columns if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")


def save_figure(fig: plt.Figure, stem: Path) -> None:
    fig.tight_layout()
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def pct(x: pd.Series) -> pd.Series:
    total = x.sum()
    if total == 0:
        return x.astype(float) * np.nan
    return 100.0 * x / total


def score_summary(
    df: pd.DataFrame,
    group_cols: list[str],
    score_cols: list[str],
) -> pd.DataFrame:
    rows = []
    for keys, g in df.groupby(group_cols, dropna=False, sort=False):
        if not isinstance(keys, tuple):
            keys = (keys,)
        base = dict(zip(group_cols, keys))

        for score in score_cols:
            values = pd.to_numeric(g[score], errors="coerce").dropna()
            if values.empty:
                continue

            rows.append(
                {
                    **base,
                    "score": score,
                    "n": int(len(values)),
                    "mean": float(values.mean()),
                    "sd": float(values.std(ddof=1)) if len(values) > 1 else np.nan,
                    "q05": float(values.quantile(0.05)),
                    "q25": float(values.quantile(0.25)),
                    "median": float(values.median()),
                    "q75": float(values.quantile(0.75)),
                    "q95": float(values.quantile(0.95)),
                }
            )

    return pd.DataFrame(rows)


def validate_prediction_contract(df: pd.DataFrame) -> dict:
    valid = df["input_status"].eq("valid")
    x = df.loc[valid].copy()

    global_score = pd.to_numeric(x["global_score"], errors="coerce")
    outside_global = int(((global_score < 0) | (global_score > 1)).sum())

    aux_outside = {}
    for col in AUX_COLUMNS:
        v = pd.to_numeric(x[col], errors="coerce")
        aux_outside[col] = int(((v < 0) | (v > 1)).sum())

    lytic = x["global_prediction"].eq("lytic")
    non = x["global_prediction"].eq("non_lytic_background")

    missing_functional_on_lytic = int(
        x.loc[lytic, "functional_prediction"].fillna("").eq("").sum()
    )
    functional_on_non_lytic = int(
        (~x.loc[non, "functional_prediction"].fillna("").eq("")).sum()
    )

    prob_cols = [
        "p_endolysin_given_lytic",
        "p_depolymerase_given_lytic",
        "p_val_given_lytic",
    ]

    probs_lytic = x.loc[lytic, prob_cols].apply(pd.to_numeric, errors="coerce")
    if len(probs_lytic):
        prob_sums = probs_lytic.sum(axis=1)
        bad_prob_sum = int((~np.isclose(prob_sums, 1.0, atol=1e-6)).sum())
        missing_probs_on_lytic = int(probs_lytic.isna().any(axis=1).sum())
    else:
        bad_prob_sum = 0
        missing_probs_on_lytic = 0

    probs_non = x.loc[non, prob_cols]
    probs_reported_on_non = int(probs_non.notna().any(axis=1).sum())

    payload = {
        "n_rows": int(len(df)),
        "n_valid": int(valid.sum()),
        "n_invalid": int((~valid).sum()),
        "global_score_outside_0_1": outside_global,
        "auxiliary_score_outside_0_1": aux_outside,
        "missing_functional_prediction_on_lytic": missing_functional_on_lytic,
        "functional_prediction_reported_on_non_lytic": functional_on_non_lytic,
        "missing_conditional_probabilities_on_lytic": missing_probs_on_lytic,
        "conditional_probability_sum_not_1_on_lytic": bad_prob_sum,
        "conditional_probabilities_reported_on_non_lytic": probs_reported_on_non,
    }

    hard_errors = (
        outside_global
        + sum(aux_outside.values())
        + missing_functional_on_lytic
        + functional_on_non_lytic
        + missing_probs_on_lytic
        + bad_prob_sum
        + probs_reported_on_non
    )
    payload["status"] = "PASS" if hard_errors == 0 else "WARN"
    return payload


def build_global_tables(valid: pd.DataFrame, tables: Path) -> dict[str, Path]:
    outputs = {}

    # Global lytic/non-lytic counts.
    global_counts = (
        valid.groupby(["source", "global_prediction"])
        .size()
        .rename("n")
        .reset_index()
    )
    global_counts["fraction_within_source"] = (
        global_counts["n"]
        / global_counts.groupby("source")["n"].transform("sum")
    )
    path = tables / "B1_global_prediction_by_source.tsv"
    global_counts.to_csv(path, sep="\t", index=False)
    outputs["global_prediction_by_source"] = path

    # Global score summary.
    global_scores = score_summary(valid, ["source"], ["global_score"])
    path = tables / "B2_global_score_summary_by_source.tsv"
    global_scores.to_csv(path, sep="\t", index=False)
    outputs["global_score_summary_by_source"] = path

    # Functional predictions, restricted to globally lytic sequences.
    lytic = valid[valid["global_prediction"] == "lytic"].copy()
    functional = (
        lytic.groupby(["source", "functional_prediction"])
        .size()
        .rename("n")
        .reset_index()
    )
    functional["fraction_among_lytic"] = (
        functional["n"]
        / functional.groupby("source")["n"].transform("sum")
    )
    path = tables / "C1_functional_prediction_by_source.tsv"
    functional.to_csv(path, sep="\t", index=False)
    outputs["functional_prediction_by_source"] = path

    # Conditional multiclass probabilities among lytic predictions.
    conditional = score_summary(
        lytic,
        ["source"],
        [
            "p_endolysin_given_lytic",
            "p_depolymerase_given_lytic",
            "p_val_given_lytic",
        ],
    )
    path = tables / "C2_conditional_probability_summary_by_source.tsv"
    conditional.to_csv(path, sep="\t", index=False)
    outputs["conditional_probability_summary_by_source"] = path

    # Auxiliary scores across all sequences.
    aux = score_summary(valid, ["source"], AUX_COLUMNS)
    path = tables / "C3_auxiliary_score_summary_by_source.tsv"
    aux.to_csv(path, sep="\t", index=False)
    outputs["auxiliary_score_summary_by_source"] = path

    # Score correlations.
    corr_rows = []
    score_cols = ["global_score", *AUX_COLUMNS]
    for source in SOURCE_ORDER:
        g = valid[valid["source"] == source][score_cols].apply(
            pd.to_numeric, errors="coerce"
        )
        corr = g.corr(method="spearman")
        for i, a in enumerate(score_cols):
            for b in score_cols[i + 1:]:
                corr_rows.append(
                    {
                        "source": source,
                        "score_a": a,
                        "score_b": b,
                        "spearman_rho": float(corr.loc[a, b]),
                        "n": int(g[[a, b]].dropna().shape[0]),
                    }
                )
    correlations = pd.DataFrame(corr_rows)
    path = tables / "C4_score_correlations_by_source.tsv"
    correlations.to_csv(path, sep="\t", index=False)
    outputs["score_correlations_by_source"] = path

    return outputs


def plot_global_prediction(valid: pd.DataFrame, figures: Path) -> Path:
    tab = (
        valid.groupby(["source", "global_prediction"])
        .size()
        .rename("n")
        .reset_index()
    )
    tab["percent"] = 100.0 * tab["n"] / tab.groupby("source")["n"].transform("sum")
    tab["source_label"] = tab["source"].map(SOURCE_LABELS)
    tab["prediction_label"] = tab["global_prediction"].replace({
        "non_lytic_background": "Non-lytic background",
        "lytic": "Lytic",
    })

    fig, ax = plt.subplots(figsize=(7.0, 4.8))
    sns.barplot(
        data=tab,
        x="source_label",
        y="percent",
        hue="prediction_label",
        order=[SOURCE_LABELS[s] for s in SOURCE_ORDER],
        hue_order=["Non-lytic background", "Lytic"],
        errorbar=None,
        ax=ax,
    )
    ax.set_xlabel("")
    ax.set_ylabel("Sequences (%)")
    ax.set_ylim(0, 100)
    ax.set_title("Global PhageLysClass predictions")
    ax.legend(title="", frameon=False)
    ax.tick_params(axis="x", rotation=15)
    sns.despine(ax=ax)

    stem = figures / "B1_global_prediction_fraction"
    save_figure(fig, stem)
    return stem


def plot_global_scores(valid: pd.DataFrame, figures: Path) -> Path:
    plot_df = valid[["source", "global_score"]].copy()
    plot_df["global_score"] = pd.to_numeric(plot_df["global_score"], errors="coerce")
    plot_df = plot_df.dropna(subset=["global_score"])
    plot_df["source_label"] = plot_df["source"].map(SOURCE_LABELS)

    fig, ax = plt.subplots(figsize=(7.0, 4.8))
    sns.boxplot(
        data=plot_df,
        x="source_label",
        y="global_score",
        order=[SOURCE_LABELS[s] for s in SOURCE_ORDER],
        showfliers=False,
        ax=ax,
    )
    ax.set_xlabel("")
    ax.set_ylabel("Global lytic score")
    ax.set_ylim(0, 1)
    ax.set_title("Global-score distributions")
    ax.tick_params(axis="x", rotation=15)
    sns.despine(ax=ax)

    stem = figures / "B2_global_score_distribution"
    save_figure(fig, stem)
    return stem


def plot_functional_composition(valid: pd.DataFrame, figures: Path) -> Path:
    lytic = valid[valid["global_prediction"] == "lytic"].copy()
    tab = (
        lytic.groupby(["source", "functional_prediction"])
        .size()
        .rename("n")
        .reset_index()
    )
    tab["percent"] = 100.0 * tab["n"] / tab.groupby("source")["n"].transform("sum")
    tab["source_label"] = tab["source"].map(SOURCE_LABELS)
    tab["function_label"] = tab["functional_prediction"].replace({
        "endolysin": "Endolysin",
        "depolymerase": "Depolymerase",
        "VAL": "VAL",
    })

    fig, ax = plt.subplots(figsize=(7.0, 4.8))
    sns.barplot(
        data=tab,
        x="source_label",
        y="percent",
        hue="function_label",
        order=[SOURCE_LABELS[s] for s in SOURCE_ORDER],
        hue_order=["Endolysin", "Depolymerase", "VAL"],
        errorbar=None,
        ax=ax,
    )
    ax.set_xlabel("")
    ax.set_ylabel("Predicted lytic sequences (%)")
    ax.set_ylim(0, 100)
    ax.set_title("Functional composition among lytic predictions")
    ax.legend(title="", frameon=False)
    ax.tick_params(axis="x", rotation=15)
    sns.despine(ax=ax)

    stem = figures / "C1_functional_composition"
    save_figure(fig, stem)
    return stem


def plot_auxiliary_scores(valid: pd.DataFrame, figures: Path) -> list[Path]:
    stems = []
    for col in AUX_COLUMNS:
        plot_df = valid[["source", col]].copy()
        plot_df[col] = pd.to_numeric(plot_df[col], errors="coerce")
        plot_df = plot_df.dropna(subset=[col])
        plot_df["source_label"] = plot_df["source"].map(SOURCE_LABELS)

        fig, ax = plt.subplots(figsize=(7.0, 4.8))
        sns.boxplot(
            data=plot_df,
            x="source_label",
            y=col,
            order=[SOURCE_LABELS[s] for s in SOURCE_ORDER],
            showfliers=False,
            ax=ax,
        )
        ax.set_xlabel("")
        ax.set_ylabel("Compatibility score")
        ax.set_ylim(0, 1)
        ax.set_title(AUX_LABELS[col])
        ax.tick_params(axis="x", rotation=15)
        sns.despine(ax=ax)

        stem = figures / f"C_auxiliary_{col}"
        save_figure(fig, stem)
        stems.append(stem)
    return stems


def build_zymctrl_analysis(
    valid: pd.DataFrame,
    tables: Path,
    figures: Path,
) -> dict:
    z = valid[valid["source"] == "zymctrl"].copy()
    outputs = {}

    if "zym_ec_number" not in z.columns:
        return {"status": "SKIPPED", "reason": "zym_ec_number not present"}

    z["zym_ec_number"] = z["zym_ec_number"].astype(str)

    rows = []
    for ec, g in z.groupby("zym_ec_number", sort=True):
        n = len(g)
        n_lytic = int(g["global_prediction"].eq("lytic").sum())
        row = {
            "ec_number": ec,
            "n": int(n),
            "n_lytic": n_lytic,
            "fraction_lytic": n_lytic / n if n else np.nan,
            "mean_global_score": float(g["global_score"].mean()),
            "median_global_score": float(g["global_score"].median()),
            "n_endolysin": int(g["functional_prediction"].eq("endolysin").sum()),
            "n_depolymerase": int(
                g["functional_prediction"].eq("depolymerase").sum()
            ),
            "n_val": int(g["functional_prediction"].eq("VAL").sum()),
        }

        for optional in [
            "zym_panel_block",
            "zym_class_relation",
            "zym_activity",
        ]:
            if optional in g.columns:
                values = sorted(g[optional].dropna().astype(str).unique())
                row[optional] = ";".join(values)

        rows.append(row)

    by_ec = pd.DataFrame(rows).sort_values("ec_number")
    path = tables / "D1_zymctrl_by_ec.tsv"
    by_ec.to_csv(path, sep="\t", index=False)
    outputs["by_ec"] = str(path)

    # EC x predicted subtype table.
    lytic = z[z["global_prediction"] == "lytic"]
    comp = (
        lytic.groupby(["zym_ec_number", "functional_prediction"])
        .size()
        .unstack(fill_value=0)
        .reindex(columns=FUNCTION_ORDER, fill_value=0)
    )
    comp_frac = comp.div(comp.sum(axis=1), axis=0)
    comp_out = comp.copy()
    for c in FUNCTION_ORDER:
        comp_out[f"fraction_{c}"] = comp_frac[c]
    comp_out = comp_out.reset_index()
    path = tables / "D2_zymctrl_functional_composition_by_ec.tsv"
    comp_out.to_csv(path, sep="\t", index=False)
    outputs["functional_by_ec"] = str(path)

    # Panel block summary when available.
    if "zym_panel_block" in z.columns:
        block_rows = []
        for block, g in z.groupby("zym_panel_block", dropna=False):
            n = len(g)
            n_lytic = int(g["global_prediction"].eq("lytic").sum())
            block_rows.append(
                {
                    "panel_block": block,
                    "n": int(n),
                    "n_lytic": n_lytic,
                    "fraction_lytic": n_lytic / n if n else np.nan,
                    "mean_global_score": float(g["global_score"].mean()),
                    "median_global_score": float(g["global_score"].median()),
                }
            )
        block_df = pd.DataFrame(block_rows)
        path = tables / "D3_zymctrl_by_panel_block.tsv"
        block_df.to_csv(path, sep="\t", index=False)
        outputs["by_panel_block"] = str(path)

    # Plot fraction lytic by EC.
    plot_df = by_ec.sort_values("fraction_lytic")
    plot_df = plot_df.copy()
    plot_df["percent_lytic"] = 100.0 * plot_df["fraction_lytic"]
    fig, ax = plt.subplots(figsize=(7.5, 5.2))
    sns.barplot(data=plot_df, y="ec_number", x="percent_lytic", order=plot_df["ec_number"].tolist(), errorbar=None, ax=ax)
    ax.set_xlabel("Predicted lytic (%)")
    ax.set_xlim(0, 100)
    ax.set_ylabel("EC number")
    ax.set_title("ZymCTRL lytic predictions by EC")
    stem = figures / "D1_zymctrl_lytic_fraction_by_ec"
    save_figure(fig, stem)
    outputs["figure_lytic_by_ec"] = str(stem)

    # Plot subtype composition among lytic by EC.
    frac = comp.div(comp.sum(axis=1), axis=0).fillna(0) * 100
    fig, ax = plt.subplots(figsize=(8.0, 5.6))
    bottom = np.zeros(len(frac))
    y = np.arange(len(frac))

    for c in FUNCTION_ORDER:
        vals = frac[c].to_numpy()
        ax.barh(
            y,
            vals,
            left=bottom,
            label="VAL" if c == "VAL" else c.capitalize(),
        )
        bottom += vals

    ax.set_yticks(y)
    ax.set_yticklabels(frac.index.astype(str))
    ax.set_xlabel("Predicted lytic sequences (%)")
    ax.set_xlim(0, 100)
    ax.set_ylabel("EC number")
    ax.set_title("ZymCTRL functional predictions by EC")
    ax.legend(frameon=False)
    stem = figures / "D2_zymctrl_functional_composition_by_ec"
    save_figure(fig, stem)
    outputs["figure_functional_by_ec"] = str(stem)

    outputs["status"] = "PASS"
    return outputs


def structural_category(name: object) -> str:
    text = "" if pd.isna(name) else str(name).lower()

    if "tailspike" in text or "tail spike" in text:
        return "Tailspike"
    if "tail fiber" in text or "tail fibre" in text:
        return "Tail fiber"
    if "baseplate" in text:
        return "Baseplate"
    if (
        "receptor" in text
        or "adsorption" in text
        or "adsorb" in text
        or "host recognition" in text
    ):
        return "Receptor-binding / adsorption"
    if "tail" in text:
        return "Other tail-associated"
    return "Other structural"


def build_ambiguous_analysis(
    valid: pd.DataFrame,
    tables: Path,
    figures: Path,
) -> dict:
    a = valid[valid["source"] == "ambiguous_structural"].copy()
    outputs = {}

    if "amb_protein_names" not in a.columns:
        return {"status": "SKIPPED", "reason": "amb_protein_names not present"}

    a["structural_category"] = a["amb_protein_names"].map(structural_category)

    rows = []
    for category, g in a.groupby("structural_category", sort=False):
        n = len(g)
        n_lytic = int(g["global_prediction"].eq("lytic").sum())

        rows.append(
            {
                "structural_category": category,
                "n": int(n),
                "n_lytic": n_lytic,
                "fraction_lytic": n_lytic / n if n else np.nan,
                "mean_global_score": float(g["global_score"].mean()),
                "median_global_score": float(g["global_score"].median()),
                "n_endolysin": int(
                    g["functional_prediction"].eq("endolysin").sum()
                ),
                "n_depolymerase": int(
                    g["functional_prediction"].eq("depolymerase").sum()
                ),
                "n_val": int(g["functional_prediction"].eq("VAL").sum()),
            }
        )

    category_df = pd.DataFrame(rows).sort_values("n", ascending=False)
    path = tables / "D4_ambiguous_structural_by_category.tsv"
    category_df.to_csv(path, sep="\t", index=False)
    outputs["by_category"] = str(path)

    # Preserve raw annotation detail as a complementary table.
    raw_rows = []
    for name, g in a.groupby("amb_protein_names", dropna=False):
        n = len(g)
        n_lytic = int(g["global_prediction"].eq("lytic").sum())
        raw_rows.append(
            {
                "protein_name": name,
                "n": int(n),
                "n_lytic": n_lytic,
                "fraction_lytic": n_lytic / n if n else np.nan,
                "n_endolysin": int(
                    g["functional_prediction"].eq("endolysin").sum()
                ),
                "n_depolymerase": int(
                    g["functional_prediction"].eq("depolymerase").sum()
                ),
                "n_val": int(g["functional_prediction"].eq("VAL").sum()),
                "median_global_score": float(g["global_score"].median()),
            }
        )

    raw_df = (
        pd.DataFrame(raw_rows)
        .sort_values(["n", "fraction_lytic"], ascending=[False, False])
    )
    path = tables / "D5_ambiguous_structural_by_protein_name.tsv"
    raw_df.to_csv(path, sep="\t", index=False)
    outputs["by_protein_name"] = str(path)

    # Category lytic fraction.
    plot_df = category_df.sort_values("fraction_lytic")
    plot_df = plot_df.copy()
    plot_df["percent_lytic"] = 100.0 * plot_df["fraction_lytic"]
    fig, ax = plt.subplots(figsize=(8.0, 5.2))
    sns.barplot(data=plot_df, y="structural_category", x="percent_lytic", order=plot_df["structural_category"].tolist(), errorbar=None, ax=ax)
    ax.set_xlabel("Predicted lytic (%)")
    ax.set_xlim(0, 100)
    ax.set_ylabel("Structural category")
    ax.set_title("Ambiguous structural challenge")
    stem = figures / "D3_ambiguous_lytic_fraction_by_category"
    save_figure(fig, stem)
    outputs["figure_lytic_by_category"] = str(stem)

    # Functional composition among lytic sequences by category.
    lytic = a[a["global_prediction"] == "lytic"]
    comp = (
        lytic.groupby(["structural_category", "functional_prediction"])
        .size()
        .unstack(fill_value=0)
        .reindex(columns=FUNCTION_ORDER, fill_value=0)
    )
    comp = comp.loc[comp.sum(axis=1).sort_values().index]
    frac = comp.div(comp.sum(axis=1), axis=0).fillna(0) * 100

    fig, ax = plt.subplots(figsize=(8.0, 5.4))
    bottom = np.zeros(len(frac))
    y = np.arange(len(frac))

    for c in FUNCTION_ORDER:
        vals = frac[c].to_numpy()
        ax.barh(
            y,
            vals,
            left=bottom,
            label="VAL" if c == "VAL" else c.capitalize(),
        )
        bottom += vals

    ax.set_yticks(y)
    ax.set_yticklabels(frac.index)
    ax.set_xlabel("Predicted lytic sequences (%)")
    ax.set_xlim(0, 100)
    ax.set_ylabel("Structural category")
    ax.set_title("Functional composition of structural lytic predictions")
    ax.legend(frameon=False)
    stem = figures / "D4_ambiguous_functional_composition_by_category"
    save_figure(fig, stem)
    outputs["figure_functional_by_category"] = str(stem)

    outputs["status"] = "PASS"
    return outputs


def build_prediction_extension_analysis(
    valid: pd.DataFrame,
    tables: Path,
    figures: Path,
) -> dict:
    p = valid[
        valid["source"] == "phagelys_prediction_extension"
    ].copy()

    if "pe_cluster_size" not in p.columns:
        return {
            "status": "SKIPPED",
            "reason": "pe_cluster_size not available from A0/BioSieve mapping",
        }

    p["pe_cluster_size"] = pd.to_numeric(
        p["pe_cluster_size"],
        errors="coerce",
    )

    if p["pe_cluster_size"].notna().sum() == 0:
        return {
            "status": "SKIPPED",
            "reason": "pe_cluster_size exists but contains no usable values",
        }

    bins = [0, 1, 5, 20, 100, np.inf]
    labels = ["1", "2–5", "6–20", "21–100", ">100"]

    p["cluster_size_bin"] = pd.cut(
        p["pe_cluster_size"],
        bins=bins,
        labels=labels,
        include_lowest=True,
    )

    rows = []
    for bin_name, g in p.groupby(
        "cluster_size_bin",
        observed=True,
        sort=False,
    ):
        n = len(g)
        n_lytic = int(g["global_prediction"].eq("lytic").sum())

        rows.append(
            {
                "cluster_size_bin": str(bin_name),
                "n_representatives": int(n),
                "n_lytic": n_lytic,
                "fraction_lytic": n_lytic / n if n else np.nan,
                "mean_global_score": float(g["global_score"].mean()),
                "median_global_score": float(g["global_score"].median()),
                "median_cluster_size": float(g["pe_cluster_size"].median()),
            }
        )

    summary = pd.DataFrame(rows)
    path = tables / "D6_prediction_extension_by_cluster_size.tsv"
    summary.to_csv(path, sep="\t", index=False)

    plot_df = summary.copy()
    plot_df["percent_lytic"] = 100.0 * plot_df["fraction_lytic"]
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    sns.barplot(data=plot_df, x="cluster_size_bin", y="percent_lytic", order=labels, errorbar=None, ax=ax)
    ax.set_xlabel("MMseqs2 cluster size")
    ax.set_ylabel("Predicted lytic (%)")
    ax.set_ylim(0, 100)
    ax.set_title("Prediction Extension by redundancy-cluster size")
    stem = figures / "D5_prediction_extension_lytic_by_cluster_size"
    save_figure(fig, stem)

    return {
        "status": "PASS",
        "table": str(path),
        "figure": str(stem),
    }


def main() -> None:
    args = parse_args()

    input_path = args.input.resolve()
    outdir = args.output_dir.resolve()
    tables = outdir / "tables"
    figures = outdir / "figures"

    if not input_path.is_file():
        raise FileNotFoundError(input_path)

    tables.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(input_path, sep="\t", compression="infer", low_memory=False)

    required = [
        "analysis_id",
        "sequence_id",
        "source",
        "sequence",
        "sequence_length",
        "input_status",
        "global_score",
        "global_prediction",
        "p_endolysin_given_lytic",
        "p_depolymerase_given_lytic",
        "p_val_given_lytic",
        "functional_prediction",
        *AUX_COLUMNS,
    ]
    require_columns(df, required)

    unknown_sources = sorted(set(df["source"].dropna()) - set(SOURCE_ORDER))
    if unknown_sources:
        raise ValueError(f"Unexpected source labels: {unknown_sources}")

    qc = validate_prediction_contract(df)
    (outdir / "prediction_contract_qc.json").write_text(
        json.dumps(qc, indent=2),
        encoding="utf-8",
    )

    valid = df[df["input_status"] == "valid"].copy()

    # Core B–C tables.
    table_outputs = build_global_tables(valid, tables)

    # Core B–C figures.
    figure_outputs = {
        "B1_global_prediction_fraction": str(
            plot_global_prediction(valid, figures)
        ),
        "B2_global_score_distribution": str(
            plot_global_scores(valid, figures)
        ),
        "C1_functional_composition": str(
            plot_functional_composition(valid, figures)
        ),
    }

    aux_figs = plot_auxiliary_scores(valid, figures)
    figure_outputs["C_auxiliary_figures"] = [str(x) for x in aux_figs]

    # Source-specific D analyses.
    zymctrl = build_zymctrl_analysis(valid, tables, figures)
    ambiguous = build_ambiguous_analysis(valid, tables, figures)
    prediction_extension = build_prediction_extension_analysis(
        valid,
        tables,
        figures,
    )

    # One compact source-level summary used later by A4/final panel builder.
    compact_rows = []
    for source in SOURCE_ORDER:
        g = valid[valid["source"] == source]
        lytic = g["global_prediction"].eq("lytic")
        compact_rows.append(
            {
                "source": source,
                "display_name": SOURCE_LABELS[source],
                "n": int(len(g)),
                "n_lytic": int(lytic.sum()),
                "fraction_lytic": float(lytic.mean()),
                "median_global_score": float(g["global_score"].median()),
                "n_endolysin": int(
                    g["functional_prediction"].eq("endolysin").sum()
                ),
                "n_depolymerase": int(
                    g["functional_prediction"].eq("depolymerase").sum()
                ),
                "n_val": int(g["functional_prediction"].eq("VAL").sum()),
            }
        )

    compact = pd.DataFrame(compact_rows)
    compact_path = tables / "A1_source_landscape_summary.tsv"
    compact.to_csv(compact_path, sep="\t", index=False)

    manifest = {
        "stage": "A1",
        "stage_name": "prediction_landscape",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "input": str(input_path),
        "n_rows_input": int(len(df)),
        "n_rows_valid": int(len(valid)),
        "prediction_contract_qc": qc,
        "analyses": {
            "global_landscape": "PASS",
            "functional_landscape": "PASS",
            "auxiliary_scores": "PASS",
            "zymctrl": zymctrl,
            "ambiguous_structural": ambiguous,
            "prediction_extension_cluster_analysis": prediction_extension,
        },
        "tables": {
            **{k: str(v) for k, v in table_outputs.items()},
            "compact_source_summary": str(compact_path),
        },
        "figures": figure_outputs,
        "notes": {
            "inference": (
                "Descriptive prediction analysis only; no models, thresholds, "
                "features, or representations are re-fitted."
            ),
            "ambiguous_structural_categories": (
                "Derived deterministically from amb_protein_names using simple "
                "keyword rules for tailspike, tail fiber, baseplate, "
                "receptor-binding/adsorption, other tail-associated, and other."
            ),
            "prediction_extension_cluster_analysis": (
                "Executed only when A0 recovered representative cluster sizes "
                "from the BioSieve/MMseqs2 mapping."
            ),
        },
    }

    (outdir / "manifest.json").write_text(
        json.dumps(manifest, indent=2),
        encoding="utf-8",
    )

    print("\nA1 — PASS")
    print("=" * 88)
    print(compact.to_string(index=False))
    print("=" * 88)
    print(f"Prediction contract: {qc['status']}")
    print(f"ZymCTRL analysis:    {zymctrl.get('status')}")
    print(f"Ambiguous analysis:  {ambiguous.get('status')}")
    print(
        "PE cluster analysis: "
        f"{prediction_extension.get('status')}"
    )
    if prediction_extension.get("status") == "SKIPPED":
        print(f"  reason: {prediction_extension.get('reason')}")
    print(f"Tables:  {tables}")
    print(f"Figures: {figures}")
    print(f"Output:  {outdir}")


if __name__ == "__main__":
    main()

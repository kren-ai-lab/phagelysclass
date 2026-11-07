#!/usr/bin/env python3
"""
A0 — Integrate external/explored sequence sources.

Run from the PhageLysClass repository root:

    python A0_integrate_external_sources.py

Default input layout
--------------------
using_model/
├── phagelys_prediction_extension_h30/
│   ├── 04_biosieve_mapping.csv
│   └── 05_phagelys_prediction_extension_H30_representatives.tsv.gz
├── phagelys_predict_H30_predictions.tsv
├── zymctrl_phagelys_panel/
│   └── zymctrl_phagelys_10000.tsv
├── zymctrl_phagelys_predictions/
│   └── zymctrl_predictions.tsv
├── ambiguous_structural_input/
│   └── ambiguous_structural_qc_pass.tsv.gz
└── ambiguous_structural_predictions.tsv

Default output
--------------
analysis_explored_seqs/
└── A0_integrate_external_sources/
    ├── external_sources_master.tsv.gz
    ├── prediction_extension_integrated.tsv.gz
    ├── zymctrl_integrated.tsv.gz
    ├── ambiguous_structural_integrated.tsv.gz
    ├── source_summary.tsv
    ├── prediction_summary.tsv
    ├── qc_summary.tsv
    └── manifest.json
"""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


PREDICTION_COLUMNS = [
    "sequence_id",
    "sequence_length",
    "input_status",
    "global_score",
    "global_prediction",
    "p_endolysin_given_lytic",
    "p_depolymerase_given_lytic",
    "p_val_given_lytic",
    "functional_prediction",
    "functional_prediction_status",
    "endolysin_compatibility_score",
    "depolymerase_compatibility_score",
    "val_compatibility_score",
    "warnings",
    "model_version",
    "representation",
]

CANONICAL_MASTER_COLUMNS = [
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
    "functional_prediction_status",
    "endolysin_compatibility_score",
    "depolymerase_compatibility_score",
    "val_compatibility_score",
    "warnings",
    "model_version",
    "representation",
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Integrate the three explored PhageLysClass sequence sources."
    )
    p.add_argument(
        "--using-model-dir",
        type=Path,
        default=Path("using_model"),
        help="Directory containing prepared sequence sets and predictions.",
    )
    p.add_argument(
        "--output-dir",
        type=Path,
        default=Path("analysis_explored_seqs/A0_integrate_external_sources"),
        help="A0 output directory.",
    )
    return p.parse_args()


def sha256_file(path: Path, block_size: int = 1024 * 1024) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(block_size)
            if not block:
                break
            h.update(block)
    return h.hexdigest()


def require_file(path: Path) -> Path:
    if not path.is_file():
        raise FileNotFoundError(f"Required file not found: {path}")
    return path


def read_table(path: Path) -> pd.DataFrame:
    suffixes = [s.lower() for s in path.suffixes]
    if ".tsv" in suffixes:
        return pd.read_csv(path, sep="\t", compression="infer")
    if ".csv" in suffixes:
        return pd.read_csv(path, compression="infer")
    if ".parquet" in suffixes:
        return pd.read_parquet(path)
    raise ValueError(f"Unsupported table format: {path}")


def normalize_sequence(series: pd.Series) -> pd.Series:
    return (
        series.fillna("")
        .astype(str)
        .str.replace(r"\s+", "", regex=True)
        .str.upper()
    )


def validate_predictions(df: pd.DataFrame, label: str) -> pd.DataFrame:
    missing = [c for c in PREDICTION_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"{label}: missing prediction columns: {missing}")

    if df["sequence_id"].isna().any():
        raise ValueError(f"{label}: missing sequence_id values.")

    df = df.copy()
    df["sequence_id"] = df["sequence_id"].astype(str)

    if df["sequence_id"].duplicated().any():
        examples = (
            df.loc[df["sequence_id"].duplicated(keep=False), "sequence_id"]
            .head(10)
            .tolist()
        )
        raise ValueError(f"{label}: duplicate prediction IDs. Examples: {examples}")

    return df


def merge_one_to_one(
    metadata: pd.DataFrame,
    predictions: pd.DataFrame,
    metadata_id: str,
    source: str,
) -> pd.DataFrame:
    if metadata_id not in metadata.columns:
        raise ValueError(f"{source}: metadata lacks ID column {metadata_id!r}")
    if "sequence" not in metadata.columns:
        raise ValueError(f"{source}: metadata lacks sequence column.")

    meta = metadata.copy()
    meta[metadata_id] = meta[metadata_id].astype(str)

    if meta[metadata_id].duplicated().any():
        examples = (
            meta.loc[meta[metadata_id].duplicated(keep=False), metadata_id]
            .head(10)
            .tolist()
        )
        raise ValueError(f"{source}: duplicate metadata IDs. Examples: {examples}")

    pred_ids = set(predictions["sequence_id"])
    meta_ids = set(meta[metadata_id])

    missing_predictions = sorted(meta_ids - pred_ids)
    extra_predictions = sorted(pred_ids - meta_ids)

    if missing_predictions or extra_predictions:
        raise RuntimeError(
            f"{source}: metadata/prediction ID mismatch. "
            f"missing_predictions={len(missing_predictions)}, "
            f"extra_predictions={len(extra_predictions)}. "
            f"Examples missing={missing_predictions[:5]}, "
            f"extra={extra_predictions[:5]}"
        )

    merged = meta.merge(
        predictions,
        left_on=metadata_id,
        right_on="sequence_id",
        how="inner",
        validate="one_to_one",
        suffixes=("_metadata", ""),
    )

    if metadata_id != "sequence_id":
        merged = merged.drop(columns=[metadata_id])

    merged.insert(1, "source", source)

    merged["sequence"] = normalize_sequence(merged["sequence"])
    merged["sequence_length"] = merged["sequence"].str.len().astype(int)

    return merged


def derive_cluster_sizes(mapping_path: Path) -> pd.DataFrame | None:
    """Best-effort extraction of representative cluster sizes."""
    if not mapping_path.is_file():
        return None

    m = read_table(mapping_path)
    if m.empty:
        return None

    lower = {str(c).lower(): c for c in m.columns}

    rep_candidates = [
        "representative_id",
        "representative",
        "rep_id",
        "cluster_representative",
        "cluster_center",
        "centroid_id",
    ]
    member_candidates = ["member_id", "member", "sequence_id", "id", "entity_id"]
    cluster_candidates = ["cluster_id", "cluster", "group_id"]

    rep_col = next((lower[x] for x in rep_candidates if x in lower), None)
    member_col = next((lower[x] for x in member_candidates if x in lower), None)
    cluster_col = next((lower[x] for x in cluster_candidates if x in lower), None)

    if rep_col is not None:
        out = (
            m.groupby(rep_col, dropna=False)
            .size()
            .rename("pe_cluster_size")
            .reset_index()
            .rename(columns={rep_col: "sequence_id"})
        )
        out["sequence_id"] = out["sequence_id"].astype(str)
        return out

    if cluster_col is not None and member_col is not None:
        flag_candidates = ["is_representative", "representative_flag", "is_rep"]
        flag_col = next((lower[x] for x in flag_candidates if x in lower), None)

        if flag_col is not None:
            sizes = (
                m.groupby(cluster_col, dropna=False)
                .size()
                .rename("pe_cluster_size")
                .reset_index()
            )
            reps = m[m[flag_col].astype(bool)][[cluster_col, member_col]].copy()
            reps = reps.rename(columns={member_col: "sequence_id"})
            reps["sequence_id"] = reps["sequence_id"].astype(str)
            return reps.merge(sizes, on=cluster_col, how="left").drop(
                columns=[cluster_col]
            )

    return None


def add_source_specific_prefixes(df: pd.DataFrame, source: str) -> pd.DataFrame:
    keep_unprefixed = set(CANONICAL_MASTER_COLUMNS)
    prefixes = {
        "phagelys_prediction_extension": "pe_",
        "zymctrl": "zym_",
        "ambiguous_structural": "amb_",
    }
    prefix = prefixes[source]

    rename = {}
    for col in df.columns:
        if col in keep_unprefixed:
            continue
        if col.startswith(prefix):
            continue
        rename[col] = f"{prefix}{col}"

    return df.rename(columns=rename)


def summarize_source(df: pd.DataFrame) -> dict:
    valid = df["input_status"].eq("valid")
    lytic = df["global_prediction"].eq("lytic")
    non = df["global_prediction"].eq("non_lytic_background")

    n_valid = int(valid.sum())
    n_lytic = int(lytic.sum())

    return {
        "source": df["source"].iloc[0],
        "n_sequences": int(len(df)),
        "n_valid": n_valid,
        "n_invalid": int((~valid).sum()),
        "n_lytic": n_lytic,
        "n_non_lytic_background": int(non.sum()),
        "fraction_lytic": n_lytic / n_valid if n_valid else np.nan,
        "n_endolysin": int(df["functional_prediction"].eq("endolysin").sum()),
        "n_depolymerase": int(df["functional_prediction"].eq("depolymerase").sum()),
        "n_val": int(df["functional_prediction"].eq("VAL").sum()),
        "mean_global_score": float(df["global_score"].mean()),
        "median_global_score": float(df["global_score"].median()),
    }


def build_prediction_summary(master: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for source, g in master.groupby("source", sort=False):
        total = len(g)
        for level, col in [
            ("global", "global_prediction"),
            ("functional", "functional_prediction"),
        ]:
            counts = g[col].fillna("").replace("", "not_reported").value_counts()
            for label, n in counts.items():
                rows.append(
                    {
                        "source": source,
                        "prediction_level": level,
                        "prediction": label,
                        "n": int(n),
                        "fraction_of_source": float(n / total) if total else np.nan,
                    }
                )
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()

    using = args.using_model_dir.resolve()
    outdir = args.output_dir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    pe_meta_path = require_file(
        using
        / "phagelys_prediction_extension_h30"
        / "05_phagelys_prediction_extension_H30_representatives.tsv.gz"
    )
    pe_pred_path = require_file(using / "phagelys_predict_H30_predictions.tsv")
    pe_mapping_path = (
        using
        / "phagelys_prediction_extension_h30"
        / "04_biosieve_mapping.csv"
    )

    zym_meta_path = require_file(
        using / "zymctrl_phagelys_panel" / "zymctrl_phagelys_10000.tsv"
    )
    zym_pred_path = require_file(
        using / "zymctrl_phagelys_predictions" / "zymctrl_predictions.tsv"
    )

    amb_meta_path = require_file(
        using
        / "ambiguous_structural_input"
        / "ambiguous_structural_qc_pass.tsv.gz"
    )
    amb_pred_path = require_file(using / "ambiguous_structural_predictions.tsv")

    pe_pred = validate_predictions(read_table(pe_pred_path), "Prediction Extension")
    zym_pred = validate_predictions(read_table(zym_pred_path), "ZymCTRL")
    amb_pred = validate_predictions(read_table(amb_pred_path), "Ambiguous Structural")

    pe_meta = read_table(pe_meta_path)
    pe = merge_one_to_one(
        pe_meta,
        pe_pred,
        metadata_id="sequence_entity_id",
        source="phagelys_prediction_extension",
    )

    cluster_sizes = derive_cluster_sizes(pe_mapping_path)
    if cluster_sizes is not None:
        pe = pe.merge(
            cluster_sizes,
            on="sequence_id",
            how="left",
            validate="one_to_one",
        )
    pe = add_source_specific_prefixes(pe, "phagelys_prediction_extension")

    zym_meta = read_table(zym_meta_path)
    zym = merge_one_to_one(
        zym_meta,
        zym_pred,
        metadata_id="sequence_id",
        source="zymctrl",
    )
    zym = add_source_specific_prefixes(zym, "zymctrl")

    amb_meta = read_table(amb_meta_path)
    amb = merge_one_to_one(
        amb_meta,
        amb_pred,
        metadata_id="entity_id",
        source="ambiguous_structural",
    )
    amb = add_source_specific_prefixes(amb, "ambiguous_structural")

    pe_out = outdir / "prediction_extension_integrated.tsv.gz"
    zym_out = outdir / "zymctrl_integrated.tsv.gz"
    amb_out = outdir / "ambiguous_structural_integrated.tsv.gz"

    pe.to_csv(pe_out, sep="\t", index=False, compression="gzip")
    zym.to_csv(zym_out, sep="\t", index=False, compression="gzip")
    amb.to_csv(amb_out, sep="\t", index=False, compression="gzip")

    all_columns = list(CANONICAL_MASTER_COLUMNS)
    for df in (pe, zym, amb):
        for col in df.columns:
            if col not in all_columns:
                all_columns.append(col)

    frames = []
    for df in (pe, zym, amb):
        x = df.copy()
        for col in all_columns:
            if col not in x.columns:
                x[col] = pd.NA
        frames.append(x[all_columns])

    master = pd.concat(frames, ignore_index=True)

    master.insert(
        0,
        "analysis_id",
        master["source"].astype(str) + "::" + master["sequence_id"].astype(str),
    )

    if master["analysis_id"].duplicated().any():
        raise RuntimeError("Duplicate analysis_id values after integration.")

    master_out = outdir / "external_sources_master.tsv.gz"
    master.to_csv(master_out, sep="\t", index=False, compression="gzip")

    summaries = pd.DataFrame(
        [
            summarize_source(pe),
            summarize_source(zym),
            summarize_source(amb),
        ]
    )

    valid_all = master["input_status"].eq("valid")
    lytic_all = master["global_prediction"].eq("lytic")

    total_row = {
        "source": "ALL",
        "n_sequences": int(len(master)),
        "n_valid": int(valid_all.sum()),
        "n_invalid": int((~valid_all).sum()),
        "n_lytic": int(lytic_all.sum()),
        "n_non_lytic_background": int(
            master["global_prediction"].eq("non_lytic_background").sum()
        ),
        "fraction_lytic": float(
            lytic_all.sum() / max(valid_all.sum(), 1)
        ),
        "n_endolysin": int(master["functional_prediction"].eq("endolysin").sum()),
        "n_depolymerase": int(
            master["functional_prediction"].eq("depolymerase").sum()
        ),
        "n_val": int(master["functional_prediction"].eq("VAL").sum()),
        "mean_global_score": float(master["global_score"].mean()),
        "median_global_score": float(master["global_score"].median()),
    }

    summaries = pd.concat(
        [summaries, pd.DataFrame([total_row])],
        ignore_index=True,
    )

    source_summary_path = outdir / "source_summary.tsv"
    summaries.to_csv(source_summary_path, sep="\t", index=False)

    prediction_summary = build_prediction_summary(master)
    prediction_summary_path = outdir / "prediction_summary.tsv"
    prediction_summary.to_csv(prediction_summary_path, sep="\t", index=False)

    qc_summary = (
        master.groupby(["source", "input_status"], dropna=False)
        .size()
        .rename("n")
        .reset_index()
    )
    qc_summary_path = outdir / "qc_summary.tsv"
    qc_summary.to_csv(qc_summary_path, sep="\t", index=False)

    inputs = {
        "prediction_extension_metadata": pe_meta_path,
        "prediction_extension_predictions": pe_pred_path,
        "prediction_extension_mapping": pe_mapping_path,
        "zymctrl_metadata": zym_meta_path,
        "zymctrl_predictions": zym_pred_path,
        "ambiguous_structural_metadata": amb_meta_path,
        "ambiguous_structural_predictions": amb_pred_path,
    }

    manifest = {
        "stage": "A0",
        "stage_name": "integrate_external_sources",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": [
            "phagelys_prediction_extension",
            "zymctrl",
            "ambiguous_structural",
        ],
        "n_sources": 3,
        "n_rows_master": int(len(master)),
        "input_files": {
            key: {
                "path": str(path.resolve()),
                "sha256": sha256_file(path) if path.is_file() else None,
            }
            for key, path in inputs.items()
        },
        "source_counts": {
            row["source"]: int(row["n_sequences"])
            for _, row in summaries[summaries["source"] != "ALL"].iterrows()
        },
        "cluster_size_added_for_prediction_extension": cluster_sizes is not None,
        "outputs": {
            "master": str(master_out),
            "prediction_extension": str(pe_out),
            "zymctrl": str(zym_out),
            "ambiguous_structural": str(amb_out),
            "source_summary": str(source_summary_path),
            "prediction_summary": str(prediction_summary_path),
            "qc_summary": str(qc_summary_path),
        },
    }

    (outdir / "manifest.json").write_text(
        json.dumps(manifest, indent=2),
        encoding="utf-8",
    )

    print("\nA0 — PASS")
    print("=" * 72)
    print(summaries.to_string(index=False))
    print("=" * 72)
    print(f"Master: {master_out}")
    print(f"Rows:   {len(master):,}")
    print(f"Cols:   {len(master.columns):,}")
    print(f"Output: {outdir}")


if __name__ == "__main__":
    main()

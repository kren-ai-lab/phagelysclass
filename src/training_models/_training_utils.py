#!/usr/bin/env python3
"""Shared utilities for PhageLysClass predictive training."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    confusion_matrix,
    f1_score,
    matthews_corrcoef,
    precision_score,
    recall_score,
    roc_auc_score,
)
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler


CANONICAL_AA = frozenset("ACDEFGHIKLMNPQRSTVWY")
MAX_LENGTH = 1024

REPRESENTATIONS = {
    "roxy_48": {
        "kind": "dense_parquet",
        "path": "data/representations/roxy_characterization/characterization_48.parquet",
        "feature_mode": "roxy48",
    },
    "onehot_1024": {
        "kind": "sparse_csr",
        "path": "data/representations/sylphy_onehot/onehot_csr.npz",
        "index": "data/representations/sylphy_onehot/row_index.tsv.gz",
    },
    "ankh2_ext1": {
        "kind": "dense_parquet",
        "path": "data/representations/sylphy_embeddings/models/elnaggarlab_ankh2-ext1/embeddings_last_mean_fp32.parquet",
        "feature_mode": "plm",
    },
    "esmc_300m": {
        "kind": "dense_parquet",
        "path": "data/representations/sylphy_embeddings/models/esmc_300m/embeddings_last_mean_fp32.parquet",
        "feature_mode": "plm",
    },
    "esm2_t6_8m": {
        "kind": "dense_parquet",
        "path": "data/representations/sylphy_embeddings/models/facebook_esm2_t6_8m_ur50d/embeddings_last_mean_fp32.parquet",
        "feature_mode": "plm",
    },
    "esm2_t12_35m": {
        "kind": "dense_parquet",
        "path": "data/representations/sylphy_embeddings/models/facebook_esm2_t12_35m_ur50d/embeddings_last_mean_fp32.parquet",
        "feature_mode": "plm",
    },
    "esm2_t30_150m": {
        "kind": "dense_parquet",
        "path": "data/representations/sylphy_embeddings/models/facebook_esm2_t30_150m_ur50d/embeddings_last_mean_fp32.parquet",
        "feature_mode": "plm",
    },
    "prot_bert": {
        "kind": "dense_parquet",
        "path": "data/representations/sylphy_embeddings/models/rostlab_prot_bert/embeddings_last_mean_fp32.parquet",
        "feature_mode": "plm",
    },
    "prot_t5_uniref50": {
        "kind": "dense_parquet",
        "path": "data/representations/sylphy_embeddings/models/rostlab_prot_t5_xl_uniref50/embeddings_last_mean_fp32.parquet",
        "feature_mode": "plm",
    },
}

ROXY_PREFIXES = (
    "aac_freq_",
    "global_basic_",
    "charge_",
    "order_disorder_",
    "structural_propensity_",
    "aaindex_",
)


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def sha256_file(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(chunk_size), b""):
            h.update(block)
    return h.hexdigest()


def read_table(path: Path) -> pd.DataFrame:
    name = path.name.lower()
    if name.endswith(".parquet"):
        return pd.read_parquet(path)
    if name.endswith(".tsv") or name.endswith(".tsv.gz"):
        return pd.read_csv(path, sep="\t", compression="infer", low_memory=False)
    if name.endswith(".csv") or name.endswith(".csv.gz"):
        return pd.read_csv(path, compression="infer", low_memory=False)
    raise ValueError(f"Unsupported table format: {path}")


def ensure_sequence_uid(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    if "sequence_uid" in df.columns:
        df["sequence_uid"] = df["sequence_uid"].astype(str)
        return df

    for candidate in ["sequence_entity_id", "entity_id"]:
        if candidate in df.columns:
            df["sequence_uid"] = df[candidate].astype(str)
            return df

    for candidate in ["sequence", "seq"]:
        if candidate in df.columns:
            df["sequence_uid"] = (
                df[candidate].astype(str).str.replace(r"\s+", "", regex=True).str.upper().map(sha256_text)
            )
            return df

    raise KeyError(
        "Could not resolve sequence_uid. Expected sequence_uid, sequence_entity_id, "
        "entity_id, sequence, or seq."
    )


def common_eligible_ids(stage10a_path: Path) -> set[str]:
    df = read_table(stage10a_path)
    df = ensure_sequence_uid(df)

    if "sequence" not in df.columns:
        if "seq" in df.columns:
            df["sequence"] = df["seq"]
        else:
            raise KeyError("Stage-10A table requires a sequence column.")

    seq = (
        df["sequence"].astype(str)
        .str.replace(r"\s+", "", regex=True)
        .str.upper()
    )
    length_ok = seq.str.len().le(MAX_LENGTH)
    canonical_ok = seq.map(lambda s: set(s).issubset(CANONICAL_AA))

    return set(df.loc[length_ok & canonical_ok, "sequence_uid"].astype(str))


def resolve_column(
    df: pd.DataFrame,
    explicit: str | None,
    candidates: Iterable[str],
    what: str,
) -> str:
    if explicit:
        if explicit not in df.columns:
            raise KeyError(f"{what} column {explicit!r} not found.")
        return explicit

    for col in candidates:
        if col in df.columns:
            return col

    raise KeyError(
        f"Could not infer {what} column. Candidates: {list(candidates)}; "
        f"available={list(df.columns)}"
    )


def normalize_outer_split(values: pd.Series) -> pd.Series:
    mapping = {
        "train": "train",
        "training": "train",
        "outer_train": "train",
        "0": "train",
        "test": "test",
        "testing": "test",
        "outer_test": "test",
        "1": "test",
    }
    out = values.astype(str).str.strip().str.lower().map(mapping)
    if out.isna().any():
        bad = sorted(values.loc[out.isna()].astype(str).unique().tolist())[:20]
        raise ValueError(f"Unrecognized outer split values: {bad}")
    return out


def load_outer_split(
    path: Path,
    split_col: str | None = None,
    label_col: str | None = None,
) -> tuple[pd.DataFrame, str]:
    df = read_table(path)
    df = ensure_sequence_uid(df)

    split_col = resolve_column(
        df,
        split_col,
        ["split", "outer_split", "partition", "set", "subset"],
        "outer split",
    )
    label_col = resolve_column(
        df,
        label_col,
        ["label", "target", "y", "class_label", "observed_label", "pu_observed_label", "class_name"],
        "label",
    )

    out = df[["sequence_uid", split_col, label_col]].copy()
    out.rename(columns={split_col: "outer_split", label_col: "raw_label"}, inplace=True)
    out["sequence_uid"] = out["sequence_uid"].astype(str)
    out["outer_split"] = normalize_outer_split(out["outer_split"])

    if out["sequence_uid"].duplicated().any():
        raise ValueError("Outer split table contains duplicate sequence_uid values.")

    return out, label_col


def encode_binary_label(
    values: pd.Series,
    positive_label: str | None = None,
) -> pd.Series:
    s = values.astype(str).str.strip()

    numeric = pd.to_numeric(s, errors="coerce")
    if numeric.notna().all() and set(numeric.astype(int).unique()).issubset({0, 1}):
        return numeric.astype(int)

    if positive_label is not None:
        return s.eq(str(positive_label)).astype(int)

    positive_tokens = {
        "1", "positive", "pos", "p", "p_core", "true", "yes",
    }
    negative_tokens = {
        "0", "negative", "neg", "n", "n_easy", "n_hard", "false", "no",
    }

    lowered = s.str.lower()
    unique = set(lowered.unique())
    if unique.issubset(positive_tokens | negative_tokens):
        return lowered.isin(positive_tokens).astype(int)

    raise ValueError(
        "Binary labels are not 0/1 and no --positive-label was supplied. "
        f"Observed values: {sorted(s.unique().tolist())[:20]}"
    )


def create_or_load_inner_split(
    *,
    outer_train: pd.DataFrame,
    y: pd.Series,
    output_path: Path,
    seed: int = 42,
    test_size: float = 0.20,
) -> pd.DataFrame:
    ids = outer_train["sequence_uid"].astype(str).reset_index(drop=True)
    y = pd.Series(y).reset_index(drop=True)

    if len(ids) != len(y):
        raise ValueError("Inner split IDs and labels have different lengths.")

    if output_path.exists():
        frozen = pd.read_csv(output_path, sep="\t", compression="infer", dtype=str)
        if set(frozen["sequence_uid"]) != set(ids):
            raise ValueError(
                f"Existing inner split {output_path} does not match current eligible outer-train IDs."
            )
        return frozen[["sequence_uid", "inner_split"]].copy()

    train_ids, val_ids = train_test_split(
        ids.to_numpy(),
        test_size=test_size,
        random_state=seed,
        stratify=y.to_numpy(),
    )

    frozen = pd.DataFrame({
        "sequence_uid": np.concatenate([train_ids, val_ids]),
        "inner_split": ["train"] * len(train_ids) + ["validation"] * len(val_ids),
    })
    frozen.sort_values("sequence_uid", inplace=True, kind="stable")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    frozen.to_csv(output_path, sep="\t", index=False, compression="gzip")
    return frozen


def feature_columns_for_dense(df: pd.DataFrame, mode: str) -> list[str]:
    if mode == "plm":
        cols = [c for c in df.columns if str(c).startswith("p_")]
    elif mode == "roxy48":
        cols = [
            c for c in df.columns
            if any(str(c).startswith(prefix) for prefix in ROXY_PREFIXES)
        ]
    else:
        raise ValueError(f"Unknown dense feature mode: {mode}")

    if not cols:
        raise ValueError(f"No features detected for dense mode={mode}.")
    return cols


def load_representation(
    *,
    project_root: Path,
    representation: str,
    requested_ids: list[str],
):
    if representation not in REPRESENTATIONS:
        raise KeyError(
            f"Unknown representation {representation!r}. "
            f"Available: {sorted(REPRESENTATIONS)}"
        )

    spec = REPRESENTATIONS[representation]

    if spec["kind"] == "dense_parquet":
        path = project_root / spec["path"]
        if not path.exists():
            raise FileNotFoundError(path)

        df = pd.read_parquet(path)
        df = ensure_sequence_uid(df)

        if df["sequence_uid"].duplicated().any():
            raise ValueError(f"Duplicate sequence_uid in representation {representation}.")

        feature_cols = feature_columns_for_dense(df, spec["feature_mode"])
        indexed = df.set_index("sequence_uid", drop=False)

        missing = [uid for uid in requested_ids if uid not in indexed.index]
        if missing:
            raise KeyError(
                f"{representation}: {len(missing):,} requested IDs are missing. "
                f"Example: {missing[:5]}"
            )

        subset = indexed.loc[requested_ids, feature_cols]
        X = subset.to_numpy(dtype=np.float32, copy=True)
        return X, feature_cols, False, path

    if spec["kind"] == "sparse_csr":
        matrix_path = project_root / spec["path"]
        index_path = project_root / spec["index"]

        if not matrix_path.exists():
            raise FileNotFoundError(matrix_path)
        if not index_path.exists():
            raise FileNotFoundError(index_path)

        matrix = sparse.load_npz(matrix_path).tocsr()
        index_df = pd.read_csv(index_path, sep="\t", compression="infer")

        if "onehot_row" not in index_df.columns or "sequence_uid" not in index_df.columns:
            raise KeyError("One-hot row index requires onehot_row and sequence_uid.")

        row_map = dict(
            zip(index_df["sequence_uid"].astype(str), index_df["onehot_row"].astype(int))
        )

        missing = [uid for uid in requested_ids if uid not in row_map]
        if missing:
            raise KeyError(
                f"onehot_1024: {len(missing):,} requested IDs are missing. Example: {missing[:5]}"
            )

        rows = np.asarray([row_map[uid] for uid in requested_ids], dtype=np.int64)
        X = matrix[rows]
        return X, None, True, matrix_path

    raise ValueError(f"Unsupported representation kind: {spec['kind']}")


def should_standardize(representation: str, algorithm: str) -> bool:
    return representation != "onehot_1024" and algorithm in {"LR", "KNN", "MLP"}


def fit_transform_scaler(X_train, X_other, *, standardize: bool):
    if not standardize:
        return X_train, X_other, None

    scaler = StandardScaler()
    X_train_t = scaler.fit_transform(X_train)
    X_other_t = scaler.transform(X_other)
    return X_train_t, X_other_t, scaler


def transform_with_scaler(X, scaler):
    return scaler.transform(X) if scaler is not None else X


def binary_metrics(y_true, y_pred, y_score) -> dict:
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred, labels=[0, 1]).ravel()
    prevalence = float(np.mean(y_true))

    out = {
        "mcc": float(matthews_corrcoef(y_true, y_pred)),
        "average_precision": float(average_precision_score(y_true, y_score)),
        "ap_baseline": prevalence,
        "ap_gain": float(average_precision_score(y_true, y_score) - prevalence),
        "f1": float(f1_score(y_true, y_pred, zero_division=0)),
        "precision": float(precision_score(y_true, y_pred, zero_division=0)),
        "recall": float(recall_score(y_true, y_pred, zero_division=0)),
        "specificity": float(tn / (tn + fp)) if (tn + fp) else np.nan,
        "balanced_accuracy": float(balanced_accuracy_score(y_true, y_pred)),
        "prevalence": prevalence,
        "predicted_positive_fraction": float(np.mean(y_pred)),
        "tp": int(tp),
        "tn": int(tn),
        "fp": int(fp),
        "fn": int(fn),
    }

    try:
        out["roc_auc"] = float(roc_auc_score(y_true, y_score))
    except Exception:
        out["roc_auc"] = np.nan

    return out


def multiclass_metrics(y_true, y_pred) -> dict:
    return {
        "mcc": float(matthews_corrcoef(y_true, y_pred)),
        "macro_f1": float(f1_score(y_true, y_pred, average="macro", zero_division=0)),
        "weighted_f1": float(f1_score(y_true, y_pred, average="weighted", zero_division=0)),
        "macro_precision": float(
            precision_score(y_true, y_pred, average="macro", zero_division=0)
        ),
        "macro_recall": float(
            recall_score(y_true, y_pred, average="macro", zero_division=0)
        ),
        "balanced_accuracy": float(balanced_accuracy_score(y_true, y_pred)),
    }


def write_json(data: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")

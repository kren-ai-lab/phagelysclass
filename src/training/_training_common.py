#!/usr/bin/env python3
"""
Shared utilities for PhageLysClass Stage 11 training.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import platform
import resource
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.base import clone
from sklearn.compose import TransformedTargetRegressor
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    brier_score_loss,
    confusion_matrix,
    f1_score,
    log_loss,
    matthews_corrcoef,
    precision_recall_fscore_support,
    precision_score,
    recall_score,
    roc_auc_score,
)
from sklearn.neighbors import KNeighborsClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

try:
    from xgboost import XGBClassifier
except Exception:
    XGBClassifier = None


SEED = 42
R_VALUES = [90, 70, 50, 30]
S_VALUES = [90, 70, 50, 30]

TASKS = {
    "global_easy": {
        "path": "global_easy",
        "type": "binary",
    },
    "global_hard": {
        "path": "global_hard",
        "type": "binary",
    },
    "global_mixed": {
        "path": "global_mixed",
        "type": "binary",
    },
    "subtype_endolysin": {
        "path": "subtype_binary/endolysin",
        "type": "binary",
    },
    "subtype_depolymerase": {
        "path": "subtype_binary/depolymerase",
        "type": "binary",
    },
    "subtype_val": {
        "path": "subtype_binary/virion_associated_lytic_enzyme",
        "type": "binary",
    },
    "multiclass_clean": {
        "path": "multiclass_clean",
        "type": "multiclass",
    },
    "pu_conservative": {
        "path": "pu_global/conservative",
        "type": "pu",
    },
    "pu_broad": {
        "path": "pu_global/broad",
        "type": "pu",
    },
}

SUPERVISED_TASKS = [
    k for k, v in TASKS.items() if v["type"] in {"binary", "multiclass"}
]
PU_TASKS = [k for k, v in TASKS.items() if v["type"] == "pu"]

ALGORITHMS = ["LR", "RF", "KNN", "XGB", "MLP"]

PARAMETER_PANELS = {
    "LR": [
        {"config_id": "LR01", "C": 0.1},
        {"config_id": "LR02", "C": 1.0},
        {"config_id": "LR03", "C": 10.0},
    ],
    "RF": [
        {
            "config_id": "RF01",
            "n_estimators": 300,
            "max_depth": None,
            "min_samples_leaf": 1,
            "max_features": "sqrt",
        },
        {
            "config_id": "RF02",
            "n_estimators": 500,
            "max_depth": None,
            "min_samples_leaf": 2,
            "max_features": "sqrt",
        },
        {
            "config_id": "RF03",
            "n_estimators": 500,
            "max_depth": 30,
            "min_samples_leaf": 1,
            "max_features": "sqrt",
        },
    ],
    "KNN": [
        {"config_id": "KNN01", "n_neighbors": 3},
        {"config_id": "KNN02", "n_neighbors": 5},
        {"config_id": "KNN03", "n_neighbors": 11},
    ],
    "XGB": [
        {
            "config_id": "XGB01",
            "n_estimators": 300,
            "max_depth": 4,
            "learning_rate": 0.10,
            "subsample": 0.8,
            "colsample_bytree": 0.8,
        },
        {
            "config_id": "XGB02",
            "n_estimators": 500,
            "max_depth": 6,
            "learning_rate": 0.05,
            "subsample": 0.8,
            "colsample_bytree": 0.8,
        },
        {
            "config_id": "XGB03",
            "n_estimators": 800,
            "max_depth": 4,
            "learning_rate": 0.03,
            "subsample": 0.8,
            "colsample_bytree": 0.8,
        },
    ],
    "MLP": [
        {
            "config_id": "MLP01",
            "hidden_layer_sizes": (256,),
            "alpha": 1e-4,
            "learning_rate_init": 1e-3,
        },
        {
            "config_id": "MLP02",
            "hidden_layer_sizes": (512, 128),
            "alpha": 1e-4,
            "learning_rate_init": 5e-4,
        },
        {
            "config_id": "MLP03",
            "hidden_layer_sizes": (256, 128, 64),
            "alpha": 1e-3,
            "learning_rate_init": 5e-4,
        },
    ],
}

PU_FIXED_CONFIG = {
    "LR": "LR02",
    "RF": "RF02",
    "XGB": "XGB02",
    "MLP": "MLP02",
}


def project_root_from_script(script_file: str) -> Path:
    # expected: <repo>/src/training/<script>.py
    return Path(script_file).resolve().parents[2]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(
        path,
        sep="\t",
        compression="infer",
        dtype=str,
        low_memory=False,
    ).fillna("")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, default=json_default) + "\n")


def json_default(obj):
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, tuple):
        return list(obj)
    raise TypeError(type(obj).__name__)


def process_peak_rss_kb() -> int:
    return int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)


def software_versions() -> dict[str, str]:
    import sklearn
    import scipy

    versions = {
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
        "scikit_learn": sklearn.__version__,
    }
    try:
        import xgboost
        versions["xgboost"] = xgboost.__version__
    except Exception:
        versions["xgboost"] = "unavailable"
    return versions


def n_jobs() -> int:
    raw = os.environ.get("PHAGELYS_N_JOBS", "-1")
    try:
        return int(raw)
    except ValueError:
        return -1


@dataclass
class RepresentationSpec:
    representation_id: str
    family: str
    format: str
    main_path: Path
    index_path: Path | None
    expected_dimension: int | None
    status: str


class RepresentationStore:
    def __init__(self, spec: RepresentationSpec):
        self.spec = spec
        self.uid_to_row: dict[str, int] = {}
        self.matrix = None
        self.feature_names: list[str] = []

        if spec.status != "READY":
            raise RuntimeError(
                f"Representation {spec.representation_id} status={spec.status}; "
                "only READY representations can train."
            )

        if spec.format == "dense_parquet":
            self._load_dense()
        elif spec.format == "sparse_csr":
            self._load_sparse()
        else:
            raise ValueError(f"Unknown representation format: {spec.format}")

    def _load_dense(self):
        path = self.spec.main_path
        if not path.is_file():
            raise FileNotFoundError(path)

        df = pd.read_parquet(path)
        if "sequence_uid" not in df.columns:
            raise KeyError(f"{path}: missing sequence_uid")

        if df["sequence_uid"].astype(str).duplicated().any():
            raise RuntimeError(f"{path}: duplicate sequence_uid values")

        if self.spec.family == "plm":
            feature_cols = [
                c for c in df.columns
                if c.startswith("p_") and c[2:].isdigit()
            ]
            feature_cols = sorted(feature_cols, key=lambda x: int(x[2:]))
        elif self.spec.family == "roxy":
            metadata = {
                "sequence_uid",
                "length",
                "is_canonical_20aa",
                "noncanonical_symbols",
            }
            feature_cols = [
                c for c in df.columns
                if c not in metadata
                and pd.api.types.is_numeric_dtype(df[c])
            ]
        else:
            raise ValueError(self.spec.family)

        if not feature_cols:
            raise RuntimeError(f"{path}: no numeric feature columns detected")

        if (
            self.spec.expected_dimension is not None
            and len(feature_cols) != self.spec.expected_dimension
        ):
            raise RuntimeError(
                f"{self.spec.representation_id}: expected "
                f"{self.spec.expected_dimension} features, found {len(feature_cols)}"
            )

        uids = df["sequence_uid"].astype(str).tolist()
        self.uid_to_row = {uid: i for i, uid in enumerate(uids)}
        self.matrix = df[feature_cols].to_numpy(dtype=np.float32, copy=True)
        self.feature_names = feature_cols

        if not np.isfinite(self.matrix).all():
            raise RuntimeError(
                f"{self.spec.representation_id}: NaN/Inf values detected"
            )

    def _load_sparse(self):
        if self.spec.index_path is None:
            raise RuntimeError("Sparse representation requires index_path")
        if not self.spec.main_path.is_file():
            raise FileNotFoundError(self.spec.main_path)
        if not self.spec.index_path.is_file():
            raise FileNotFoundError(self.spec.index_path)

        matrix = sparse.load_npz(self.spec.main_path).tocsr()
        idx = read_tsv(self.spec.index_path)

        if "sequence_uid" not in idx.columns:
            raise KeyError(f"{self.spec.index_path}: missing sequence_uid")

        row_col = None
        for candidate in ("onehot_row", "row", "row_index"):
            if candidate in idx.columns:
                row_col = candidate
                break
        if row_col is None:
            raise KeyError(
                f"{self.spec.index_path}: no onehot row-index column found"
            )

        idx[row_col] = idx[row_col].astype(int)
        if len(idx) != matrix.shape[0]:
            raise RuntimeError(
                f"One-hot row index count {len(idx)} != matrix rows {matrix.shape[0]}"
            )
        if idx["sequence_uid"].duplicated().any():
            raise RuntimeError("One-hot index contains duplicate sequence_uid")

        self.uid_to_row = dict(
            zip(idx["sequence_uid"].astype(str), idx[row_col].astype(int))
        )
        self.matrix = matrix
        self.feature_names = [f"onehot_{i}" for i in range(matrix.shape[1])]

    @property
    def dimension(self) -> int:
        return int(self.matrix.shape[1])

    @property
    def sparse(self) -> bool:
        return sparse.issparse(self.matrix)

    def require_uids(self, uids: list[str]) -> None:
        missing = [uid for uid in uids if uid not in self.uid_to_row]
        if missing:
            preview = ", ".join(missing[:5])
            raise RuntimeError(
                f"{self.spec.representation_id}: {len(missing):,} eligible "
                f"sequence_uids are missing. First: {preview}. "
                "Representation incompleteness cannot redefine the benchmark."
            )

    def get(self, uids: list[str]):
        rows = np.fromiter(
            (self.uid_to_row[uid] for uid in uids),
            dtype=np.int64,
            count=len(uids),
        )
        return self.matrix[rows]


def load_representation_registry(
    registry_path: Path,
    root: Path,
) -> dict[str, RepresentationSpec]:
    if not registry_path.is_file():
        raise FileNotFoundError(registry_path)

    df = pd.read_csv(registry_path, sep="\t", dtype=str).fillna("")
    required = {
        "representation_id",
        "family",
        "format",
        "main_path",
        "status",
    }
    missing = required - set(df.columns)
    if missing:
        raise KeyError(f"Registry missing columns: {sorted(missing)}")

    specs = {}
    for row in df.to_dict("records"):
        dim = None
        if row.get("expected_dimension", ""):
            dim = int(row["expected_dimension"])
        main = Path(row["main_path"])
        if not main.is_absolute():
            main = root / main
        idx = None
        if row.get("index_path", ""):
            idx = Path(row["index_path"])
            if not idx.is_absolute():
                idx = root / idx
        spec = RepresentationSpec(
            representation_id=row["representation_id"],
            family=row["family"],
            format=row["format"],
            main_path=main.resolve(),
            index_path=idx.resolve() if idx else None,
            expected_dimension=dim,
            status=row["status"],
        )
        specs[spec.representation_id] = spec
    return specs


def supervised_target(df: pd.DataFrame, task_type: str):
    if task_type == "binary":
        if "label" not in df.columns:
            raise KeyError("Binary task requires label column")
        raw = df["label"].astype(str)
        values = set(raw.unique())
        if values <= {"0", "1"}:
            return raw.astype(int).to_numpy(), {0: "negative", 1: "positive"}
        # Conservative support for text labels.
        mapping = {}
        for value in sorted(values):
            low = value.strip().lower()
            if low in {"positive", "pos", "p", "true"}:
                mapping[value] = 1
            elif low in {"negative", "neg", "n", "false"}:
                mapping[value] = 0
        if set(mapping.values()) == {0, 1} and len(mapping) == len(values):
            return raw.map(mapping).astype(int).to_numpy(), {
                0: "negative",
                1: "positive",
            }
        raise RuntimeError(f"Unsupported binary labels: {sorted(values)}")

    if task_type == "multiclass":
        candidate = None
        for c in (
            "class_name",
            "canonical_target_class",
            "label",
        ):
            if c in df.columns and df[c].nunique() >= 2:
                candidate = c
                break
        if candidate is None:
            raise RuntimeError("No multiclass target column found")
        names = df[candidate].astype(str)
        classes = sorted(names.unique())
        mapping = {name: i for i, name in enumerate(classes)}
        inverse = {i: name for name, i in mapping.items()}
        return names.map(mapping).astype(int).to_numpy(), inverse

    raise ValueError(task_type)


def pu_target(df: pd.DataFrame):
    candidate = None
    for c in ("pnu_role", "training_role", "label", "class_name"):
        if c in df.columns and df[c].nunique() >= 2:
            values = {str(x).strip().lower() for x in df[c].unique()}
            if (
                {"p", "u"} <= values
                or {"positive", "unlabeled"} <= values
                or values <= {"0", "1"}
            ):
                candidate = c
                break
    if candidate is None:
        raise RuntimeError("No P/U role column found")

    raw = df[candidate].astype(str)
    values = sorted(raw.unique())
    mapping = {}
    for value in values:
        low = value.strip().lower()
        if low in {"p", "positive", "1", "labeled_positive"}:
            mapping[value] = 1
        elif low in {"u", "unlabeled", "0", "unknown"}:
            mapping[value] = 0
    if set(mapping.values()) != {0, 1}:
        raise RuntimeError(
            f"Cannot map PU roles from {candidate}: {values}"
        )
    return raw.map(mapping).astype(int).to_numpy(), candidate


def get_config(algorithm: str, config_id: str) -> dict[str, Any]:
    for cfg in PARAMETER_PANELS[algorithm]:
        if cfg["config_id"] == config_id:
            return cfg
    raise KeyError((algorithm, config_id))


def build_model(
    algorithm: str,
    config: dict[str, Any],
    *,
    task_type: str,
    representation_family: str,
    n_classes: int,
):
    cfg = dict(config)
    cfg.pop("config_id", None)

    if algorithm == "LR":
        estimator = LogisticRegression(
            C=cfg["C"],
            l1_ratio=0.0,
            solver="lbfgs",
            max_iter=2000,
            class_weight="balanced",
            random_state=SEED,
        )
    elif algorithm == "RF":
        estimator = RandomForestClassifier(
            n_estimators=cfg["n_estimators"],
            max_depth=cfg["max_depth"],
            min_samples_leaf=cfg["min_samples_leaf"],
            max_features=cfg["max_features"],
            class_weight="balanced_subsample",
            random_state=SEED,
            n_jobs=n_jobs(),
        )
    elif algorithm == "KNN":
        estimator = KNeighborsClassifier(
            n_neighbors=cfg["n_neighbors"],
            weights="distance",
            metric="cosine",
            algorithm="brute",
            n_jobs=n_jobs(),
        )
    elif algorithm == "XGB":
        if XGBClassifier is None:
            raise RuntimeError("xgboost is not installed")
        common = dict(
            n_estimators=cfg["n_estimators"],
            max_depth=cfg["max_depth"],
            learning_rate=cfg["learning_rate"],
            subsample=cfg["subsample"],
            colsample_bytree=cfg["colsample_bytree"],
            tree_method="hist",
            random_state=SEED,
            n_jobs=n_jobs(),
            verbosity=0,
        )
        if task_type == "multiclass":
            estimator = XGBClassifier(
                objective="multi:softprob",
                num_class=n_classes,
                eval_metric="mlogloss",
                **common,
            )
        else:
            estimator = XGBClassifier(
                objective="binary:logistic",
                eval_metric="logloss",
                **common,
            )
    elif algorithm == "MLP":
        estimator = MLPClassifier(
            hidden_layer_sizes=cfg["hidden_layer_sizes"],
            alpha=cfg["alpha"],
            learning_rate_init=cfg["learning_rate_init"],
            activation="relu",
            solver="adam",
            batch_size=256,
            max_iter=300,
            early_stopping=False,
            tol=1e-4,
            n_iter_no_change=20,
            random_state=SEED,
        )
    else:
        raise ValueError(algorithm)

    use_scaler = (
        representation_family in {"roxy", "plm"}
        and algorithm in {"LR", "KNN", "MLP"}
    )

    if use_scaler:
        return Pipeline(
            [
                ("scaler", StandardScaler()),
                ("model", estimator),
            ]
        )
    return estimator


def fit_model(model, X, y, algorithm: str):
    if algorithm == "XGB":
        # Class balancing is derived strictly from the current fitting partition.
        from sklearn.utils.class_weight import compute_sample_weight

        sw = compute_sample_weight(class_weight="balanced", y=y)
        model.fit(X, y, sample_weight=sw)
    else:
        model.fit(X, y)
    return model


def positive_scores(model, X) -> np.ndarray:
    prob = model.predict_proba(X)
    classes = np.asarray(model.classes_)
    idx = np.where(classes == 1)[0]
    if len(idx) != 1:
        raise RuntimeError(f"Cannot find positive class in model.classes_={classes}")
    return np.asarray(prob[:, idx[0]], dtype=float)


def multiclass_scores(model, X) -> np.ndarray:
    return np.asarray(model.predict_proba(X), dtype=float)


def binary_metrics(y_true, scores, threshold: float = 0.5) -> dict[str, Any]:
    y_true = np.asarray(y_true, dtype=int)
    scores = np.asarray(scores, dtype=float)
    y_pred = (scores >= threshold).astype(int)

    tn, fp, fn, tp = confusion_matrix(
        y_true, y_pred, labels=[0, 1]
    ).ravel()

    prevalence = float(np.mean(y_true == 1))
    out = {
        "mcc": float(matthews_corrcoef(y_true, y_pred)),
        "ap": float(average_precision_score(y_true, scores)),
        "ap_baseline": prevalence,
        "ap_gain": float(average_precision_score(y_true, scores) - prevalence),
        "roc_auc": float(roc_auc_score(y_true, scores))
        if len(np.unique(y_true)) == 2
        else np.nan,
        "f1": float(f1_score(y_true, y_pred, zero_division=0)),
        "precision": float(
            precision_score(y_true, y_pred, zero_division=0)
        ),
        "recall": float(recall_score(y_true, y_pred, zero_division=0)),
        "specificity": float(tn / (tn + fp)) if (tn + fp) else np.nan,
        "balanced_accuracy": float(
            balanced_accuracy_score(y_true, y_pred)
        ),
        "brier": float(brier_score_loss(y_true, scores)),
        "log_loss": float(
            log_loss(
                y_true,
                np.column_stack([1.0 - scores, scores]),
                labels=[0, 1],
            )
        ),
        "threshold": float(threshold),
        "prevalence": prevalence,
        "predicted_positive_fraction": float(np.mean(y_pred == 1)),
        "tp": int(tp),
        "tn": int(tn),
        "fp": int(fp),
        "fn": int(fn),
        "n": int(len(y_true)),
        "n_positive": int(np.sum(y_true == 1)),
        "n_negative": int(np.sum(y_true == 0)),
    }
    return out


def multiclass_metrics(
    y_true,
    prob,
    class_names: dict[int, str],
) -> tuple[dict[str, Any], pd.DataFrame, np.ndarray]:
    y_true = np.asarray(y_true, dtype=int)
    prob = np.asarray(prob, dtype=float)
    labels = sorted(class_names)
    y_pred = np.argmax(prob, axis=1).astype(int)

    pr, rc, f1, support = precision_recall_fscore_support(
        y_true,
        y_pred,
        labels=labels,
        zero_division=0,
    )
    per_class = pd.DataFrame(
        {
            "class_id": labels,
            "class_name": [class_names[i] for i in labels],
            "precision": pr,
            "recall": rc,
            "f1": f1,
            "support": support,
        }
    )

    summary = {
        "mcc": float(matthews_corrcoef(y_true, y_pred)),
        "macro_f1": float(
            f1_score(y_true, y_pred, average="macro", zero_division=0)
        ),
        "balanced_accuracy": float(
            balanced_accuracy_score(y_true, y_pred)
        ),
        "macro_precision": float(
            precision_score(
                y_true, y_pred, average="macro", zero_division=0
            )
        ),
        "macro_recall": float(
            recall_score(
                y_true, y_pred, average="macro", zero_division=0
            )
        ),
        "weighted_f1": float(
            f1_score(
                y_true, y_pred, average="weighted", zero_division=0
            )
        ),
        "log_loss": float(log_loss(y_true, prob, labels=labels)),
        "n": int(len(y_true)),
    }
    cm = confusion_matrix(y_true, y_pred, labels=labels)
    return summary, per_class, cm


def choose_best_config(
    metrics_df: pd.DataFrame,
    task_type: str,
) -> str:
    if task_type == "binary":
        ranked = metrics_df.sort_values(
            ["mcc", "ap", "config_order"],
            ascending=[False, False, True],
            kind="stable",
        )
    else:
        ranked = metrics_df.sort_values(
            ["mcc", "macro_f1", "config_order"],
            ascending=[False, False, True],
            kind="stable",
        )
    return str(ranked.iloc[0]["config_id"])


def make_binary_prediction_table(
    metadata: pd.DataFrame,
    y_true,
    scores,
    *,
    experiment_id: str,
    config_id: str | None = None,
) -> pd.DataFrame:
    out = metadata[
        ["entity_id", "sequence_uid", "cluster_id"]
    ].reset_index(drop=True).copy()
    out.insert(0, "experiment_id", experiment_id)
    if config_id is not None:
        out.insert(1, "config_id", config_id)
    out["y_true"] = np.asarray(y_true, dtype=int)
    out["score_positive"] = np.asarray(scores, dtype=float)
    out["y_pred"] = (out["score_positive"] >= 0.5).astype(int)
    return out


def make_multiclass_prediction_table(
    metadata: pd.DataFrame,
    y_true,
    prob,
    class_names: dict[int, str],
    *,
    experiment_id: str,
    config_id: str | None = None,
) -> pd.DataFrame:
    out = metadata[
        ["entity_id", "sequence_uid", "cluster_id"]
    ].reset_index(drop=True).copy()
    out.insert(0, "experiment_id", experiment_id)
    if config_id is not None:
        out.insert(1, "config_id", config_id)
    y_true = np.asarray(y_true, dtype=int)
    prob = np.asarray(prob, dtype=float)
    y_pred = np.argmax(prob, axis=1)
    out["y_true"] = y_true
    out["y_true_name"] = [class_names[int(x)] for x in y_true]
    out["y_pred"] = y_pred
    out["y_pred_name"] = [class_names[int(x)] for x in y_pred]
    for idx in sorted(class_names):
        slug = (
            class_names[idx]
            .lower()
            .replace(" ", "_")
            .replace("-", "_")
            .replace("/", "_")
        )
        out[f"p_{slug}"] = prob[:, idx]
    return out


def pu_diagnostics(y_observed, scores) -> dict[str, Any]:
    y_observed = np.asarray(y_observed, dtype=int)
    scores = np.asarray(scores, dtype=float)
    p = scores[y_observed == 1]
    u = scores[y_observed == 0]

    if len(p) == 0 or len(u) == 0:
        raise RuntimeError("PU diagnostics require both P and U")

    q90, q95, q99 = np.quantile(u, [0.90, 0.95, 0.99])
    return {
        "n": int(len(scores)),
        "n_p": int(len(p)),
        "n_u": int(len(u)),
        "median_p_score": float(np.median(p)),
        "median_u_score": float(np.median(u)),
        "u_q90": float(q90),
        "u_q95": float(q95),
        "u_q99": float(q99),
        "p_recall_above_u_q90": float(np.mean(p >= q90)),
        "p_recall_above_u_q95": float(np.mean(p >= q95)),
        "p_recall_above_u_q99": float(np.mean(p >= q99)),
        "fraction_u_ge_0_5": float(np.mean(u >= 0.5)),
        "fraction_p_ge_0_5": float(np.mean(p >= 0.5)),
    }

#!/usr/bin/env python3
"""
PhageLysClass inference script.

Single-file inference entry point for the frozen PhageLysClass v1 production bundle.

Pipeline
--------
FASTA -> Sylphy ESM2-T30 mean-pooled embedding -> global lytic score
      -> gated multiclass functional assignment
      -> optional subtype compatibility scores
      -> TSV/CSV predictions + JSON run metadata

The same PhageLysClassPredictor class can later be imported by a web API.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from collections import Counter
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import joblib
import numpy as np
import pandas as pd


CANONICAL_AA = frozenset("ACDEFGHIKLMNPQRSTVWY")
DEFAULT_MAX_LENGTH = 1024
DEFAULT_THRESHOLD = 0.5

CORE_GLOBAL_TASK = "global_mixed"
CORE_MULTICLASS_TASK = "multiclass_clean"
AUX_TASKS = [
    "subtype_endolysin",
    "subtype_depolymerase",
    "subtype_val",
]

BIOLOGICAL_CLASS_ORDER = [
    "endolysin",
    "depolymerase",
    "virion_associated_lytic_enzyme",
]

BIOLOGICAL_DISPLAY = {
    "endolysin": "endolysin",
    "depolymerase": "depolymerase",
    "virion_associated_lytic_enzyme": "VAL",
}


@dataclass
class SequenceRecord:
    sequence_id: str
    sequence: str


def parse_fasta(path: Path) -> list[SequenceRecord]:
    """Minimal FASTA parser with duplicate-ID protection."""
    records: list[SequenceRecord] = []
    current_id: str | None = None
    chunks: list[str] = []

    with path.open("r", encoding="utf-8") as handle:
        for line_no, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line:
                continue

            if line.startswith(">"):
                if current_id is not None:
                    records.append(
                        SequenceRecord(
                            sequence_id=current_id,
                            sequence="".join(chunks).replace(" ", "").upper(),
                        )
                    )
                header = line[1:].strip()
                if not header:
                    raise ValueError(f"FASTA line {line_no}: empty header.")
                current_id = header.split()[0]
                chunks = []
            else:
                if current_id is None:
                    raise ValueError(
                        f"FASTA line {line_no}: sequence found before first header."
                    )
                chunks.append(line)

    if current_id is not None:
        records.append(
            SequenceRecord(
                sequence_id=current_id,
                sequence="".join(chunks).replace(" ", "").upper(),
            )
        )

    if not records:
        raise ValueError(f"No FASTA records found in {path}")

    ids = [r.sequence_id for r in records]
    id_counts = Counter(ids)
    duplicates = sorted(x for x, count in id_counts.items() if count > 1)
    if duplicates:
        raise ValueError(
            "Duplicate FASTA identifiers are not allowed. "
            f"Examples: {duplicates[:10]}"
        )

    return records


def validate_sequence(sequence: str, max_length: int) -> list[str]:
    warnings: list[str] = []

    if not sequence:
        warnings.append("empty_sequence")
        return warnings

    invalid = sorted(set(sequence) - CANONICAL_AA)
    if invalid:
        warnings.append("unsupported_residue:" + ",".join(invalid))

    if len(sequence) > max_length:
        warnings.append(f"sequence_too_long:{len(sequence)}>{max_length}")

    return warnings


def _extract_estimator_classes(model) -> np.ndarray:
    if not hasattr(model, "classes_"):
        raise RuntimeError(
            f"Serialized estimator {type(model).__name__} has no classes_ attribute."
        )
    return np.asarray(model.classes_)


def _binary_positive_probability(model, X: np.ndarray) -> np.ndarray:
    proba = np.asarray(model.predict_proba(X), dtype=float)
    classes = _extract_estimator_classes(model)

    matches = np.where(classes.astype(str) == "1")[0]
    if len(matches) != 1:
        try:
            matches = np.where(classes.astype(int) == 1)[0]
        except Exception:
            matches = np.array([], dtype=int)

    if len(matches) != 1:
        raise RuntimeError(
            f"Could not identify positive class=1 in estimator classes_: {classes.tolist()}"
        )

    return proba[:, int(matches[0])]


def _multiclass_probabilities(
    model,
    X: np.ndarray,
    encoded_class_order: list[int],
) -> np.ndarray:
    proba = np.asarray(model.predict_proba(X), dtype=float)
    classes = _extract_estimator_classes(model)

    positions = []
    for encoded in encoded_class_order:
        matches = np.where(classes.astype(str) == str(encoded))[0]
        if len(matches) != 1:
            try:
                matches = np.where(classes.astype(int) == int(encoded))[0]
            except Exception:
                matches = np.array([], dtype=int)

        if len(matches) != 1:
            raise RuntimeError(
                f"Could not align encoded class {encoded} in classes_={classes.tolist()}"
            )
        positions.append(int(matches[0]))

    ordered = proba[:, positions]
    row_sums = ordered.sum(axis=1, keepdims=True)
    if (row_sums <= 0).any():
        raise RuntimeError("Multiclass estimator returned non-positive probability sum.")
    return ordered / row_sums


class PhageLysClassPredictor:
    """Reusable predictor; the CLI below is only a thin wrapper."""

    def __init__(
        self,
        model_root: Path,
        device: str = "auto",
        batch_size: int = 16,
    ) -> None:
        self.model_root = model_root.resolve()
        self.batch_size = int(batch_size)

        if self.batch_size < 1:
            raise ValueError("batch_size must be >= 1")

        manifest_path = self.model_root / "manifest.json"
        if not manifest_path.is_file():
            raise FileNotFoundError(f"Missing model manifest: {manifest_path}")

        self.manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        self.model_version = str(self.manifest.get("model_version", self.model_root.name))

        if self.manifest.get("architecture_id") != "C_2_core_plus_3_aux":
            raise RuntimeError(
                "This script expects the frozen '2 core + 3 auxiliary' architecture."
            )

        self.representation_metadata = self._load_representation_metadata()
        self.max_length = int(
            self.representation_metadata
            .get("eligibility", {})
            .get("max_sequence_length", DEFAULT_MAX_LENGTH)
        )

        self.checkpoint = str(
            self.representation_metadata.get(
                "checkpoint",
                "facebook/esm2_t30_150m_ur50d",
            )
        )

        self.device = self._resolve_device(device)
        self.tasks = self._load_tasks()

        global_ensemble = self.tasks[CORE_GLOBAL_TASK]["ensemble"]
        self.global_threshold = float(
            global_ensemble.get("threshold", DEFAULT_THRESHOLD)
            if global_ensemble is not None
            else DEFAULT_THRESHOLD
        )


    def _load_representation_metadata(self) -> dict:
        rep_dir = self.model_root / "representation"
        json_files = sorted(rep_dir.glob("*.json"))
        if len(json_files) != 1:
            raise RuntimeError(
                f"Expected exactly one representation metadata JSON in {rep_dir}; "
                f"found {len(json_files)}."
            )
        metadata = json.loads(json_files[0].read_text(encoding="utf-8"))
        if int(metadata.get("dimension", -1)) != 640:
            raise RuntimeError(
                f"Expected 640-dimensional ESM2-T30 representation; "
                f"observed {metadata.get('dimension')}."
            )
        return metadata

    @staticmethod
    def _resolve_device(device: str) -> str:
        import torch

        device = device.lower()
        if device == "auto":
            return "cuda" if torch.cuda.is_available() else "cpu"
        if device == "cuda" and not torch.cuda.is_available():
            raise RuntimeError("--device cuda requested but CUDA is not available.")
        if device not in {"cpu", "cuda"}:
            raise ValueError("--device must be one of: auto, cpu, cuda")
        return device

    def _load_tasks(self) -> dict:
        task_payload = {}

        for task in [CORE_GLOBAL_TASK, CORE_MULTICLASS_TASK, *AUX_TASKS]:
            task_dir = self.model_root / "tasks" / task
            metadata_path = task_dir / "metadata.json"

            if not metadata_path.is_file():
                raise FileNotFoundError(metadata_path)

            metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
            member_records = metadata.get("members", [])
            if not member_records:
                raise RuntimeError(f"{task}: no serialized member models found in metadata.")

            loaded_members = []
            for member in member_records:
                model_path = task_dir / member["model_file"]
                if not model_path.is_file():
                    raise FileNotFoundError(model_path)
                loaded_members.append(
                    {
                        "algorithm": member["algorithm"],
                        "strategy_id": member["strategy_id"],
                        "model": joblib.load(model_path),
                        "model_path": model_path,
                    }
                )

            ensemble_path = task_dir / "ensemble.json"
            ensemble = (
                json.loads(ensemble_path.read_text(encoding="utf-8"))
                if ensemble_path.is_file()
                else None
            )

            task_payload[task] = {
                "metadata": metadata,
                "ensemble": ensemble,
                "members": loaded_members,
            }

        return task_payload

    @staticmethod
    def _canonical_sylphy_checkpoint(checkpoint: str) -> str:
        """
        Normalize the frozen ESM2-T30 checkpoint to the canonical Hugging Face ID
        used by Sylphy.
        """
        low = str(checkpoint).strip().lower()
        if low == "facebook/esm2_t30_150m_ur50d":
            return "facebook/esm2_t30_150M_UR50D"
        return str(checkpoint).strip()

    @staticmethod
    def _extract_sylphy_matrix(coded_dataset, n_expected: int) -> np.ndarray:
        """
        Robustly recover the sequence-level embedding matrix from Sylphy output.

        Supports:
        - 640 scalar embedding columns;
        - one vector/list column containing 640 values per row.

        Sequence/id metadata columns are ignored automatically.
        """
        columns = list(getattr(coded_dataset, "columns", []))
        if not columns:
            raise RuntimeError("Sylphy coded_dataset exposes no columns.")

        # First try common scalar-column schemas.
        preferred_schemas = [
            [f"p_{i}" for i in range(640)],
            [f"embedding_{i}" for i in range(640)],
            [f"emb_{i}" for i in range(640)],
        ]
        for cols in preferred_schemas:
            if all(c in columns for c in cols):
                arr = coded_dataset.select(cols).to_numpy()
                arr = np.asarray(arr, dtype=np.float32)
                if arr.shape == (n_expected, 640):
                    return arr

        # Next detect exactly 640 scalar numeric columns by conversion.
        numeric_cols = []
        for col in columns:
            if str(col).lower() in {
                "sequence", "sequence_id", "id", "entity_id",
                "name", "header", "label",
            }:
                continue
            try:
                values = np.asarray(coded_dataset[col].to_numpy())
                if values.ndim != 1 or len(values) != n_expected:
                    continue
                np.asarray(values, dtype=np.float32)
                numeric_cols.append(col)
            except Exception:
                continue

        if len(numeric_cols) == 640:
            arr = coded_dataset.select(numeric_cols).to_numpy()
            arr = np.asarray(arr, dtype=np.float32)
            if arr.shape == (n_expected, 640):
                return arr

        # Finally, support a single vector/list embedding column.
        for col in columns:
            try:
                values = coded_dataset[col].to_list()
            except Exception:
                continue
            if len(values) != n_expected or not values:
                continue
            first = values[0]
            if isinstance(first, (list, tuple, np.ndarray)) and len(first) == 640:
                arr = np.vstack(values).astype(np.float32, copy=False)
                if arr.shape == (n_expected, 640):
                    return arr

        raise RuntimeError(
            "Could not identify a 640-dimensional Sylphy embedding matrix. "
            f"Observed columns ({len(columns)}): {columns[:30]}"
        )

    def embed_sequences(self, sequences: list[str]) -> tuple[np.ndarray, float]:
        """
        Generate the frozen production representation through Sylphy.

        Contract:
        - model: facebook/esm2_t30_150M_UR50D
        - final hidden representation
        - sequence mean pooling
        - FP32
        - 640 dimensions

        This intentionally uses the same representation backend as the benchmark
        instead of reimplementing tokenization/pooling with Transformers.
        """
        try:
            import polars as pl
            from sylphy.embedding_extractor import create_embedding
        except ImportError as exc:
            raise RuntimeError(
                "PhageLysClass inference requires Sylphy with embedding support. "
                "Install it with: pip install 'sylphy[embeddings]'"
            ) from exc

        checkpoint = self._canonical_sylphy_checkpoint(self.checkpoint)
        dataset = pl.DataFrame({"sequence": sequences})

        start = time.perf_counter()

        embedder = create_embedding(
            model_name=checkpoint,
            dataset=dataset,
            column_seq="sequence",
            name_device=self.device,
            precision="fp32",
        )
        embedder.run_process(
            batch_size=self.batch_size,
            pool="mean",
        )

        coded = embedder.coded_dataset
        X = self._extract_sylphy_matrix(coded, n_expected=len(sequences))

        if X.shape != (len(sequences), 640):
            raise RuntimeError(
                f"Sylphy returned embedding shape {X.shape}; "
                f"expected ({len(sequences)}, 640)."
            )
        if not np.isfinite(X).all():
            raise RuntimeError("Sylphy produced non-finite embedding values.")

        elapsed = time.perf_counter() - start
        return X.astype(np.float32, copy=False), elapsed

    def _predict_binary_task(self, task: str, X: np.ndarray) -> np.ndarray:
        payload = self.tasks[task]
        member_scores = np.column_stack(
            [
                _binary_positive_probability(member["model"], X)
                for member in payload["members"]
            ]
        )

        ensemble = payload["ensemble"]
        if ensemble is None:
            if member_scores.shape[1] != 1:
                raise RuntimeError(f"{task}: missing ensemble recipe for multi-member task.")
            return member_scores[:, 0]

        if ensemble.get("type") != "unweighted_soft_voting":
            raise RuntimeError(
                f"{task}: unsupported ensemble type {ensemble.get('type')!r}"
            )

        return member_scores.mean(axis=1)

    def _predict_multiclass_task(self, task: str, X: np.ndarray) -> np.ndarray:
        payload = self.tasks[task]
        ensemble = payload["ensemble"]

        if ensemble is None:
            encoded_order = [0, 1, 2]
        else:
            if ensemble.get("type") != "unweighted_soft_voting":
                raise RuntimeError(
                    f"{task}: unsupported ensemble type {ensemble.get('type')!r}"
                )
            encoded_order = [int(x) for x in ensemble.get("encoded_class_order", [0, 1, 2])]

        member_probs = np.stack(
            [
                _multiclass_probabilities(member["model"], X, encoded_order)
                for member in payload["members"]
            ],
            axis=0,
        )

        if ensemble is None:
            if member_probs.shape[0] != 1:
                raise RuntimeError(f"{task}: missing ensemble recipe for multi-member task.")
            proba = member_probs[0]
        else:
            proba = member_probs.mean(axis=0)

        proba = proba / proba.sum(axis=1, keepdims=True)
        return proba

    def predict_valid_embeddings(
        self,
        sequence_ids: list[str],
        sequences: list[str],
        X: np.ndarray,
        auxiliary: bool,
    ) -> pd.DataFrame:
        if len(sequence_ids) != len(sequences) or len(sequence_ids) != len(X):
            raise ValueError("sequence_ids, sequences and embeddings must have the same length.")

        global_score = self._predict_binary_task(CORE_GLOBAL_TASK, X)
        global_positive = global_score >= self.global_threshold

        multiclass = self._predict_multiclass_task(CORE_MULTICLASS_TASK, X)
        multiclass_idx = np.argmax(multiclass, axis=1)

        rows = []
        for i, sequence_id in enumerate(sequence_ids):
            global_prediction = "lytic" if global_positive[i] else "non_lytic_background"

            if global_positive[i]:
                functional_prediction = BIOLOGICAL_DISPLAY[
                    BIOLOGICAL_CLASS_ORDER[int(multiclass_idx[i])]
                ]
                functional_status = "reported"
            else:
                functional_prediction = ""
                functional_status = "not_reported_global_non_lytic"

            warnings_list: list[str] = []
            if (
                global_positive[i]
                and BIOLOGICAL_CLASS_ORDER[int(multiclass_idx[i])]
                == "virion_associated_lytic_enzyme"
            ):
                warnings_list.append("VAL_evidence_limited")

            # Multiclass probabilities are conditional on the sequence being
            # evaluated inside the lytic-class space. They are therefore
            # reported only when the global core supports a lytic prediction.
            if global_positive[i]:
                p_endolysin_given_lytic = float(multiclass[i, 0])
                p_depolymerase_given_lytic = float(multiclass[i, 1])
                p_val_given_lytic = float(multiclass[i, 2])
            else:
                p_endolysin_given_lytic = np.nan
                p_depolymerase_given_lytic = np.nan
                p_val_given_lytic = np.nan

            rows.append(
                {
                    "sequence_id": sequence_id,
                    "sequence_length": len(sequences[i]),
                    "input_status": "valid",
                    "global_score": float(global_score[i]),
                    "global_prediction": global_prediction,
                    "p_endolysin_given_lytic": p_endolysin_given_lytic,
                    "p_depolymerase_given_lytic": p_depolymerase_given_lytic,
                    "p_val_given_lytic": p_val_given_lytic,
                    "functional_prediction": functional_prediction,
                    "functional_prediction_status": functional_status,
                    "warnings": ";".join(warnings_list),
                }
            )

        result = pd.DataFrame(rows)

        if auxiliary:
            aux_scores = {
                "endolysin_compatibility_score": self._predict_binary_task(
                    "subtype_endolysin", X
                ),
                "depolymerase_compatibility_score": self._predict_binary_task(
                    "subtype_depolymerase", X
                ),
                "val_compatibility_score": self._predict_binary_task(
                    "subtype_val", X
                ),
            }
            for column, values in aux_scores.items():
                result[column] = values.astype(float)

        return result

    def predict_fasta(
        self,
        fasta_path: Path,
        auxiliary: bool = False,
    ) -> tuple[pd.DataFrame, dict]:
        records = parse_fasta(fasta_path)

        valid_ids: list[str] = []
        valid_sequences: list[str] = []
        invalid_rows: list[dict] = []

        for record in records:
            seq_warnings = validate_sequence(record.sequence, self.max_length)

            hard_invalid = any(
                warning == "empty_sequence"
                or warning.startswith("unsupported_residue:")
                or warning.startswith("sequence_too_long:")
                for warning in seq_warnings
            )

            if hard_invalid:
                invalid_rows.append(
                    {
                        "sequence_id": record.sequence_id,
                        "sequence_length": len(record.sequence),
                        "input_status": "invalid",
                        "global_score": np.nan,
                        "global_prediction": "",
                        "p_endolysin_given_lytic": np.nan,
                        "p_depolymerase_given_lytic": np.nan,
                        "p_val_given_lytic": np.nan,
                        "functional_prediction": "",
                        "functional_prediction_status": "not_reported_invalid_input",
                        "warnings": ";".join(seq_warnings),
                    }
                )
            else:
                valid_ids.append(record.sequence_id)
                valid_sequences.append(record.sequence)

        total_start = time.perf_counter()

        embedding_seconds = 0.0
        inference_seconds = 0.0

        if valid_sequences:
            X, embedding_seconds = self.embed_sequences(valid_sequences)

            inference_start = time.perf_counter()
            valid_predictions = self.predict_valid_embeddings(
                valid_ids,
                valid_sequences,
                X,
                auxiliary=auxiliary,
            )
            inference_seconds = time.perf_counter() - inference_start
        else:
            valid_predictions = pd.DataFrame()

        if invalid_rows:
            invalid_df = pd.DataFrame(invalid_rows)
            if auxiliary:
                for col in [
                    "endolysin_compatibility_score",
                    "depolymerase_compatibility_score",
                    "val_compatibility_score",
                ]:
                    invalid_df[col] = np.nan
        else:
            invalid_df = pd.DataFrame()

        # Restore original FASTA order.
        parts = [df for df in [valid_predictions, invalid_df] if len(df)]
        predictions = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()

        order = {record.sequence_id: i for i, record in enumerate(records)}
        predictions["_input_order"] = predictions["sequence_id"].map(order)
        predictions = (
            predictions
            .sort_values("_input_order")
            .drop(columns="_input_order")
            .reset_index(drop=True)
        )

        predictions["model_version"] = self.model_version
        predictions["representation"] = self.representation_metadata.get(
            "representation_id",
            "esm2_t30_150m",
        )

        total_seconds = time.perf_counter() - total_start

        run_metadata = {
            "schema_version": "phagelysclass-inference-run-v1",
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "input_fasta": str(fasta_path.resolve()),
            "model_root": str(self.model_root),
            "model_version": self.model_version,
            "architecture_id": self.manifest.get("architecture_id"),
            "representation": self.representation_metadata.get("representation_id"),
            "embedding_backend": "sylphy",
            "checkpoint": self._canonical_sylphy_checkpoint(self.checkpoint),
            "embedding_pool": "mean",
            "embedding_precision": "fp32",
            "embedding_dimension": 640,
            "device": self.device,
            "batch_size": self.batch_size,
            "auxiliary": bool(auxiliary),
            "global_threshold": float(self.global_threshold),
            "n_sequences": len(records),
            "n_valid": len(valid_sequences),
            "n_invalid": len(invalid_rows),
            "runtime_seconds": {
                "embedding": float(embedding_seconds),
                "model_inference": float(inference_seconds),
                "total": float(total_seconds),
            },
        }

        return predictions, run_metadata


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Run the frozen PhageLysClass inference pipeline on protein sequences in FASTA format."
        )
    )
    parser.add_argument(
        "--input",
        required=True,
        type=Path,
        help="Input protein FASTA file.",
    )
    parser.add_argument(
        "--models",
        required=True,
        type=Path,
        help="Path to the exported production bundle, e.g. models/phagelysclass_v1.",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Output prediction table (.tsv or .csv).",
    )
    parser.add_argument(
        "--metadata-output",
        type=Path,
        default=None,
        help=(
            "Optional JSON run metadata path. "
            "Default: <output>.metadata.json"
        ),
    )
    parser.add_argument(
        "--auxiliary",
        action="store_true",
        help="Also emit Endolysin/DPO/VAL binary compatibility scores.",
    )
    parser.add_argument(
        "--device",
        choices=["auto", "cpu", "cuda"],
        default="auto",
        help="Embedding device. Default: auto.",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=16,
        help="ESM2 embedding batch size. Default: 16.",
    )
    return parser


def main() -> int:
    parser = build_arg_parser()
    args = parser.parse_args()

    if not args.input.is_file():
        parser.error(f"Input FASTA does not exist: {args.input}")
    if not args.models.is_dir():
        parser.error(f"Model bundle does not exist: {args.models}")

    if args.output.suffix.lower() not in {".tsv", ".csv"}:
        parser.error("--output must end in .tsv or .csv")

    metadata_output = args.metadata_output
    if metadata_output is None:
        metadata_output = Path(str(args.output) + ".metadata.json")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    metadata_output.parent.mkdir(parents=True, exist_ok=True)

    predictor = PhageLysClassPredictor(
        model_root=args.models,
        device=args.device,
        batch_size=args.batch_size,
    )

    predictions, metadata = predictor.predict_fasta(
        fasta_path=args.input,
        auxiliary=args.auxiliary,
    )

    sep = "\t" if args.output.suffix.lower() == ".tsv" else ","
    predictions.to_csv(args.output, sep=sep, index=False)

    metadata_output.write_text(
        json.dumps(metadata, indent=2),
        encoding="utf-8",
    )

    print(f"Predictions: {args.output}")
    print(f"Metadata:    {metadata_output}")
    print(
        "Sequences:   "
        f"{metadata['n_sequences']} "
        f"(valid={metadata['n_valid']}, invalid={metadata['n_invalid']})"
    )
    print(
        "Runtime:     "
        f"{metadata['runtime_seconds']['total']:.3f} s "
        f"(embedding={metadata['runtime_seconds']['embedding']:.3f} s, "
        f"inference={metadata['runtime_seconds']['model_inference']:.3f} s)"
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""
10d_generate_plm_embeddings_sylphy.py

PhageLysClass Stage 10D — Frozen PLM representations with Sylphy.

Purpose
-------
Generate one pooled FP32 embedding matrix per frozen protein language model
for the exact-sequence universe materialized in Stage 10A.

Frozen eligibility contract
---------------------------
A sequence is PLM-eligible only when:
  1. length <= 1024 aa
  2. every residue belongs to the canonical 20-aa alphabet:
     ACDEFGHIKLMNPQRSTVWY

This is intended to be the SAME representation universe used for the frozen
one-hot baseline. Therefore Stage 10C should be run with:

    --max-length 1024 --canonical-only

Ineligible sequences are never deleted from the Stage-10A universe. They are
recorded explicitly in the eligibility table and excluded only from the
representation matrices.

Frozen representation contract
------------------------------
- 11 predefined PLMs inherited from the validated PhageLysData workflow
- Sylphy CLI
- max length: 1024 aa
- canonical 20-aa sequences only
- final/last hidden layer
- mean layer aggregation
- mean sequence pooling
- FP32
- CUDA by default
- model-level resumability
- one Parquet embedding file per model
- sequence_uid is the authoritative join key
- NO CSV embedding export

Sylphy itself writes directly to Parquet. A temporary raw Parquet is then
stream-normalized to a final file with:
    sequence_uid, length, p_0, p_1, ...

Default outputs
---------------
data/representations/sylphy_embeddings/
├── eligibility.tsv.gz
├── eligible_sequences.parquet
├── filtered_sequences.tsv.gz
├── model_registry.tsv
├── run_status.tsv
├── manifest.json
└── models/
    ├── facebook_esm2_t6_8m_ur50d/
    │   ├── embeddings_last_mean_fp32.parquet
    │   ├── metadata.json
    │   └── run.log
    └── ...

The script is read-only with respect to Stage 10A and all task/split datasets.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import logging
import math
import os
import re
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

try:
    import pyarrow as pa
    import pyarrow.compute as pc
    import pyarrow.parquet as pq
except Exception as exc:
    raise RuntimeError(
        "Stage 10D requires pyarrow for direct/streaming Parquet handling. "
        "Install with conda/pip before running."
    ) from exc


SCRIPT_NAME = "10d_generate_plm_embeddings_sylphy.py"
SCRIPT_VERSION = "10d.v1.0.0"
CONTRACT_VERSION = "phagelysclass_stage10d_plm_v1.0.0"

MAX_LENGTH = 1024
CANONICAL_AA_STRING = "ACDEFGHIKLMNPQRSTVWY"
CANONICAL_AA = frozenset(CANONICAL_AA_STRING)

DEVICE = "cuda"
PRECISION = "fp32"
LAYERS = "last"
LAYER_AGG = "mean"
POOL = "mean"

GPU_INDEX = 0
POST_MODEL_SLEEP = 5
GPU_CHECK_INTERVAL = 2
GPU_RELEASE_TIMEOUT = 90
GPU_MEMORY_TOLERANCE_MB = 512

# Frozen model set carried over from the previously validated PhageLysData
# embedding layer supplied by the user.
FROZEN_MODELS = [
    "elnaggarlab/ankh2-ext1",
    "esmc_300m",
    "facebook/esm2_t12_35m_ur50d",
    "facebook/esm2_t30_150m_ur50d",
    "rostlab/prot_t5_xl_uniref50",
]

# Conservative initial batches for sequences up to 1024 aa.
# Sylphy's CUDA OOM backoff remains enabled and may reduce these automatically.
MODEL_BATCH_SIZES = {
    "elnaggarlab/ankh2-ext1": 1,
    "esmc_300m": 2,
    "facebook/esm2_t12_35m_ur50d": 16,
    "facebook/esm2_t30_150m_ur50d": 4,
    "rostlab/prot_t5_xl_uniref50": 1,
}

DEFAULT_BATCH_SIZE = 1
POSTPROCESS_BATCH_ROWS = 4096

LOGGER = logging.getLogger("phagelysclass.stage10d")


# -----------------------------------------------------------------------------
# General helpers
# -----------------------------------------------------------------------------

def project_root_from_script(script_file: str) -> Path:
    return Path(script_file).resolve().parent.parents[1]


def now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def package_version(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except Exception:
        return "unknown"


def sha256_file(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def slugify(value: str) -> str:
    value = value.lower()
    value = re.sub(r"[/\s]+", "_", value)
    value = re.sub(r"[^a-z0-9._-]", "_", value)
    value = re.sub(r"_+", "_", value)
    return value.strip("_")


def configure_logging(level: str, log_file: Path) -> None:
    log_file.parent.mkdir(parents=True, exist_ok=True)
    LOGGER.setLevel(getattr(logging, level.upper()))
    LOGGER.handlers.clear()

    formatter = logging.Formatter(
        "%(asctime)s | %(levelname)s | %(message)s"
    )

    stream = logging.StreamHandler(sys.stdout)
    stream.setFormatter(formatter)
    LOGGER.addHandler(stream)

    file_handler = logging.FileHandler(
        log_file,
        mode="a",
        encoding="utf-8",
    )
    file_handler.setFormatter(formatter)
    LOGGER.addHandler(file_handler)


def write_json(data: dict[str, Any], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(
            data,
            indent=2,
            sort_keys=True,
            ensure_ascii=False,
        ) + "\n",
        encoding="utf-8",
    )


def normalize_sequence(sequence: str) -> str:
    return "".join(str(sequence).split()).upper()


def noncanonical_residues(sequence: str) -> str:
    return "".join(
        sorted(set(sequence).difference(CANONICAL_AA))
    )


# -----------------------------------------------------------------------------
# Eligibility
# -----------------------------------------------------------------------------

def load_sequence_universe(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(path)

    df = pd.read_csv(
        path,
        sep="\t",
        compression="infer",
        dtype=str,
        low_memory=False,
    ).fillna("")

    required = {"sequence_uid", "sequence"}
    missing = required - set(df.columns)
    if missing:
        raise KeyError(
            f"Stage-10A universe missing columns: {sorted(missing)}"
        )

    if df["sequence_uid"].duplicated().any():
        raise RuntimeError(
            "Duplicate sequence_uid values in Stage-10A universe."
        )

    df["sequence"] = df["sequence"].map(normalize_sequence)
    df["length"] = df["sequence"].str.len().astype(int)

    # Stable model input order, independent of any accidental file reordering.
    df = df.sort_values(
        "sequence_uid",
        kind="stable",
    ).reset_index(drop=True)

    return df


def build_eligibility(universe: pd.DataFrame) -> pd.DataFrame:
    out = universe[
        ["sequence_uid", "sequence", "length"]
    ].copy()

    out["length_eligible"] = out["length"].le(MAX_LENGTH)
    out["noncanonical_residues"] = out["sequence"].map(
        noncanonical_residues
    )
    out["canonical_amino_acids_only"] = out[
        "noncanonical_residues"
    ].eq("")
    out["representation_eligible"] = (
        out["length_eligible"]
        & out["canonical_amino_acids_only"]
    )

    def reason(row: pd.Series) -> str:
        reasons = []
        if not bool(row["length_eligible"]):
            reasons.append("sequence_length_gt_1024")
        if not bool(row["canonical_amino_acids_only"]):
            reasons.append("contains_noncanonical_amino_acid")
        return ";".join(reasons)

    out["representation_filter_reason"] = out.apply(
        reason,
        axis=1,
    )

    return out


def materialize_eligibility(
    eligibility: pd.DataFrame,
    output_dir: Path,
) -> tuple[pd.DataFrame, Path, Path, Path]:
    eligible = (
        eligibility.loc[eligibility["representation_eligible"]]
        .copy()
        .reset_index(drop=True)
    )
    filtered = (
        eligibility.loc[~eligibility["representation_eligible"]]
        .copy()
        .reset_index(drop=True)
    )

    eligibility_path = output_dir / "eligibility.tsv.gz"
    filtered_path = output_dir / "filtered_sequences.tsv.gz"
    eligible_path = output_dir / "eligible_sequences.parquet"

    eligibility.to_csv(
        eligibility_path,
        sep="\t",
        index=False,
        compression="gzip",
    )
    filtered.to_csv(
        filtered_path,
        sep="\t",
        index=False,
        compression="gzip",
    )

    # This is the reusable Sylphy input. Keeping sequence_uid here makes the
    # frozen representation universe independently inspectable, although the
    # CLI output is normalized by row order after extraction.
    eligible[
        ["sequence_uid", "sequence", "length"]
    ].to_parquet(
        eligible_path,
        index=False,
        compression="zstd",
    )

    return eligible, eligibility_path, filtered_path, eligible_path


# -----------------------------------------------------------------------------
# Sylphy execution
# -----------------------------------------------------------------------------

def sylphy_command(
    *,
    model: str,
    batch_size: int,
    input_parquet: Path,
    raw_output: Path,
    device: str,
    precision: str,
) -> list[str]:
    return [
        "sylphy",
        "embed",
        "--input",
        str(input_parquet),
        "--seq-col",
        "sequence",
        "--output",
        str(raw_output),
        "--model",
        model,
        "--device",
        device,
        "--precision",
        precision,
        "--batch-size",
        str(batch_size),
        "--max-length",
        str(MAX_LENGTH),
        "--layers",
        LAYERS,
        "--layer-agg",
        LAYER_AGG,
        "--pool",
        POOL,
        "--oom-backoff",
    ]


def run_sylphy(
    command: list[str],
    log_file: Path,
) -> tuple[int, float]:
    start = time.monotonic()

    with log_file.open("w", encoding="utf-8") as handle:
        handle.write(
            "# PhageLysClass Stage 10D — Sylphy execution\n"
        )
        handle.write("# command: " + " ".join(command) + "\n\n")
        handle.flush()

        result = subprocess.run(
            command,
            stdout=handle,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        )

    return result.returncode, time.monotonic() - start


# -----------------------------------------------------------------------------
# Parquet normalization
# -----------------------------------------------------------------------------

NON_EMBEDDING_COLUMN_NAMES = {
    "sequence",
    "seq",
    "sequence_uid",
    "id",
    "length",
    "index",
}


def is_numeric_arrow_type(data_type: pa.DataType) -> bool:
    return (
        pa.types.is_integer(data_type)
        or pa.types.is_floating(data_type)
    )


def detect_embedding_columns(
    parquet_file: pq.ParquetFile,
) -> list[str]:
    schema = parquet_file.schema_arrow
    cols = []

    for field in schema:
        if field.name.lower() in NON_EMBEDDING_COLUMN_NAMES:
            continue
        if not is_numeric_arrow_type(field.type):
            continue
        cols.append(field.name)

    if not cols:
        raise RuntimeError(
            "No numeric embedding columns detected in Sylphy Parquet output."
        )

    return cols


def _scalar_bool(value: Any) -> bool:
    try:
        return bool(value.as_py())
    except Exception:
        return bool(value)


def validate_finite_arrow(array: pa.Array, column_name: str) -> None:
    if pa.types.is_floating(array.type):
        nan_any = pc.any(pc.is_nan(array))
        inf_any = pc.any(pc.is_inf(array))

        if _scalar_bool(nan_any) or _scalar_bool(inf_any):
            raise RuntimeError(
                f"Invalid NaN/Inf values detected in {column_name}."
            )


def postprocess_raw_parquet(
    *,
    raw_output: Path,
    final_output: Path,
    eligible: pd.DataFrame,
) -> dict[str, Any]:
    pf = pq.ParquetFile(raw_output)

    expected_rows = len(eligible)
    observed_rows = pf.metadata.num_rows

    if observed_rows != expected_rows:
        raise RuntimeError(
            f"Row-count mismatch: eligible={expected_rows:,}, "
            f"Sylphy={observed_rows:,}"
        )

    embedding_columns = detect_embedding_columns(pf)
    embedding_dimensions = len(embedding_columns)

    schema_names = set(pf.schema_arrow.names)
    sequence_col = None
    for candidate in ["sequence", "seq"]:
        if candidate in schema_names:
            sequence_col = candidate
            break

    read_columns = list(embedding_columns)
    if sequence_col is not None:
        read_columns.append(sequence_col)

    out_fields = [
        pa.field("sequence_uid", pa.string()),
        pa.field("length", pa.int32()),
    ] + [
        pa.field(f"p_{i}", pa.float32())
        for i in range(embedding_dimensions)
    ]
    out_schema = pa.schema(out_fields)

    tmp_output = final_output.with_suffix(".tmp.parquet")
    tmp_output.unlink(missing_ok=True)

    writer = pq.ParquetWriter(
        tmp_output,
        schema=out_schema,
        compression="zstd",
        use_dictionary=["sequence_uid"],
        write_statistics=True,
    )

    offset = 0

    try:
        for batch in pf.iter_batches(
            batch_size=POSTPROCESS_BATCH_ROWS,
            columns=read_columns,
        ):
            table = pa.Table.from_batches([batch])
            n = table.num_rows

            expected_block = eligible.iloc[
                offset: offset + n
            ]

            if len(expected_block) != n:
                raise RuntimeError(
                    "Internal row alignment error during Parquet normalization."
                )

            if sequence_col is not None:
                observed_sequences = (
                    table[sequence_col]
                    .combine_chunks()
                    .to_pylist()
                )
                expected_sequences = expected_block[
                    "sequence"
                ].astype(str).tolist()

                mismatch_count = sum(
                    a != b
                    for a, b in zip(
                        observed_sequences,
                        expected_sequences,
                    )
                )
                if mismatch_count:
                    raise RuntimeError(
                        f"Sylphy output sequence order/content mismatch in "
                        f"stream block: {mismatch_count:,} rows."
                    )

            arrays: list[pa.Array] = [
                pa.array(
                    expected_block["sequence_uid"].astype(str).tolist(),
                    type=pa.string(),
                ),
                pa.array(
                    expected_block["length"].astype(int).tolist(),
                    type=pa.int32(),
                ),
            ]

            for col in embedding_columns:
                arr = table[col].combine_chunks()
                validate_finite_arrow(arr, col)
                arr32 = pc.cast(
                    arr,
                    pa.float32(),
                    safe=False,
                )
                arrays.append(arr32)

            out_table = pa.Table.from_arrays(
                arrays,
                schema=out_schema,
            )

            writer.write_table(
                out_table,
                row_group_size=n,
            )

            offset += n

    finally:
        writer.close()

    if offset != expected_rows:
        tmp_output.unlink(missing_ok=True)
        raise RuntimeError(
            f"Streamed {offset:,} rows; expected {expected_rows:,}."
        )

    tmp_output.replace(final_output)

    verify = pq.ParquetFile(final_output)

    if verify.metadata.num_rows != expected_rows:
        raise RuntimeError(
            "Final Parquet row count changed after serialization."
        )

    verify_schema = verify.schema_arrow

    if verify_schema.field("sequence_uid").type != pa.string():
        raise RuntimeError(
            "Final sequence_uid dtype is not string."
        )

    feature_names = [
        name
        for name in verify_schema.names
        if name.startswith("p_")
    ]

    if len(feature_names) != embedding_dimensions:
        raise RuntimeError(
            "Final embedding dimension mismatch."
        )

    non_fp32 = [
        name
        for name in feature_names
        if verify_schema.field(name).type != pa.float32()
    ]

    if non_fp32:
        raise RuntimeError(
            f"{len(non_fp32):,} embedding columns are not FP32."
        )

    return {
        "row_count": int(expected_rows),
        "embedding_dimensions": int(embedding_dimensions),
        "parquet_size_bytes": int(final_output.stat().st_size),
        "parquet_sha256": sha256_file(final_output),
        "embedding_dtype": "float32",
        "sequence_column_verified": sequence_col is not None,
    }


# -----------------------------------------------------------------------------
# Model completion / status
# -----------------------------------------------------------------------------

def model_paths(
    models_dir: Path,
    model: str,
) -> dict[str, Path]:
    directory = models_dir / slugify(model)

    return {
        "dir": directory,
        "raw": directory / ".sylphy_raw.parquet",
        "parquet": directory / "embeddings_last_mean_fp32.parquet",
        "metadata": directory / "metadata.json",
        "log": directory / "run.log",
    }


def valid_completed_model(
    *,
    model: str,
    model_dir: Path,
    expected_rows: int,
    device: str,
    precision: str,
) -> tuple[bool, dict[str, Any] | None]:
    metadata_path = model_dir / "metadata.json"
    parquet_path = model_dir / "embeddings_last_mean_fp32.parquet"

    if not metadata_path.is_file() or not parquet_path.is_file():
        return False, None

    try:
        md = json.loads(
            metadata_path.read_text(encoding="utf-8")
        )

        expected_contract = {
            "status": "completed",
            "model": model,
            "eligible_sequence_count": expected_rows,
            "max_length": MAX_LENGTH,
            "canonical_amino_acids_only": True,
            "precision": precision,
            "layers": LAYERS,
            "layer_agg": LAYER_AGG,
            "pool": POOL,
        }

        for key, expected in expected_contract.items():
            if md.get(key) != expected:
                return False, md

        pf = pq.ParquetFile(parquet_path)

        if pf.metadata.num_rows != expected_rows:
            return False, md

        schema = pf.schema_arrow

        if "sequence_uid" not in schema.names:
            return False, md

        features = [
            name for name in schema.names
            if name.startswith("p_")
        ]

        if len(features) != int(md.get("embedding_dimensions", -1)):
            return False, md

        if any(
            schema.field(name).type != pa.float32()
            for name in features
        ):
            return False, md

        if parquet_path.stat().st_size != int(
            md.get("parquet_size_bytes", -1)
        ):
            return False, md

        return True, md

    except Exception:
        return False, None


def collect_status(
    *,
    models_dir: Path,
    expected_rows: int,
    device: str,
    precision: str,
) -> pd.DataFrame:
    rows = []

    for model in FROZEN_MODELS:
        paths = model_paths(models_dir, model)

        valid, metadata = valid_completed_model(
            model=model,
            model_dir=paths["dir"],
            expected_rows=expected_rows,
            device=device,
            precision=precision,
        )

        if valid:
            status = "PASS"
        else:
            status = str(
                (metadata or {}).get("status", "PENDING")
            ).upper()

        rows.append(
            {
                "model": model,
                "model_slug": slugify(model),
                "requested_batch_size": MODEL_BATCH_SIZES.get(
                    model,
                    DEFAULT_BATCH_SIZE,
                ),
                "eligible_sequences": expected_rows,
                "embedding_dimensions": (
                    (metadata or {}).get(
                        "embedding_dimensions",
                        "",
                    )
                ),
                "precision": precision,
                "layers": LAYERS,
                "layer_agg": LAYER_AGG,
                "pool": POOL,
                "max_length": MAX_LENGTH,
                "parquet_output": str(paths["parquet"]),
                "status": status,
            }
        )

    return pd.DataFrame(rows)


# -----------------------------------------------------------------------------
# GPU helpers
# -----------------------------------------------------------------------------

def get_gpu_memory_used_mb() -> int | None:
    try:
        result = subprocess.run(
            [
                "nvidia-smi",
                "-i",
                str(GPU_INDEX),
                "--query-gpu=memory.used",
                "--format=csv,noheader,nounits",
            ],
            capture_output=True,
            text=True,
            check=True,
        )
        text = result.stdout.strip()
        return int(text.splitlines()[0]) if text else None
    except Exception:
        return None


def wait_for_gpu_release(
    baseline_memory_mb: int | None,
) -> None:
    time.sleep(POST_MODEL_SLEEP)

    if baseline_memory_mb is None:
        return

    target = baseline_memory_mb + GPU_MEMORY_TOLERANCE_MB
    deadline = time.monotonic() + GPU_RELEASE_TIMEOUT

    while time.monotonic() < deadline:
        current = get_gpu_memory_used_mb()

        if current is None or current <= target:
            return

        time.sleep(GPU_CHECK_INTERVAL)

    current = get_gpu_memory_used_mb()

    LOGGER.warning(
        "GPU memory remained above expected baseline after %s s "
        "| baseline=%s MB | current=%s MB",
        GPU_RELEASE_TIMEOUT,
        baseline_memory_mb,
        current,
    )


# -----------------------------------------------------------------------------
# CLI / main
# -----------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    root = project_root_from_script(__file__)

    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        "--input",
        default=str(
            root
            / "data/representations/sequence_universe/unique_sequences.tsv.gz"
        ),
    )
    parser.add_argument(
        "--output-dir",
        default=str(
            root / "data/representations/sylphy_embeddings"
        ),
    )
    parser.add_argument(
        "--model",
        action="append",
        choices=FROZEN_MODELS,
        help=(
            "Run only selected frozen model(s). Repeatable. "
            "If omitted, all 11 models are run sequentially."
        ),
    )
    parser.add_argument(
        "--force-model",
        action="append",
        choices=FROZEN_MODELS,
        default=[],
        help="Recompute model even if a valid completed Parquet exists.",
    )
    parser.add_argument(
        "--device",
        default=DEVICE,
        choices=["cuda", "cpu"],
    )
    parser.add_argument(
        "--precision",
        default=PRECISION,
        choices=["fp32", "fp16", "bf16"],
        help=(
            "Frozen paper default is fp32. Use another value only for "
            "explicit exploratory runs."
        ),
    )
    parser.add_argument(
        "--preflight-only",
        action="store_true",
    )
    parser.add_argument(
        "--require-all",
        action="store_true",
        help="Fail at closure unless all 11 frozen models are complete.",
    )
    parser.add_argument(
        "--keep-raw",
        action="store_true",
        help="Keep Sylphy's temporary raw Parquet after normalization.",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
    )

    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if args.precision != PRECISION:
        LOGGER.warning(
            "Non-frozen precision requested: %s. Paper contract is fp32.",
            args.precision,
        )

    input_path = Path(args.input).expanduser().resolve()
    output_dir = Path(args.output_dir).expanduser().resolve()
    models_dir = output_dir / "models"

    output_dir.mkdir(parents=True, exist_ok=True)
    models_dir.mkdir(parents=True, exist_ok=True)

    configure_logging(
        args.log_level,
        output_dir / "stage10d.log",
    )

    LOGGER.info(
        "Stage 10D: load Stage-10A exact-sequence universe."
    )

    universe = load_sequence_universe(input_path)
    eligibility = build_eligibility(universe)

    eligible, eligibility_path, filtered_path, eligible_path = (
        materialize_eligibility(
            eligibility,
            output_dir,
        )
    )

    n_total = len(universe)
    n_eligible = len(eligible)
    n_filtered = n_total - n_eligible

    length_filtered = int(
        (~eligibility["length_eligible"]).sum()
    )
    noncanonical_filtered = int(
        (~eligibility["canonical_amino_acids_only"]).sum()
    )

    LOGGER.info(
        "Eligibility | total=%s | eligible=%s | filtered=%s "
        "| >1024=%s | noncanonical=%s",
        f"{n_total:,}",
        f"{n_eligible:,}",
        f"{n_filtered:,}",
        f"{length_filtered:,}",
        f"{noncanonical_filtered:,}",
    )

    # Freeze registry table.
    registry = pd.DataFrame(
        [
            {
                "model": model,
                "model_slug": slugify(model),
                "requested_batch_size": MODEL_BATCH_SIZES.get(
                    model,
                    DEFAULT_BATCH_SIZE,
                ),
                "max_length": MAX_LENGTH,
                "canonical_only": True,
                "layers": LAYERS,
                "layer_agg": LAYER_AGG,
                "pool": POOL,
                "precision": PRECISION,
            }
            for model in FROZEN_MODELS
        ]
    )

    registry_path = output_dir / "model_registry.tsv"
    registry.to_csv(
        registry_path,
        sep="\t",
        index=False,
    )

    sylphy_exe = shutil.which("sylphy")

    if sylphy_exe is None:
        raise RuntimeError("'sylphy' executable not found in PATH.")

    if args.device == "cuda" and shutil.which("nvidia-smi") is None:
        raise RuntimeError(
            "CUDA requested but nvidia-smi is unavailable."
        )

    runtime = {
        "python": sys.version.split()[0],
        "sylphy": package_version("sylphy"),
        "pandas": pd.__version__,
        "numpy": np.__version__,
        "pyarrow": package_version("pyarrow"),
        "device": args.device,
        "precision": args.precision,
        "sylphy_executable": sylphy_exe,
    }

    if args.preflight_only:
        print("=" * 96)
        print("Stage 10D preflight PASS")
        print("=" * 96)
        print(f"Total sequences      : {n_total:,}")
        print(f"Eligible             : {n_eligible:,}")
        print(f"Filtered             : {n_filtered:,}")
        print(f"Length >1024         : {length_filtered:,}")
        print(f"Noncanonical         : {noncanonical_filtered:,}")
        print(f"Frozen PLMs          : {len(FROZEN_MODELS)}")
        print(f"Input Parquet        : {eligible_path}")
        return

    selected_models = (
        args.model
        if args.model
        else list(FROZEN_MODELS)
    )
    force_models = set(args.force_model or [])

    for idx, model in enumerate(selected_models, start=1):
        paths = model_paths(models_dir, model)
        paths["dir"].mkdir(parents=True, exist_ok=True)

        valid, metadata = valid_completed_model(
            model=model,
            model_dir=paths["dir"],
            expected_rows=n_eligible,
            device=args.device,
            precision=args.precision,
        )

        if valid and model not in force_models:
            LOGGER.info(
                "[%s/%s] %s | valid completed output found; skip.",
                idx,
                len(selected_models),
                model,
            )
            continue

        batch_size = MODEL_BATCH_SIZES.get(
            model,
            DEFAULT_BATCH_SIZE,
        )

        paths["raw"].unlink(missing_ok=True)
        paths["parquet"].unlink(missing_ok=True)

        command = sylphy_command(
            model=model,
            batch_size=batch_size,
            input_parquet=eligible_path,
            raw_output=paths["raw"],
            device=args.device,
            precision=args.precision,
        )

        LOGGER.info(
            "[%s/%s] %s | batch=%s | eligible=%s | start.",
            idx,
            len(selected_models),
            model,
            batch_size,
            f"{n_eligible:,}",
        )

        baseline_gpu = get_gpu_memory_used_mb()
        started_at = now_iso()

        try:
            return_code, runtime_seconds = run_sylphy(
                command,
                paths["log"],
            )

            if return_code != 0:
                write_json(
                    {
                        "script": SCRIPT_NAME,
                        "script_version": SCRIPT_VERSION,
                        "contract_version": CONTRACT_VERSION,
                        "status": "failed",
                        "failure_stage": "sylphy_execution",
                        "return_code": return_code,
                        "model": model,
                        "requested_batch_size": batch_size,
                        "eligible_sequence_count": n_eligible,
                        "max_length": MAX_LENGTH,
                        "canonical_amino_acids_only": True,
                        "precision": args.precision,
                        "layers": LAYERS,
                        "layer_agg": LAYER_AGG,
                        "pool": POOL,
                        "started_at": started_at,
                        "finished_at": now_iso(),
                        "runtime_seconds": runtime_seconds,
                        "command": command,
                    },
                    paths["metadata"],
                )

                LOGGER.error(
                    "[%s/%s] %s | Sylphy FAIL; see %s",
                    idx,
                    len(selected_models),
                    model,
                    paths["log"],
                )
                continue

            if not paths["raw"].is_file():
                raise RuntimeError(
                    "Sylphy returned success but raw Parquet was not created."
                )

            LOGGER.info(
                "[%s/%s] %s | normalize Parquet / enforce FP32.",
                idx,
                len(selected_models),
                model,
            )

            stats = postprocess_raw_parquet(
                raw_output=paths["raw"],
                final_output=paths["parquet"],
                eligible=eligible,
            )

            metadata = {
                "script": SCRIPT_NAME,
                "script_version": SCRIPT_VERSION,
                "contract_version": CONTRACT_VERSION,
                "status": "completed",
                "model": model,
                "model_slug": slugify(model),
                "requested_batch_size": batch_size,
                "frozen_sequence_count": n_total,
                "eligible_sequence_count": n_eligible,
                "filtered_sequence_count": n_filtered,
                "max_length": MAX_LENGTH,
                "canonical_alphabet": CANONICAL_AA_STRING,
                "canonical_amino_acids_only": True,
                "precision": args.precision,
                "layers": LAYERS,
                "layer_agg": LAYER_AGG,
                "pool": POOL,
                "device": args.device,
                "started_at": started_at,
                "finished_at": now_iso(),
                "runtime_seconds": runtime_seconds,
                "command": command,
                "parquet_output": str(paths["parquet"]),
                **stats,
            }

            write_json(
                metadata,
                paths["metadata"],
            )

            if not args.keep_raw:
                paths["raw"].unlink(missing_ok=True)

            LOGGER.info(
                "[%s/%s] %s | PASS | dim=%s | rows=%s | %.2f h",
                idx,
                len(selected_models),
                model,
                stats["embedding_dimensions"],
                f"{stats['row_count']:,}",
                runtime_seconds / 3600.0,
            )

        except Exception as exc:
            paths["parquet"].unlink(missing_ok=True)

            if not args.keep_raw:
                paths["raw"].unlink(missing_ok=True)

            write_json(
                {
                    "script": SCRIPT_NAME,
                    "script_version": SCRIPT_VERSION,
                    "contract_version": CONTRACT_VERSION,
                    "status": "failed",
                    "failure_stage": "postprocess_or_validation",
                    "error": repr(exc),
                    "model": model,
                    "requested_batch_size": batch_size,
                    "eligible_sequence_count": n_eligible,
                    "max_length": MAX_LENGTH,
                    "canonical_amino_acids_only": True,
                    "precision": args.precision,
                    "layers": LAYERS,
                    "layer_agg": LAYER_AGG,
                    "pool": POOL,
                    "started_at": started_at,
                    "finished_at": now_iso(),
                    "command": command,
                },
                paths["metadata"],
            )

            LOGGER.exception(
                "[%s/%s] %s | FAIL.",
                idx,
                len(selected_models),
                model,
            )

        finally:
            wait_for_gpu_release(baseline_gpu)

            status = collect_status(
                models_dir=models_dir,
                expected_rows=n_eligible,
                device=args.device,
                precision=args.precision,
            )
            status.to_csv(
                output_dir / "run_status.tsv",
                sep="\t",
                index=False,
            )

    status = collect_status(
        models_dir=models_dir,
        expected_rows=n_eligible,
        device=args.device,
        precision=args.precision,
    )

    status_path = output_dir / "run_status.tsv"
    status.to_csv(
        status_path,
        sep="\t",
        index=False,
    )

    completed_models = int(
        status["status"].eq("PASS").sum()
    )

    reason_counts = (
        eligibility.loc[
            ~eligibility["representation_eligible"],
            "representation_filter_reason",
        ]
        .value_counts()
        .to_dict()
    )

    manifest = {
        "stage": "10D",
        "name": "frozen_sylphy_plm_embeddings",
        "contract_version": CONTRACT_VERSION,
        "input": str(input_path),
        "input_sha256": sha256_file(input_path),
        "eligibility": {
            "maximum_sequence_length": MAX_LENGTH,
            "canonical_amino_acids_only": True,
            "canonical_alphabet": CANONICAL_AA_STRING,
            "frozen_unique_sequences": n_total,
            "eligible_sequences": n_eligible,
            "filtered_sequences": n_filtered,
            "filtered_length_gt_1024": length_filtered,
            "filtered_noncanonical": noncanonical_filtered,
            "filter_reason_counts": {
                str(k): int(v)
                for k, v in reason_counts.items()
            },
        },
        "representation": {
            "frozen_models": FROZEN_MODELS,
            "frozen_model_count": len(FROZEN_MODELS),
            "completed_model_count": completed_models,
            "max_length": MAX_LENGTH,
            "precision": PRECISION,
            "layers": LAYERS,
            "layer_aggregation": LAYER_AGG,
            "pooling": POOL,
            "authoritative_join_key": "sequence_uid",
            "output_format": "parquet",
            "embedding_dtype": "float32",
            "csv_embedding_export": False,
        },
        "runtime": runtime,
        "outputs": {
            "eligibility": str(eligibility_path),
            "eligible_sequences": str(eligible_path),
            "filtered_sequences": str(filtered_path),
            "model_registry": str(registry_path),
            "run_status": str(status_path),
            "models_directory": str(models_dir),
        },
        "cross_representation_contract": {
            "onehot_should_use_max_length_1024": True,
            "onehot_should_use_canonical_only": True,
            "recommended_stage10c_cli": (
                "--max-length 1024 --canonical-only"
            ),
        },
    }

    manifest_path = output_dir / "manifest.json"
    write_json(
        manifest,
        manifest_path,
    )

    print()
    print("=" * 96)
    print("PhageLysClass Stage 10D — PLM representation summary")
    print("=" * 96)
    print(f"Frozen sequences       : {n_total:,}")
    print(f"Eligible sequences     : {n_eligible:,}")
    print(f"Filtered sequences     : {n_filtered:,}")
    print(f"  length >1024         : {length_filtered:,}")
    print(f"  noncanonical         : {noncanonical_filtered:,}")
    print(f"Frozen models complete : {completed_models}/{len(FROZEN_MODELS)}")
    print(f"Max length             : {MAX_LENGTH}")
    print(f"Precision              : {PRECISION}")
    print(f"Layer / pooling        : {LAYERS} / {POOL}")
    print(f"Output                  : Parquet only")
    print(f"Run status              : {status_path}")
    print(f"Manifest                : {manifest_path}")

    if args.require_all and completed_models != len(FROZEN_MODELS):
        missing = status.loc[
            status["status"] != "PASS",
            "model",
        ].tolist()

        raise RuntimeError(
            "Stage 10D is incomplete. Missing/failed frozen models: "
            + ", ".join(missing)
        )


if __name__ == "__main__":
    main()

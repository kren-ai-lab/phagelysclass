#!/usr/bin/env python3
"""
10c_encode_unique_sequences_onehot_sylphy.py

Stage 10C — One-hot encoding of the frozen exact-sequence universe with Sylphy.

The exact sequence universe is encoded once and indexed by sequence_uid.
Because full-length protein one-hot matrices are extremely wide, Sylphy's
dense uint8 output is converted immediately to SciPy CSR sparse matrices.

Default behavior:
- use Sylphy's extended residue alphabet;
- set max_length to the longest frozen sequence;
- therefore retain all sequences supported by Sylphy's alphabet;
- process in resumable chunks;
- persist sparse CSR, never a dense master matrix.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import math
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

from sylphy.constants import residues
from sylphy.sequence_encoder.encoder_base import encode_onehot


def project_root_from_script(script_file: str) -> Path:
    return Path(script_file).resolve().parent.parents[1]


def parse_args() -> argparse.Namespace:
    root = project_root_from_script(__file__)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        default=str(root / "data/representations/sequence_universe/unique_sequences.tsv.gz"),
    )
    parser.add_argument(
        "--output-dir",
        default=str(root / "data/representations/sylphy_onehot"),
    )
    parser.add_argument(
        "--max-length",
        type=int,
        default=0,
        help="0 = use longest eligible input sequence; set 1024 explicitly for a capped baseline.",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=512,
        help="Sequences per dense Sylphy batch before immediate CSR conversion.",
    )
    parser.add_argument("--canonical-only", action="store_true")
    parser.add_argument("--max-chunks", type=int, default=None)
    parser.add_argument("--no-resume", action="store_true")
    parser.add_argument("--no-finalize", action="store_true")
    parser.add_argument("--keep-chunks", action="store_true")
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Fail if any frozen sequence is excluded.",
    )
    return parser.parse_args()


def package_version(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except Exception:
        return "unknown"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def normalize_sequence(seq: str) -> str:
    return "".join(str(seq).split()).upper()


def build_alphabet(canonical_only: bool) -> tuple[str, ...]:
    alphabet = tuple(residues(extended=not canonical_only))
    if not alphabet or len(alphabet) != len(set(alphabet)):
        raise RuntimeError("Invalid Sylphy residue alphabet.")
    return alphabet


def classify_eligibility(
    sequence: str,
    alphabet_set: set[str],
    max_length: int,
) -> tuple[bool, str, str]:
    unsupported = "".join(sorted(set(sequence) - alphabet_set))
    if unsupported:
        return False, "unsupported_alphabet", unsupported
    if len(sequence) > max_length:
        return False, "longer_than_max_length", ""
    return True, "", ""


def validate_chunk(
    path: Path,
    expected_shape: tuple[int, int],
    expected_nnz: int,
) -> bool:
    if not path.is_file():
        return False
    try:
        matrix = sparse.load_npz(path)
        return (
            sparse.isspmatrix_csr(matrix)
            and matrix.shape == expected_shape
            and matrix.dtype == np.uint8
            and int(matrix.nnz) == int(expected_nnz)
        )
    except Exception:
        return False


def sylphy_encode_chunk(
    sequences: list[str],
    alphabet: tuple[str, ...],
    max_length: int,
) -> sparse.csr_matrix:
    dense = encode_onehot(
        sequences=sequences,
        alphabet=alphabet,
        max_length=max_length,
    )
    expected_shape = (len(sequences), max_length * len(alphabet))
    if dense.shape != expected_shape:
        raise RuntimeError(
            f"Sylphy returned {dense.shape}, expected {expected_shape}."
        )
    if dense.dtype != np.uint8:
        dense = dense.astype(np.uint8, copy=False)

    matrix = sparse.csr_matrix(dense, dtype=np.uint8)
    matrix.sum_duplicates()
    matrix.sort_indices()
    del dense
    return matrix


def main() -> None:
    args = parse_args()

    input_path = Path(args.input).expanduser().resolve()
    output_dir = Path(args.output_dir).expanduser().resolve()
    chunk_dir = output_dir / "chunks"

    if not input_path.is_file():
        raise FileNotFoundError(input_path)
    if args.chunk_size <= 0:
        raise ValueError("--chunk-size must be > 0")

    output_dir.mkdir(parents=True, exist_ok=True)
    chunk_dir.mkdir(parents=True, exist_ok=True)

    universe = pd.read_csv(
        input_path,
        sep="\t",
        compression="infer",
        dtype=str,
        low_memory=False,
    ).fillna("")

    required = {"sequence_uid", "sequence"}
    missing = required - set(universe.columns)
    if missing:
        raise KeyError(f"Missing required columns: {sorted(missing)}")
    if universe["sequence_uid"].duplicated().any():
        raise RuntimeError("Duplicate sequence_uid values in Stage 10A universe.")

    universe["sequence"] = universe["sequence"].map(normalize_sequence)
    universe["length"] = universe["sequence"].str.len().astype(int)

    alphabet = build_alphabet(args.canonical_only)
    alphabet_set = set(alphabet)

    observed_max = int(universe["length"].max())
    max_length = int(args.max_length) if args.max_length > 0 else observed_max

    eligibility = universe["sequence"].map(
        lambda seq: classify_eligibility(seq, alphabet_set, max_length)
    )
    universe["eligible"] = eligibility.map(lambda x: x[0])
    universe["exclusion_reason"] = eligibility.map(lambda x: x[1])
    universe["unsupported_symbols"] = eligibility.map(lambda x: x[2])

    eligible = universe.loc[universe["eligible"]].copy().reset_index(drop=True)
    excluded = universe.loc[~universe["eligible"]].copy()

    eligible["onehot_row"] = np.arange(len(eligible), dtype=np.int64)

    n_total = len(universe)
    n_eligible = len(eligible)
    n_excluded = len(excluded)
    n_channels = len(alphabet)
    n_features = max_length * n_channels
    n_chunks = math.ceil(n_eligible / args.chunk_size)

    print("=" * 96)
    print("PhageLysClass — Stage 10C Sylphy one-hot")
    print("=" * 96)
    print(f"Input sequences    : {n_total:,}")
    print(f"Eligible           : {n_eligible:,}")
    print(f"Excluded           : {n_excluded:,}")
    print(f"Sylphy version     : {package_version('sylphy')}")
    print(f"Alphabet mode      : {'canonical' if args.canonical_only else 'extended'}")
    print(f"Alphabet ({n_channels})        : {''.join(alphabet)}")
    print(f"max_length         : {max_length:,}")
    print(f"Feature dimension  : {n_features:,}")
    print(f"Chunk size         : {args.chunk_size:,}")
    print(f"Chunks             : {n_chunks:,}")

    excluded_path = output_dir / "excluded_sequences.tsv.gz"
    excluded[
        ["sequence_uid", "sequence", "length", "exclusion_reason", "unsupported_symbols"]
    ].to_csv(excluded_path, sep="\t", index=False, compression="gzip")

    row_index_path = output_dir / "row_index.tsv.gz"
    eligible[["onehot_row", "sequence_uid", "length"]].to_csv(
        row_index_path, sep="\t", index=False, compression="gzip"
    )

    if n_excluded:
        print("\nExclusions:")
        for reason, count in excluded["exclusion_reason"].value_counts().items():
            print(f"  {reason:28s}: {count:,}")

    n_chunks_to_run = n_chunks if args.max_chunks is None else min(
        n_chunks, max(0, args.max_chunks)
    )

    status_rows = []

    for chunk_idx in range(n_chunks_to_run):
        start = chunk_idx * args.chunk_size
        stop = min(start + args.chunk_size, n_eligible)
        block = eligible.iloc[start:stop]
        sequences = block["sequence"].tolist()
        expected_nnz = int(sum(len(seq) for seq in sequences))
        expected_shape = (len(block), n_features)

        chunk_path = chunk_dir / f"part_{chunk_idx:06d}.npz"

        if (
            not args.no_resume
            and validate_chunk(chunk_path, expected_shape, expected_nnz)
        ):
            print(
                f"[{chunk_idx + 1:>5}/{n_chunks_to_run:<5}] "
                f"reuse {chunk_path.name} rows={len(block):,} nnz={expected_nnz:,}"
            )
            status_rows.append({
                "chunk": chunk_idx,
                "start_row": start,
                "stop_row_exclusive": stop,
                "n_sequences": len(block),
                "n_features": n_features,
                "nnz": expected_nnz,
                "status": "reused",
                "path": str(chunk_path),
            })
            continue

        print(
            f"[{chunk_idx + 1:>5}/{n_chunks_to_run:<5}] "
            f"encode rows {start:,}:{stop:,} ({len(block):,})"
        )

        matrix = sylphy_encode_chunk(sequences, alphabet, max_length)

        if int(matrix.nnz) != expected_nnz:
            raise RuntimeError(
                f"Chunk {chunk_idx}: nnz={matrix.nnz:,}, expected={expected_nnz:,}."
            )

        tmp_path = chunk_path.with_suffix(".tmp.npz")
        sparse.save_npz(tmp_path, matrix, compressed=True)
        tmp_path.replace(chunk_path)

        status_rows.append({
            "chunk": chunk_idx,
            "start_row": start,
            "stop_row_exclusive": stop,
            "n_sequences": len(block),
            "n_features": n_features,
            "nnz": int(matrix.nnz),
            "status": "computed",
            "path": str(chunk_path),
        })
        del matrix

    status_path = output_dir / "chunk_status.tsv"
    pd.DataFrame(status_rows).to_csv(status_path, sep="\t", index=False)

    complete_run = n_chunks_to_run == n_chunks
    final_path = output_dir / "onehot_csr.npz"
    final_shape = None
    final_nnz = None

    if complete_run and not args.no_finalize:
        expected_paths = [
            chunk_dir / f"part_{idx:06d}.npz"
            for idx in range(n_chunks)
        ]
        missing_paths = [p for p in expected_paths if not p.is_file()]
        if missing_paths:
            raise RuntimeError(
                f"Cannot finalize: {len(missing_paths):,} chunk files missing."
            )

        print("\n[finalize] Loading sparse chunks...")
        matrices = [sparse.load_npz(path) for path in expected_paths]

        print("[finalize] Concatenating CSR matrices...")
        final_matrix = sparse.vstack(matrices, format="csr", dtype=np.uint8)
        del matrices

        final_shape = final_matrix.shape
        final_nnz = int(final_matrix.nnz)
        expected_total_nnz = int(eligible["length"].sum())

        if final_shape != (n_eligible, n_features):
            raise RuntimeError(
                f"Final shape {final_shape} != {(n_eligible, n_features)}"
            )
        if final_nnz != expected_total_nnz:
            raise RuntimeError(
                f"Final nnz={final_nnz:,} != expected={expected_total_nnz:,}"
            )

        print(f"[finalize] Saving {final_path.name}...")
        tmp_final = final_path.with_suffix(".tmp.npz")
        sparse.save_npz(tmp_final, final_matrix, compressed=True)
        tmp_final.replace(final_path)
        del final_matrix

        check = sparse.load_npz(final_path)
        if check.shape != final_shape or int(check.nnz) != final_nnz:
            raise RuntimeError("Final CSR failed post-save integrity validation.")
        del check

        if not args.keep_chunks:
            shutil.rmtree(chunk_dir, ignore_errors=True)

    density = (
        final_nnz / (n_eligible * n_features)
        if final_nnz is not None and n_eligible and n_features
        else None
    )

    manifest = {
        "stage": "10C",
        "name": "sylphy_onehot_exact_sequence_universe",
        "input": str(input_path),
        "input_sha256": sha256_file(input_path),
        "sylphy_version": package_version("sylphy"),
        "scipy_version": package_version("scipy"),
        "encoding": {
            "implementation": "sylphy.sequence_encoder.encoder_base.encode_onehot",
            "storage": "scipy_csr_uint8",
            "alphabet_mode": "canonical" if args.canonical_only else "extended",
            "alphabet": list(alphabet),
            "n_channels": n_channels,
            "max_length": max_length,
            "n_features": n_features,
            "padding": "all-zero",
        },
        "counts": {
            "frozen_unique_sequences": n_total,
            "eligible_sequences": n_eligible,
            "excluded_sequences": n_excluded,
            "total_chunks": n_chunks,
            "chunks_executed_or_reused": n_chunks_to_run,
            "complete_run": complete_run,
            "finalized": bool(complete_run and not args.no_finalize),
            "final_shape": list(final_shape) if final_shape is not None else None,
            "final_nnz": final_nnz,
            "final_density": density,
        },
        "exclusion_counts": (
            excluded["exclusion_reason"].value_counts().to_dict()
            if n_excluded else {}
        ),
        "outputs": {
            "final_onehot_csr": str(final_path) if final_shape is not None else None,
            "row_index": str(row_index_path),
            "excluded_sequences": str(excluded_path),
            "chunk_status": str(status_path),
            "chunk_directory": str(chunk_dir) if chunk_dir.exists() else None,
        },
        "design_decisions": {
            "computed_once_per_exact_sequence": True,
            "task_specific": False,
            "homology_regime_specific": False,
            "dense_master_matrix_persisted": False,
            "sparse_storage_required": True,
            "default_max_length_policy": (
                "observed_maximum" if args.max_length <= 0 else "user_fixed"
            ),
        },
    }

    manifest_path = output_dir / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    print("\n" + "=" * 96)
    print("Stage 10C summary")
    print("=" * 96)
    print(f"Eligible sequences : {n_eligible:,}/{n_total:,}")
    print(f"max_length         : {max_length:,}")
    print(f"Feature dimension  : {n_features:,}")
    print(f"Complete run       : {complete_run}")
    if final_shape is not None:
        print(f"Final shape        : {final_shape}")
        print(f"Final nnz          : {final_nnz:,}")
        print(f"Final density      : {density:.10f}")
        print(f"Final CSR          : {final_path}")
    print(f"Row index          : {row_index_path}")
    print(f"Excluded           : {excluded_path}")
    print(f"Manifest           : {manifest_path}")

    if args.strict and n_excluded:
        raise SystemExit(
            f"STRICT FAIL — {n_excluded:,} frozen sequences were excluded."
        )

    if complete_run and not args.no_finalize:
        print("PASS — Stage 10C Sylphy one-hot representation materialized.")
    elif complete_run:
        print("PASS — all one-hot chunks materialized; final CSR intentionally skipped.")
    else:
        print("PARTIAL PASS — requested one-hot chunks completed.")


if __name__ == "__main__":
    main()

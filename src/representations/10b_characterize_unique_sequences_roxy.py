#!/usr/bin/env python3
"""
10b_characterize_unique_sequences_roxy.py

Stage 10B — Compact Roxy characterization of the frozen exact-sequence universe.

Design
------
Compute a deliberately compact, interpretable descriptor panel once per exact
unique sequence. The panel is capped below 50 Roxy-derived descriptors and is
intended for dataset characterization / EDA, not as the main predictive
representation.

Selected panel:
    - 20 amino-acid frequencies
    - 18 physicochemical / acid-base / thermodynamic-proxy descriptors
    - 6 order / structural-propensity descriptors
    - 4 AAIndex residue-property descriptors

Total:
    48 Roxy descriptors

Important:
    AAIndex entries are treated here as published residue-property scales.
    They are NOT labelled as genuine evolutionary-profile descriptors.
    True evolutionary information (e.g. PSSM/MSA conservation) would require
    a separate workflow and is intentionally not fabricated from sequence-only
    descriptors.

Input
-----
data/representations/sequence_universe/unique_sequences.tsv.gz

Default outputs
---------------
data/representations/roxy_characterization/
├── chunks/
│   ├── part_000000.parquet
│   ├── part_000001.parquet
│   └── ...
├── characterization_48.parquet
├── descriptor_schema.tsv
├── chunk_status.tsv
├── failed_sequences.tsv
└── manifest.json

The computation is chunked and resumable. Existing valid chunk files are
reused by default.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import math
import shutil
from pathlib import Path
from typing import Any

import pandas as pd
import polars as pl

from roxy import DESCRIPTOR_REGISTRY


AA_ORDER = list("ACDEFGHIKLMNPQRSTVWY")

AAC_FEATURES = [f"aac_freq_{aa}" for aa in AA_ORDER]

GLOBAL_BASIC_FEATURES = [
    "global_basic_molecular_weight",
    "global_basic_aliphatic_index",
    "global_basic_hydropathy_mean",
    "global_basic_hydropathy_std",
    "global_basic_polarity_mean",
    "global_basic_polarity_std",
    "global_basic_aromatic_fraction",
    "global_basic_hydrophobic_fraction",
    "global_basic_charged_fraction",
    "global_basic_boman_index",
    "global_basic_net_charge_ph7",
    "global_basic_fcr",
    "global_basic_ncpr",
    "global_basic_donors_per_residue",
    "global_basic_acceptors_per_residue",
    "global_basic_shannon_entropy",
]

CHARGE_FEATURES = [
    "charge_net_charge_ph5p0",
    "charge_net_charge_ph9p0",
]

ORDER_DISORDER_FEATURES = [
    "order_disorder_disorder_fraction",
    "order_disorder_order_fraction",
    "order_disorder_flexibility_fraction",
]

STRUCTURAL_FEATURES = [
    "structural_propensity_helix_mean",
    "structural_propensity_sheet_mean",
    "structural_propensity_turn_mean",
]

AAINDEX_FEATURES = [
    "aaindex_ANDN920101_mean",
    "aaindex_ARGP820101_mean",
    "aaindex_KYTJ820101_mean",
    "aaindex_FAUJ880104_mean",
]

SELECTED_FEATURES = (
    AAC_FEATURES
    + GLOBAL_BASIC_FEATURES
    + CHARGE_FEATURES
    + ORDER_DISORDER_FEATURES
    + STRUCTURAL_FEATURES
    + AAINDEX_FEATURES
)

EXPECTED_DESCRIPTOR_COUNT = 48

DESCRIPTOR_FAMILIES = {
    "aac": AAC_FEATURES,
    "global_basic": GLOBAL_BASIC_FEATURES,
    "charge": CHARGE_FEATURES,
    "order_disorder": ORDER_DISORDER_FEATURES,
    "structural_propensity": STRUCTURAL_FEATURES,
    "aaindex": AAINDEX_FEATURES,
}

SCHEMA_ROWS = []

for feature in AAC_FEATURES:
    SCHEMA_ROWS.append(
        {
            "feature": feature,
            "category": "amino_acid_frequency",
            "roxy_descriptor": "aac",
            "interpretation": "Relative amino-acid frequency",
        }
    )

INTERPRETATIONS = {
    "global_basic_molecular_weight": "Approximate molecular weight",
    "global_basic_aliphatic_index": "Aliphatic side-chain content / thermostability proxy",
    "global_basic_hydropathy_mean": "Mean sequence hydropathy",
    "global_basic_hydropathy_std": "Hydropathy heterogeneity",
    "global_basic_polarity_mean": "Mean residue polarity",
    "global_basic_polarity_std": "Polarity heterogeneity",
    "global_basic_aromatic_fraction": "Fraction of aromatic residues",
    "global_basic_hydrophobic_fraction": "Fraction of hydrophobic residues",
    "global_basic_charged_fraction": "Fraction of charged residues",
    "global_basic_boman_index": "Protein-binding / interaction propensity proxy",
    "global_basic_net_charge_ph7": "Net charge near neutral pH",
    "global_basic_fcr": "Fraction of charged residues (FCR)",
    "global_basic_ncpr": "Net charge per residue (NCPR)",
    "global_basic_donors_per_residue": "Hydrogen-bond donors per residue",
    "global_basic_acceptors_per_residue": "Hydrogen-bond acceptors per residue",
    "global_basic_shannon_entropy": "Sequence compositional entropy",
    "charge_net_charge_ph5p0": "Net charge at pH 5",
    "charge_net_charge_ph9p0": "Net charge at pH 9",
    "order_disorder_disorder_fraction": "Disorder-promoting residue fraction",
    "order_disorder_order_fraction": "Order-promoting residue fraction",
    "order_disorder_flexibility_fraction": "Flexibility-promoting residue fraction",
    "structural_propensity_helix_mean": "Mean alpha-helix propensity",
    "structural_propensity_sheet_mean": "Mean beta-sheet propensity",
    "structural_propensity_turn_mean": "Mean turn propensity",
    "aaindex_ANDN920101_mean": "Mean ANDN920101 residue-property scale",
    "aaindex_ARGP820101_mean": "Mean ARGP820101 hydrophobicity scale",
    "aaindex_KYTJ820101_mean": "Mean KYTJ820101 hydropathy scale",
    "aaindex_FAUJ880104_mean": "Mean FAUJ880104 side-chain STERIMOL length scale",
}

for family, features in [
    ("global_basic", GLOBAL_BASIC_FEATURES),
    ("charge", CHARGE_FEATURES),
]:
    for feature in features:
        SCHEMA_ROWS.append(
            {
                "feature": feature,
                "category": "physicochemical_thermodynamic_proxy",
                "roxy_descriptor": family,
                "interpretation": INTERPRETATIONS[feature],
            }
        )

for family, features in [
    ("order_disorder", ORDER_DISORDER_FEATURES),
    ("structural_propensity", STRUCTURAL_FEATURES),
]:
    for feature in features:
        SCHEMA_ROWS.append(
            {
                "feature": feature,
                "category": "order_structural_propensity",
                "roxy_descriptor": family,
                "interpretation": INTERPRETATIONS[feature],
            }
        )

for feature in AAINDEX_FEATURES:
    SCHEMA_ROWS.append(
        {
            "feature": feature,
            "category": "aaindex_residue_property",
            "roxy_descriptor": "aaindex",
            "interpretation": INTERPRETATIONS[feature],
        }
    )


def project_root_from_script(script_file: str) -> Path:
    script_dir = Path(script_file).resolve().parent
    return script_dir.parents[1]


def parse_args() -> argparse.Namespace:
    root = project_root_from_script(__file__)

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        default=str(
            root
            / "data/representations/sequence_universe/unique_sequences.tsv.gz"
        ),
        help="Stage 10A exact-sequence universe.",
    )
    parser.add_argument(
        "--output-dir",
        default=str(
            root / "data/representations/roxy_characterization"
        ),
        help="Output directory.",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=5000,
        help="Sequences per Roxy batch (default: 5000).",
    )
    parser.add_argument(
        "--no-resume",
        action="store_true",
        help="Recompute chunks even if valid outputs already exist.",
    )
    parser.add_argument(
        "--max-chunks",
        type=int,
        default=None,
        help="Optional debugging limit.",
    )
    parser.add_argument(
        "--keep-chunks",
        action="store_true",
        help=(
            "Keep chunk Parquet files after final aggregation. "
            "Recommended until the full run is validated."
        ),
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Fail if any individual sequence cannot be characterized.",
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


def roxy_to_polars(obj: Any) -> pl.DataFrame:
    if isinstance(obj, pl.DataFrame):
        return obj

    if isinstance(obj, pd.DataFrame):
        return pl.from_pandas(obj)

    if hasattr(obj, "to_dict"):
        try:
            payload = obj.to_dict(as_series=False)
            return pl.DataFrame(payload)
        except TypeError:
            payload = obj.to_dict(orient="list")
            return pl.DataFrame(payload)

    raise TypeError(
        f"Unsupported Roxy output type: {type(obj)!r}"
    )


def compute_family(
    family: str,
    sequences: list[str],
    ids: list[str],
) -> pl.DataFrame:
    descriptor = DESCRIPTOR_REGISTRY[family]()
    result = descriptor.compute(sequences, ids=ids)
    frame = roxy_to_polars(result)

    if "id" not in frame.columns:
        raise RuntimeError(
            f"Roxy descriptor '{family}' did not return an 'id' column."
        )

    expected = DESCRIPTOR_FAMILIES[family]
    missing = [col for col in expected if col not in frame.columns]
    if missing:
        raise RuntimeError(
            f"Descriptor '{family}' missing selected columns: {missing}"
        )

    return frame.select(["id", *expected])


def compute_batch(
    sequences: list[str],
    ids: list[str],
) -> pl.DataFrame:
    frames = []

    for family in DESCRIPTOR_FAMILIES:
        frame = compute_family(family, sequences, ids)
        frames.append(frame)

    merged = frames[0]
    for frame in frames[1:]:
        merged = merged.join(
            frame,
            on="id",
            how="inner",
            validate="1:1",
        )

    if merged.height != len(ids):
        raise RuntimeError(
            f"Roxy batch returned {merged.height:,} rows for "
            f"{len(ids):,} input sequences."
        )

    missing_features = [
        col for col in SELECTED_FEATURES if col not in merged.columns
    ]
    if missing_features:
        raise RuntimeError(
            f"Selected feature columns missing after merge: {missing_features}"
        )

    return merged.select(["id", *SELECTED_FEATURES])


def compute_with_fallback(
    sequences: list[str],
    ids: list[str],
    failures: list[dict[str, str]],
) -> pl.DataFrame:
    """
    Compute a batch. If Roxy fails on a heterogeneous batch, recursively bisect
    until the problematic individual sequence is isolated.

    Failed singleton rows are retained with null descriptors so sequence count
    remains auditable.
    """
    try:
        return compute_batch(sequences, ids)
    except Exception as exc:
        if len(ids) == 1:
            failures.append(
                {
                    "sequence_uid": ids[0],
                    "error": repr(exc),
                }
            )
            row = {"id": ids[0]}
            row.update({feature: None for feature in SELECTED_FEATURES})
            return pl.DataFrame([row])

        midpoint = len(ids) // 2

        left = compute_with_fallback(
            sequences[:midpoint],
            ids[:midpoint],
            failures,
        )
        right = compute_with_fallback(
            sequences[midpoint:],
            ids[midpoint:],
            failures,
        )
        return pl.concat([left, right], how="vertical_relaxed")


def validate_chunk(path: Path, expected_rows: int) -> bool:
    if not path.is_file():
        return False

    try:
        schema = pl.read_parquet_schema(path)
        required = {"sequence_uid", *SELECTED_FEATURES}
        if not required.issubset(set(schema)):
            return False

        n_rows = (
            pl.scan_parquet(path)
            .select(pl.len().alias("n"))
            .collect()
            .item()
        )
        return int(n_rows) == int(expected_rows)
    except Exception:
        return False


def write_descriptor_schema(path: Path) -> None:
    schema = pd.DataFrame(SCHEMA_ROWS)
    schema.insert(0, "feature_order", range(1, len(schema) + 1))
    schema.to_csv(path, sep="\t", index=False)


def main() -> None:
    args = parse_args()

    if len(SELECTED_FEATURES) != EXPECTED_DESCRIPTOR_COUNT:
        raise RuntimeError(
            f"Internal feature panel has {len(SELECTED_FEATURES)} columns; "
            f"expected {EXPECTED_DESCRIPTOR_COUNT}."
        )

    input_path = Path(args.input).expanduser().resolve()
    output_dir = Path(args.output_dir).expanduser().resolve()
    chunk_dir = output_dir / "chunks"

    output_dir.mkdir(parents=True, exist_ok=True)
    chunk_dir.mkdir(parents=True, exist_ok=True)

    if not input_path.is_file():
        raise FileNotFoundError(input_path)

    print("=" * 100)
    print("PhageLysClass — Stage 10B compact Roxy characterization")
    print("=" * 100)
    print(f"Input       : {input_path}")
    print(f"Output dir  : {output_dir}")
    print(f"Chunk size  : {args.chunk_size:,}")
    print(f"Descriptors : {len(SELECTED_FEATURES)}")
    print(f"Roxy        : {package_version('roxy')}")
    print()

    # The full 442k-row index is small enough to keep in memory while all heavy
    # numerical computation remains chunked.
    universe = pd.read_csv(
        input_path,
        sep="\t",
        compression="infer",
        dtype={
            "sequence_uid": str,
            "sequence": str,
        },
        low_memory=False,
    ).fillna("")

    required = {"sequence_uid", "sequence", "length"}
    missing = required - set(universe.columns)
    if missing:
        raise KeyError(
            f"Sequence universe missing required columns: {sorted(missing)}"
        )

    if universe["sequence_uid"].duplicated().any():
        raise RuntimeError("Duplicate sequence_uid values in Stage 10A input.")

    n_total = len(universe)
    n_chunks = math.ceil(n_total / args.chunk_size)

    if args.max_chunks is not None:
        n_chunks_to_run = min(n_chunks, args.max_chunks)
    else:
        n_chunks_to_run = n_chunks

    print(f"Unique sequences : {n_total:,}")
    print(f"Total chunks     : {n_chunks:,}")
    if args.max_chunks is not None:
        print(f"Debug run chunks : {n_chunks_to_run:,}")
    print()

    descriptor_schema_path = output_dir / "descriptor_schema.tsv"
    write_descriptor_schema(descriptor_schema_path)

    all_failures: list[dict[str, str]] = []
    status_rows = []

    for chunk_idx in range(n_chunks_to_run):
        start = chunk_idx * args.chunk_size
        stop = min(start + args.chunk_size, n_total)

        chunk = universe.iloc[start:stop].copy()
        expected_rows = len(chunk)

        chunk_path = chunk_dir / f"part_{chunk_idx:06d}.parquet"

        if (
            not args.no_resume
            and validate_chunk(chunk_path, expected_rows)
        ):
            print(
                f"[{chunk_idx + 1:>5}/{n_chunks_to_run:<5}] "
                f"reuse {chunk_path.name} ({expected_rows:,})"
            )
            status_rows.append(
                {
                    "chunk": chunk_idx,
                    "start_row": start,
                    "stop_row_exclusive": stop,
                    "n_sequences": expected_rows,
                    "status": "reused",
                    "n_failures": 0,
                    "path": str(chunk_path),
                }
            )
            continue

        ids = chunk["sequence_uid"].astype(str).tolist()
        sequences = (
            chunk["sequence"]
            .astype(str)
            .str.replace(r"\s+", "", regex=True)
            .str.upper()
            .tolist()
        )

        chunk_failures: list[dict[str, str]] = []

        print(
            f"[{chunk_idx + 1:>5}/{n_chunks_to_run:<5}] "
            f"compute rows {start:,}:{stop:,} "
            f"({expected_rows:,} sequences)"
        )

        features = compute_with_fallback(
            sequences,
            ids,
            chunk_failures,
        )

        metadata = pl.DataFrame(
            {
                "sequence_uid": ids,
                "length": chunk["length"].astype(int).tolist(),
                "is_canonical_20aa": (
                    chunk["is_canonical_20aa"].tolist()
                    if "is_canonical_20aa" in chunk.columns
                    else [None] * expected_rows
                ),
                "noncanonical_symbols": (
                    chunk["noncanonical_symbols"].astype(str).tolist()
                    if "noncanonical_symbols" in chunk.columns
                    else [""] * expected_rows
                ),
            }
        )

        features = features.rename({"id": "sequence_uid"})
        out = metadata.join(
            features,
            on="sequence_uid",
            how="left",
            validate="1:1",
        )

        if out.height != expected_rows:
            raise RuntimeError(
                f"Chunk {chunk_idx}: output rows {out.height:,} != "
                f"input rows {expected_rows:,}"
            )

        tmp_path = chunk_path.with_suffix(".tmp.parquet")
        out.write_parquet(
            tmp_path,
            compression="zstd",
            statistics=True,
        )
        tmp_path.replace(chunk_path)

        all_failures.extend(chunk_failures)

        status_rows.append(
            {
                "chunk": chunk_idx,
                "start_row": start,
                "stop_row_exclusive": stop,
                "n_sequences": expected_rows,
                "status": "computed",
                "n_failures": len(chunk_failures),
                "path": str(chunk_path),
            }
        )

        if chunk_failures:
            print(
                f"    WARN: {len(chunk_failures):,} sequence(s) "
                "required singleton fallback and failed."
            )

    status_df = pd.DataFrame(status_rows)
    status_path = output_dir / "chunk_status.tsv"
    status_df.to_csv(status_path, sep="\t", index=False)

    failures_path = output_dir / "failed_sequences.tsv"
    pd.DataFrame(
        all_failures,
        columns=["sequence_uid", "error"],
    ).to_csv(failures_path, sep="\t", index=False)

    # A --max-chunks run is deliberately partial and should not materialize the
    # final dataset as if it were complete.
    complete_run = n_chunks_to_run == n_chunks

    final_path = output_dir / "characterization_48.parquet"

    if complete_run:
        expected_parts = [
            chunk_dir / f"part_{i:06d}.parquet"
            for i in range(n_chunks)
        ]
        missing_parts = [p for p in expected_parts if not p.is_file()]
        if missing_parts:
            raise RuntimeError(
                f"Cannot finalize: {len(missing_parts):,} chunk files missing."
            )

        print()
        print("[finalize] Concatenating chunk Parquet files...")

        # Polars lazy streaming avoids loading the full 442k × 48 matrix into
        # Python object memory during finalization.
        lazy = pl.scan_parquet(
            str(chunk_dir / "part_*.parquet")
        )

        try:
            lazy.sink_parquet(
                final_path,
                compression="zstd",
                statistics=True,
                maintain_order=True,
            )
        except TypeError:
            # Compatibility with Polars versions lacking maintain_order here.
            lazy.sink_parquet(
                final_path,
                compression="zstd",
                statistics=True,
            )

        final_rows = (
            pl.scan_parquet(final_path)
            .select(pl.len().alias("n"))
            .collect()
            .item()
        )

        if int(final_rows) != n_total:
            raise RuntimeError(
                f"Final output rows {final_rows:,} != input rows {n_total:,}."
            )

        final_null_counts = (
            pl.scan_parquet(final_path)
            .select(
                [
                    pl.col(col).null_count().alias(col)
                    for col in SELECTED_FEATURES
                ]
            )
            .collect()
            .row(0, named=True)
        )

        features_with_nulls = {
            key: int(value)
            for key, value in final_null_counts.items()
            if int(value) > 0
        }
    else:
        final_rows = None
        features_with_nulls = {}

    manifest = {
        "stage": "10B",
        "name": "compact_roxy_characterization",
        "input": str(input_path),
        "input_sha256": sha256_file(input_path),
        "roxy_version": package_version("roxy"),
        "polars_version": package_version("polars"),
        "descriptor_count": len(SELECTED_FEATURES),
        "selected_features": SELECTED_FEATURES,
        "descriptor_families": {
            key: list(value)
            for key, value in DESCRIPTOR_FAMILIES.items()
        },
        "categories": {
            "amino_acid_frequency": len(AAC_FEATURES),
            "physicochemical_thermodynamic_proxy": (
                len(GLOBAL_BASIC_FEATURES) + len(CHARGE_FEATURES)
            ),
            "order_structural_propensity": (
                len(ORDER_DISORDER_FEATURES) + len(STRUCTURAL_FEATURES)
            ),
            "aaindex_residue_property": len(AAINDEX_FEATURES),
        },
        "counts": {
            "input_sequences": n_total,
            "chunk_size": args.chunk_size,
            "total_chunks": n_chunks,
            "chunks_executed_or_reused": n_chunks_to_run,
            "failed_sequences_this_run": len(all_failures),
            "final_rows": int(final_rows) if final_rows is not None else None,
        },
        "validation": {
            "complete_run": complete_run,
            "features_with_null_values": features_with_nulls,
        },
        "design_decisions": {
            "computed_once_per_exact_sequence": True,
            "homology_regime_specific": False,
            "task_specific": False,
            "descriptor_cap": 50,
            "selected_descriptor_count": len(SELECTED_FEATURES),
            "true_evolutionary_profile_features_included": False,
            "evolutionary_note": (
                "The probed Roxy registry provides AAIndex and sequence-order "
                "descriptors but no genuine MSA/PSSM/conservation profile. "
                "AAIndex features are therefore reported as residue-property "
                "scales, not mislabeled as evolutionary evidence."
            ),
        },
        "outputs": {
            "final_characterization": (
                str(final_path) if complete_run else None
            ),
            "chunk_directory": str(chunk_dir),
            "descriptor_schema": str(descriptor_schema_path),
            "chunk_status": str(status_path),
            "failed_sequences": str(failures_path),
        },
    }

    manifest_path = output_dir / "manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2),
        encoding="utf-8",
    )

    print()
    print("=" * 100)
    print("Stage 10B summary")
    print("=" * 100)
    print(f"Selected descriptors : {len(SELECTED_FEATURES)}")
    print(f"Failed sequences     : {len(all_failures):,}")
    print(f"Complete run         : {complete_run}")

    if complete_run:
        print(f"Final rows           : {int(final_rows):,}")
        print(f"Final file           : {final_path}")

    print(f"Schema               : {descriptor_schema_path}")
    print(f"Failures             : {failures_path}")
    print(f"Manifest             : {manifest_path}")

    if args.strict and all_failures:
        raise SystemExit(
            f"STRICT FAIL — {len(all_failures):,} sequences failed Roxy."
        )

    if complete_run:
        print("PASS — Stage 10B Roxy characterization materialized.")
    else:
        print(
            "PARTIAL PASS — requested subset of chunks completed; "
            "final characterization was not materialized."
        )


if __name__ == "__main__":
    main()

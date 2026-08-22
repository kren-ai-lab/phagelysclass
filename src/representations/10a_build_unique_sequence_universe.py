#!/usr/bin/env python3
"""
10a_build_unique_sequence_universe.py

Stage 10A — Build the exact-sequence universe used by downstream
characterization and representation workflows.

Design
------
Representations/descriptors are computed ONCE per exact unique amino-acid
sequence and are later joined back to task-specific datasets by sequence UID.

The default union uses the smallest set of frozen Stage-06 datasets that spans
the complete PhageLysClass paper universe:

    global_mixed
        -> P_core + N_easy + N_hard

    pu_global/broad
        -> P_core + U_explicit + U_other_unresolved

    challenge_sets/uniprot_lytic_related
        -> challenge_lytic_related

    challenge_sets/uniprot_ambiguous_structural
        -> challenge_ambiguous_structural

P_core is shared by global_mixed and PU broad and is deduplicated exactly.

No homology reduction is applied here. H90/H70/H50/H30 are NOT used because
the purpose of this stage is to create the full exact-sequence universe from
which all later task/regime-specific matrices can be indexed.

Expected frozen unique-sequence count:
    442,329

Default outputs
---------------
data/representations/sequence_universe/
├── unique_sequences.tsv.gz
├── sequence_state_summary.tsv
├── source_contributions.tsv
├── validation_summary.tsv
└── manifest.json
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import pandas as pd


EXPECTED_UNIQUE_COUNT = 442_329

MINIMAL_SOURCES = [
    {
        "source_id": "global_mixed",
        "path": "global_mixed/dataset.tsv.gz",
        "role": "P_core + N_easy + N_hard",
    },
    {
        "source_id": "pu_broad",
        "path": "pu_global/broad/dataset.tsv.gz",
        "role": "P_core + U_explicit + U_other_unresolved",
    },
    {
        "source_id": "challenge_lytic_related",
        "path": "challenge_sets/uniprot_lytic_related/dataset.tsv.gz",
        "role": "challenge_lytic_related",
    },
    {
        "source_id": "challenge_ambiguous_structural",
        "path": "challenge_sets/uniprot_ambiguous_structural/dataset.tsv.gz",
        "role": "challenge_ambiguous_structural",
    },
]


CANONICAL_AA = frozenset("ACDEFGHIKLMNPQRSTVWY")


def project_root_from_script(script_file: str) -> Path:
    script_dir = Path(script_file).resolve().parent
    return script_dir.parents[1]


def parse_args() -> argparse.Namespace:
    root = project_root_from_script(__file__)

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--tasks-root",
        default=str(root / "data/tasks"),
        help="Frozen Stage-06 task root.",
    )
    parser.add_argument(
        "--output-dir",
        default=str(root / "data/representations/sequence_universe"),
        help="Output directory.",
    )
    parser.add_argument(
        "--expected-unique-count",
        type=int,
        default=EXPECTED_UNIQUE_COUNT,
        help=(
            "Expected exact-sequence count from the frozen Stage-05 universe. "
            "Use <=0 to disable this count check."
        ),
    )
    parser.add_argument(
        "--audit-all-stage06",
        action="store_true",
        help=(
            "Additionally scan every Stage-06 dataset.tsv.gz and verify that "
            "the minimal source union contains every exact sequence used by "
            "any task/challenge dataset."
        ),
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Convert validation warnings into errors.",
    )
    return parser.parse_args()


def read_dataset(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(path)

    df = pd.read_csv(
        path,
        sep="\t",
        compression="infer",
        dtype=str,
        low_memory=False,
    ).fillna("")

    required = {"entity_id", "sequence"}
    missing = required - set(df.columns)
    if missing:
        raise KeyError(
            f"{path} is missing required columns: {sorted(missing)}"
        )

    df["entity_id"] = df["entity_id"].astype(str).str.strip()
    df["sequence"] = (
        df["sequence"]
        .astype(str)
        .str.replace(r"\s+", "", regex=True)
        .str.upper()
        .str.strip()
    )

    if (df["entity_id"] == "").any():
        raise ValueError(f"Empty entity_id values found in {path}")
    if (df["sequence"] == "").any():
        raise ValueError(f"Empty sequences found in {path}")

    return df


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def sequence_uid(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def canonical_flag(sequence: str) -> bool:
    return set(sequence).issubset(CANONICAL_AA)


def noncanonical_symbols(sequence: str) -> str:
    symbols = sorted(set(sequence) - CANONICAL_AA)
    return "".join(symbols)


def coalesce_consistent(
    block: pd.DataFrame,
    column: str,
) -> tuple[str, bool]:
    """
    Return one non-empty value if all non-empty values agree.
    """
    if column not in block.columns:
        return "", True

    values = sorted(
        {
            str(x).strip()
            for x in block[column].tolist()
            if str(x).strip()
        }
    )
    if len(values) == 0:
        return "", True
    if len(values) == 1:
        return values[0], True
    return "|".join(values), False


def main() -> None:
    args = parse_args()

    tasks_root = Path(args.tasks_root).expanduser().resolve()
    output_dir = Path(args.output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 96)
    print("PhageLysClass — Stage 10A exact-sequence universe")
    print("=" * 96)
    print(f"Stage-06 root : {tasks_root}")
    print(f"Output dir    : {output_dir}")
    print()

    source_frames = []
    source_rows = []
    source_hashes = {}

    for source in MINIMAL_SOURCES:
        path = tasks_root / source["path"]
        print(f"[load] {source['source_id']}: {path}")

        df = read_dataset(path)
        source_hashes[str(path)] = sha256_file(path)

        before_exact_duplicates = len(df)
        within_unique = int(df["sequence"].nunique())

        source_rows.append(
            {
                "source_id": source["source_id"],
                "role": source["role"],
                "path": str(path),
                "rows": len(df),
                "exact_unique_sequences": within_unique,
                "within_source_exact_duplicates": (
                    before_exact_duplicates - within_unique
                ),
            }
        )

        keep_cols = [
            col
            for col in [
                "entity_id",
                "sequence",
                "sequence_state",
                "canonical_target_class",
            ]
            if col in df.columns
        ]

        block = df[keep_cols].copy()
        block["source_id"] = source["source_id"]
        source_frames.append(block)

    combined = pd.concat(source_frames, ignore_index=True)

    print()
    print(f"Combined source rows : {len(combined):,}")
    print(f"Exact unique seqs    : {combined['sequence'].nunique():,}")

    # Stable exact-sequence UID.
    combined["sequence_uid"] = combined["sequence"].map(sequence_uid)

    # A SHA-256 collision is practically impossible, but verify the contract.
    uid_to_seq_n = combined.groupby("sequence_uid")["sequence"].nunique()
    hash_collisions = int((uid_to_seq_n > 1).sum())
    if hash_collisions:
        raise RuntimeError(
            f"Detected {hash_collisions:,} sequence SHA-256 collisions."
        )

    # Build one row per exact sequence.
    records = []
    inconsistent_state_count = 0
    inconsistent_class_count = 0

    for uid, block in combined.groupby("sequence_uid", sort=False):
        seq = block["sequence"].iloc[0]

        sequence_state, state_consistent = coalesce_consistent(
            block, "sequence_state"
        )
        canonical_class, class_consistent = coalesce_consistent(
            block, "canonical_target_class"
        )

        if not state_consistent:
            inconsistent_state_count += 1
        if not class_consistent:
            inconsistent_class_count += 1

        entity_ids = sorted(set(block["entity_id"].astype(str)))
        source_ids = sorted(set(block["source_id"].astype(str)))

        records.append(
            {
                "sequence_uid": uid,
                "sequence": seq,
                "length": len(seq),
                "is_canonical_20aa": canonical_flag(seq),
                "noncanonical_symbols": noncanonical_symbols(seq),
                "sequence_state": sequence_state,
                "canonical_target_class": canonical_class,
                "representative_entity_id": entity_ids[0],
                "n_entity_ids": len(entity_ids),
                "n_minimal_sources": len(source_ids),
                "minimal_sources": ";".join(source_ids),
                "sequence_state_consistent": state_consistent,
                "canonical_target_class_consistent": class_consistent,
                # Informative eligibility flags only. No sequence is filtered.
                "length_le_1024": len(seq) <= 1024,
                "length_le_2048": len(seq) <= 2048,
            }
        )

    universe = pd.DataFrame(records)

    # Stable output order independent of input row order.
    universe = universe.sort_values("sequence_uid").reset_index(drop=True)

    unique_count = len(universe)
    duplicate_rows_collapsed = len(combined) - unique_count

    validation_rows: list[dict[str, Any]] = [
        {
            "check": "exact_unique_sequences",
            "value": unique_count,
            "expected": (
                args.expected_unique_count
                if args.expected_unique_count > 0
                else ""
            ),
            "status": (
                "PASS"
                if args.expected_unique_count <= 0
                or unique_count == args.expected_unique_count
                else "WARN"
            ),
        },
        {
            "check": "sha256_sequence_collisions",
            "value": hash_collisions,
            "expected": 0,
            "status": "PASS" if hash_collisions == 0 else "FAIL",
        },
        {
            "check": "inconsistent_sequence_state",
            "value": inconsistent_state_count,
            "expected": 0,
            "status": (
                "PASS" if inconsistent_state_count == 0 else "WARN"
            ),
        },
        {
            "check": "inconsistent_canonical_target_class",
            "value": inconsistent_class_count,
            "expected": 0,
            "status": (
                "PASS" if inconsistent_class_count == 0 else "WARN"
            ),
        },
        {
            "check": "empty_sequences",
            "value": int((universe["sequence"] == "").sum()),
            "expected": 0,
            "status": "PASS",
        },
    ]

    # Optional exhaustive Stage-06 audit.
    audit_unique_count = None
    missing_from_minimal = None

    if args.audit_all_stage06:
        print()
        print("[audit] Scanning all Stage-06 dataset.tsv.gz files...")
        all_paths = sorted(tasks_root.rglob("dataset.tsv.gz"))
        all_sequences = set()

        for path in all_paths:
            df = read_dataset(path)
            all_sequences.update(df["sequence"].astype(str))

        minimal_sequences = set(universe["sequence"].astype(str))
        missing = all_sequences - minimal_sequences

        audit_unique_count = len(all_sequences)
        missing_from_minimal = len(missing)

        validation_rows.append(
            {
                "check": "all_stage06_unique_sequences",
                "value": audit_unique_count,
                "expected": unique_count,
                "status": (
                    "PASS" if audit_unique_count == unique_count else "WARN"
                ),
            }
        )
        validation_rows.append(
            {
                "check": "stage06_sequences_missing_from_minimal_union",
                "value": missing_from_minimal,
                "expected": 0,
                "status": "PASS" if missing_from_minimal == 0 else "FAIL",
            }
        )

    # Summary tables.
    state_summary = (
        universe.groupby("sequence_state", dropna=False)
        .size()
        .reset_index(name="count")
        .sort_values("count", ascending=False)
    )
    state_summary["fraction"] = state_summary["count"] / unique_count

    source_df = pd.DataFrame(source_rows)
    validation_df = pd.DataFrame(validation_rows)

    # Additional high-level sequence characteristics useful for downstream
    # computational planning.
    sequence_characteristics = {
        "n_unique_sequences": unique_count,
        "combined_source_rows": len(combined),
        "duplicate_source_rows_collapsed": duplicate_rows_collapsed,
        "canonical_20aa": int(universe["is_canonical_20aa"].sum()),
        "noncanonical": int((~universe["is_canonical_20aa"]).sum()),
        "length_le_1024": int(universe["length_le_1024"].sum()),
        "length_gt_1024": int((~universe["length_le_1024"]).sum()),
        "length_le_2048": int(universe["length_le_2048"].sum()),
        "length_gt_2048": int((~universe["length_le_2048"]).sum()),
        "min_length": int(universe["length"].min()),
        "median_length": float(universe["length"].median()),
        "mean_length": float(universe["length"].mean()),
        "max_length": int(universe["length"].max()),
    }

    universe_path = output_dir / "unique_sequences.tsv.gz"
    state_path = output_dir / "sequence_state_summary.tsv"
    source_path = output_dir / "source_contributions.tsv"
    validation_path = output_dir / "validation_summary.tsv"
    manifest_path = output_dir / "manifest.json"

    universe.to_csv(
        universe_path,
        sep="\t",
        index=False,
        compression="gzip",
    )
    state_summary.to_csv(state_path, sep="\t", index=False)
    source_df.to_csv(source_path, sep="\t", index=False)
    validation_df.to_csv(validation_path, sep="\t", index=False)

    manifest = {
        "stage": "10A",
        "name": "exact_sequence_representation_universe",
        "purpose": (
            "Compute descriptors/encodings/embeddings once per exact unique "
            "sequence and later join them to task/regime datasets."
        ),
        "input_root": str(tasks_root),
        "minimal_sources": MINIMAL_SOURCES,
        "input_sha256": source_hashes,
        "expected_unique_count": (
            args.expected_unique_count
            if args.expected_unique_count > 0
            else None
        ),
        "counts": sequence_characteristics,
        "validation": validation_rows,
        "audit_all_stage06": {
            "enabled": bool(args.audit_all_stage06),
            "all_stage06_unique_sequences": audit_unique_count,
            "missing_from_minimal_union": missing_from_minimal,
        },
        "decisions": {
            "deduplication_key": "exact amino-acid sequence",
            "stable_sequence_id": "SHA-256(sequence)",
            "homology_reduction_applied": False,
            "task_specific_redundancy_regimes_used": False,
            "sequences_filtered_by_length": False,
            "sequences_filtered_by_alphabet": False,
        },
        "outputs": {
            "unique_sequences": str(universe_path),
            "sequence_state_summary": str(state_path),
            "source_contributions": str(source_path),
            "validation_summary": str(validation_path),
        },
    }

    manifest_path.write_text(
        json.dumps(manifest, indent=2),
        encoding="utf-8",
    )

    print()
    print("Sequence universe")
    print("-" * 96)
    for key, value in sequence_characteristics.items():
        if isinstance(value, float):
            print(f"{key:32s}: {value:.2f}")
        else:
            print(f"{key:32s}: {value:,}")

    print()
    print("Validation")
    print("-" * 96)
    print(validation_df.to_string(index=False))

    statuses = set(validation_df["status"].astype(str))
    if "FAIL" in statuses:
        raise SystemExit("FAIL — sequence universe validation failed.")
    if args.strict and "WARN" in statuses:
        raise SystemExit(
            "STRICT FAIL — one or more sequence-universe checks returned WARN."
        )

    print()
    print(f"Unique sequences : {universe_path}")
    print(f"State summary    : {state_path}")
    print(f"Source summary   : {source_path}")
    print(f"Validation       : {validation_path}")
    print(f"Manifest         : {manifest_path}")
    print()
    print("PASS — Stage 10A exact-sequence universe materialized.")


if __name__ == "__main__":
    main()

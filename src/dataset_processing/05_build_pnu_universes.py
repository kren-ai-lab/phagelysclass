#!/usr/bin/env python3
"""
05_build_pnu_universes.py

Freeze exact-sequence-level P/N/U candidate universes for PhageLysClass.

Inputs
------
1. PhageLysData Core raw extraction from script 01.
2. UniProt phage background classified by script 03.

Principles
----------
- P comes only from PhageLysData Core.
- UniProt exact matches to Core are decontamination overlaps, not new P records.
- UniProt accessions are collapsed to exact amino-acid sequence entities here.
- No homology threshold is applied at this stage.
- Explicit target-like or ambiguous annotations block a sequence from N/U.
- N_easy and N_hard remain candidate negative strata.
- U is unlabeled, never treated as negative ground truth.
- Task construction, balancing, and train/validation/test splitting happen later.

Sequence-level precedence for duplicated UniProt accessions
-----------------------------------------------------------
known_p_overlap
    > challenge_lytic_related
    > challenge_ambiguous_structural
    > N_hard
    > N_easy
    > U_explicit
    > U_other_unresolved

This precedence is conservative: a target-like or ambiguous annotation prevents
a duplicate sequence from entering a negative pool.

Outputs
-------
data/processed/pnu_universes/
├── pnu_sequence_master.tsv.gz
├── p_core.tsv.gz
├── n_easy.tsv.gz
├── n_hard.tsv.gz
├── u_explicit.tsv.gz
├── u_other_unresolved.tsv.gz
├── challenge_lytic_related.tsv.gz
├── challenge_ambiguous_structural.tsv.gz
├── background_known_p_overlap.tsv.gz
├── pnu_universe_summary.tsv
├── phagelysdata_overlap_summary.tsv
└── pnu_manifest.json
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import pandas as pd


PNU_VERSION = "v1"

DEFAULT_CORE = "data/raw/phagelysdata/phagelysdata_core_raw.tsv.gz"
DEFAULT_BACKGROUND = (
    "data/processed/uniprot_background/"
    "uniprot_background_classified.tsv.gz"
)

UNIPROT_USECOLS = [
    "Entry",
    "Protein names",
    "Reviewed",
    "Protein existence",
    "Annotation",
    "Sequence",
    "retrieval_scope",
    "sequence_norm",
    "sequence_sha256",
    "evidence_state",
    "primary_rule_id",
    "matched_rule_ids",
    "matched_rule_groups",
    "has_rule_conflict",
    "phagelysdata_core_match",
]

STATE_PRECEDENCE = [
    "known_p_overlap",
    "lytic_related",
    "ambiguous_structural",
    "n_candidate_hard",
    "n_candidate_easy",
    "u_uncharacterized",
    "u_domain_unknown",
    "u_low_information",
    "u_other_unresolved",
]

U_EXPLICIT_STATES = {
    "u_uncharacterized",
    "u_domain_unknown",
    "u_low_information",
}


def parse_args() -> argparse.Namespace:
    script_dir = Path(__file__).resolve().parent
    project_root = script_dir.parents[1]

    parser = argparse.ArgumentParser(
        description=(
            "Build exact-sequence-level P/N/U candidate universes from "
            "PhageLysData Core and the classified UniProt phage background."
        )
    )
    parser.add_argument(
        "--phagelysdata-core",
        default=str(project_root / DEFAULT_CORE),
        help="PhageLysData Core raw TSV.GZ from script 01.",
    )
    parser.add_argument(
        "--classified-background",
        default=str(project_root / DEFAULT_BACKGROUND),
        help="Classified UniProt background TSV.GZ from script 03.",
    )
    parser.add_argument(
        "--output-dir",
        default=str(project_root / "data/processed/pnu_universes"),
        help="Output directory.",
    )
    return parser.parse_args()


def normalize_sequence(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return re.sub(r"\s+", "", str(value)).upper()


def sequence_sha256(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("utf-8")).hexdigest()


def join_unique(series: pd.Series) -> str:
    values = sorted(
        {
            str(value).strip()
            for value in series
            if pd.notna(value) and str(value).strip()
        }
    )
    return ";".join(values)


def state_set(serialized: str) -> set[str]:
    if not serialized:
        return set()
    return {x for x in serialized.split(";") if x}


def assign_sequence_state(observed_states: str) -> tuple[str, str, str]:
    """
    Return:
        pnu_role
        sequence_state
        stratum
    """
    states = state_set(observed_states)

    if "known_p_overlap" in states:
        return "P_overlap", "known_p_overlap", "known_p_overlap"

    if "lytic_related" in states:
        return "challenge", "challenge_lytic_related", "lytic_related"

    if "ambiguous_structural" in states:
        return (
            "challenge",
            "challenge_ambiguous_structural",
            "ambiguous_structural",
        )

    if "n_candidate_hard" in states:
        return "N", "N_hard", "hard"

    if "n_candidate_easy" in states:
        return "N", "N_easy", "easy"

    if states & U_EXPLICIT_STATES:
        explicit = sorted(states & U_EXPLICIT_STATES)
        return "U", "U_explicit", ";".join(explicit)

    if "u_other_unresolved" in states:
        return "U", "U_other_unresolved", "other_unresolved"

    raise ValueError(
        "Could not assign sequence-level state from observed states: "
        + observed_states
    )


def load_core(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"PhageLysData Core not found: {path}")

    core = pd.read_csv(
        path,
        sep="\t",
        compression="infer",
        dtype=str,
    ).fillna("")

    required = {
        "sequence_entity_id",
        "sequence",
        "canonical_target_class",
        "canonical_class_basis",
        "assigned_tier",
    }
    missing = required - set(core.columns)
    if missing:
        raise KeyError(
            "PhageLysData Core missing columns: "
            + ", ".join(sorted(missing))
        )

    core["sequence"] = core["sequence"].map(normalize_sequence)

    if (core["sequence"] == "").any():
        raise ValueError("PhageLysData Core contains empty sequences.")

    if core["sequence"].duplicated().any():
        raise ValueError(
            "PhageLysData Core is expected to contain one exact-sequence "
            "entity per sequence."
        )

    core["sequence_sha256"] = core["sequence"].map(sequence_sha256)
    core["entity_id"] = core["sequence_entity_id"]
    core["source_universe"] = "phagelysdata_core"
    core["pnu_role"] = "P"
    core["sequence_state"] = "P_core"
    core["stratum"] = core["canonical_target_class"]
    core["n_source_accessions"] = 0
    core["source_accessions"] = ""
    core["protein_names"] = ""
    core["observed_evidence_states"] = ""
    core["observed_primary_rules"] = ""
    core["observed_rule_groups"] = ""
    core["retrieval_scopes"] = ""
    core["any_reviewed"] = ""
    core["max_annotation_score"] = ""
    core["has_cross_accession_state_conflict"] = False

    return core


def load_background(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"Classified background not found: {path}")

    header = pd.read_csv(
        path,
        sep="\t",
        compression="infer",
        nrows=0,
    ).columns.tolist()

    missing = set(UNIPROT_USECOLS) - set(header)
    if missing:
        raise KeyError(
            "Classified UniProt background missing columns: "
            + ", ".join(sorted(missing))
        )

    df = pd.read_csv(
        path,
        sep="\t",
        compression="infer",
        dtype=str,
        usecols=UNIPROT_USECOLS,
        low_memory=False,
    ).fillna("")

    df["sequence_norm"] = df["sequence_norm"].map(normalize_sequence)

    if (df["sequence_norm"] == "").any():
        raise ValueError(
            "Classified UniProt background contains empty normalized sequences."
        )

    return df


def parse_annotation_score(value: object) -> float | None:
    if value is None or pd.isna(value):
        return None

    text = str(value).strip()
    if not text:
        return None

    match = re.search(r"([0-9]+(?:\.[0-9]+)?)", text)
    return float(match.group(1)) if match else None


def collapse_background(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()

    df["annotation_score_numeric"] = df["Annotation"].map(
        parse_annotation_score
    )
    df["reviewed_bool"] = (
        df["Reviewed"]
        .astype(str)
        .str.strip()
        .str.lower()
        .eq("reviewed")
    )

    grouped = (
        df.groupby("sequence_norm", sort=False)
        .agg(
            sequence_sha256=("sequence_sha256", "first"),
            n_source_accessions=("Entry", "nunique"),
            source_accessions=("Entry", join_unique),
            protein_names=("Protein names", join_unique),
            observed_evidence_states=("evidence_state", join_unique),
            observed_primary_rules=("primary_rule_id", join_unique),
            observed_rule_groups=("matched_rule_groups", join_unique),
            retrieval_scopes=("retrieval_scope", join_unique),
            phagelysdata_core_matches=("phagelysdata_core_match", join_unique),
            any_reviewed=("reviewed_bool", "max"),
            max_annotation_score=("annotation_score_numeric", "max"),
        )
        .reset_index()
        .rename(columns={"sequence_norm": "sequence"})
    )

    assignments = grouped["observed_evidence_states"].map(
        assign_sequence_state
    )

    grouped[
        ["pnu_role", "sequence_state", "stratum"]
    ] = pd.DataFrame(assignments.tolist(), index=grouped.index)

    grouped["has_cross_accession_state_conflict"] = (
        grouped["observed_evidence_states"]
        .map(lambda x: len(state_set(x)) > 1)
    )

    grouped["entity_id"] = (
        "UNIPROTSEQ_" + grouped["sequence_sha256"].str[:24]
    )
    grouped["source_universe"] = "uniprot_background"

    return grouped


def build_master(
    core: pd.DataFrame,
    background_entities: pd.DataFrame,
) -> pd.DataFrame:
    # Known P overlaps are represented by the PhageLysData Core entity only.
    background_nonoverlap = background_entities.loc[
        background_entities["sequence_state"].ne("known_p_overlap")
    ].copy()

    common_columns = [
        "entity_id",
        "sequence",
        "sequence_sha256",
        "source_universe",
        "pnu_role",
        "sequence_state",
        "stratum",
        "canonical_target_class",
        "canonical_class_basis",
        "assigned_tier",
        "n_source_accessions",
        "source_accessions",
        "protein_names",
        "observed_evidence_states",
        "observed_primary_rules",
        "observed_rule_groups",
        "retrieval_scopes",
        "any_reviewed",
        "max_annotation_score",
        "has_cross_accession_state_conflict",
    ]

    core_part = core.copy()
    background_part = background_nonoverlap.copy()

    for col in common_columns:
        if col not in core_part:
            core_part[col] = ""
        if col not in background_part:
            background_part[col] = ""

    master = pd.concat(
        [
            core_part[common_columns],
            background_part[common_columns],
        ],
        ignore_index=True,
    )

    if master["sequence"].duplicated().any():
        duplicated = int(master["sequence"].duplicated().sum())
        raise RuntimeError(
            f"P/N/U master contains {duplicated:,} duplicated exact sequences."
        )

    return master


def write_subset(
    master: pd.DataFrame,
    path: Path,
    mask: pd.Series,
) -> int:
    subset = master.loc[mask].copy()
    subset.to_csv(
        path,
        sep="\t",
        index=False,
        compression="gzip",
    )
    return len(subset)


def main() -> None:
    args = parse_args()

    core_path = Path(args.phagelysdata_core).expanduser().resolve()
    background_path = Path(
        args.classified_background
    ).expanduser().resolve()
    output_dir = Path(args.output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 92)
    print("PhageLysClass — build exact-sequence P/N/U universes")
    print("=" * 92)
    print(f"PhageLysData Core : {core_path}")
    print(f"UniProt background: {background_path}")
    print(f"Output dir        : {output_dir}")
    print()

    core = load_core(core_path)
    background = load_background(background_path)

    print(f"Core exact-sequence entities        : {len(core):,}")
    print(f"UniProt accession-level records     : {len(background):,}")

    background_entities = collapse_background(background)

    print(
        f"UniProt exact-sequence entities     : "
        f"{len(background_entities):,}"
    )

    overlap_entities = background_entities.loc[
        background_entities["sequence_state"].eq("known_p_overlap")
    ].copy()

    core_sequences = set(core["sequence"])
    overlap_core_sequences = set(overlap_entities["sequence"])
    overlap_core_sequences &= core_sequences

    print(
        f"Unique Core sequences seen in UniProt: "
        f"{len(overlap_core_sequences):,}/{len(core):,} "
        f"({100*len(overlap_core_sequences)/len(core):.2f}%)"
    )

    master = build_master(core, background_entities)

    master_path = output_dir / "pnu_sequence_master.tsv.gz"
    master.to_csv(
        master_path,
        sep="\t",
        index=False,
        compression="gzip",
    )

    subset_counts = {}

    subset_counts["P_core"] = write_subset(
        master,
        output_dir / "p_core.tsv.gz",
        master["pnu_role"].eq("P"),
    )
    subset_counts["N_easy"] = write_subset(
        master,
        output_dir / "n_easy.tsv.gz",
        master["sequence_state"].eq("N_easy"),
    )
    subset_counts["N_hard"] = write_subset(
        master,
        output_dir / "n_hard.tsv.gz",
        master["sequence_state"].eq("N_hard"),
    )
    subset_counts["U_explicit"] = write_subset(
        master,
        output_dir / "u_explicit.tsv.gz",
        master["sequence_state"].eq("U_explicit"),
    )
    subset_counts["U_other_unresolved"] = write_subset(
        master,
        output_dir / "u_other_unresolved.tsv.gz",
        master["sequence_state"].eq("U_other_unresolved"),
    )
    subset_counts["challenge_lytic_related"] = write_subset(
        master,
        output_dir / "challenge_lytic_related.tsv.gz",
        master["sequence_state"].eq("challenge_lytic_related"),
    )
    subset_counts["challenge_ambiguous_structural"] = write_subset(
        master,
        output_dir / "challenge_ambiguous_structural.tsv.gz",
        master["sequence_state"].eq(
            "challenge_ambiguous_structural"
        ),
    )

    overlap_path = output_dir / "background_known_p_overlap.tsv.gz"
    overlap_entities.to_csv(
        overlap_path,
        sep="\t",
        index=False,
        compression="gzip",
    )

    summary_rows = []

    for state, n in (
        master["sequence_state"]
        .value_counts()
        .sort_index()
        .items()
    ):
        summary_rows.append(
            {
                "section": "sequence_state",
                "name": state,
                "count": int(n),
                "fraction_of_master": n / len(master),
            }
        )

    for role, n in (
        master["pnu_role"]
        .value_counts()
        .sort_index()
        .items()
    ):
        summary_rows.append(
            {
                "section": "pnu_role",
                "name": role,
                "count": int(n),
                "fraction_of_master": n / len(master),
            }
        )

    summary_rows.extend(
        [
            {
                "section": "audit",
                "name": "uniprot_accession_rows",
                "count": int(len(background)),
                "fraction_of_master": "",
            },
            {
                "section": "audit",
                "name": "uniprot_unique_exact_sequences",
                "count": int(len(background_entities)),
                "fraction_of_master": "",
            },
            {
                "section": "audit",
                "name": "background_exact_sequences_overlapping_core",
                "count": int(len(overlap_entities)),
                "fraction_of_master": "",
            },
            {
                "section": "audit",
                "name": "unique_phagelysdata_core_sequences_seen_in_uniprot",
                "count": int(len(overlap_core_sequences)),
                "fraction_of_master": len(overlap_core_sequences) / len(core),
            },
        ]
    )

    summary = pd.DataFrame(summary_rows)
    summary_path = output_dir / "pnu_universe_summary.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)

    # Correct overlap accounting by unique PhageLysData sequence entities.
    core_overlap = core.copy()
    core_overlap["present_in_uniprot_background"] = core_overlap[
        "sequence"
    ].isin(overlap_core_sequences)

    overlap_summary = (
        core_overlap.groupby(
            ["canonical_target_class", "present_in_uniprot_background"],
            dropna=False,
        )
        .size()
        .rename("count")
        .reset_index()
    )

    totals = (
        core_overlap.groupby("canonical_target_class")
        .size()
        .rename("class_total")
        .reset_index()
    )

    overlap_summary = overlap_summary.merge(
        totals,
        on="canonical_target_class",
        how="left",
    )
    overlap_summary["fraction_within_class"] = (
        overlap_summary["count"] / overlap_summary["class_total"]
    )

    overlap_summary_path = (
        output_dir / "phagelysdata_overlap_summary.tsv"
    )
    overlap_summary.to_csv(
        overlap_summary_path,
        sep="\t",
        index=False,
    )

    manifest = {
        "pnu_version": PNU_VERSION,
        "phagelysdata_core": str(core_path),
        "classified_uniprot_background": str(background_path),
        "n_core_entities": int(len(core)),
        "n_uniprot_accession_rows": int(len(background)),
        "n_uniprot_exact_sequence_entities": int(
            len(background_entities)
        ),
        "n_unique_core_sequences_seen_in_uniprot": int(
            len(overlap_core_sequences)
        ),
        "n_master_exact_sequence_entities": int(len(master)),
        "state_precedence": STATE_PRECEDENCE,
        "u_explicit_states": sorted(U_EXPLICIT_STATES),
        "subset_counts": subset_counts,
        "design_note": (
            "Exact-sequence P/N/U universe freeze only. No homology filtering, "
            "balancing, train/test splitting, or task-specific labels applied."
        ),
    }

    manifest_path = output_dir / "pnu_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2),
        encoding="utf-8",
    )

    print()
    print("Sequence-level universe counts")
    print("-" * 92)

    for state, n in (
        master["sequence_state"]
        .value_counts()
        .sort_index()
        .items()
    ):
        print(f"{state:34s} {n:>10,}")

    print()
    print(f"Master                 : {master_path}")
    print(f"Summary                : {summary_path}")
    print(f"Core overlap summary   : {overlap_summary_path}")
    print(f"Manifest               : {manifest_path}")
    print()
    print(
        "PASS — exact-sequence P/N/U universes frozen; no homology, "
        "balancing, splitting, or task construction applied."
    )


if __name__ == "__main__":
    main()

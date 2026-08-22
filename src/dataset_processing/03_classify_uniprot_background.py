#!/usr/bin/env python3
"""
03_classify_uniprot_background.py

Classify the broad raw UniProtKB phage-protein background into evidence-state
candidate pools for PhageLysClass.

Rule-set note (v2): lytic identity is not inferred from free-text Function [CC]
mentions, and generic protease/peptidase terms are not clean-negative evidence.

IMPORTANT
---------
This script does NOT create final ML datasets and does NOT declare ground-truth
negatives. It creates auditable evidence-state candidates:

    known_p_overlap
        Exact amino-acid sequence already present in PhageLysData Core.

    lytic_related
        UniProt annotation contains direct/near-direct lytic terminology.

    ambiguous_structural
        Tail/adsorption-associated proteins that may include VAL/DPO functions.

    n_candidate_hard
        Lysis-adjacent non-enzymatic machinery (e.g. holin/spanin).

    n_candidate_easy
        Explicitly annotated non-lytic phage machinery/functions.

    u_uncharacterized
        Explicit uncharacterized/hypothetical/unknown naming.

    u_domain_unknown
        DUF/UPF-style unresolved domain/family naming.

    u_low_information
        Generic ORF/predicted/putative low-information naming.

    u_other_unresolved
        Remaining records not safely assigned to another state.

The rule registry is kept outside the code:
    config/rules/uniprot_background_evidence_rules.tsv

Precedence
----------
1. exact sequence overlap with PhageLysData Core
2. highest-priority enabled local evidence rule
3. u_other_unresolved

Rule design guardrails
----------------------
- Free-text Function [CC] may support context, but it is not used by the direct
  lytic rules to infer protein identity.
- Generic "protease" / "peptidase" terminology is not sufficient evidence for a
  clean negative.

No sequence is discarded because of homology. Exact identity to PhageLysData
Core is marked only to prevent known positives from later entering U/N pools.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable

import pandas as pd


DEFAULT_BACKGROUND = "data/raw/uniprot_background/uniprot_background_raw.tsv.gz"
DEFAULT_PHAGELYSDATA = "data/raw/phagelysdata/phagelysdata_core_raw.tsv.gz"

REQUIRED_BACKGROUND_COLUMNS = {
    "Entry",
    "Sequence",
}

REQUIRED_RULE_COLUMNS = {
    "rule_id",
    "evidence_state",
    "rule_group",
    "priority",
    "fields",
    "pattern",
    "rationale",
    "enabled",
}

FIELD_NAME_MAP = {
    # stable internal name -> possible UniProt TSV header
    "accession": "Entry",
    "entry_name": "Entry Name",
    "protein_name": "Protein names",
    "gene_names": "Gene Names",
    "organism_name": "Organism",
    "organism_id": "Organism (ID)",
    "lineage": "Taxonomic lineage",
    "lineage_ids": "Taxonomic lineage (Ids)",
    "virus_hosts": "Virus hosts",
    "reviewed": "Reviewed",
    "protein_existence": "Protein existence",
    "annotation_score": "Annotation",
    "fragment": "Fragment",
    "length": "Length",
    "sequence": "Sequence",
    "cc_function": "Function [CC]",
    "cc_catalytic_activity": "Catalytic activity",
    "ec": "EC number",
    "go_id": "Gene Ontology IDs",
    "keyword": "Keywords",
    "date_created": "Date of creation",
    "date_modified": "Date of last modification",
    "retrieval_query_id": "retrieval_query_id",
    "retrieval_scope": "retrieval_scope",
}

CANONICAL_AA = set("ACDEFGHIKLMNPQRSTVWY")


def parse_args() -> argparse.Namespace:
    script_dir = Path(__file__).resolve().parent
    project_root = script_dir.parents[1]

    parser = argparse.ArgumentParser(
        description=(
            "Classify the raw UniProt phage background into auditable "
            "evidence-state candidate pools."
        )
    )
    parser.add_argument(
        "--background",
        default=str(project_root / DEFAULT_BACKGROUND),
        help="Raw UniProt background TSV.GZ from script 02.",
    )
    parser.add_argument(
        "--phagelysdata-core",
        default=str(project_root / DEFAULT_PHAGELYSDATA),
        help="Raw PhageLysData Core TSV.GZ from script 01.",
    )
    parser.add_argument(
        "--rules",
        default=str(
            project_root
            / "config"
            / "rules"
            / "uniprot_background_evidence_rules.tsv"
        ),
        help="Evidence-rule registry.",
    )
    parser.add_argument(
        "--output-dir",
        default=str(project_root / "data/processed/uniprot_background"),
        help="Output directory.",
    )
    parser.add_argument(
        "--chunksize",
        type=int,
        default=50000,
        help="Rows processed per pandas chunk (default: 50,000).",
    )
    parser.add_argument(
        "--write-state-files",
        action="store_true",
        help=(
            "Also write one TSV.GZ per evidence state. The master classified "
            "file is always written."
        ),
    )
    return parser.parse_args()


def truthy(value: object) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def normalize_sequence(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return re.sub(r"\s+", "", str(value)).upper()


def sequence_sha256(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("utf-8")).hexdigest()


def load_positive_sequences(path: Path) -> tuple[set[str], dict[str, str]]:
    if not path.is_file():
        raise FileNotFoundError(f"PhageLysData Core not found: {path}")

    core = pd.read_csv(
        path,
        sep="\t",
        compression="infer",
        dtype=str,
        usecols=["sequence_entity_id", "sequence", "canonical_target_class"],
    ).fillna("")

    core["sequence_norm"] = core["sequence"].map(normalize_sequence)

    if (core["sequence_norm"] == "").any():
        raise ValueError("PhageLysData Core contains empty sequences.")

    if core["sequence_norm"].duplicated().any():
        # Exact-sequence entities should be unique. Fail instead of silently
        # choosing one label.
        duplicates = int(core["sequence_norm"].duplicated().sum())
        raise ValueError(
            f"PhageLysData Core contains {duplicates:,} duplicated exact sequences."
        )

    sequence_set = set(core["sequence_norm"])

    # Map exact sequence -> compact provenance string.
    mapping = {
        row.sequence_norm: (
            f"{row.sequence_entity_id}|{row.canonical_target_class}"
        )
        for row in core.itertuples(index=False)
    }

    return sequence_set, mapping


def load_rules(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"Evidence-rule registry not found: {path}")

    rules = pd.read_csv(path, sep="\t", dtype=str).fillna("")

    missing = REQUIRED_RULE_COLUMNS - set(rules.columns)
    if missing:
        raise KeyError(
            "Rule registry missing columns: " + ", ".join(sorted(missing))
        )

    rules = rules.loc[rules["enabled"].map(truthy)].copy()

    if rules.empty:
        raise ValueError("No enabled evidence rules.")

    if rules["rule_id"].duplicated().any():
        duplicated = sorted(
            rules.loc[rules["rule_id"].duplicated(), "rule_id"].unique()
        )
        raise ValueError(
            "Duplicate rule_id values: " + ", ".join(duplicated)
        )

    rules["priority"] = pd.to_numeric(rules["priority"], errors="raise").astype(int)
    rules = rules.sort_values(
        ["priority", "rule_id"],
        ascending=[False, True],
        kind="stable",
    ).reset_index(drop=True)

    compiled = []
    for row in rules.itertuples(index=False):
        field_keys = [x.strip() for x in row.fields.split(";") if x.strip()]
        unknown_fields = sorted(set(field_keys) - set(FIELD_NAME_MAP))
        if unknown_fields:
            raise ValueError(
                f"Rule {row.rule_id} references unknown fields: "
                + ", ".join(unknown_fields)
            )

        try:
            regex = re.compile(row.pattern, flags=re.IGNORECASE)
        except re.error as exc:
            raise ValueError(
                f"Invalid regex in rule {row.rule_id}: {exc}"
            ) from exc

        compiled.append(
            {
                "rule_id": row.rule_id,
                "evidence_state": row.evidence_state,
                "rule_group": row.rule_group,
                "priority": int(row.priority),
                "field_keys": field_keys,
                "regex": regex,
                "pattern": row.pattern,
                "rationale": row.rationale,
            }
        )

    return rules, compiled


def validate_background_header(path: Path) -> list[str]:
    if not path.is_file():
        raise FileNotFoundError(f"Raw UniProt background not found: {path}")

    header = pd.read_csv(
        path,
        sep="\t",
        compression="infer",
        nrows=0,
    ).columns.tolist()

    missing = REQUIRED_BACKGROUND_COLUMNS - set(header)
    if missing:
        raise KeyError(
            "Raw UniProt background missing required columns: "
            + ", ".join(sorted(missing))
        )

    return header


def field_value(record: pd.Series, field_key: str) -> str:
    column = FIELD_NAME_MAP[field_key]
    value = record.get(column, "")
    if value is None or pd.isna(value):
        return ""
    return str(value)


def evaluate_rules(
    record: pd.Series,
    compiled_rules: list[dict],
) -> tuple[str, str, str, str]:
    """
    Return:
        primary_state
        primary_rule_id
        matched_rule_ids
        matched_rule_groups
    """
    matches: list[dict] = []

    for rule in compiled_rules:
        text_parts = [
            field_value(record, key)
            for key in rule["field_keys"]
        ]
        text = " | ".join(part for part in text_parts if part)

        if not text:
            continue

        if rule["regex"].search(text):
            matches.append(rule)

    if not matches:
        return (
            "u_other_unresolved",
            "",
            "",
            "",
        )

    # compiled_rules are already priority-descending.
    primary = matches[0]

    return (
        primary["evidence_state"],
        primary["rule_id"],
        ";".join(rule["rule_id"] for rule in matches),
        ";".join(sorted({rule["rule_group"] for rule in matches})),
    )


def classify_chunk(
    chunk: pd.DataFrame,
    positive_sequences: set[str],
    positive_mapping: dict[str, str],
    compiled_rules: list[dict],
) -> pd.DataFrame:
    chunk = chunk.fillna("").copy()

    chunk["sequence_norm"] = chunk["Sequence"].map(normalize_sequence)
    chunk["sequence_sha256"] = chunk["sequence_norm"].map(sequence_sha256)

    chunk["sequence_length_observed"] = chunk["sequence_norm"].str.len()
    chunk["sequence_is_canonical20"] = chunk["sequence_norm"].map(
        lambda seq: bool(seq) and set(seq).issubset(CANONICAL_AA)
    )

    # Annotation-derived state first. Exact PhageLysData overlap then takes
    # precedence but the underlying rule matches are preserved for audit.
    rule_results = chunk.apply(
        lambda row: evaluate_rules(row, compiled_rules),
        axis=1,
        result_type="expand",
    )
    rule_results.columns = [
        "local_evidence_state",
        "primary_rule_id",
        "matched_rule_ids",
        "matched_rule_groups",
    ]

    chunk = pd.concat(
        [chunk.reset_index(drop=True), rule_results.reset_index(drop=True)],
        axis=1,
    )

    chunk["is_phagelysdata_core_exact"] = chunk["sequence_norm"].isin(
        positive_sequences
    )

    chunk["phagelysdata_core_match"] = chunk["sequence_norm"].map(
        positive_mapping
    ).fillna("")

    chunk["evidence_state"] = chunk["local_evidence_state"]
    chunk.loc[
        chunk["is_phagelysdata_core_exact"],
        "evidence_state",
    ] = "known_p_overlap"

    # Explicitly expose whether multiple local rule groups matched.
    chunk["n_matched_rules"] = chunk["matched_rule_ids"].map(
        lambda x: 0 if not x else len(x.split(";"))
    )
    chunk["n_matched_rule_groups"] = chunk["matched_rule_groups"].map(
        lambda x: 0 if not x else len(x.split(";"))
    )
    chunk["has_rule_conflict"] = chunk["n_matched_rule_groups"] > 1

    return chunk


def write_header_if_needed(
    df: pd.DataFrame,
    path: Path,
    first_write: bool,
) -> None:
    mode = "wt" if first_write else "at"
    with gzip.open(path, mode, encoding="utf-8", newline="") as handle:
        df.to_csv(
            handle,
            sep="\t",
            index=False,
            header=first_write,
            lineterminator="\n",
        )


def main() -> None:
    args = parse_args()

    background_path = Path(args.background).expanduser().resolve()
    phagelys_path = Path(args.phagelysdata_core).expanduser().resolve()
    rules_path = Path(args.rules).expanduser().resolve()
    output_dir = Path(args.output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    header = validate_background_header(background_path)
    positive_sequences, positive_mapping = load_positive_sequences(
        phagelys_path
    )
    rules_df, compiled_rules = load_rules(rules_path)

    print("=" * 92)
    print("PhageLysClass — classify UniProt phage background")
    print("=" * 92)
    print(f"Background       : {background_path}")
    print(f"PhageLysData Core: {phagelys_path}")
    print(f"Rules            : {rules_path}")
    print(f"Enabled rules    : {len(compiled_rules):,}")
    print(f"Core sequences   : {len(positive_sequences):,}")
    print(f"Chunksize        : {args.chunksize:,}")
    print(f"Output dir       : {output_dir}")
    print()

    master_path = output_dir / "uniprot_background_classified.tsv.gz"

    if master_path.exists():
        master_path.unlink()

    state_paths: dict[str, Path] = {}
    state_first_write: dict[str, bool] = {}

    counts = Counter()
    local_counts = Counter()
    rule_counts = Counter()
    scope_counts = Counter()
    canonical_counts = Counter()
    total_rows = 0
    exact_overlap_rows = 0

    for chunk_index, chunk in enumerate(
        pd.read_csv(
            background_path,
            sep="\t",
            compression="infer",
            dtype=str,
            chunksize=args.chunksize,
            low_memory=False,
        ),
        start=1,
    ):
        classified = classify_chunk(
            chunk=chunk,
            positive_sequences=positive_sequences,
            positive_mapping=positive_mapping,
            compiled_rules=compiled_rules,
        )

        first_master = chunk_index == 1
        write_header_if_needed(
            classified,
            master_path,
            first_master,
        )

        if args.write_state_files:
            for state, subset in classified.groupby(
                "evidence_state",
                sort=False,
            ):
                safe_state = re.sub(r"[^A-Za-z0-9_.-]+", "_", state)
                path = state_paths.setdefault(
                    state,
                    output_dir / f"uniprot_background_{safe_state}.tsv.gz",
                )
                first = state_first_write.get(state, True)
                write_header_if_needed(subset, path, first)
                state_first_write[state] = False

        total_rows += len(classified)
        exact_overlap_rows += int(
            classified["is_phagelysdata_core_exact"].sum()
        )

        counts.update(classified["evidence_state"].value_counts().to_dict())
        local_counts.update(
            classified["local_evidence_state"].value_counts().to_dict()
        )
        canonical_counts.update(
            classified["sequence_is_canonical20"]
            .astype(str)
            .value_counts()
            .to_dict()
        )

        if "retrieval_scope" in classified.columns:
            scope_counts.update(
                classified["retrieval_scope"].value_counts().to_dict()
            )

        for rule_ids in classified["matched_rule_ids"]:
            if not rule_ids:
                continue
            for rule_id in rule_ids.split(";"):
                rule_counts[rule_id] += 1

        print(
            f"[CHUNK {chunk_index:>3}] "
            f"processed={total_rows:,} | "
            f"P-overlap={exact_overlap_rows:,}"
        )

    if total_rows == 0:
        raise ValueError("No records were processed from UniProt background.")

    summary_rows = []

    for state, n in sorted(counts.items()):
        summary_rows.append(
            {
                "section": "final_evidence_state",
                "name": state,
                "count": int(n),
                "fraction": n / total_rows,
            }
        )

    for state, n in sorted(local_counts.items()):
        summary_rows.append(
            {
                "section": "local_evidence_state_before_P_overlap",
                "name": state,
                "count": int(n),
                "fraction": n / total_rows,
            }
        )

    for rule_id, n in sorted(rule_counts.items()):
        summary_rows.append(
            {
                "section": "rule_match_count",
                "name": rule_id,
                "count": int(n),
                "fraction": n / total_rows,
            }
        )

    for scope, n in sorted(scope_counts.items()):
        summary_rows.append(
            {
                "section": "retrieval_scope",
                "name": scope,
                "count": int(n),
                "fraction": n / total_rows,
            }
        )

    for value, n in sorted(canonical_counts.items()):
        summary_rows.append(
            {
                "section": "sequence_is_canonical20",
                "name": value,
                "count": int(n),
                "fraction": n / total_rows,
            }
        )

    summary = pd.DataFrame(summary_rows)
    summary_path = output_dir / "uniprot_background_classification_summary.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)

    # Freeze the active rule registry next to outputs for provenance.
    frozen_rules_path = output_dir / "uniprot_background_evidence_rules_frozen.tsv"
    rules_df.to_csv(frozen_rules_path, sep="\t", index=False)

    manifest = {
        "background": str(background_path),
        "phagelysdata_core": str(phagelys_path),
        "rules": str(rules_path),
        "n_enabled_rules": len(compiled_rules),
        "n_phagelysdata_core_sequences": len(positive_sequences),
        "n_rows_processed": total_rows,
        "n_exact_phagelysdata_core_overlaps": exact_overlap_rows,
        "master_output": str(master_path),
        "summary_output": str(summary_path),
        "write_state_files": bool(args.write_state_files),
        "design_note": (
            "Evidence-state candidates only. Final P/N/U and task-specific "
            "datasets are constructed downstream."
        ),
    }

    manifest_path = output_dir / "classification_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2),
        encoding="utf-8",
    )

    print()
    print("=" * 92)
    print("Classification complete")
    print("=" * 92)
    print(f"Rows processed          : {total_rows:,}")
    print(f"Exact Core P overlaps   : {exact_overlap_rows:,}")
    print()
    print("Evidence-state counts")
    print("-" * 92)
    for state, n in counts.most_common():
        print(f"{state:32s} {n:>10,}  ({100*n/total_rows:6.2f}%)")
    print()
    print(f"Master classified table : {master_path}")
    print(f"Summary                 : {summary_path}")
    print(f"Frozen rules            : {frozen_rules_path}")
    print(f"Manifest                : {manifest_path}")
    print()
    print(
        "PASS — evidence-state classification only; no final negative labels, "
        "homology filtering, balancing, or task construction applied."
    )


if __name__ == "__main__":
    main()

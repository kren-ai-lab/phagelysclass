#!/usr/bin/env python3
"""
04_audit_uniprot_background.py

Small audit of the evidence-state classification produced by
03_classify_uniprot_background.py.

The audit is descriptive only. It does NOT change evidence states, define final
negatives, remove sequences, balance data, or create task-specific datasets.

Main questions
--------------
1. What composes `lytic_related`?
2. Which rules dominate n_candidate_easy / n_candidate_hard /
   ambiguous_structural?
3. How do Reviewed status, annotation score, protein existence, and fragment
   status vary by evidence state?
4. Which PhageLysData canonical classes appear among exact Core overlaps?
5. Do local rules conflict with one another?
6. Can we inspect a reproducible random sample from every evidence state?

Outputs
-------
<output-dir>/
├── audit_state_overview.tsv
├── audit_primary_rules_by_state.tsv
├── audit_rule_group_combinations.tsv
├── audit_metadata_distributions.tsv
├── audit_annotation_score_by_state.tsv
├── audit_known_p_overlap_classes.tsv
├── audit_lytic_related_protein_names.tsv
├── audit_lytic_related_rules.tsv
├── audit_rule_conflicts_by_state.tsv
├── audit_random_samples.tsv
└── audit_report.md
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd


DEFAULT_CLASSIFIED = (
    "data/processed/uniprot_background/"
    "uniprot_background_classified.tsv.gz"
)

STATE_ORDER = [
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

# Internal audit names -> columns written by UniProt / script 03.
COLUMN_CANDIDATES = {
    "accession": ["Entry", "accession"],
    "protein_name": ["Protein names", "protein_name"],
    "organism": ["Organism", "organism_name"],
    "reviewed": ["Reviewed", "reviewed"],
    "protein_existence": ["Protein existence", "protein_existence"],
    "annotation_score": ["Annotation", "annotation_score"],
    "fragment": ["Fragment", "fragment"],
    "cc_function": ["Function [CC]", "cc_function"],
    "cc_catalytic_activity": [
        "Catalytic activity",
        "cc_catalytic_activity",
    ],
    "keywords": ["Keywords", "keyword"],
    "retrieval_scope": ["retrieval_scope", "retrieval_scopes"],
    "evidence_state": ["evidence_state"],
    "local_evidence_state": ["local_evidence_state"],
    "primary_rule_id": ["primary_rule_id"],
    "matched_rule_ids": ["matched_rule_ids"],
    "matched_rule_groups": ["matched_rule_groups"],
    "has_rule_conflict": ["has_rule_conflict"],
    "n_matched_rules": ["n_matched_rules"],
    "n_matched_rule_groups": ["n_matched_rule_groups"],
    "is_phagelysdata_core_exact": ["is_phagelysdata_core_exact"],
    "phagelysdata_core_match": ["phagelysdata_core_match"],
}


def parse_args() -> argparse.Namespace:
    script_dir = Path(__file__).resolve().parent
    project_root = script_dir.parents[1]

    parser = argparse.ArgumentParser(
        description="Audit the classified UniProt phage background."
    )
    parser.add_argument(
        "--classified",
        default=str(project_root / DEFAULT_CLASSIFIED),
        help="Classified TSV.GZ from script 03.",
    )
    parser.add_argument(
        "--output-dir",
        default=str(
            project_root
            / "data"
            / "processed"
            / "uniprot_background"
            / "audit"
        ),
        help="Audit output directory.",
    )
    parser.add_argument(
        "--sample-per-state",
        type=int,
        default=75,
        help="Random audit records sampled per evidence state.",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=50,
        help="Number of top protein names / combinations to report.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for reproducible audit samples.",
    )
    return parser.parse_args()


def resolve_columns(header: list[str]) -> dict[str, str]:
    resolved = {}

    for internal, candidates in COLUMN_CANDIDATES.items():
        for candidate in candidates:
            if candidate in header:
                resolved[internal] = candidate
                break

    required = {
        "evidence_state",
        "primary_rule_id",
        "matched_rule_ids",
        "matched_rule_groups",
        "phagelysdata_core_match",
    }
    missing = required - set(resolved)
    if missing:
        raise KeyError(
            "Classified table is missing required audit columns: "
            + ", ".join(sorted(missing))
        )

    return resolved


def normalize_bool(series: pd.Series) -> pd.Series:
    return (
        series.fillna("")
        .astype(str)
        .str.strip()
        .str.lower()
        .isin({"true", "1", "yes", "y"})
    )


def parse_annotation_score(value: object) -> float | None:
    """
    Parse common UniProt annotation-score representations:
        "5 out of 5"
        "5"
        "5.0"
    """
    if value is None or pd.isna(value):
        return None

    text = str(value).strip()
    if not text:
        return None

    match = re.search(r"([0-9]+(?:\.[0-9]+)?)", text)
    if not match:
        return None

    return float(match.group(1))


def ordered_state_frame(
    df: pd.DataFrame,
    state_col: str,
) -> pd.DataFrame:
    states = list(df[state_col].dropna().astype(str).unique())
    rank = {state: i for i, state in enumerate(STATE_ORDER)}

    return df.assign(
        _state_rank=df[state_col].map(lambda x: rank.get(str(x), 999))
    ).sort_values(
        ["_state_rank", state_col],
        kind="stable",
    ).drop(columns="_state_rank")


def build_state_overview(
    df: pd.DataFrame,
    cols: dict[str, str],
) -> pd.DataFrame:
    state = cols["evidence_state"]
    total = len(df)

    out = (
        df.groupby(state, dropna=False)
        .size()
        .rename("count")
        .reset_index()
        .rename(columns={state: "evidence_state"})
    )
    out["fraction"] = out["count"] / total
    out["percent"] = 100.0 * out["fraction"]

    return ordered_state_frame(out, "evidence_state")


def build_primary_rule_table(
    df: pd.DataFrame,
    cols: dict[str, str],
) -> pd.DataFrame:
    state = cols["evidence_state"]
    rule = cols["primary_rule_id"]

    work = df[[state, rule]].copy()
    work[rule] = work[rule].replace("", "NO_PRIMARY_RULE")

    out = (
        work.groupby([state, rule], dropna=False)
        .size()
        .rename("count")
        .reset_index()
        .rename(
            columns={
                state: "evidence_state",
                rule: "primary_rule_id",
            }
        )
    )

    totals = out.groupby("evidence_state")["count"].transform("sum")
    out["fraction_within_state"] = out["count"] / totals

    out = ordered_state_frame(out, "evidence_state")
    return out.sort_values(
        ["_dummy"] if False else ["evidence_state", "count"],
        ascending=[True, False],
        kind="stable",
    )


def build_rule_group_combinations(
    df: pd.DataFrame,
    cols: dict[str, str],
    top_n: int,
) -> pd.DataFrame:
    state = cols["evidence_state"]
    groups = cols["matched_rule_groups"]

    work = df[[state, groups]].copy()
    work[groups] = work[groups].replace("", "NO_RULE_GROUP")

    out = (
        work.groupby([state, groups], dropna=False)
        .size()
        .rename("count")
        .reset_index()
        .rename(
            columns={
                state: "evidence_state",
                groups: "matched_rule_groups",
            }
        )
    )

    totals = out.groupby("evidence_state")["count"].transform("sum")
    out["fraction_within_state"] = out["count"] / totals

    out = (
        out.sort_values(
            ["evidence_state", "count"],
            ascending=[True, False],
            kind="stable",
        )
        .groupby("evidence_state", group_keys=False)
        .head(top_n)
        .reset_index(drop=True)
    )

    return ordered_state_frame(out, "evidence_state")


def build_metadata_distributions(
    df: pd.DataFrame,
    cols: dict[str, str],
) -> pd.DataFrame:
    state_col = cols["evidence_state"]
    outputs = []

    for internal_name in [
        "reviewed",
        "protein_existence",
        "fragment",
        "retrieval_scope",
    ]:
        if internal_name not in cols:
            continue

        source_col = cols[internal_name]

        work = df[[state_col, source_col]].copy()
        work[source_col] = (
            work[source_col]
            .fillna("")
            .astype(str)
            .str.strip()
            .replace("", "NA")
        )

        counts = (
            work.groupby([state_col, source_col], dropna=False)
            .size()
            .rename("count")
            .reset_index()
            .rename(
                columns={
                    state_col: "evidence_state",
                    source_col: "value",
                }
            )
        )

        totals = counts.groupby("evidence_state")["count"].transform("sum")
        counts["fraction_within_state"] = counts["count"] / totals
        counts.insert(1, "metadata_field", internal_name)

        outputs.append(counts)

    if not outputs:
        return pd.DataFrame(
            columns=[
                "evidence_state",
                "metadata_field",
                "value",
                "count",
                "fraction_within_state",
            ]
        )

    result = pd.concat(outputs, ignore_index=True)
    return ordered_state_frame(result, "evidence_state")


def build_annotation_score_summary(
    df: pd.DataFrame,
    cols: dict[str, str],
) -> pd.DataFrame:
    if "annotation_score" not in cols:
        return pd.DataFrame()

    state = cols["evidence_state"]
    annotation = cols["annotation_score"]

    work = df[[state, annotation]].copy()
    work["annotation_score_numeric"] = work[annotation].map(
        parse_annotation_score
    )

    summary = (
        work.groupby(state, dropna=False)["annotation_score_numeric"]
        .agg(
            n_nonmissing="count",
            mean="mean",
            median="median",
            min="min",
            max="max",
        )
        .reset_index()
        .rename(columns={state: "evidence_state"})
    )

    total_by_state = work.groupby(state).size()
    summary["n_total"] = summary["evidence_state"].map(total_by_state)
    summary["fraction_nonmissing"] = (
        summary["n_nonmissing"] / summary["n_total"]
    )

    return ordered_state_frame(summary, "evidence_state")


def build_known_p_overlap_classes(
    df: pd.DataFrame,
    cols: dict[str, str],
) -> pd.DataFrame:
    state = cols["evidence_state"]
    match_col = cols["phagelysdata_core_match"]

    overlap = df.loc[
        df[state].eq("known_p_overlap"),
        [match_col],
    ].copy()

    if overlap.empty:
        return pd.DataFrame(
            columns=["canonical_target_class", "count", "fraction"]
        )

    # Script 03 stores: sequence_entity_id|canonical_target_class
    overlap["canonical_target_class"] = (
        overlap[match_col]
        .fillna("")
        .astype(str)
        .str.rsplit("|", n=1)
        .str[-1]
    )

    out = (
        overlap["canonical_target_class"]
        .value_counts(dropna=False)
        .rename_axis("canonical_target_class")
        .reset_index(name="count")
    )
    out["fraction"] = out["count"] / len(overlap)

    return out


def build_lytic_name_table(
    df: pd.DataFrame,
    cols: dict[str, str],
    top_n: int,
) -> pd.DataFrame:
    if "protein_name" not in cols:
        return pd.DataFrame()

    state = cols["evidence_state"]
    protein = cols["protein_name"]

    subset = df.loc[
        df[state].eq("lytic_related"),
        [protein],
    ].copy()

    subset[protein] = (
        subset[protein]
        .fillna("")
        .astype(str)
        .str.strip()
        .replace("", "NA")
    )

    out = (
        subset[protein]
        .value_counts(dropna=False)
        .head(top_n)
        .rename_axis("protein_name")
        .reset_index(name="count")
    )

    if len(subset):
        out["fraction_within_lytic_related"] = out["count"] / len(subset)

    return out


def build_lytic_rule_table(
    df: pd.DataFrame,
    cols: dict[str, str],
) -> pd.DataFrame:
    state = cols["evidence_state"]
    rule = cols["primary_rule_id"]

    subset = df.loc[
        df[state].eq("lytic_related"),
        [rule],
    ].copy()

    subset[rule] = subset[rule].replace("", "NO_PRIMARY_RULE")

    out = (
        subset[rule]
        .value_counts(dropna=False)
        .rename_axis("primary_rule_id")
        .reset_index(name="count")
    )

    if len(subset):
        out["fraction_within_lytic_related"] = out["count"] / len(subset)

    return out


def build_conflict_table(
    df: pd.DataFrame,
    cols: dict[str, str],
) -> pd.DataFrame:
    if "has_rule_conflict" not in cols:
        return pd.DataFrame()

    state = cols["evidence_state"]
    conflict = cols["has_rule_conflict"]

    work = df[[state, conflict]].copy()
    work["has_rule_conflict_bool"] = normalize_bool(work[conflict])

    out = (
        work.groupby(state, dropna=False)["has_rule_conflict_bool"]
        .agg(
            n_total="size",
            n_conflict="sum",
        )
        .reset_index()
        .rename(columns={state: "evidence_state"})
    )
    out["fraction_conflict"] = out["n_conflict"] / out["n_total"]

    return ordered_state_frame(out, "evidence_state")


def build_random_samples(
    df: pd.DataFrame,
    cols: dict[str, str],
    sample_per_state: int,
    seed: int,
) -> pd.DataFrame:
    state = cols["evidence_state"]

    preferred = [
        "accession",
        "protein_name",
        "organism",
        "reviewed",
        "protein_existence",
        "annotation_score",
        "fragment",
        "cc_function",
        "cc_catalytic_activity",
        "keywords",
        "retrieval_scope",
        "evidence_state",
        "local_evidence_state",
        "primary_rule_id",
        "matched_rule_ids",
        "matched_rule_groups",
        "has_rule_conflict",
        "phagelysdata_core_match",
    ]

    selected_internal = [x for x in preferred if x in cols]
    selected_columns = [cols[x] for x in selected_internal]

    samples = []

    for evidence_state, subset in df.groupby(state, sort=False):
        n = min(len(subset), sample_per_state)
        if n == 0:
            continue

        sampled = subset.sample(
            n=n,
            random_state=seed,
        )[selected_columns].copy()

        rename = {
            cols[internal]: internal
            for internal in selected_internal
        }
        sampled = sampled.rename(columns=rename)
        samples.append(sampled)

    if not samples:
        return pd.DataFrame()

    result = pd.concat(samples, ignore_index=True)

    if "evidence_state" in result.columns:
        result = ordered_state_frame(result, "evidence_state")

    return result


def write_report(
    output_path: Path,
    overview: pd.DataFrame,
    p_classes: pd.DataFrame,
    lytic_rules: pd.DataFrame,
    conflicts: pd.DataFrame,
    annotation_summary: pd.DataFrame,
) -> None:
    total = int(overview["count"].sum())

    lines = [
        "# UniProt background evidence-state audit",
        "",
        f"Total classified records: **{total:,}**.",
        "",
        "## Evidence-state overview",
        "",
        "| Evidence state | Count | Percent |",
        "|---|---:|---:|",
    ]

    for row in overview.itertuples(index=False):
        lines.append(
            f"| {row.evidence_state} | {int(row.count):,} | "
            f"{float(row.percent):.2f}% |"
        )

    lines.extend(
        [
            "",
            "## Known PhageLysData Core overlaps",
            "",
        ]
    )

    if p_classes.empty:
        lines.append("No exact Core overlaps were detected.")
    else:
        lines.extend(
            [
                "| Canonical target class | Count | Fraction |",
                "|---|---:|---:|",
            ]
        )
        for row in p_classes.itertuples(index=False):
            lines.append(
                f"| {row.canonical_target_class} | {int(row.count):,} | "
                f"{float(row.fraction):.3f} |"
            )

    lines.extend(
        [
            "",
            "## Lytic-related primary rules",
            "",
        ]
    )

    if lytic_rules.empty:
        lines.append("No lytic-related records were detected.")
    else:
        lines.extend(
            [
                "| Primary rule | Count | Fraction |",
                "|---|---:|---:|",
            ]
        )
        for row in lytic_rules.itertuples(index=False):
            lines.append(
                f"| {row.primary_rule_id} | {int(row.count):,} | "
                f"{float(row.fraction_within_lytic_related):.3f} |"
            )

    if not conflicts.empty:
        lines.extend(
            [
                "",
                "## Rule conflicts",
                "",
                "| Evidence state | Conflicts | Total | Fraction |",
                "|---|---:|---:|---:|",
            ]
        )
        for row in conflicts.itertuples(index=False):
            lines.append(
                f"| {row.evidence_state} | {int(row.n_conflict):,} | "
                f"{int(row.n_total):,} | {float(row.fraction_conflict):.3f} |"
            )

    if not annotation_summary.empty:
        lines.extend(
            [
                "",
                "## Annotation score",
                "",
                "| Evidence state | N with score | Mean | Median |",
                "|---|---:|---:|---:|",
            ]
        )
        for row in annotation_summary.itertuples(index=False):
            mean = "NA" if pd.isna(row.mean) else f"{row.mean:.2f}"
            median = "NA" if pd.isna(row.median) else f"{row.median:.2f}"
            lines.append(
                f"| {row.evidence_state} | {int(row.n_nonmissing):,} | "
                f"{mean} | {median} |"
            )

    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "- This audit is descriptive only; candidate states are not final ML labels.",
            "- `lytic_related`, `ambiguous_structural`, and rule-conflict records should be manually inspected before freezing task definitions.",
            "- `n_candidate_easy` and `n_candidate_hard` remain candidate negative pools, not biological ground truth negatives.",
            "- `u_other_unresolved` should not be converted to negative by default.",
            "",
        ]
    )

    output_path.write_text(
        "\n".join(lines),
        encoding="utf-8",
    )


def main() -> None:
    args = parse_args()

    classified_path = Path(args.classified).expanduser().resolve()
    output_dir = Path(args.output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    if not classified_path.is_file():
        raise FileNotFoundError(
            f"Classified background not found: {classified_path}"
        )

    header = pd.read_csv(
        classified_path,
        sep="\t",
        compression="infer",
        nrows=0,
    ).columns.tolist()

    cols = resolve_columns(header)

    # Audit intentionally avoids loading Sequence to keep memory usage modest.
    usecols = sorted(set(cols.values()))

    print("=" * 92)
    print("PhageLysClass — UniProt background audit")
    print("=" * 92)
    print(f"Classified table : {classified_path}")
    print(f"Output dir       : {output_dir}")
    print(f"Sample/state     : {args.sample_per_state}")
    print(f"Random seed      : {args.seed}")
    print()

    df = pd.read_csv(
        classified_path,
        sep="\t",
        compression="infer",
        dtype=str,
        usecols=usecols,
        low_memory=False,
    ).fillna("")

    overview = build_state_overview(df, cols)
    primary_rules = build_primary_rule_table(df, cols)
    rule_groups = build_rule_group_combinations(
        df,
        cols,
        top_n=args.top_n,
    )
    metadata = build_metadata_distributions(df, cols)
    annotation_summary = build_annotation_score_summary(df, cols)
    p_classes = build_known_p_overlap_classes(df, cols)
    lytic_names = build_lytic_name_table(
        df,
        cols,
        top_n=args.top_n,
    )
    lytic_rules = build_lytic_rule_table(df, cols)
    conflicts = build_conflict_table(df, cols)
    samples = build_random_samples(
        df,
        cols,
        sample_per_state=args.sample_per_state,
        seed=args.seed,
    )

    outputs = {
        "state_overview": output_dir / "audit_state_overview.tsv",
        "primary_rules": output_dir / "audit_primary_rules_by_state.tsv",
        "rule_groups": output_dir / "audit_rule_group_combinations.tsv",
        "metadata": output_dir / "audit_metadata_distributions.tsv",
        "annotation": output_dir / "audit_annotation_score_by_state.tsv",
        "p_classes": output_dir / "audit_known_p_overlap_classes.tsv",
        "lytic_names": output_dir / "audit_lytic_related_protein_names.tsv",
        "lytic_rules": output_dir / "audit_lytic_related_rules.tsv",
        "conflicts": output_dir / "audit_rule_conflicts_by_state.tsv",
        "samples": output_dir / "audit_random_samples.tsv",
        "report": output_dir / "audit_report.md",
    }

    overview.to_csv(outputs["state_overview"], sep="\t", index=False)
    primary_rules.to_csv(outputs["primary_rules"], sep="\t", index=False)
    rule_groups.to_csv(outputs["rule_groups"], sep="\t", index=False)
    metadata.to_csv(outputs["metadata"], sep="\t", index=False)
    annotation_summary.to_csv(outputs["annotation"], sep="\t", index=False)
    p_classes.to_csv(outputs["p_classes"], sep="\t", index=False)
    lytic_names.to_csv(outputs["lytic_names"], sep="\t", index=False)
    lytic_rules.to_csv(outputs["lytic_rules"], sep="\t", index=False)
    conflicts.to_csv(outputs["conflicts"], sep="\t", index=False)
    samples.to_csv(outputs["samples"], sep="\t", index=False)

    write_report(
        output_path=outputs["report"],
        overview=overview,
        p_classes=p_classes,
        lytic_rules=lytic_rules,
        conflicts=conflicts,
        annotation_summary=annotation_summary,
    )

    print("Evidence-state overview")
    print("-" * 92)
    for row in overview.itertuples(index=False):
        print(
            f"{row.evidence_state:32s} "
            f"{int(row.count):>10,} "
            f"({float(row.percent):6.2f}%)"
        )

    print()
    print("Audit outputs")
    print("-" * 92)
    for path in outputs.values():
        print(path)

    print()
    print(
        "PASS — descriptive audit only; no evidence states or ML labels modified."
    )


if __name__ == "__main__":
    main()

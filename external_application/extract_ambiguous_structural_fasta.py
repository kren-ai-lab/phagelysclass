#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


CANONICAL_AA = set("ACDEFGHIKLMNPQRSTVWY")
MAX_LENGTH = 1024


def parse_args():
    p = argparse.ArgumentParser(
        description="Extract the uniprot_ambiguous_structural challenge set to a QC-filtered FASTA."
    )
    p.add_argument("--input", required=True, type=Path)
    p.add_argument("--output-dir", type=Path, default=Path("ambiguous_structural_input"))
    return p.parse_args()


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.input, sep="\t", compression="infer")

    required = {"entity_id", "sequence", "challenge_set", "training_role"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    # Keep only the intended frozen challenge set.
    df = df[
        (df["challenge_set"] == "uniprot_ambiguous_structural")
        & (df["training_role"] == "challenge_only")
    ].copy()

    df["sequence"] = (
        df["sequence"]
        .fillna("")
        .astype(str)
        .str.replace(r"\s+", "", regex=True)
        .str.upper()
    )
    df["sequence_length_qc"] = df["sequence"].str.len()

    def qc_reason(seq):
        if not seq:
            return "empty_sequence"
        invalid = sorted(set(seq) - CANONICAL_AA)
        if invalid:
            return "noncanonical:" + "".join(invalid)
        if len(seq) > MAX_LENGTH:
            return f"length_gt_{MAX_LENGTH}"
        return ""

    df["qc_reason"] = df["sequence"].map(qc_reason)

    passed = df[df["qc_reason"] == ""].copy()
    rejected = df[df["qc_reason"] != ""].copy()

    # Exact duplicate sequences: keep first deterministically.
    dup_mask = passed.duplicated(subset=["sequence"], keep="first")
    duplicates = passed[dup_mask].copy()
    duplicates["qc_reason"] = "exact_sequence_duplicate"
    rejected = pd.concat([rejected, duplicates], ignore_index=True)
    passed = passed[~dup_mask].copy()

    passed = passed.sort_values("entity_id").reset_index(drop=True)

    passed.to_csv(
        args.output_dir / "ambiguous_structural_qc_pass.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    rejected.to_csv(
        args.output_dir / "ambiguous_structural_qc_rejected.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    fasta = args.output_dir / "ambiguous_structural_qc_pass.fasta"
    with fasta.open("w") as handle:
        for row in passed.itertuples(index=False):
            handle.write(f">{row.entity_id}\n")
            seq = row.sequence
            for i in range(0, len(seq), 80):
                handle.write(seq[i:i+80] + "\n")

    print(f"Input challenge rows: {len(df):,}")
    print(f"QC-pass sequences:    {len(passed):,}")
    print(f"Rejected:             {len(rejected):,}")
    print(f"FASTA:                {fasta}")


if __name__ == "__main__":
    main()

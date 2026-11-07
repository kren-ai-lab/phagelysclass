#!/usr/bin/env python3
"""
Prepare the PhageLysData prediction-extension sequence set for PhageLysClass use.

Workflow
--------
1. Load phagelysdata_prediction_extension.parquet
2. Recover prediction-extension sequences only
3. QC:
   - non-empty ID and sequence
   - canonical 20-aa alphabet
   - length <= 1024 aa
   - exact-sequence deduplication
4. Apply one fixed BioSieve/MMseqs2 redundancy reduction:
   - minimum sequence identity = 0.30
   - coverage = 0.80
5. Recover full metadata for MMseqs2 representatives
6. Export representative TSV + FASTA + mapping + reports

Expected source columns
-----------------------
sequence_entity_id
sequence_sha256
sequence
sequence_length
entity_scope
prediction_extension_member
canonical_target_class
canonical_class_basis
assigned_tier
source_observation_count
source_family_count
source_families_json
valid_target_assertion_count
prediction_target_assertion_count
nonprediction_target_assertion_count
prediction_target_classes_json
prediction_target_sources_json
phage_association_status
publication_count

Requirements
------------
pip install pandas pyarrow biosieve
MMseqs2 must be available as `mmseqs` in PATH.

Example
-------
python prepare_phagelys_prediction_extension.py \
    --input phagelysdata_prediction_extension.parquet \
    --output-dir phagelys_prediction_extension_h30
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path

import pandas as pd


# ---------------------------------------------------------------------
# Frozen test configuration
# ---------------------------------------------------------------------

IDENTITY_THRESHOLD = 0.30
COVERAGE = 0.80

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")
MAX_LENGTH = 1024

REQUIRED_COLUMNS = {
    "sequence_entity_id",
    "sequence_sha256",
    "sequence",
    "sequence_length",
    "entity_scope",
    "prediction_extension_member",
}


# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Prepare and redundancy-reduce PhageLysData prediction-extension sequences."
    )
    p.add_argument(
        "--input",
        type=Path,
        default=Path("phagelysdata_prediction_extension.parquet"),
        help="PhageLysData prediction-extension parquet.",
    )
    p.add_argument(
        "--output-dir",
        type=Path,
        default=Path("phagelys_prediction_extension_h30"),
        help="Output directory.",
    )
    p.add_argument(
        "--keep-existing",
        action="store_true",
        help="Allow re-use of an existing output directory.",
    )
    p.add_argument(
        "--identity-param",
        default=None,
        help=(
            "Optional explicit BioSieve MMseqs2 identity parameter name. "
            "If omitted, the script probes the installed BioSieve version."
        ),
    )
    p.add_argument(
        "--coverage-param",
        default=None,
        help=(
            "Optional explicit BioSieve MMseqs2 coverage parameter name. "
            "If omitted, the script probes the installed BioSieve version."
        ),
    )
    return p.parse_args()


def run(cmd: list[str], *, check: bool = True, capture: bool = False) -> subprocess.CompletedProcess:
    print("$", " ".join(map(str, cmd)), flush=True)
    return subprocess.run(
        cmd,
        check=check,
        text=True,
        capture_output=capture,
    )


def command_version(cmd: list[str]) -> str | None:
    try:
        r = subprocess.run(cmd, text=True, capture_output=True, check=False)
        txt = (r.stdout or r.stderr or "").strip()
        return txt.splitlines()[0] if txt else None
    except Exception:
        return None


def sha256_file(path: Path, block_size: int = 1024 * 1024) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(block_size)
            if not block:
                break
            h.update(block)
    return h.hexdigest()


def normalize_sequence(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return "".join(str(value).split()).upper()


def qc_reason(seq_id: object, seq: str) -> str | None:
    if seq_id is None or pd.isna(seq_id) or not str(seq_id).strip():
        return "missing_sequence_entity_id"

    if not seq:
        return "empty_sequence"

    invalid = sorted(set(seq) - VALID_AA)
    if invalid:
        return "noncanonical:" + "".join(invalid)

    if len(seq) > MAX_LENGTH:
        return "length_gt_1024"

    return None


def write_fasta(df: pd.DataFrame, path: Path) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for row in df.itertuples(index=False):
            seq_id = str(row.sequence_entity_id)
            seq = str(row.sequence)
            handle.write(f">{seq_id}\n")
            for i in range(0, len(seq), 80):
                handle.write(seq[i:i + 80] + "\n")


def preflight() -> dict[str, str | None]:
    missing = []

    biosieve = shutil.which("biosieve")
    if biosieve is None:
        missing.append("biosieve")

    mmseqs = shutil.which("mmseqs")
    if mmseqs is None:
        missing.append("mmseqs")

    if missing:
        raise RuntimeError(
            "Missing required executable(s) in PATH: " + ", ".join(missing)
        )

    return {
        "biosieve": biosieve,
        "mmseqs": mmseqs,
        "biosieve_version": command_version(["biosieve", "--version"]),
        "mmseqs_version": command_version(["mmseqs", "version"]),
    }


def biosieve_command(
    input_csv: Path,
    output_csv: Path,
    mapping_csv: Path,
    report_json: Path,
    identity_param: str,
    coverage_param: str,
) -> list[str]:
    return [
        "biosieve",
        "reduce",
        "-i",
        str(input_csv),
        "-o",
        str(output_csv),
        "--strategy",
        "mmseqs2",
        "--mapping-output",
        str(mapping_csv),
        "--report-output",
        str(report_json),
        "--set",
        f"mmseqs2.{identity_param}={IDENTITY_THRESHOLD}",
        "--set",
        f"mmseqs2.{coverage_param}={COVERAGE}",
    ]


def probe_biosieve_mmseqs_params(
    outdir: Path,
    identity_param: str | None,
    coverage_param: str | None,
) -> tuple[str, str]:
    """
    Probe parameter names on a tiny dataset so an incompatible BioSieve version
    fails before the large reduction starts.
    """
    if identity_param and coverage_param:
        return identity_param, coverage_param

    identity_candidates = (
        [identity_param]
        if identity_param
        else ["min_seq_id", "threshold", "identity_threshold", "min_identity"]
    )
    coverage_candidates = (
        [coverage_param]
        if coverage_param
        else ["coverage", "min_coverage", "cov"]
    )

    probe = pd.DataFrame(
        {
            "id": ["probe_1", "probe_2", "probe_3"],
            "sequence": [
                "ACDEFGHIKLMNPQRSTVWYACDEFGHIK",
                "ACDEFGHIKLMNPQRSTVWYACDEFGHIK",
                "MNPQRSTVWYACDEFGHIKLMNPQRSTVW",
            ],
            "label": ["prediction_extension"] * 3,
        }
    )

    probe_root = outdir / "_biosieve_probe"
    if probe_root.exists():
        shutil.rmtree(probe_root)
    probe_root.mkdir(parents=True)

    probe_input = probe_root / "probe.csv"
    probe.to_csv(probe_input, index=False)

    failures: list[dict[str, str]] = []

    for ikey in identity_candidates:
        for ckey in coverage_candidates:
            out_csv = probe_root / f"out_{ikey}_{ckey}.csv"
            map_csv = probe_root / f"map_{ikey}_{ckey}.csv"
            report_json = probe_root / f"report_{ikey}_{ckey}.json"

            cmd = biosieve_command(
                probe_input,
                out_csv,
                map_csv,
                report_json,
                ikey,
                ckey,
            )

            r = subprocess.run(cmd, text=True, capture_output=True, check=False)

            if r.returncode == 0 and out_csv.exists():
                shutil.rmtree(probe_root)
                return ikey, ckey

            failures.append(
                {
                    "identity_param": ikey,
                    "coverage_param": ckey,
                    "stderr": (r.stderr or "")[-2000:],
                    "stdout": (r.stdout or "")[-2000:],
                }
            )

    failure_path = outdir / "biosieve_parameter_probe_failures.json"
    failure_path.write_text(json.dumps(failures, indent=2), encoding="utf-8")

    raise RuntimeError(
        "Could not identify the installed BioSieve MMseqs2 parameter names. "
        f"Probe details were written to {failure_path}. "
        "Re-run with --identity-param and --coverage-param if needed."
    )


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def main() -> None:
    args = parse_args()

    source = args.input.resolve()
    outdir = args.output_dir.resolve()

    if not source.exists():
        raise FileNotFoundError(f"Input not found: {source}")

    if outdir.exists() and any(outdir.iterdir()) and not args.keep_existing:
        raise FileExistsError(
            f"Output directory already exists and is not empty: {outdir}\n"
            "Use a new directory or pass --keep-existing."
        )

    outdir.mkdir(parents=True, exist_ok=True)

    software = preflight()

    # -------------------------------------------------------------
    # 1. Recover prediction-extension sequences
    # -------------------------------------------------------------

    print("\n[1/4] Loading PhageLysData prediction extension", flush=True)
    df = pd.read_parquet(source)

    missing_cols = sorted(REQUIRED_COLUMNS - set(df.columns))
    if missing_cols:
        raise ValueError(
            "Input does not match the expected PhageLysData prediction-extension schema. "
            f"Missing columns: {missing_cols}"
        )

    n_source = len(df)

    extension_mask = df["prediction_extension_member"].fillna(False).astype(bool)

    # Keep the explicit prediction-only scope when present.
    # Rows marked as extension members remain the primary criterion.
    if "entity_scope" in df.columns:
        scope_mask = df["entity_scope"].astype("string").eq("prediction_only_target")
        extension_mask &= scope_mask

    ext = df.loc[extension_mask].copy()
    n_extension = len(ext)

    if n_extension == 0:
        raise RuntimeError("No prediction-extension rows were recovered.")

    ext_path = outdir / "01_prediction_extension_extracted.tsv.gz"
    ext.to_csv(ext_path, sep="\t", index=False, compression="gzip")

    # -------------------------------------------------------------
    # 2. QC
    # -------------------------------------------------------------

    print("[2/4] QC and exact-sequence deduplication", flush=True)

    ext["sequence_original"] = ext["sequence"]
    ext["sequence"] = ext["sequence"].map(normalize_sequence)
    ext["sequence_length_qc"] = ext["sequence"].str.len()

    ext["qc_reason"] = [
        qc_reason(seq_id, seq)
        for seq_id, seq in zip(ext["sequence_entity_id"], ext["sequence"])
    ]

    rejected = ext.loc[ext["qc_reason"].notna()].copy()
    passed = ext.loc[ext["qc_reason"].isna()].copy()

    # Deterministic exact deduplication.
    passed = passed.sort_values("sequence_entity_id", kind="stable").reset_index(drop=True)
    duplicate_mask = passed.duplicated(subset=["sequence"], keep="first")

    exact_duplicates = passed.loc[duplicate_mask].copy()
    exact_duplicates["qc_reason"] = "exact_sequence_duplicate"

    passed = passed.loc[~duplicate_mask].copy().reset_index(drop=True)

    if not exact_duplicates.empty:
        rejected = pd.concat([rejected, exact_duplicates], ignore_index=True)

    qc_pass_path = outdir / "02_qc_pass.tsv.gz"
    qc_fail_path = outdir / "02_qc_rejected.tsv.gz"

    passed.to_csv(qc_pass_path, sep="\t", index=False, compression="gzip")
    rejected.to_csv(qc_fail_path, sep="\t", index=False, compression="gzip")

    if passed.empty:
        raise RuntimeError("No sequences remained after QC.")

    # BioSieve input: use canonical column names.
    biosieve_input = pd.DataFrame(
        {
            "id": passed["sequence_entity_id"].astype(str),
            "sequence": passed["sequence"].astype(str),
            "label": "prediction_extension",
        }
    )

    biosieve_input_path = outdir / "03_biosieve_input.csv"
    biosieve_input.to_csv(biosieve_input_path, index=False)

    # -------------------------------------------------------------
    # 3. BioSieve/MMseqs2 H30
    # -------------------------------------------------------------

    print(
        f"[3/4] BioSieve/MMseqs2 reduction: "
        f"min identity={IDENTITY_THRESHOLD:.2f}, coverage={COVERAGE:.2f}",
        flush=True,
    )

    identity_param, coverage_param = probe_biosieve_mmseqs_params(
        outdir,
        args.identity_param,
        args.coverage_param,
    )

    print(
        f"Detected BioSieve parameters: "
        f"mmseqs2.{identity_param}, mmseqs2.{coverage_param}",
        flush=True,
    )

    reduced_csv = outdir / "04_biosieve_representatives.csv"
    mapping_csv = outdir / "04_biosieve_mapping.csv"
    report_json = outdir / "04_biosieve_report.json"

    cmd = biosieve_command(
        biosieve_input_path,
        reduced_csv,
        mapping_csv,
        report_json,
        identity_param,
        coverage_param,
    )

    run(cmd)

    if not reduced_csv.exists():
        raise RuntimeError("BioSieve completed without producing the reduced dataset.")

    reduced = pd.read_csv(reduced_csv)

    if "id" not in reduced.columns:
        raise ValueError(
            f"BioSieve reduced output does not contain an 'id' column: {list(reduced.columns)}"
        )

    representative_ids = set(reduced["id"].astype(str))

    # -------------------------------------------------------------
    # 4. Recover full PhageLysData metadata
    # -------------------------------------------------------------

    print("[4/4] Recovering representative metadata and writing final files", flush=True)

    reps = passed[
        passed["sequence_entity_id"].astype(str).isin(representative_ids)
    ].copy()

    # Keep BioSieve representative order.
    order = {seq_id: i for i, seq_id in enumerate(reduced["id"].astype(str))}
    reps["_biosieve_order"] = reps["sequence_entity_id"].astype(str).map(order)
    reps = reps.sort_values("_biosieve_order").drop(columns="_biosieve_order").reset_index(drop=True)

    final_tsv = outdir / "05_phagelys_prediction_extension_H30_representatives.tsv.gz"
    final_fasta = outdir / "05_phagelys_prediction_extension_H30_representatives.fasta"

    reps.to_csv(final_tsv, sep="\t", index=False, compression="gzip")
    write_fasta(reps, final_fasta)

    rejection_counts = Counter(rejected["qc_reason"].astype(str)) if not rejected.empty else Counter()

    summary = {
        "source_file": str(source),
        "source_sha256": sha256_file(source),
        "source_rows": int(n_source),
        "prediction_extension_rows": int(n_extension),
        "qc_pass_unique_rows": int(len(passed)),
        "qc_rejected_rows": int(len(rejected)),
        "qc_rejection_counts": dict(rejection_counts),
        "mmseqs_identity_threshold": IDENTITY_THRESHOLD,
        "mmseqs_coverage": COVERAGE,
        "biosieve_identity_param": identity_param,
        "biosieve_coverage_param": coverage_param,
        "representatives": int(len(reps)),
        "removed_by_mmseqs2": int(len(passed) - len(reps)),
        "retained_fraction_after_mmseqs2": float(len(reps) / len(passed)),
        "max_sequence_length": MAX_LENGTH,
        "canonical_alphabet": "".join(sorted(VALID_AA)),
        "software": software,
        "biosieve_command": cmd,
        "outputs": {
            "prediction_extension_extracted": str(ext_path),
            "qc_pass": str(qc_pass_path),
            "qc_rejected": str(qc_fail_path),
            "biosieve_input": str(biosieve_input_path),
            "biosieve_representatives": str(reduced_csv),
            "biosieve_mapping": str(mapping_csv),
            "biosieve_report": str(report_json),
            "final_representatives_tsv": str(final_tsv),
            "final_representatives_fasta": str(final_fasta),
        },
    }

    summary_path = outdir / "05_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print("\nDONE")
    print(f"Source rows:                 {n_source:,}")
    print(f"Prediction-extension rows:   {n_extension:,}")
    print(f"QC-pass unique sequences:    {len(passed):,}")
    print(f"H30 representatives:         {len(reps):,}")
    print(f"Retained fraction:           {len(reps)/len(passed):.4f}")
    print(f"Final FASTA:                 {final_fasta}")
    print(f"Final metadata:              {final_tsv}")
    print(f"Summary:                     {summary_path}")


if __name__ == "__main__":
    main()

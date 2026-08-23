#!/usr/bin/env python3
"""
08_build_factorial_homology_splits.py

Materialize the full PhageLysClass redundancy × split factorial:

    R90/R70/R50/R30 × S90/S70/S50/S30

for all nine trainable Stage-06/07 tasks.

Design
------
- Input datasets are the frozen Stage-07 redundancy-controlled datasets.
- Redundancy threshold R and split threshold S are independent factors.
- Homology clusters are built across the WHOLE task, never separately by class.
- MMseqs2 criterion:
      min_seq_id = S / 100
      coverage   = 0.80
      cov_mode   = 0
      cluster_mode = 0
- BioSieve homology_aware splitting uses the precomputed MMseqs2 clusters.
- Outer split target: train ~= 70%, test ~= 30%; no validation partition.
- Challenge sets are untouched.
- No representation or model information enters this stage.

Canonical outputs
-----------------
data/splits/<task_path>/<R>_<S>/
├── train.tsv.gz
├── test.tsv.gz
├── dataset_with_split.tsv.gz
├── clusters.tsv.gz
├── biosieve_report/...
├── split_summary.tsv
└── manifest.json

At the task level:
data/splits/<task_path>/factorial_split_overview.tsv
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any

import pandas as pd


R_THRESHOLDS = [90, 70, 50, 30]
S_THRESHOLDS = [90, 70, 50, 30]

COVERAGE = 0.80
COV_MODE = 0
CLUSTER_MODE = 0
TEST_SIZE = 0.30
SEED = 42

TASKS = [
    {
        "task_id": "global_easy",
        "task_path": "global_easy",
        "primary_group": "class_name",
        "secondary_group": "sequence_state",
    },
    {
        "task_id": "global_hard",
        "task_path": "global_hard",
        "primary_group": "class_name",
        "secondary_group": "sequence_state",
    },
    {
        "task_id": "global_mixed",
        "task_path": "global_mixed",
        "primary_group": "class_name",
        "secondary_group": "sequence_state",
    },
    {
        "task_id": "subtype_endolysin",
        "task_path": "subtype_binary/endolysin",
        "primary_group": "class_name",
        "secondary_group": "sequence_state",
    },
    {
        "task_id": "subtype_depolymerase",
        "task_path": "subtype_binary/depolymerase",
        "primary_group": "class_name",
        "secondary_group": "sequence_state",
    },
    {
        "task_id": "subtype_val",
        "task_path": "subtype_binary/virion_associated_lytic_enzyme",
        "primary_group": "class_name",
        "secondary_group": "sequence_state",
    },
    {
        "task_id": "multiclass_clean",
        "task_path": "multiclass_clean",
        "primary_group": "class_name",
        "secondary_group": "canonical_target_class",
    },
    {
        "task_id": "pu_conservative",
        "task_path": "pu_global/conservative",
        "primary_group": "pu_role",
        "secondary_group": "sequence_state",
    },
    {
        "task_id": "pu_broad",
        "task_path": "pu_global/broad",
        "primary_group": "pu_role",
        "secondary_group": "sequence_state",
    },
]


def project_root_from_script(script_file: str) -> Path:
    script_dir = Path(script_file).resolve().parent
    return script_dir.parents[1]


def parse_args() -> argparse.Namespace:
    root = project_root_from_script(__file__)

    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--redundancy-root",
        default=str(root / "data/redundancy"),
    )
    p.add_argument(
        "--splits-root",
        default=str(root / "data/splits"),
    )
    p.add_argument(
        "--threads",
        type=int,
        default=int(os.environ.get("MMSEQS_THREADS", "16")),
    )
    p.add_argument(
        "--tmp-dir",
        default=os.environ.get("MMSEQS_TMPDIR"),
        help="Optional parent directory for temporary MMseqs2 work.",
    )
    p.add_argument(
        "--only-task",
        action="append",
        choices=[t["task_id"] for t in TASKS],
        help="Restrict execution to one or more task IDs.",
    )
    p.add_argument(
        "--only-r",
        action="append",
        type=int,
        choices=R_THRESHOLDS,
        help="Restrict redundancy threshold(s).",
    )
    p.add_argument(
        "--only-s",
        action="append",
        type=int,
        choices=S_THRESHOLDS,
        help="Restrict split threshold(s).",
    )
    p.add_argument("--overwrite", action="store_true")
    return p.parse_args()


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def package_version(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except Exception:
        return "unknown"


def command_version(command: list[str]) -> str:
    try:
        result = subprocess.run(
            command,
            capture_output=True,
            text=True,
            check=True,
        )
        return (result.stdout.strip() or result.stderr.strip()).splitlines()[0]
    except Exception:
        return "unknown"


def ensure_tools() -> None:
    for executable in ["mmseqs", "biosieve"]:
        if shutil.which(executable) is None:
            raise RuntimeError(
                f"'{executable}' not found in PATH. "
                "Activate/install the environment containing BioSieve and MMseqs2."
            )


def read_dataset(path: Path) -> pd.DataFrame:
    return pd.read_csv(
        path,
        sep="\t",
        compression="infer",
        dtype=str,
        low_memory=False,
    ).fillna("")


def resolve_sequence_column(df: pd.DataFrame) -> str:
    for col in ["sequence", "sequence_norm", "seq"]:
        if col in df.columns:
            return col
    raise KeyError(
        "No sequence column found. Expected one of: sequence, sequence_norm, seq."
    )


def resolve_entity_column(df: pd.DataFrame) -> str:
    for col in ["entity_id", "sequence_uid", "sequence_sha256"]:
        if col in df.columns:
            return col

    # Stage 07 should normally preserve entity_id. For robustness, derive a
    # stable exact-sequence SHA256 identifier only when no upstream ID exists.
    seq_col = resolve_sequence_column(df)
    df["sequence_uid"] = (
        df[seq_col]
        .astype(str)
        .str.replace(r"\s+", "", regex=True)
        .str.upper()
        .map(lambda s: hashlib.sha256(s.encode("utf-8")).hexdigest())
    )
    return "sequence_uid"


def validate_input(df: pd.DataFrame) -> tuple[str, str]:
    if df.empty:
        raise ValueError("Stage-07 dataset is empty.")

    seq_col = resolve_sequence_column(df)
    entity_col = resolve_entity_column(df)

    if df[seq_col].astype(str).eq("").any():
        raise ValueError("Empty sequences detected.")
    if df[entity_col].astype(str).eq("").any():
        raise ValueError(f"Empty {entity_col} values detected.")
    if df[entity_col].duplicated().any():
        raise ValueError(f"Duplicate {entity_col} values detected.")

    return seq_col, entity_col


def write_fasta(df: pd.DataFrame, seq_col: str, path: Path) -> None:
    # Avoid itertuples() here because pandas renames columns beginning with
    # an underscore (e.g. "_split_id"), so getattr(row, "_split_id") fails.
    with path.open("w", encoding="utf-8") as handle:
        for split_id, sequence in zip(df["_split_id"], df[seq_col]):
            sequence = str(sequence).replace(" ", "").upper()
            handle.write(f">{split_id}\n{sequence}\n")


def run_logged(command: list[str], log_handle) -> None:
    print("  $ " + " ".join(command), flush=True)
    completed = subprocess.run(
        command,
        stdout=log_handle,
        stderr=subprocess.STDOUT,
        text=True,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            f"Command failed ({completed.returncode}): {' '.join(command)}"
        )


def build_mmseqs_clusters(
    df: pd.DataFrame,
    *,
    seq_col: str,
    identity: float,
    threads: int,
    tmp_parent: Path | None,
    output_dir: Path,
) -> pd.DataFrame:
    raw_membership = output_dir / "mmseqs_clusters_raw.tsv"
    mmseqs_log = output_dir / "mmseqs.log"

    with tempfile.TemporaryDirectory(
        prefix=f"phagelys_split_s{int(identity * 100)}_",
        dir=str(tmp_parent) if tmp_parent else None,
    ) as td:
        td = Path(td)
        fasta = td / "input.fasta"
        seqdb = td / "sequence_db"
        cludb = td / "cluster_db"
        mm_tmp = td / "tmp"
        mm_tmp.mkdir(parents=True, exist_ok=True)

        write_fasta(df, seq_col, fasta)

        with mmseqs_log.open("w", encoding="utf-8") as log:
            run_logged(
                ["mmseqs", "createdb", str(fasta), str(seqdb)],
                log,
            )
            run_logged(
                [
                    "mmseqs",
                    "cluster",
                    str(seqdb),
                    str(cludb),
                    str(mm_tmp),
                    "--min-seq-id",
                    str(identity),
                    "-c",
                    str(COVERAGE),
                    "--cov-mode",
                    str(COV_MODE),
                    "--cluster-mode",
                    str(CLUSTER_MODE),
                    "--threads",
                    str(threads),
                ],
                log,
            )
            run_logged(
                [
                    "mmseqs",
                    "createtsv",
                    str(seqdb),
                    str(seqdb),
                    str(cludb),
                    str(raw_membership),
                ],
                log,
            )

    clusters = pd.read_csv(
        raw_membership,
        sep="\t",
        header=None,
        names=["cluster_representative", "_split_id"],
        dtype=str,
    )

    if clusters.empty:
        raise RuntimeError("MMseqs2 produced an empty cluster-membership table.")

    expected = set(df["_split_id"].astype(str))
    observed = set(clusters["_split_id"].astype(str))

    missing = expected - observed
    extra = observed - expected

    if missing or extra:
        raise RuntimeError(
            f"MMseqs2 membership mismatch: missing={len(missing)}, extra={len(extra)}"
        )

    if clusters["_split_id"].duplicated().any():
        raise RuntimeError("A sequence belongs to more than one MMseqs2 cluster.")

    # Use representative ID as a stable cluster label.
    clusters["cluster_id"] = clusters["cluster_representative"].astype(str)

    return clusters[["_split_id", "cluster_id", "cluster_representative"]]


def write_biosieve_params(clusters_csv: Path, params_path: Path) -> None:
    # BioSieve expects the input identifier column to be named "id".
    text = f"""homology_aware:
  mode: precomputed
  clusters_path: {json.dumps(str(clusters_csv))}
  member_col: id
  cluster_col: cluster_id
  test_size: {TEST_SIZE}
"""
    params_path.write_text(text, encoding="utf-8")


def discover_biosieve_assignment(
    biosieve_dir: Path,
    expected_ids: set[str],
) -> pd.DataFrame:
    """
    Recover train/test IDs from BioSieve outputs without depending on a single
    filename convention. BioSieve uses "id" as its canonical identifier column;
    internally we normalize it back to "_split_id".
    """
    csv_files = sorted(biosieve_dir.rglob("*.csv"))

    def resolve_id_column(frame: pd.DataFrame) -> str | None:
        if "id" in frame.columns:
            return "id"
        if "_split_id" in frame.columns:
            return "_split_id"
        return None

    # 1) Separate train/test CSV files.
    train_candidates = [
        p for p in csv_files
        if "train" in p.stem.lower() and "report" not in p.stem.lower()
    ]
    test_candidates = [
        p for p in csv_files
        if "test" in p.stem.lower() and "report" not in p.stem.lower()
    ]

    for train_path in train_candidates:
        train = pd.read_csv(train_path, dtype=str).fillna("")
        train_id_col = resolve_id_column(train)
        if train_id_col is None:
            continue

        for test_path in test_candidates:
            test = pd.read_csv(test_path, dtype=str).fillna("")
            test_id_col = resolve_id_column(test)
            if test_id_col is None:
                continue

            train_ids = set(train[train_id_col].astype(str))
            test_ids = set(test[test_id_col].astype(str))

            if (
                train_ids.isdisjoint(test_ids)
                and train_ids | test_ids == expected_ids
            ):
                return pd.concat(
                    [
                        pd.DataFrame({
                            "_split_id": sorted(train_ids),
                            "split": "train",
                        }),
                        pd.DataFrame({
                            "_split_id": sorted(test_ids),
                            "split": "test",
                        }),
                    ],
                    ignore_index=True,
                )

    # 2) One combined file containing a split/partition column.
    for path in csv_files:
        frame = pd.read_csv(path, dtype=str).fillna("")
        id_col = resolve_id_column(frame)
        if id_col is None:
            continue

        split_col = next(
            (
                c for c in ["split", "partition", "set", "subset"]
                if c in frame.columns
            ),
            None,
        )
        if split_col is None:
            continue

        normalized = frame[split_col].astype(str).str.strip().str.lower()
        mapping = {
            "train": "train",
            "training": "train",
            "test": "test",
            "testing": "test",
        }
        split = normalized.map(mapping)

        if split.isna().any():
            continue

        candidate = pd.DataFrame({
            "_split_id": frame[id_col].astype(str),
            "split": split,
        })

        if (
            set(candidate["_split_id"]) == expected_ids
            and not candidate["_split_id"].duplicated().any()
            and set(candidate["split"]) == {"train", "test"}
        ):
            return candidate

    raise RuntimeError(
        f"Could not recover a complete train/test assignment from BioSieve outputs under {biosieve_dir}"
    )


def enforce_cluster_integrity(
    assignment: pd.DataFrame,
    clusters: pd.DataFrame,
    *,
    output_dir: Path,
    target_test_fraction: float = TEST_SIZE,
) -> pd.DataFrame:
    """
    BioSieve 0.1.2 can return the requested row-level test size while splitting
    a small number of precomputed homology clusters across train/test.

    This deterministic post-processing step repairs ONLY those crossing
    clusters. Each affected cluster is assigned wholly to train or test using
    the choice that keeps the global test size closest to the requested target.
    Ties are resolved by the majority BioSieve assignment, then to train.

    The biological entities and cluster definitions are never changed.
    """
    merged = assignment.merge(
        clusters[["_split_id", "cluster_id"]],
        on="_split_id",
        how="left",
        validate="1:1",
    )

    if merged["cluster_id"].isna().any():
        raise RuntimeError("Missing cluster IDs during cluster-integrity enforcement.")

    per_cluster = (
        merged
        .groupby(["cluster_id", "split"])
        .size()
        .unstack(fill_value=0)
    )

    for col in ["train", "test"]:
        if col not in per_cluster.columns:
            per_cluster[col] = 0

    shared = per_cluster.loc[
        (per_cluster["train"] > 0) & (per_cluster["test"] > 0)
    ].copy()

    repair_path = output_dir / "biosieve_cluster_integrity_repair.tsv"

    if shared.empty:
        pd.DataFrame(
            columns=[
                "cluster_id",
                "cluster_size",
                "initial_train",
                "initial_test",
                "final_split",
                "rows_reassigned",
            ]
        ).to_csv(repair_path, sep="\t", index=False)
        return assignment

    repaired = assignment.copy()
    target_test_n = target_test_fraction * len(repaired)
    current_test_n = int(repaired["split"].eq("test").sum())

    # Largest crossing clusters first; stable cluster_id ordering breaks ties.
    shared = shared.assign(
        cluster_size=shared["train"] + shared["test"]
    ).sort_values(
        ["cluster_size", "cluster_id"],
        ascending=[False, True],
        kind="stable",
    )

    repair_rows = []

    for cluster_id, row in shared.iterrows():
        initial_train = int(row["train"])
        initial_test = int(row["test"])
        cluster_size = initial_train + initial_test

        # Remove this cluster's current test contribution before comparing the
        # two legal whole-cluster assignments.
        test_without_cluster = current_test_n - initial_test
        test_if_train = test_without_cluster
        test_if_test = test_without_cluster + cluster_size

        distance_train = abs(test_if_train - target_test_n)
        distance_test = abs(test_if_test - target_test_n)

        if distance_test < distance_train:
            final_split = "test"
        elif distance_train < distance_test:
            final_split = "train"
        else:
            # Preserve as much as possible of BioSieve's initial proposal.
            if initial_test > initial_train:
                final_split = "test"
            else:
                final_split = "train"

        member_ids = set(
            merged.loc[
                merged["cluster_id"].eq(cluster_id),
                "_split_id",
            ].astype(str)
        )
        mask = repaired["_split_id"].astype(str).isin(member_ids)

        old = repaired.loc[mask, "split"].astype(str)
        rows_reassigned = int((old != final_split).sum())
        repaired.loc[mask, "split"] = final_split

        current_test_n = (
            test_if_test if final_split == "test" else test_if_train
        )

        repair_rows.append({
            "cluster_id": cluster_id,
            "cluster_size": cluster_size,
            "initial_train": initial_train,
            "initial_test": initial_test,
            "final_split": final_split,
            "rows_reassigned": rows_reassigned,
        })

    repair_df = pd.DataFrame(repair_rows)
    repair_df.to_csv(repair_path, sep="\t", index=False)

    # Hard closure audit.
    check = repaired.merge(
        clusters[["_split_id", "cluster_id"]],
        on="_split_id",
        how="left",
        validate="1:1",
    )
    crossing_after = (
        check.groupby("cluster_id")["split"].nunique().gt(1).sum()
    )
    if crossing_after:
        raise RuntimeError(
            f"Cluster-integrity repair failed: {crossing_after} clusters still cross train/test."
        )

    print(
        f"  BioSieve cluster-integrity repair: "
        f"{len(repair_df)} crossing clusters repaired; "
        f"{int(repair_df['rows_reassigned'].sum())} rows reassigned; "
        f"final test={int(repaired['split'].eq('test').sum())}/{len(repaired)} "
        f"({repaired['split'].eq('test').mean():.4f})"
    )

    return repaired

def run_biosieve_split(
    df: pd.DataFrame,
    clusters: pd.DataFrame,
    *,
    output_dir: Path,
) -> pd.DataFrame:
    biosieve_dir = output_dir / "biosieve_report"
    biosieve_dir.mkdir(parents=True, exist_ok=True)

    input_csv = output_dir / "_biosieve_input.csv"
    clusters_csv = output_dir / "_biosieve_clusters.tsv"
    params_yaml = output_dir / "_biosieve_params.yaml"

    # BioSieve requires an input column literally named "id".
    # Keep the biological/internal columns untouched and add "id" only to the
    # temporary interchange files used by BioSieve.
    biosieve_input = df.copy()
    biosieve_input["id"] = biosieve_input["_split_id"]
    biosieve_input.to_csv(input_csv, index=False)

    biosieve_clusters = (
        clusters[["_split_id", "cluster_id"]]
        .rename(columns={"_split_id": "id"})
    )
    # BioSieve 0.1.2 reads precomputed cluster files as a two-column
    # tab-separated table. Writing comma-separated data makes pandas see one
    # column and triggers: "2 column names provided for a DataFrame of width 1".
    biosieve_clusters.to_csv(
        clusters_csv,
        sep="\t",
        index=False,
    )

    write_biosieve_params(clusters_csv, params_yaml)

    command = [
        "biosieve",
        "split",
        "-i",
        str(input_csv),
        "-o",
        str(biosieve_dir),
        "--strategy",
        "homology_aware",
        "--params",
        str(params_yaml),
    ]

    print("  $ " + " ".join(command), flush=True)
    completed = subprocess.run(command)
    if completed.returncode != 0:
        raise RuntimeError(
            f"BioSieve split failed with exit code {completed.returncode}."
        )

    assignment = discover_biosieve_assignment(
        biosieve_dir,
        expected_ids=set(df["_split_id"].astype(str)),
    )

    # BioSieve 0.1.2 may satisfy the requested row count by splitting a small
    # number of precomputed clusters. Repair those crossings deterministically
    # before any canonical split is materialized.
    assignment = enforce_cluster_integrity(
        assignment,
        clusters,
        output_dir=output_dir,
        target_test_fraction=TEST_SIZE,
    )

    # Clean only temporary interchange files; persistent BioSieve outputs stay.
    input_csv.unlink(missing_ok=True)
    clusters_csv.unlink(missing_ok=True)
    params_yaml.unlink(missing_ok=True)

    return assignment


def composition(df: pd.DataFrame, column: str) -> str:
    if column not in df.columns:
        return ""
    counts = df[column].astype(str).value_counts(dropna=False)
    return "; ".join(
        f"{k if str(k).strip() else 'NA'}={int(v)}"
        for k, v in counts.items()
    )


def split_complete(output_dir: Path, input_sha256: str) -> bool:
    manifest_path = output_dir / "manifest.json"
    required = [
        output_dir / "train.tsv.gz",
        output_dir / "test.tsv.gz",
        output_dir / "dataset_with_split.tsv.gz",
        output_dir / "clusters.tsv.gz",
        output_dir / "split_summary.tsv",
        manifest_path,
    ]
    if not all(p.is_file() for p in required):
        return False

    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        return (
            manifest.get("status") == "PASS"
            and manifest.get("input_sha256") == input_sha256
        )
    except Exception:
        return False


def run_one(
    *,
    task: dict[str, Any],
    r: int,
    s: int,
    redundancy_root: Path,
    splits_root: Path,
    threads: int,
    tmp_parent: Path | None,
    overwrite: bool,
    mmseqs_version: str,
    biosieve_version: str,
) -> dict[str, Any]:
    input_path = (
        redundancy_root
        / task["task_path"]
        / f"H{r}"
        / "dataset.tsv.gz"
    )
    if not input_path.is_file():
        raise FileNotFoundError(input_path)

    input_sha = sha256_file(input_path)
    rs_id = f"R{r}_S{s}"
    output_dir = splits_root / task["task_path"] / rs_id
    output_dir.mkdir(parents=True, exist_ok=True)

    if not overwrite and split_complete(output_dir, input_sha):
        print(f"SKIP {task['task_id']} {rs_id} — already PASS")
        summary = pd.read_csv(output_dir / "split_summary.tsv", sep="\t").iloc[0]
        return summary.to_dict()

    # Remove old output products for this cell only.
    for name in [
        "train.tsv.gz",
        "test.tsv.gz",
        "dataset_with_split.tsv.gz",
        "clusters.tsv.gz",
        "split_summary.tsv",
        "manifest.json",
        "mmseqs_clusters_raw.tsv",
        "mmseqs.log",
        "biosieve_cluster_integrity_repair.tsv",
    ]:
        (output_dir / name).unlink(missing_ok=True)

    biosieve_dir = output_dir / "biosieve_report"
    if biosieve_dir.exists():
        shutil.rmtree(biosieve_dir)

    df = read_dataset(input_path)
    seq_col, entity_col = validate_input(df)

    # Stable internal split IDs protect against spaces/special characters in
    # biological IDs and allow transparent mapping back to the source dataset.
    df = df.reset_index(drop=True)
    df["_split_id"] = [
        f"SEQ_{i:09d}" for i in range(1, len(df) + 1)
    ]

    clusters = build_mmseqs_clusters(
        df,
        seq_col=seq_col,
        identity=s / 100.0,
        threads=threads,
        tmp_parent=tmp_parent,
        output_dir=output_dir,
    )

    assignment = run_biosieve_split(
        df,
        clusters,
        output_dir=output_dir,
    )

    work = (
        df
        .merge(clusters, on="_split_id", how="left", validate="1:1")
        .merge(assignment, on="_split_id", how="left", validate="1:1")
    )

    if work["split"].isna().any():
        raise RuntimeError("Unassigned sequences detected.")
    if set(work["split"]) != {"train", "test"}:
        raise RuntimeError(
            f"Expected train/test only, observed: {sorted(work['split'].unique())}"
        )

    train = work.loc[work["split"].eq("train")].copy()
    test = work.loc[work["split"].eq("test")].copy()

    # Critical audits.
    entity_overlap = (
        set(train[entity_col].astype(str))
        & set(test[entity_col].astype(str))
    )
    shared_clusters = (
        set(train["cluster_id"].astype(str))
        & set(test["cluster_id"].astype(str))
    )

    if entity_overlap:
        raise RuntimeError(
            f"Exact-entity leakage: {len(entity_overlap)} shared IDs."
        )
    if shared_clusters:
        raise RuntimeError(
            f"Homology-cluster leakage: {len(shared_clusters)} clusters cross train/test."
        )
    if len(train) + len(test) != len(df):
        raise RuntimeError("train + test count does not equal input count.")

    # Drop private interchange ID from canonical datasets.
    canonical_cols = [
        c for c in work.columns
        if c != "_split_id"
    ]
    work[canonical_cols].to_csv(
        output_dir / "dataset_with_split.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    train[canonical_cols].to_csv(
        output_dir / "train.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    test[canonical_cols].to_csv(
        output_dir / "test.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    cluster_export = (
        work[
            [entity_col, "cluster_id", "cluster_representative"]
        ]
        .copy()
        .rename(columns={entity_col: "entity_id"})
    )
    cluster_export.to_csv(
        output_dir / "clusters.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    test_fraction = len(test) / len(work) if len(work) else 0.0
    warnings = []
    if abs(test_fraction - TEST_SIZE) > 0.05:
        warnings.append(
            f"test_fraction={test_fraction:.4f} differs from target={TEST_SIZE:.2f} by >0.05"
        )

    repair_path = output_dir / "biosieve_cluster_integrity_repair.tsv"
    if repair_path.is_file():
        repair_df = pd.read_csv(repair_path, sep="\t")
        n_biosieve_crossing_clusters = int(len(repair_df))
        n_rows_reassigned_for_integrity = (
            int(repair_df["rows_reassigned"].sum())
            if len(repair_df) and "rows_reassigned" in repair_df.columns
            else 0
        )
    else:
        n_biosieve_crossing_clusters = 0
        n_rows_reassigned_for_integrity = 0

    summary = {
        "task_id": task["task_id"],
        "task_path": task["task_path"],
        "redundancy_threshold": r,
        "split_threshold": s,
        "rs_id": rs_id,
        "n_total": len(work),
        "n_train": len(train),
        "n_test": len(test),
        "actual_test_fraction": test_fraction,
        "n_homology_clusters": work["cluster_id"].nunique(),
        "largest_cluster_size": int(work["cluster_id"].value_counts().max()),
        "shared_train_test_clusters": len(shared_clusters),
        "train_test_entity_overlap": len(entity_overlap),
        "biosieve_crossing_clusters_before_repair": n_biosieve_crossing_clusters,
        "rows_reassigned_for_cluster_integrity": n_rows_reassigned_for_integrity,
        "primary_group_column": task["primary_group"],
        "train_primary_composition": composition(train, task["primary_group"]),
        "test_primary_composition": composition(test, task["primary_group"]),
        "secondary_group_column": task["secondary_group"],
        "train_secondary_composition": composition(train, task["secondary_group"]),
        "test_secondary_composition": composition(test, task["secondary_group"]),
        "n_warnings": len(warnings),
        "status": "PASS",
    }

    pd.DataFrame([summary]).to_csv(
        output_dir / "split_summary.tsv",
        sep="\t",
        index=False,
    )

    manifest = {
        "stage": "08-factorial",
        "status": "PASS",
        "task_id": task["task_id"],
        "task_path": task["task_path"],
        "redundancy_threshold": r,
        "split_threshold": s,
        "rs_id": rs_id,
        "input_dataset": str(input_path),
        "input_sha256": input_sha,
        "entity_column": entity_col,
        "sequence_column": seq_col,
        "split_design": {
            "outer_train_target": 1.0 - TEST_SIZE,
            "outer_test_target": TEST_SIZE,
            "validation_size": 0.0,
            "whole_task_clustering": True,
            "challenge_sets_untouched": True,
            "representation_independent": True,
        },
        "mmseqs2": {
            "version": mmseqs_version,
            "min_seq_id": s / 100.0,
            "coverage": COVERAGE,
            "cov_mode": COV_MODE,
            "cluster_mode": CLUSTER_MODE,
            "threads": threads,
        },
        "biosieve": {
            "version": biosieve_version,
            "strategy": "homology_aware",
            "mode": "precomputed",
            "test_size": TEST_SIZE,
            "cluster_integrity_postprocessing": {
                "enabled": True,
                "reason": (
                    "BioSieve 0.1.2 may split a small number of precomputed "
                    "clusters to satisfy the requested row-level test size."
                ),
                "crossing_clusters_before_repair": n_biosieve_crossing_clusters,
                "rows_reassigned": n_rows_reassigned_for_integrity,
                "rule": (
                    "Assign each crossing cluster wholly to the partition that "
                    "minimizes deviation from target test size; majority "
                    "BioSieve assignment breaks exact ties."
                ),
            },
        },
        "audit": {
            "n_total": len(work),
            "n_train": len(train),
            "n_test": len(test),
            "actual_test_fraction": test_fraction,
            "n_homology_clusters": work["cluster_id"].nunique(),
            "shared_train_test_clusters": len(shared_clusters),
            "train_test_entity_overlap": len(entity_overlap),
            "warnings": warnings,
        },
        "outputs": {
            "train": str(output_dir / "train.tsv.gz"),
            "test": str(output_dir / "test.tsv.gz"),
            "dataset_with_split": str(output_dir / "dataset_with_split.tsv.gz"),
            "clusters": str(output_dir / "clusters.tsv.gz"),
            "summary": str(output_dir / "split_summary.tsv"),
        },
    }
    (output_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n",
        encoding="utf-8",
    )

    print(
        f"PASS {task['task_id']} {rs_id} | "
        f"n={len(work):,} train={len(train):,} test={len(test):,} "
        f"clusters={work['cluster_id'].nunique():,}"
    )
    return summary


def main() -> None:
    args = parse_args()
    ensure_tools()

    redundancy_root = Path(args.redundancy_root).expanduser().resolve()
    splits_root = Path(args.splits_root).expanduser().resolve()
    splits_root.mkdir(parents=True, exist_ok=True)

    tmp_parent = (
        Path(args.tmp_dir).expanduser().resolve()
        if args.tmp_dir else None
    )
    if tmp_parent:
        tmp_parent.mkdir(parents=True, exist_ok=True)

    selected_tasks = [
        t for t in TASKS
        if not args.only_task or t["task_id"] in args.only_task
    ]
    r_values = args.only_r or R_THRESHOLDS
    s_values = args.only_s or S_THRESHOLDS

    mmseqs_version = command_version(["mmseqs", "version"])
    biosieve_version = package_version("biosieve")

    print("=" * 100)
    print("PhageLysClass — Stage 08 factorial homology-aware splitting")
    print("=" * 100)
    print(f"Redundancy root : {redundancy_root}")
    print(f"Splits root     : {splits_root}")
    print(f"Tasks           : {len(selected_tasks)}")
    print(f"R thresholds    : {r_values}")
    print(f"S thresholds    : {s_values}")
    print(f"MMseqs2         : {mmseqs_version}")
    print(f"BioSieve        : {biosieve_version}")
    print(f"Threads         : {args.threads}")
    print()

    all_rows = []

    for task in selected_tasks:
        task_rows = []

        for r in r_values:
            for s in s_values:
                row = run_one(
                    task=task,
                    r=r,
                    s=s,
                    redundancy_root=redundancy_root,
                    splits_root=splits_root,
                    threads=args.threads,
                    tmp_parent=tmp_parent,
                    overwrite=args.overwrite,
                    mmseqs_version=mmseqs_version,
                    biosieve_version=biosieve_version,
                )
                task_rows.append(row)
                all_rows.append(row)

        task_out = splits_root / task["task_path"]
        pd.DataFrame(task_rows).sort_values(
            ["redundancy_threshold", "split_threshold"],
            ascending=[False, False],
        ).to_csv(
            task_out / "factorial_split_overview.tsv",
            sep="\t",
            index=False,
        )

    global_overview = pd.DataFrame(all_rows)
    if not global_overview.empty:
        global_overview.to_csv(
            splits_root / "factorial_split_overview.tsv",
            sep="\t",
            index=False,
        )

    expected = len(selected_tasks) * len(r_values) * len(s_values)
    passed = int((global_overview["status"] == "PASS").sum()) if not global_overview.empty else 0

    print()
    print("=" * 100)
    print(f"Stage 08 factorial cells PASS: {passed}/{expected}")
    print("=" * 100)

    if passed != expected:
        raise SystemExit(2)

    print("PASS — requested factorial homology-aware splits are complete.")


if __name__ == "__main__":
    main()

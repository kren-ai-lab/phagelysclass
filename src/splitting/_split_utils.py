#!/usr/bin/env python3
"""Shared utilities for Stage 08 homology-aware train/test splitting."""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import os
import re
import shutil
import subprocess
from pathlib import Path

import pandas as pd


DEFAULT_THRESHOLDS = [90, 70, 50, 30]


def project_root_from_script(script_file: str) -> Path:
    script_dir = Path(script_file).resolve().parent
    return script_dir.parents[1]


def package_version(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except Exception:
        return "unknown"


def command_version(executable: str, args: list[str]) -> str:
    try:
        proc = subprocess.run(
            [executable, *args],
            capture_output=True,
            text=True,
            check=False,
        )
        text = (proc.stdout or proc.stderr or "").strip()
        return text.splitlines()[0] if text else "unknown"
    except Exception:
        return "unknown"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def add_split_args(
    parser: argparse.ArgumentParser,
    *,
    script_file: str,
    default_input_root: str | None = None,
    default_output_root: str | None = None,
) -> argparse.ArgumentParser:
    root = project_root_from_script(script_file)

    if default_input_root is not None:
        parser.add_argument(
            "--input-root",
            default=str(root / default_input_root),
            help="Stage 07 redundancy root for this task.",
        )

    if default_output_root is not None:
        parser.add_argument(
            "--output-root",
            default=str(root / default_output_root),
            help="Stage 08 output root.",
        )

    parser.add_argument(
        "--thresholds",
        nargs="+",
        type=int,
        default=DEFAULT_THRESHOLDS,
        help="Redundancy/homology regimes to split (default: 90 70 50 30).",
    )
    parser.add_argument(
        "--test-size",
        type=float,
        default=0.30,
        help="Target external test fraction (default: 0.30).",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Requested BioSieve split random state (default: 42).",
    )
    parser.add_argument(
        "--coverage",
        type=float,
        default=0.80,
        help="MMseqs2 alignment coverage used to build split clusters.",
    )
    parser.add_argument(
        "--cov-mode",
        type=int,
        default=0,
        help="MMseqs2 coverage mode (default: 0, query and target).",
    )
    parser.add_argument(
        "--cluster-mode",
        type=int,
        default=1,
        help=(
            "MMseqs2 cluster mode. Default 1 = connected components, chosen "
            "for leakage-aware cluster construction."
        ),
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=max(1, (os.cpu_count() or 2) // 2),
    )
    parser.add_argument(
        "--reuse-existing",
        action="store_true",
        help="Reuse a completed threshold split if manifest.json exists.",
    )
    parser.add_argument(
        "--keep-work",
        action="store_true",
        help="Keep FASTA/MMseqs2/BioSieve temporary work files.",
    )
    return parser


def preflight() -> dict[str, str]:
    biosieve = shutil.which("biosieve")
    mmseqs = shutil.which("mmseqs")

    if biosieve is None:
        raise RuntimeError("BioSieve CLI not found in PATH.")
    if mmseqs is None:
        raise RuntimeError("MMseqs2 binary ('mmseqs') not found in PATH.")

    return {
        "biosieve_executable": biosieve,
        "biosieve_version": package_version("biosieve"),
        "mmseqs_executable": mmseqs,
        "mmseqs_version": command_version(mmseqs, ["version"]),
    }


def load_dataset(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"Stage 07 dataset not found: {path}")

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
            "Dataset missing columns: " + ", ".join(sorted(missing))
        )

    if df["entity_id"].duplicated().any():
        raise ValueError("Duplicated entity_id values in input dataset.")
    if df["sequence"].duplicated().any():
        raise ValueError("Exact duplicate sequences remain in Stage 07 dataset.")

    return df


def write_fasta(df: pd.DataFrame, path: Path) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for row in df[["entity_id", "sequence"]].itertuples(index=False):
            entity_id = str(row.entity_id).strip()
            sequence = re.sub(r"\s+", "", str(row.sequence)).upper()
            if not entity_id or not sequence:
                raise ValueError("Empty entity_id or sequence while writing FASTA.")
            handle.write(f">{entity_id}\n{sequence}\n")


def build_connected_components(
    *,
    df: pd.DataFrame,
    threshold_percent: int,
    coverage: float,
    cov_mode: int,
    cluster_mode: int,
    threads: int,
    mmseqs_executable: str,
    work_dir: Path,
) -> pd.DataFrame:
    """
    Build full-task MMseqs2 clusters.

    Connected components are used by default so any detected homology edge at
    the requested criterion belongs to one indivisible split group.
    """
    work_dir.mkdir(parents=True, exist_ok=True)

    fasta_path = work_dir / "dataset.fasta"
    prefix = work_dir / "homology"
    tmp_dir = work_dir / "mmseqs_tmp"

    write_fasta(df, fasta_path)

    cmd = [
        mmseqs_executable,
        "easy-cluster",
        str(fasta_path),
        str(prefix),
        str(tmp_dir),
        "--min-seq-id", str(threshold_percent / 100.0),
        "-c", str(coverage),
        "--cov-mode", str(cov_mode),
        "--cluster-mode", str(cluster_mode),
        "--single-step-clustering", "1",
        "--threads", str(threads),
    ]

    print(
        f"    MMseqs2 components: identity={threshold_percent}% "
        f"coverage={coverage:.2f} cov_mode={cov_mode} "
        f"cluster_mode={cluster_mode}"
    )
    proc = subprocess.run(cmd, check=False)
    if proc.returncode != 0:
        raise RuntimeError(
            f"MMseqs2 clustering failed at H{threshold_percent}."
        )

    cluster_tsv = Path(str(prefix) + "_cluster.tsv")
    if not cluster_tsv.is_file():
        raise FileNotFoundError(
            f"Expected MMseqs2 cluster table not found: {cluster_tsv}"
        )

    raw = pd.read_csv(
        cluster_tsv,
        sep="\t",
        header=None,
        names=["representative_id", "id"],
        dtype=str,
    )

    # MMseqs2 normally emits every member, including representatives.
    # Still, explicitly add missing singletons for a complete contract.
    assigned = set(raw["id"].astype(str))
    all_ids = set(df["entity_id"].astype(str))
    missing = sorted(all_ids - assigned)

    if missing:
        extra = pd.DataFrame(
            {
                "representative_id": missing,
                "id": missing,
            }
        )
        raw = pd.concat([raw, extra], ignore_index=True)

    unknown = set(raw["id"].astype(str)) - all_ids
    if unknown:
        raise RuntimeError(
            f"MMseqs2 cluster output contains {len(unknown):,} unknown IDs."
        )

    if raw["id"].duplicated().any():
        # For a valid hard partition each member must belong to one component.
        dup = int(raw["id"].duplicated().sum())
        raise RuntimeError(
            f"MMseqs2 component table has {dup:,} duplicated members."
        )

    clusters = raw.rename(
        columns={"representative_id": "cluster_id"}
    )[["id", "cluster_id"]].copy()

    return clusters


def write_biosieve_params(
    *,
    path: Path,
    clusters_path: Path,
    test_size: float,
    seed: int,
) -> None:
    """
    Write BioSieve HomologyAwareSplitter parameters.

    JSON is valid YAML, so the params file is serialized with json.dumps()
    instead of hand-built YAML text. This avoids escaping/newline/path parsing
    issues while preserving the BioSieve `--params` interface.
    """
    payload = {
        "homology_aware": {
            "mode": "precomputed",
            "clusters_path": str(clusters_path),
            "clusters_format": "csv",
            "member_col": "id",
            "cluster_col": "cluster_id",
            "test_size": float(test_size),
            "val_size": 0.0,
            "seed": int(seed),
        }
    }

    path.write_text(
        json.dumps(payload, indent=2) + "\n",
        encoding="utf-8",
    )

def _extract_ids_from_split_dir(split_dir: Path) -> dict[str, set[str]]:
    """
    Discover BioSieve outputs without assuming one exact file naming scheme.
    """
    result: dict[str, set[str]] = {}

    csv_files = sorted(split_dir.rglob("*.csv"))

    # First preference: separate train/test/val files.
    for path in csv_files:
        stem = path.stem.lower()
        split_name = None
        if "train" in stem:
            split_name = "train"
        elif "test" in stem:
            split_name = "test"
        elif "val" in stem or "valid" in stem:
            split_name = "val"

        if split_name is None:
            continue

        try:
            df = pd.read_csv(path, dtype=str).fillna("")
        except Exception:
            continue

        id_col = "id" if "id" in df.columns else (
            "entity_id" if "entity_id" in df.columns else None
        )
        if id_col:
            result.setdefault(split_name, set()).update(
                df[id_col].astype(str)
            )

    if "train" in result and "test" in result:
        return result

    # Fallback: one table with a split/partition column.
    for path in csv_files:
        try:
            df = pd.read_csv(path, dtype=str).fillna("")
        except Exception:
            continue

        id_col = "id" if "id" in df.columns else (
            "entity_id" if "entity_id" in df.columns else None
        )
        split_col = next(
            (
                col for col in ["split", "partition", "set"]
                if col in df.columns
            ),
            None,
        )
        if id_col and split_col:
            for value, block in df.groupby(split_col):
                name = str(value).strip().lower()
                if name.startswith("train"):
                    key = "train"
                elif name.startswith("test"):
                    key = "test"
                elif name.startswith("val"):
                    key = "val"
                else:
                    continue
                result.setdefault(key, set()).update(
                    block[id_col].astype(str)
                )

    return result


def run_biosieve_homology_split(
    *,
    df: pd.DataFrame,
    clusters: pd.DataFrame,
    test_size: float,
    seed: int,
    biosieve_executable: str,
    work_dir: Path,
) -> tuple[set[str], set[str], dict]:
    work_dir.mkdir(parents=True, exist_ok=True)

    adapter_path = work_dir / "dataset.csv"
    clusters_path = work_dir / "clusters.csv"
    params_path = work_dir / "params.yaml"
    split_dir = work_dir / "splits"

    # BioSieve only needs stable IDs/sequences for the split. Original task
    # metadata are restored after split membership is resolved.
    df[["entity_id", "sequence"]].rename(
        columns={"entity_id": "id"}
    ).to_csv(adapter_path, index=False)

    # Headered CSV expected by the explicit `clusters_format: csv` contract:
    # id,cluster_id
    clusters[["id", "cluster_id"]].to_csv(
        clusters_path,
        index=False,
    )

    write_biosieve_params(
        path=params_path,
        clusters_path=clusters_path.resolve(),
        test_size=test_size,
        seed=seed,
    )

    # JSON is valid YAML; ensure the file is syntactically valid and has the
    # exact strategy block expected before invoking BioSieve.
    try:
        parsed_params = json.loads(params_path.read_text(encoding="utf-8"))
    except Exception as exc:
        raise RuntimeError(
            f"Generated BioSieve params file is invalid: {params_path}"
        ) from exc

    if "homology_aware" not in parsed_params:
        raise RuntimeError(
            "Generated BioSieve params file lacks 'homology_aware' block."
        )

    if split_dir.exists():
        shutil.rmtree(split_dir)

    cmd = [
        biosieve_executable,
        "split",
        "-i", str(adapter_path),
        "-o", str(split_dir),
        "--strategy", "homology_aware",
        "--params", str(params_path),
    ]

    print(
        f"    BioSieve homology_aware split "
        f"(test_size={test_size:.2f}, seed={seed}, val_size=0)"
    )

    proc = subprocess.run(cmd, check=False)
    if proc.returncode != 0:
        raise RuntimeError(
            "BioSieve homology-aware split failed. The current Stage-08 "
            "adapter writes a headered CSV precomputed cluster table and "
            "passes clusters_format=csv, seed=<seed>, and val_size=0."
        )

    discovered = _extract_ids_from_split_dir(split_dir)
    train_ids = discovered.get("train", set())
    test_ids = discovered.get("test", set())
    val_ids = discovered.get("val", set())

    if not train_ids or not test_ids:
        raise RuntimeError(
            f"Could not discover non-empty train/test outputs in {split_dir}."
        )

    if val_ids:
        raise RuntimeError(
            "BioSieve produced a non-empty validation split, but Stage 08 "
            "is frozen as train/test only."
        )

    report = {
        "biosieve_split_directory": str(split_dir),
        "seed": seed,
        "clusters_format": "csv",
        "val_size": 0.0,
    }

    # Include any JSON output emitted by BioSieve.
    json_payloads = {}
    for json_path in split_dir.rglob("*.json"):
        try:
            json_payloads[json_path.name] = json.loads(
                json_path.read_text(encoding="utf-8")
            )
        except Exception:
            pass
    if json_payloads:
        report["biosieve_json_outputs"] = json_payloads

    return train_ids, test_ids, report

def summarize_split(
    *,
    full: pd.DataFrame,
    train: pd.DataFrame,
    test: pd.DataFrame,
    summary_cols: list[str],
    target_test_size: float,
) -> pd.DataFrame:
    rows = []

    def add_row(section, group, total, train_n, test_n):
        rows.append(
            {
                "section": section,
                "group": group,
                "total": int(total),
                "train": int(train_n),
                "test": int(test_n),
                "test_fraction": (
                    float(test_n / total) if total else 0.0
                ),
                "target_test_fraction": target_test_size,
                "absolute_test_fraction_deviation": (
                    abs((test_n / total) - target_test_size)
                    if total else 0.0
                ),
            }
        )

    add_row("overall", "ALL", len(full), len(train), len(test))

    for col in summary_cols:
        if col not in full.columns:
            continue
        values = list(full[col].astype(str).unique())
        for value in values:
            total_n = int((full[col].astype(str) == value).sum())
            train_n = int((train[col].astype(str) == value).sum())
            test_n = int((test[col].astype(str) == value).sum())
            add_row(col, value, total_n, train_n, test_n)

    return pd.DataFrame(rows)


def split_one_threshold(
    *,
    task_name: str,
    input_path: Path,
    output_dir: Path,
    threshold_percent: int,
    summary_cols: list[str],
    test_size: float,
    seed: int,
    coverage: float,
    cov_mode: int,
    cluster_mode: int,
    threads: int,
    reuse_existing: bool,
    keep_work: bool,
) -> dict:
    info = preflight()

    input_path = input_path.expanduser().resolve()
    output_dir = output_dir.expanduser().resolve()

    manifest_path = output_dir / "manifest.json"
    if reuse_existing and manifest_path.is_file():
        print(f"[H{threshold_percent}] reuse existing: {output_dir}")
        return json.loads(manifest_path.read_text(encoding="utf-8"))

    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    df = load_dataset(input_path)
    work_dir = output_dir / "work"

    clusters = build_connected_components(
        df=df,
        threshold_percent=threshold_percent,
        coverage=coverage,
        cov_mode=cov_mode,
        cluster_mode=cluster_mode,
        threads=threads,
        mmseqs_executable=info["mmseqs_executable"],
        work_dir=work_dir / "mmseqs",
    )

    train_ids, test_ids, biosieve_report = run_biosieve_homology_split(
        df=df,
        clusters=clusters,
        test_size=test_size,
        seed=seed,
        biosieve_executable=info["biosieve_executable"],
        work_dir=work_dir / "biosieve",
    )

    all_ids = set(df["entity_id"].astype(str))
    if train_ids & test_ids:
        raise RuntimeError("Train/test overlap detected.")
    if train_ids | test_ids != all_ids:
        missing = all_ids - (train_ids | test_ids)
        extra = (train_ids | test_ids) - all_ids
        raise RuntimeError(
            f"Split membership mismatch: missing={len(missing):,}, "
            f"extra={len(extra):,}."
        )

    cluster_map = clusters.set_index("id")["cluster_id"].to_dict()

    train_clusters = {cluster_map[x] for x in train_ids}
    test_clusters = {cluster_map[x] for x in test_ids}
    shared_clusters = train_clusters & test_clusters

    if shared_clusters:
        raise RuntimeError(
            f"Homology leakage: {len(shared_clusters):,} MMseqs2 components "
            "appear in both train and test."
        )

    train = df.loc[df["entity_id"].isin(train_ids)].copy()
    test = df.loc[df["entity_id"].isin(test_ids)].copy()

    train["outer_split"] = "train"
    test["outer_split"] = "test"
    train["homology_cluster_id"] = train["entity_id"].map(cluster_map)
    test["homology_cluster_id"] = test["entity_id"].map(cluster_map)

    train.to_csv(
        output_dir / "train.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    test.to_csv(
        output_dir / "test.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    assignments = pd.concat(
        [
            train[["entity_id", "homology_cluster_id", "outer_split"]],
            test[["entity_id", "homology_cluster_id", "outer_split"]],
        ],
        ignore_index=True,
    )
    assignments.to_csv(
        output_dir / "split_assignments.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    clusters.rename(
        columns={"id": "entity_id", "cluster_id": "homology_cluster_id"}
    ).to_csv(
        output_dir / "homology_clusters.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    summary = summarize_split(
        full=df,
        train=train,
        test=test,
        summary_cols=summary_cols,
        target_test_size=test_size,
    )
    summary.to_csv(
        output_dir / "split_summary.tsv",
        sep="\t",
        index=False,
    )

    cluster_sizes = clusters["cluster_id"].value_counts()
    cluster_summary = pd.DataFrame(
        [
            {
                "n_sequences": len(df),
                "n_clusters": int(cluster_sizes.size),
                "n_singleton_clusters": int((cluster_sizes == 1).sum()),
                "largest_cluster_size": int(cluster_sizes.max()),
                "train_clusters": len(train_clusters),
                "test_clusters": len(test_clusters),
                "shared_train_test_clusters": len(shared_clusters),
            }
        ]
    )
    cluster_summary.to_csv(
        output_dir / "cluster_summary.tsv",
        sep="\t",
        index=False,
    )

    # Warnings are descriptive only; they do not alter the split.
    warning_rows = summary.loc[
        (summary["section"] != "overall")
        & (
            (summary["absolute_test_fraction_deviation"] > 0.15)
            | (summary["test"] < 5)
        )
    ]
    warnings = []
    for row in warning_rows.itertuples(index=False):
        warnings.append(
            f"{row.section}:{row.group} test={row.test}/{row.total} "
            f"({100*row.test_fraction:.1f}%)"
        )

    manifest = {
        "stage": "08",
        "task_name": task_name,
        "redundancy_regime": f"H{threshold_percent}",
        "input": str(input_path),
        "input_sha256": sha256_file(input_path),
        "outer_split": {
            "train_target_fraction": 1.0 - test_size,
            "test_target_fraction": test_size,
            "validation_materialized": False,
            "purpose": (
                "External homology-aware train/test benchmark split. "
                "Any internal random holdout used during model fitting is "
                "created later from the training partition only."
            ),
        },
        "homology_clustering": {
            "engine": "MMseqs2 easy-cluster",
            "min_seq_id": threshold_percent / 100.0,
            "coverage": coverage,
            "cov_mode": cov_mode,
            "cluster_mode": cluster_mode,
            "cluster_mode_interpretation": (
                "connected components" if cluster_mode == 1 else "custom"
            ),
        },
        "split_engine": "BioSieve homology_aware (precomputed clusters)",
        "seed_requested": seed,
        "biosieve": {
            "version": info["biosieve_version"],
            "executable": info["biosieve_executable"],
            **biosieve_report,
        },
        "mmseqs2": {
            "version": info["mmseqs_version"],
            "executable": info["mmseqs_executable"],
        },
        "counts": {
            "total": len(df),
            "train": len(train),
            "test": len(test),
            "actual_test_fraction": len(test) / len(df),
            "n_homology_clusters": int(cluster_sizes.size),
            "shared_train_test_clusters": len(shared_clusters),
        },
        "warnings": warnings,
        "design_note": (
            "Stage 08 freezes only the external train/test partition. "
            "No class balancing, feature computation, hyperparameter tuning, "
            "or internal training holdout is performed here."
        ),
    }
    manifest_path.write_text(
        json.dumps(manifest, indent=2),
        encoding="utf-8",
    )

    if not keep_work:
        shutil.rmtree(work_dir, ignore_errors=True)

    print(
        f"[H{threshold_percent}] total={len(df):,} "
        f"train={len(train):,} test={len(test):,} "
        f"test={100*len(test)/len(df):.2f}% "
        f"clusters={cluster_sizes.size:,} leakage_clusters=0"
    )
    if warnings:
        print("    WARN:", " | ".join(warnings))

    return manifest


def split_task_thresholds(
    *,
    task_name: str,
    input_root: Path,
    output_root: Path,
    thresholds: list[int],
    summary_cols: list[str],
    test_size: float,
    seed: int,
    coverage: float,
    cov_mode: int,
    cluster_mode: int,
    threads: int,
    reuse_existing: bool,
    keep_work: bool,
) -> pd.DataFrame:
    input_root = input_root.expanduser().resolve()
    output_root = output_root.expanduser().resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    if not 0.0 < test_size < 1.0:
        raise ValueError("--test-size must be between 0 and 1.")

    rows = []

    print("=" * 92)
    print(f"PhageLysClass — Stage 08 split: {task_name}")
    print("=" * 92)
    print(f"Input root : {input_root}")
    print(f"Output root: {output_root}")
    print(f"Outer split: {100*(1-test_size):.0f}/{100*test_size:.0f} train/test")
    print(f"Thresholds : {', '.join('H'+str(x) for x in thresholds)}")
    print()

    for threshold in thresholds:
        input_path = input_root / f"H{threshold}" / "dataset.tsv.gz"
        out_dir = output_root / f"H{threshold}"

        manifest = split_one_threshold(
            task_name=task_name,
            input_path=input_path,
            output_dir=out_dir,
            threshold_percent=threshold,
            summary_cols=summary_cols,
            test_size=test_size,
            seed=seed,
            coverage=coverage,
            cov_mode=cov_mode,
            cluster_mode=cluster_mode,
            threads=threads,
            reuse_existing=reuse_existing,
            keep_work=keep_work,
        )

        counts = manifest["counts"]
        rows.append(
            {
                "threshold": threshold,
                "threshold_label": f"H{threshold}",
                "total": counts["total"],
                "train": counts["train"],
                "test": counts["test"],
                "actual_test_fraction": counts["actual_test_fraction"],
                "n_homology_clusters": counts["n_homology_clusters"],
                "shared_train_test_clusters": counts[
                    "shared_train_test_clusters"
                ],
                "n_warnings": len(manifest.get("warnings", [])),
            }
        )

    overview = pd.DataFrame(rows).sort_values(
        "threshold",
        ascending=False,
    )
    overview.to_csv(
        output_root / "split_overview.tsv",
        sep="\t",
        index=False,
    )

    print()
    print(f"Overview: {output_root / 'split_overview.tsv'}")
    print(
        "PASS — external homology-aware train/test splits materialized; "
        "no internal validation split created."
    )
    return overview

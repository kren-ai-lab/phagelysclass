#!/usr/bin/env python3
"""
02_build_uniprot_background.py

Build the broad UniProtKB phage-protein background for PhageLysClass using
native UniProt REST cursor pagination.

Design
------
This stage retrieves a broad phage-protein universe by taxonomy/context only.
It deliberately does NOT assign U/N states or use functional terms to define
the retrieval universe.

Large UniProt queries are retrieved from /uniprotkb/search with cursor
pagination (default page size: 500), retries are applied per page, and each
page is checkpointed locally. This avoids relying on one long /stream request.

This script does NOT:
    - define positives or negatives,
    - classify U/N states,
    - remove PhageLysData overlaps,
    - collapse exact-sequence entities,
    - perform homology filtering,
    - balance classes,
    - create task-specific datasets.

Expected registry
-----------------
config/queries/uniprot_background_queries_taxonomic.tsv

Required columns:
    query_id, scope, query, rationale, enabled

Outputs
-------
data/raw/uniprot_background/
├── pages/
│   ├── PHAGE_CAUDOVIRICETES/
│   │   ├── checkpoint.json
│   │   ├── page_000001.tsv.gz
│   │   └── ...
│   └── PHAGE_NAMED_NONCAUDOVIRICETES/
├── queries/
│   ├── PHAGE_CAUDOVIRICETES.tsv.gz
│   └── PHAGE_NAMED_NONCAUDOVIRICETES.tsv.gz
├── uniprot_background_raw.tsv.gz
├── uniprot_background_query_manifest.tsv
└── uniprot_background_summary.tsv

The final union is accession-level but exact-sequence duplicates are preserved.
The two default registry scopes are biologically disjoint by construction
(Caudoviricetes vs named phages outside Caudoviricetes).
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import re
import shutil
import sys
import time
from pathlib import Path
from typing import Optional

import pandas as pd
import requests


SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
DEFAULT_EXPECTED_RELEASE = "2026_02"
DEFAULT_PAGE_SIZE = 500

RETURN_FIELDS = [
    "accession",
    "id",
    "protein_name",
    "gene_names",
    "organism_name",
    "organism_id",
    "lineage",
    "lineage_ids",
    "virus_hosts",
    "reviewed",
    "protein_existence",
    "annotation_score",
    "fragment",
    "length",
    "sequence",
    "cc_function",
    "cc_catalytic_activity",
    "ec",
    "go_id",
    "keyword",
    "date_created",
    "date_modified",
]

REQUIRED_REGISTRY_COLUMNS = {
    "query_id",
    "scope",
    "query",
    "rationale",
    "enabled",
}


def parse_args() -> argparse.Namespace:
    script_dir = Path(__file__).resolve().parent
    project_root = script_dir.parents[1]

    parser = argparse.ArgumentParser(
        description=(
            "Retrieve a broad taxonomy-defined UniProtKB phage background "
            "using native cursor pagination."
        )
    )
    parser.add_argument(
        "--query-registry",
        default=str(
            project_root
            / "config"
            / "queries"
            / "uniprot_background_queries_taxonomic.tsv"
        ),
        help="Taxonomic UniProt query registry.",
    )
    parser.add_argument(
        "--output-dir",
        default="data/raw/uniprot_background",
        help="Output directory.",
    )
    parser.add_argument(
        "--mode",
        choices=["count", "download"],
        default="count",
        help="count = query-size audit; download = paginated retrieval.",
    )
    parser.add_argument(
        "--page-size",
        type=int,
        default=DEFAULT_PAGE_SIZE,
        help="UniProt search page size (default: 500; max recommended: 500).",
    )
    parser.add_argument(
        "--expected-release",
        default=DEFAULT_EXPECTED_RELEASE,
        help="Expected UniProtKB release.",
    )
    parser.add_argument(
        "--allow-release-mismatch",
        action="store_true",
        help="Continue if the live UniProtKB release differs from expected.",
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=120.0,
        help="Timeout per page request in seconds.",
    )
    parser.add_argument(
        "--retries",
        type=int,
        default=8,
        help="Maximum attempts per page.",
    )
    parser.add_argument(
        "--sleep",
        type=float,
        default=2.0,
        help="Base exponential-backoff sleep in seconds.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help=(
            "Restart each query from page 1, deleting its previous pages and "
            "checkpoint."
        ),
    )
    parser.add_argument(
        "--cleanup-pages",
        action="store_true",
        help="Delete page files after a query has been merged successfully.",
    )
    return parser.parse_args()


def truthy(value: object) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def load_registry(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"Query registry not found: {path}")

    registry = pd.read_csv(path, sep="\t", dtype=str).fillna("")

    missing = REQUIRED_REGISTRY_COLUMNS - set(registry.columns)
    if missing:
        raise KeyError(
            "Query registry missing columns: " + ", ".join(sorted(missing))
        )

    registry = registry.loc[registry["enabled"].map(truthy)].copy()

    if registry.empty:
        raise ValueError("No enabled queries in registry.")

    if registry["query_id"].duplicated().any():
        duplicates = sorted(
            registry.loc[
                registry["query_id"].duplicated(), "query_id"
            ].unique()
        )
        raise ValueError(
            "Duplicate query_id values: " + ", ".join(duplicates)
        )

    return registry.reset_index(drop=True)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def release_metadata(response: requests.Response) -> dict[str, str]:
    return {
        "uniprot_release": response.headers.get("X-UniProt-Release", ""),
        "uniprot_release_date": response.headers.get(
            "X-UniProt-Release-Date", ""
        ),
    }


def enforce_release(
    metadata: dict[str, str],
    expected: str,
    allow_mismatch: bool,
) -> None:
    observed = metadata.get("uniprot_release", "")

    if not observed:
        print(
            "WARNING: X-UniProt-Release header not exposed; "
            "release could not be verified.",
            file=sys.stderr,
        )
        return

    if expected and observed != expected:
        message = (
            f"UniProt release mismatch: expected {expected}, "
            f"observed {observed}."
        )
        if allow_mismatch:
            print("WARNING: " + message, file=sys.stderr)
        else:
            raise RuntimeError(
                message
                + " Use --allow-release-mismatch only if intentional."
            )


def get_next_link(link_header: Optional[str]) -> Optional[str]:
    """
    Extract UniProt's rel="next" cursor URL from the Link header.

    IMPORTANT:
    Do not split the Link header on commas. UniProt places the requested
    `fields` list inside the URL itself, and that list is comma-separated.
    Splitting on commas therefore corrupts the URL and makes pagination appear
    to end after page 1.
    """
    if not link_header:
        return None

    match = re.search(
        r'<([^>]+)>\s*;\s*rel=["\']?next["\']?',
        link_header,
        flags=re.IGNORECASE,
    )
    return match.group(1) if match else None


def request_page(
    session: requests.Session,
    *,
    url: str,
    params: Optional[dict[str, object]],
    args: argparse.Namespace,
) -> requests.Response:
    """
    Retrieve one UniProt page with retries.

    The entire page body is materialized by requests before returning. With a
    page size of 500, a transient connection failure affects only that page.
    """
    last_error: Exception | None = None

    for attempt in range(1, args.retries + 1):
        try:
            response = session.get(
                url,
                params=params,
                timeout=args.timeout,
                headers={
                    "User-Agent": (
                        "PhageLysClass/0.1 "
                        "(native paginated UniProt acquisition)"
                    ),
                    "Accept": "text/plain; format=tsv",
                },
            )

            if response.status_code == 429 or 500 <= response.status_code < 600:
                wait = args.sleep * (2 ** (attempt - 1))
                retry_after = response.headers.get("Retry-After")
                if retry_after:
                    try:
                        wait = max(wait, float(retry_after))
                    except ValueError:
                        pass

                print(
                    f"           transient HTTP {response.status_code}; "
                    f"retry {attempt}/{args.retries} in {wait:.1f}s",
                    file=sys.stderr,
                )
                response.close()
                time.sleep(wait)
                continue

            response.raise_for_status()

            # Force body download here so connection errors are caught by the
            # retry loop for this specific page.
            _ = response.content
            return response

        except (
            requests.RequestException,
            ConnectionError,
            TimeoutError,
        ) as exc:
            last_error = exc

            if attempt == args.retries:
                break

            wait = args.sleep * (2 ** (attempt - 1))
            print(
                f"           page attempt {attempt}/{args.retries} failed: "
                f"{exc}. Retrying in {wait:.1f}s",
                file=sys.stderr,
            )
            time.sleep(wait)

    raise RuntimeError(
        f"UniProt page failed after {args.retries} attempts: {last_error}"
    )


def validate_page_tsv(content: bytes) -> tuple[str, int]:
    """
    Validate one TSV page and return (header, data_rows).
    """
    text = content.decode("utf-8")

    if not text.strip():
        raise RuntimeError("UniProt returned an empty TSV page.")

    lines = text.splitlines()
    header = lines[0]

    if "\t" not in header:
        raise RuntimeError(
            "UniProt response does not look like a TSV page. "
            f"First line: {header[:200]!r}"
        )

    n_rows = max(0, len(lines) - 1)
    return header, n_rows


def write_page(path: Path, content: bytes) -> tuple[str, int]:
    header, n_rows = validate_page_tsv(content)

    path.parent.mkdir(parents=True, exist_ok=True)

    with gzip.open(path, "wb") as handle:
        handle.write(content)

    return header, n_rows


def load_checkpoint(path: Path) -> Optional[dict]:
    if not path.is_file():
        return None
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def save_checkpoint(path: Path, checkpoint: dict) -> None:
    tmp = path.with_suffix(".json.tmp")
    with tmp.open("w", encoding="utf-8") as handle:
        json.dump(checkpoint, handle, indent=2, sort_keys=True)
    tmp.replace(path)


def reset_query_state(page_dir: Path, query_output: Path) -> None:
    if page_dir.exists():
        shutil.rmtree(page_dir)
    if query_output.exists():
        query_output.unlink()


def count_query(
    session: requests.Session,
    query: str,
    args: argparse.Namespace,
) -> tuple[int | None, dict[str, str]]:
    response = request_page(
        session,
        url=SEARCH_URL,
        params={
            "query": query,
            "format": "tsv",
            "fields": "accession",
            "size": 1,
        },
        args=args,
    )

    try:
        metadata = release_metadata(response)
        enforce_release(
            metadata,
            args.expected_release,
            args.allow_release_mismatch,
        )

        total = response.headers.get("X-Total-Results")
        return (int(total) if total is not None else None), metadata
    finally:
        response.close()


def merge_pages(
    page_dir: Path,
    destination: Path,
) -> tuple[int, int]:
    """
    Merge page TSV.GZ files into one query-level TSV.GZ.

    Returns:
        n_pages, n_rows
    """
    pages = sorted(page_dir.glob("page_*.tsv.gz"))

    if not pages:
        raise FileNotFoundError(f"No page files found in {page_dir}")

    destination.parent.mkdir(parents=True, exist_ok=True)

    expected_header: Optional[str] = None
    total_rows = 0

    with gzip.open(destination, "wt", encoding="utf-8", newline="") as out:
        for page_idx, page in enumerate(pages, start=1):
            with gzip.open(page, "rt", encoding="utf-8", newline="") as inp:
                header = inp.readline().rstrip("\n")

                if expected_header is None:
                    expected_header = header
                    out.write(header + "\n")
                elif header != expected_header:
                    raise RuntimeError(
                        f"Header mismatch while merging {page}"
                    )

                for line in inp:
                    out.write(line)
                    total_rows += 1

    return len(pages), total_rows


def download_query_paginated(
    session: requests.Session,
    row: object,
    output_dir: Path,
    args: argparse.Namespace,
) -> dict[str, object]:
    page_dir = output_dir / "pages" / row.query_id
    query_output = output_dir / "queries" / f"{row.query_id}.tsv.gz"
    checkpoint_path = page_dir / "checkpoint.json"

    if args.force:
        reset_query_state(page_dir, query_output)

    page_dir.mkdir(parents=True, exist_ok=True)
    query_output.parent.mkdir(parents=True, exist_ok=True)

    checkpoint = load_checkpoint(checkpoint_path)

    if checkpoint:
        # Refuse to resume if the biological retrieval contract changed.
        for key, expected in {
            "query_id": row.query_id,
            "query": row.query,
            "page_size": args.page_size,
            "fields": RETURN_FIELDS,
        }.items():
            if checkpoint.get(key) != expected:
                raise RuntimeError(
                    f"Checkpoint mismatch for {row.query_id}: {key}. "
                    "Use --force to restart this query."
                )

        # A previous buggy version could mark page 1 as complete when Link
        # parsing failed. Never trust `complete=True` unless downloaded rows
        # equal UniProt's reported total.
        rows_ckpt = int(checkpoint.get("rows_downloaded", 0))
        total_ckpt = checkpoint.get("total_results_reported")
        checkpoint_is_consistent = (
            checkpoint.get("complete")
            and query_output.is_file()
            and total_ckpt is not None
            and rows_ckpt == int(total_ckpt)
        )

        if checkpoint_is_consistent:
            print(
                f"[RESUME] {row.query_id}: already complete; "
                f"reusing {query_output}"
            )
            return {
                "query_id": row.query_id,
                "scope": row.scope,
                "query": row.query,
                "rationale": row.rationale,
                "n_rows": rows_ckpt,
                "n_pages": int(checkpoint["pages_downloaded"]),
                "total_results_reported": total_ckpt,
                "uniprot_release": checkpoint.get("uniprot_release", ""),
                "uniprot_release_date": checkpoint.get(
                    "uniprot_release_date", ""
                ),
                "file": str(query_output),
                "sha256": sha256(query_output),
                "resumed": True,
            }

        if checkpoint.get("complete") and not checkpoint_is_consistent:
            print(
                f"[REPAIR] {row.query_id}: inconsistent old checkpoint "
                f"({rows_ckpt:,}/{int(total_ckpt) if total_ckpt is not None else 'NA'} rows). "
                "Restarting this query from page 1."
            )
            reset_query_state(page_dir, query_output)
            page_dir.mkdir(parents=True, exist_ok=True)
            checkpoint = None

    if checkpoint:
        next_url = checkpoint.get("next_url")
        page_index = int(checkpoint.get("next_page_index", 1))
        rows_downloaded = int(checkpoint.get("rows_downloaded", 0))
        pages_downloaded = int(checkpoint.get("pages_downloaded", 0))
        total_results = checkpoint.get("total_results_reported")
        metadata = {
            "uniprot_release": checkpoint.get("uniprot_release", ""),
            "uniprot_release_date": checkpoint.get(
                "uniprot_release_date", ""
            ),
        }

        print(
            f"[RESUME] {row.query_id}: page {page_index:,}, "
            f"{rows_downloaded:,} rows already saved"
        )

        # A non-complete checkpoint should always have a next URL.
        if not next_url:
            raise RuntimeError(
                f"Incomplete checkpoint for {row.query_id} has no next_url. "
                "Use --force to restart."
            )

        current_url = next_url
        current_params = None

    else:
        page_index = 1
        rows_downloaded = 0
        pages_downloaded = 0
        total_results = None
        metadata = {
            "uniprot_release": "",
            "uniprot_release_date": "",
        }

        current_url = SEARCH_URL
        current_params = {
            "query": row.query,
            "format": "tsv",
            "fields": ",".join(RETURN_FIELDS),
            "size": args.page_size,
        }

        print(f"[DOWNLOAD] {row.query_id}")
        print(f"           scope: {row.scope}")
        print(f"           query: {row.query}")
        print(f"           page size: {args.page_size}")

    while True:
        response = request_page(
            session,
            url=current_url,
            params=current_params,
            args=args,
        )

        try:
            page_metadata = release_metadata(response)
            enforce_release(
                page_metadata,
                args.expected_release,
                args.allow_release_mismatch,
            )

            if metadata["uniprot_release"] == "":
                metadata = page_metadata
            elif (
                page_metadata.get("uniprot_release")
                and page_metadata["uniprot_release"]
                != metadata["uniprot_release"]
            ):
                raise RuntimeError(
                    "UniProt release changed during paginated retrieval."
                )

            if total_results is None:
                total_header = response.headers.get("X-Total-Results")
                if total_header is not None:
                    total_results = int(total_header)

            next_url = get_next_link(response.headers.get("Link"))

            page_path = page_dir / f"page_{page_index:06d}.tsv.gz"

            # Page-level atomicity: write to temporary file, validate, then move.
            tmp_page = page_path.with_suffix(".tsv.gz.tmp")
            header, n_rows = write_page(tmp_page, response.content)
            tmp_page.replace(page_path)

            rows_downloaded += n_rows
            pages_downloaded += 1

            if total_results:
                pct = 100.0 * rows_downloaded / total_results
                print(
                    f"           page {page_index:>5,} | "
                    f"{rows_downloaded:>9,}/{total_results:,} "
                    f"({pct:6.2f}%)"
                )
            else:
                print(
                    f"           page {page_index:>5,} | "
                    f"{rows_downloaded:>9,} rows"
                )

            pagination_finished = next_url is None

            # `complete` is reserved for a fully merged + row-count-validated
            # query. Finishing cursor traversal alone is not enough.
            save_checkpoint(
                checkpoint_path,
                {
                    "query_id": row.query_id,
                    "scope": row.scope,
                    "query": row.query,
                    "fields": RETURN_FIELDS,
                    "page_size": args.page_size,
                    "expected_release": args.expected_release,
                    "uniprot_release": metadata["uniprot_release"],
                    "uniprot_release_date": metadata[
                        "uniprot_release_date"
                    ],
                    "total_results_reported": total_results,
                    "rows_downloaded": rows_downloaded,
                    "pages_downloaded": pages_downloaded,
                    "next_page_index": page_index + 1,
                    "next_url": next_url,
                    "pagination_finished": pagination_finished,
                    "complete": False,
                },
            )

            if pagination_finished:
                break

            current_url = next_url
            current_params = None
            page_index += 1

        finally:
            response.close()

    if total_results is not None and rows_downloaded < total_results:
        raise RuntimeError(
            f"Pagination stopped early for {row.query_id}: "
            f"{rows_downloaded:,}/{total_results:,} rows. "
            "The UniProt Link header did not provide a usable next cursor."
        )

    n_pages_merged, n_rows_merged = merge_pages(
        page_dir=page_dir,
        destination=query_output,
    )

    if n_rows_merged != rows_downloaded:
        raise RuntimeError(
            f"Row-count mismatch for {row.query_id}: "
            f"checkpoint={rows_downloaded:,}, merged={n_rows_merged:,}"
        )

    if total_results is not None and n_rows_merged != total_results:
        raise RuntimeError(
            f"UniProt total mismatch for {row.query_id}: "
            f"reported={total_results:,}, downloaded={n_rows_merged:,}"
        )

    # Mark checkpoint complete only after successful merged-file validation.
    checkpoint = load_checkpoint(checkpoint_path) or {}
    checkpoint["complete"] = True
    checkpoint["merged_file"] = str(query_output)
    checkpoint["merged_sha256"] = sha256(query_output)
    save_checkpoint(checkpoint_path, checkpoint)

    if args.cleanup_pages:
        for page in page_dir.glob("page_*.tsv.gz"):
            page.unlink()

    return {
        "query_id": row.query_id,
        "scope": row.scope,
        "query": row.query,
        "rationale": row.rationale,
        "n_rows": n_rows_merged,
        "n_pages": n_pages_merged,
        "total_results_reported": total_results,
        "uniprot_release": metadata["uniprot_release"],
        "uniprot_release_date": metadata["uniprot_release_date"],
        "file": str(query_output),
        "sha256": sha256(query_output),
        "resumed": False,
    }


def build_final_union(
    registry: pd.DataFrame,
    query_dir: Path,
    destination: Path,
) -> dict[str, int]:
    """
    Stream-merge query-level TSV.GZ files into the final raw background.

    Two provenance columns are appended:
        retrieval_query_id
        retrieval_scope

    Accessions and sequences are NOT collapsed here.
    """
    destination.parent.mkdir(parents=True, exist_ok=True)

    global_accessions: set[str] = set()
    sequence_hashes: set[bytes] = set()
    duplicate_accessions = 0
    total_rows = 0

    canonical_header: Optional[list[str]] = None
    accession_idx: Optional[int] = None
    sequence_idx: Optional[int] = None

    with gzip.open(destination, "wt", encoding="utf-8", newline="") as out_handle:
        writer = csv.writer(
            out_handle,
            delimiter="\t",
            lineterminator="\n",
            quoting=csv.QUOTE_MINIMAL,
        )

        for row in registry.itertuples(index=False):
            path = query_dir / f"{row.query_id}.tsv.gz"

            if not path.is_file():
                raise FileNotFoundError(f"Missing query-level file: {path}")

            with gzip.open(path, "rt", encoding="utf-8", newline="") as in_handle:
                reader = csv.reader(in_handle, delimiter="\t")
                header = next(reader)

                if canonical_header is None:
                    canonical_header = header
                    try:
                        accession_idx = header.index("Entry")
                        sequence_idx = header.index("Sequence")
                    except ValueError as exc:
                        raise RuntimeError(
                            "Expected UniProt TSV columns 'Entry' and "
                            "'Sequence' were not found. Header:\n"
                            + "\t".join(header)
                        ) from exc

                    writer.writerow(
                        canonical_header
                        + ["retrieval_query_id", "retrieval_scope"]
                    )
                elif header != canonical_header:
                    raise RuntimeError(
                        f"Header mismatch in {path}"
                    )

                assert accession_idx is not None
                assert sequence_idx is not None

                for record in reader:
                    if not record:
                        continue

                    accession = record[accession_idx]
                    sequence = record[sequence_idx]

                    if accession in global_accessions:
                        duplicate_accessions += 1
                    else:
                        global_accessions.add(accession)

                    # Keep only a fixed-size digest for exact-sequence counting.
                    if sequence:
                        sequence_hashes.add(
                            hashlib.sha1(
                                sequence.encode("utf-8")
                            ).digest()
                        )

                    writer.writerow(
                        record + [row.query_id, row.scope]
                    )
                    total_rows += 1

    return {
        "rows_union": total_rows,
        "unique_accessions": len(global_accessions),
        "duplicate_accessions_across_scopes": duplicate_accessions,
        "unique_exact_sequences_before_consolidation": len(sequence_hashes),
    }


def run_count(
    registry: pd.DataFrame,
    output_dir: Path,
    args: argparse.Namespace,
) -> None:
    records = []

    with requests.Session() as session:
        for row in registry.itertuples(index=False):
            print(f"[COUNT] {row.query_id}")
            print(f"        query: {row.query}")

            total, metadata = count_query(
                session=session,
                query=row.query,
                args=args,
            )

            print(
                f"        results: "
                f"{total:,}" if total is not None else "        results: NA"
            )

            records.append(
                {
                    "query_id": row.query_id,
                    "scope": row.scope,
                    "n_results": total,
                    "uniprot_release": metadata["uniprot_release"],
                    "uniprot_release_date": metadata[
                        "uniprot_release_date"
                    ],
                    "query": row.query,
                    "rationale": row.rationale,
                }
            )

    counts = pd.DataFrame(records)
    path = output_dir / "uniprot_background_query_counts.tsv"
    counts.to_csv(path, sep="\t", index=False)

    print()
    print(f"Query audit written: {path}")


def run_download(
    registry: pd.DataFrame,
    output_dir: Path,
    args: argparse.Namespace,
) -> None:
    manifest_rows = []

    with requests.Session() as session:
        # Be explicit about persistent connection reuse across small pages.
        session.headers.update(
            {
                "Accept-Encoding": "gzip, deflate",
            }
        )

        for row in registry.itertuples(index=False):
            result = download_query_paginated(
                session=session,
                row=row,
                output_dir=output_dir,
                args=args,
            )
            manifest_rows.append(result)
            print(
                f"           COMPLETE: {result['n_rows']:,} rows "
                f"in {result['n_pages']:,} pages\n"
            )

    manifest = pd.DataFrame(manifest_rows)
    manifest_path = output_dir / "uniprot_background_query_manifest.tsv"
    manifest.to_csv(manifest_path, sep="\t", index=False)

    raw_path = output_dir / "uniprot_background_raw.tsv.gz"
    union_stats = build_final_union(
        registry=registry,
        query_dir=output_dir / "queries",
        destination=raw_path,
    )

    summary_rows = [
        {
            "metric": "query_rows_sum",
            "value": int(manifest["n_rows"].sum()),
        },
        *[
            {"metric": key, "value": int(value)}
            for key, value in union_stats.items()
        ],
    ]

    summary = pd.DataFrame(summary_rows)
    summary_path = output_dir / "uniprot_background_summary.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)

    if union_stats["duplicate_accessions_across_scopes"] != 0:
        print(
            "WARNING: duplicate accessions were observed across retrieval "
            "scopes. They are preserved in the raw file for audit.",
            file=sys.stderr,
        )

    print("=" * 88)
    print("PhageLysClass — paginated UniProt background complete")
    print("=" * 88)
    print(f"Rows                  : {union_stats['rows_union']:,}")
    print(f"Unique accessions     : {union_stats['unique_accessions']:,}")
    print(
        "Unique exact sequences: "
        f"{union_stats['unique_exact_sequences_before_consolidation']:,}"
    )
    print(f"Canonical raw dataset : {raw_path}")
    print(f"Manifest              : {manifest_path}")
    print(f"Summary               : {summary_path}")
    print()
    print(
        "PASS — taxonomy-defined raw acquisition only. No U/N assignment, "
        "PhageLysData decontamination, sequence consolidation, homology "
        "filtering, balancing, or task construction applied."
    )


def main() -> None:
    args = parse_args()

    if not 1 <= args.page_size <= 500:
        raise ValueError("--page-size must be between 1 and 500.")

    registry_path = Path(args.query_registry).expanduser().resolve()
    output_dir = Path(args.output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    registry = load_registry(registry_path)

    print("=" * 88)
    print("PhageLysClass — broad UniProtKB phage background")
    print("=" * 88)
    print(f"Mode             : {args.mode}")
    print(
        f"Acquisition      : "
        f"{'native cursor pagination' if args.mode == 'download' else 'count audit'}"
    )
    print(f"Query registry   : {registry_path}")
    print(f"Enabled queries  : {len(registry):,}")
    print(f"Expected release : {args.expected_release}")
    print(f"Page size        : {args.page_size}")
    print(f"Output dir       : {output_dir}")
    print()

    for row in registry.itertuples(index=False):
        print(f"  - {row.query_id}: {row.query}")
    print()

    if args.mode == "count":
        run_count(
            registry=registry,
            output_dir=output_dir,
            args=args,
        )
    else:
        run_download(
            registry=registry,
            output_dir=output_dir,
            args=args,
        )


if __name__ == "__main__":
    main()

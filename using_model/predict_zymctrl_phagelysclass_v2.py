#!/usr/bin/env python3
from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Run frozen PhageLysClass inference on ZymCTRL-generated sequences."
    )
    p.add_argument("--input-fasta", required=True, type=Path)
    p.add_argument("--input-metadata", required=True, type=Path)
    p.add_argument("--models", required=True, type=Path)
    p.add_argument(
        "--predictor",
        required=True,
        type=Path,
        help="Path to phagelysclass_predict_large_fasta.py",
    )
    p.add_argument(
        "--output-dir",
        type=Path,
        default=Path("zymctrl_phagelys_predictions"),
    )
    p.add_argument("--device", choices=["auto", "cpu", "cuda"], default="auto")
    p.add_argument("--batch-size", type=int, default=2)
    return p.parse_args()


def load_predictor_module(path: Path):
    spec = importlib.util.spec_from_file_location("phagelysclass_predictor_module", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not import predictor script: {path}")

    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    if not hasattr(module, "PhageLysClassPredictor"):
        raise RuntimeError(f"{path} does not expose PhageLysClassPredictor.")

    return module


def main() -> None:
    args = parse_args()

    fasta = args.input_fasta.resolve()
    metadata_path = args.input_metadata.resolve()
    models = args.models.resolve()
    predictor_path = args.predictor.resolve()
    outdir = args.output_dir.resolve()

    for path in (fasta, metadata_path, predictor_path):
        if not path.exists():
            raise FileNotFoundError(path)

    if not models.is_dir():
        raise FileNotFoundError(models)

    outdir.mkdir(parents=True, exist_ok=True)

    raw_predictions_path = outdir / "zymctrl_phagelysclass_predictions.tsv"
    merged_path = outdir / "zymctrl_phagelysclass_predictions_with_metadata.tsv"
    summary_ec_path = outdir / "zymctrl_phagelysclass_summary_by_ec.tsv"
    summary_overall_path = outdir / "zymctrl_phagelysclass_summary_overall.tsv"
    run_metadata_path = outdir / "zymctrl_phagelysclass_run_metadata.json"

    module = load_predictor_module(predictor_path)

    predictor = module.PhageLysClassPredictor(
        model_root=models,
        device=args.device,
        batch_size=args.batch_size,
    )

    predictions, run_metadata = predictor.predict_fasta(
        fasta_path=fasta,
        auxiliary=True,
    )

    predictions.to_csv(raw_predictions_path, sep="\t", index=False)

    meta = pd.read_csv(metadata_path, sep="\t")

    if "sequence_id" not in meta.columns:
        raise ValueError("ZymCTRL metadata must contain a sequence_id column.")

    if meta["sequence_id"].duplicated().any():
        raise ValueError("Duplicate sequence_id values found in ZymCTRL metadata.")

    if predictions["sequence_id"].duplicated().any():
        raise ValueError("Duplicate sequence_id values found in prediction output.")

    merged = meta.merge(
        predictions,
        on="sequence_id",
        how="left",
        validate="one_to_one",
        suffixes=("_zymctrl", "_prediction"),
    )

    missing_predictions = int(merged["global_prediction"].isna().sum())
    if missing_predictions:
        raise RuntimeError(
            f"{missing_predictions} ZymCTRL sequences did not receive a prediction."
        )

    merged.to_csv(merged_path, sep="\t", index=False)

    n_valid = int((merged["input_status"] == "valid").sum())
    n_lytic = int((merged["global_prediction"] == "lytic").sum())

    overall = {
        "n_sequences": int(len(merged)),
        "n_valid": n_valid,
        "n_invalid": int((merged["input_status"] == "invalid").sum()),
        "n_lytic": n_lytic,
        "n_non_lytic_background": int(
            (merged["global_prediction"] == "non_lytic_background").sum()
        ),
        "fraction_lytic": n_lytic / n_valid if n_valid else float("nan"),
        "n_endolysin": int(
            (merged["functional_prediction"] == "endolysin").sum()
        ),
        "n_depolymerase": int(
            (merged["functional_prediction"] == "depolymerase").sum()
        ),
        "n_val": int(
            (merged["functional_prediction"] == "VAL").sum()
        ),
    }

    overall["fraction_endolysin_among_lytic"] = (
        overall["n_endolysin"] / n_lytic if n_lytic else 0.0
    )
    overall["fraction_depolymerase_among_lytic"] = (
        overall["n_depolymerase"] / n_lytic if n_lytic else 0.0
    )
    overall["fraction_val_among_lytic"] = (
        overall["n_val"] / n_lytic if n_lytic else 0.0
    )

    pd.DataFrame([overall]).to_csv(
        summary_overall_path,
        sep="\t",
        index=False,
    )

    if "ec_number" in merged.columns:
        rows = []

        for ec, group in merged.groupby("ec_number", dropna=False):
            n_valid_ec = int((group["input_status"] == "valid").sum())
            n_lytic_ec = int((group["global_prediction"] == "lytic").sum())

            row = {
                "ec_number": ec,
                "n_sequences": int(len(group)),
                "n_valid": n_valid_ec,
                "n_lytic": n_lytic_ec,
                "n_non_lytic_background": int(
                    (group["global_prediction"] == "non_lytic_background").sum()
                ),
                "fraction_lytic": (
                    n_lytic_ec / n_valid_ec if n_valid_ec else float("nan")
                ),
                "mean_global_score": group["global_score"].mean(),
                "median_global_score": group["global_score"].median(),
                "n_endolysin": int(
                    (group["functional_prediction"] == "endolysin").sum()
                ),
                "n_depolymerase": int(
                    (group["functional_prediction"] == "depolymerase").sum()
                ),
                "n_val": int(
                    (group["functional_prediction"] == "VAL").sum()
                ),
            }

            if "panel_block" in group.columns:
                vals = group["panel_block"].dropna().astype(str).unique()
                row["panel_block"] = ";".join(sorted(vals))

            if "activity" in group.columns:
                vals = group["activity"].dropna().astype(str).unique()
                row["activity"] = ";".join(sorted(vals))

            if "perplexity" in group.columns:
                row["mean_perplexity"] = group["perplexity"].mean()
                row["median_perplexity"] = group["perplexity"].median()

            rows.append(row)

        pd.DataFrame(rows).sort_values("ec_number").to_csv(
            summary_ec_path,
            sep="\t",
            index=False,
        )

    payload = {
        "analysis": "ZymCTRL -> PhageLysClass",
        "input_fasta": str(fasta),
        "input_metadata": str(metadata_path),
        "model_bundle": str(models),
        "predictor_script": str(predictor_path),
        "device": args.device,
        "batch_size": args.batch_size,
        "auxiliary": True,
        "predictor_run_metadata": run_metadata,
        "n_sequences_merged": int(len(merged)),
        "outputs": {
            "raw_predictions": str(raw_predictions_path),
            "merged_predictions": str(merged_path),
            "summary_overall": str(summary_overall_path),
            "summary_by_ec": str(summary_ec_path),
        },
    }

    run_metadata_path.write_text(
        json.dumps(payload, indent=2),
        encoding="utf-8",
    )

    print("\nDONE")
    print(f"Sequences:      {len(merged):,}")
    print(f"Lytic:          {overall['n_lytic']:,}")
    print(f"Fraction lytic: {overall['fraction_lytic']:.4f}")
    print(f"Endolysin:      {overall['n_endolysin']:,}")
    print(f"Depolymerase:   {overall['n_depolymerase']:,}")
    print(f"VAL:            {overall['n_val']:,}")
    print()
    print(f"Merged output:  {merged_path}")
    print(f"Summary by EC:  {summary_ec_path}")
    print(f"Metadata:       {run_metadata_path}")


if __name__ == "__main__":
    main()

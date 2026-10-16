#!/usr/bin/env python3
"""
Generate the fixed PhageLys ZymCTRL EC panel.

Design
------
10 EC conditions x 1,000 accepted unique sequences = 10,000 sequences.

The panel deliberately uses catalytic activities, not biological localization:
EC numbers cannot distinguish an endolysin from a virion-associated lysin (VAL)
when both catalyse the same reaction. Therefore the five peptidoglycan-lytic ECs
are interpreted jointly as Endolysin/VAL-compatible catalytic space, while the
five surface-polysaccharide ECs represent Depolymerase-compatible catalytic space.

No perplexity threshold is used for acceptance. Perplexity is recorded as
metadata so the generated cohort is not cherry-picked before PhageLysClass
inference.

Outputs
-------
<output-dir>/zymctrl_phagelys_10000.fasta
<output-dir>/zymctrl_phagelys_10000.tsv
<output-dir>/zymctrl_phagelys_10000.summary.tsv
<output-dir>/zymctrl_phagelys_10000.config.json

Example
-------
python generate_zymctrl_phagelys_panel.py --device cuda --output-dir zymctrl_panel

Dependencies
------------
pip install torch transformers pandas tqdm
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import random
from collections import Counter
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import pandas as pd
import torch
from tqdm import tqdm
from transformers import AutoTokenizer, GPT2LMHeadModel

MODEL_ID = "AI4PD/ZymCTRL"
VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")
SPECIAL_TOKENS = ["<start>", "<end>", "<|endoftext|>", "<pad>", " ", "<sep>", "\n", "\r", "\t"]

PANEL = [
    {
        "ec_number": "3.2.1.17",
        "panel_block": "peptidoglycan_lytic",
        "class_relation": "Endolysin_or_VAL",
        "activity": "lysozyme / N-acetylmuramidase",
    },
    {
        "ec_number": "3.2.1.52",
        "panel_block": "peptidoglycan_lytic",
        "class_relation": "Endolysin_or_VAL",
        "activity": "beta-N-acetylhexosaminidase / glucosaminidase",
    },
    {
        "ec_number": "3.5.1.28",
        "panel_block": "peptidoglycan_lytic",
        "class_relation": "Endolysin_or_VAL",
        "activity": "N-acetylmuramoyl-L-alanine amidase",
    },
    {
        "ec_number": "4.2.2.29",
        "panel_block": "peptidoglycan_lytic",
        "class_relation": "Endolysin_or_VAL",
        "activity": "peptidoglycan lytic transglycosylase",
    },
    {
        "ec_number": "3.4.24.75",
        "panel_block": "peptidoglycan_lytic",
        "class_relation": "Endolysin_or_VAL",
        "activity": "glycyl-glycine endopeptidase / M23 peptidase",
    },
    {
        "ec_number": "3.2.1.87",
        "panel_block": "surface_depolymerase",
        "class_relation": "Depolymerase",
        "activity": "capsular-polysaccharide endo-1,3-alpha-galactosidase",
    },
    {
        "ec_number": "3.2.1.129",
        "panel_block": "surface_depolymerase",
        "class_relation": "Depolymerase",
        "activity": "endo-alpha-sialidase / endosialidase",
    },
    {
        "ec_number": "4.2.2.1",
        "panel_block": "surface_depolymerase",
        "class_relation": "Depolymerase",
        "activity": "hyaluronate lyase",
    },
    {
        "ec_number": "4.2.2.3",
        "panel_block": "surface_depolymerase",
        "class_relation": "Depolymerase",
        "activity": "mannuronate-specific alginate lyase",
    },
    {
        "ec_number": "3.2.1.65",
        "panel_block": "surface_depolymerase",
        "class_relation": "Depolymerase",
        "activity": "levanase",
    },
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Generate the fixed 10-EC PhageLys ZymCTRL panel.")
    p.add_argument("--output-dir", type=Path, default=Path("zymctrl_phagelys_panel"))
    p.add_argument("--model", default=MODEL_ID, help="Hugging Face model ID or local ZymCTRL path.")
    p.add_argument("--per-ec", type=int, default=1000, help="Accepted unique sequences per EC.")
    p.add_argument("--batch-size", type=int, default=10)
    p.add_argument("--device", choices=["auto", "cuda", "cpu"], default="auto")
    p.add_argument("--max-length", type=int, default=1024)
    p.add_argument("--top-k", type=int, default=9)
    p.add_argument("--repetition-penalty", type=float, default=1.2)
    p.add_argument("--temperature", type=float, default=1.0)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument(
        "--max-batches-factor",
        type=float,
        default=5.0,
        help="Safety factor above the theoretical number of batches required per EC.",
    )
    return p.parse_args()


def select_device(name: str) -> torch.device:
    if name == "auto":
        return torch.device("cuda" if torch.cuda.is_available() else "cpu")
    if name == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA requested but no CUDA device is available.")
    return torch.device(name)


def decode_sequence(decoded: str) -> Tuple[str | None, str | None]:
    if "<sep>" not in decoded:
        return None, "missing_sep"

    seq = decoded.split("<sep>", 1)[1]
    for token in SPECIAL_TOKENS:
        seq = seq.replace(token, "")
    seq = seq.strip().upper()

    if not seq:
        return None, "empty_sequence"

    invalid = sorted(set(seq) - VALID_AA)
    if invalid:
        return None, "noncanonical:" + "".join(invalid)

    return seq, None


def calculate_perplexity(token_ids: torch.Tensor, model: GPT2LMHeadModel) -> float:
    """Match the ZymCTRL model-card definition: exp(model loss) on the full generation."""
    with torch.inference_mode():
        outputs = model(token_ids.unsqueeze(0), labels=token_ids.unsqueeze(0))
    return float(math.exp(float(outputs.loss.detach().cpu())))


def wrap_fasta(seq: str, width: int = 80) -> str:
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def main() -> None:
    args = parse_args()
    if args.per_ec <= 0:
        raise ValueError("--per-ec must be > 0")
    if args.batch_size <= 0:
        raise ValueError("--batch-size must be > 0")

    device = select_device(args.device)
    outdir = args.output_dir
    outdir.mkdir(parents=True, exist_ok=True)

    total_target = len(PANEL) * args.per_ec
    prefix = outdir / f"zymctrl_phagelys_{total_target}"

    fasta_path = Path(str(prefix) + ".fasta")
    tsv_path = Path(str(prefix) + ".tsv")
    summary_path = Path(str(prefix) + ".summary.tsv")
    config_path = Path(str(prefix) + ".config.json")

    for path in (fasta_path, tsv_path, summary_path, config_path):
        if path.exists():
            raise FileExistsError(f"Refusing to overwrite existing output: {path}")

    random.seed(args.seed)
    torch.manual_seed(args.seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(args.seed)

    print(f"Loading tokenizer/model: {args.model}")
    tokenizer = AutoTokenizer.from_pretrained(args.model)
    model = GPT2LMHeadModel.from_pretrained(args.model).to(device)
    model.eval()

    # ZymCTRL's published generation recipe uses eos_token_id=1 and pad_token_id=0.
    eos_token_id = 1
    pad_token_id = 0

    config = {
        "model": args.model,
        "panel": PANEL,
        "per_ec": args.per_ec,
        "total_requested": total_target,
        "batch_size": args.batch_size,
        "max_length": args.max_length,
        "top_k": args.top_k,
        "repetition_penalty": args.repetition_penalty,
        "temperature": args.temperature,
        "seed": args.seed,
        "device": str(device),
        "eos_token_id": eos_token_id,
        "pad_token_id": pad_token_id,
        "acceptance_policy": (
            "completed before max_length + canonical amino acids + globally unique sequence; "
            "perplexity recorded but not used for filtering"
        ),
    }
    config_path.write_text(json.dumps(config, indent=2), encoding="utf-8")

    fields = [
        "sequence_id",
        "ec_number",
        "panel_block",
        "class_relation",
        "activity",
        "sequence",
        "length",
        "perplexity",
        "batch_id",
        "seed",
        "sha256",
    ]

    seen_hashes = set()
    global_index = 0
    global_batch = 0
    summaries: List[Dict[str, object]] = []

    with fasta_path.open("w", encoding="utf-8") as fasta_handle, \
         tsv_path.open("w", encoding="utf-8", newline="") as tsv_handle:

        writer = csv.DictWriter(tsv_handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()

        progress = tqdm(total=total_target, desc="Accepted", unit="seq")

        for entry in PANEL:
            ec = entry["ec_number"]
            target = args.per_ec
            accepted = 0
            attempted = 0
            rejected = Counter()

            prompt = tokenizer.encode(ec, return_tensors="pt").to(device)

            theoretical_batches = math.ceil(target / args.batch_size)
            max_batches = max(
                theoretical_batches + 10,
                math.ceil(theoretical_batches * args.max_batches_factor),
            )
            ec_batches = 0

            print(f"\n{ec} | {entry['activity']} | target={target}")

            while accepted < target:
                if ec_batches >= max_batches:
                    raise RuntimeError(
                        f"Safety limit reached for EC {ec}: {accepted}/{target} accepted. "
                        f"Rejected: {dict(rejected)}. Increase --max-batches-factor if needed."
                    )

                remaining = target - accepted
                # Generate full batches when possible. For the last batch, only ask for what is needed.
                n_return = min(args.batch_size, remaining)

                batch_seed = args.seed + global_batch
                torch.manual_seed(batch_seed)
                if torch.cuda.is_available():
                    torch.cuda.manual_seed_all(batch_seed)

                with torch.inference_mode():
                    outputs = model.generate(
                        prompt,
                        top_k=args.top_k,
                        repetition_penalty=args.repetition_penalty,
                        temperature=args.temperature,
                        max_length=args.max_length,
                        eos_token_id=eos_token_id,
                        pad_token_id=pad_token_id,
                        do_sample=True,
                        num_return_sequences=n_return,
                    )

                attempted += int(outputs.shape[0])
                current_batch = global_batch
                global_batch += 1
                ec_batches += 1

                for output in outputs:
                    if accepted >= target:
                        break

                    # Official ZymCTRL recipe rejects max-length-truncated outputs.
                    # In a padded generation batch, completed samples end in PAD=0.
                    # For single/longest rows, EOS=1 can be the last token instead.
                    token_list = output.tolist()
                    if eos_token_id not in token_list:
                        rejected["truncated_no_eos"] += 1
                        continue

                    eos_pos = token_list.index(eos_token_id)
                    trimmed = output[: eos_pos + 1]

                    decoded = tokenizer.decode(trimmed, skip_special_tokens=False)
                    seq, reason = decode_sequence(decoded)
                    if seq is None:
                        rejected[reason] += 1
                        continue

                    digest = hashlib.sha256(seq.encode("ascii")).hexdigest()
                    if digest in seen_hashes:
                        rejected["duplicate_sequence"] += 1
                        continue

                    ppl = calculate_perplexity(trimmed.to(device), model)

                    global_index += 1
                    accepted += 1
                    seen_hashes.add(digest)
                    seq_id = f"ZYMCTRL_{global_index:08d}"

                    row = {
                        "sequence_id": seq_id,
                        "ec_number": ec,
                        "panel_block": entry["panel_block"],
                        "class_relation": entry["class_relation"],
                        "activity": entry["activity"],
                        "sequence": seq,
                        "length": len(seq),
                        "perplexity": f"{ppl:.6f}",
                        "batch_id": current_batch,
                        "seed": batch_seed,
                        "sha256": digest,
                    }
                    writer.writerow(row)

                    fasta_handle.write(
                        f">{seq_id} ec={ec} block={entry['panel_block']} "
                        f"class_relation={entry['class_relation']} ppl={ppl:.6f}\n"
                        f"{wrap_fasta(seq)}\n"
                    )

                    progress.update(1)

                tsv_handle.flush()
                fasta_handle.flush()

            summaries.append(
                {
                    "ec_number": ec,
                    "panel_block": entry["panel_block"],
                    "class_relation": entry["class_relation"],
                    "activity": entry["activity"],
                    "target": target,
                    "accepted": accepted,
                    "attempted": attempted,
                    "batches": ec_batches,
                    "rejected_total": sum(rejected.values()),
                    "rejection_reasons": json.dumps(dict(rejected), sort_keys=True),
                }
            )

        progress.close()

    pd.DataFrame(summaries).to_csv(summary_path, sep="\t", index=False)

    print("\nGeneration complete")
    print(f"Accepted unique sequences: {global_index}")
    print(f"FASTA:   {fasta_path}")
    print(f"TSV:     {tsv_path}")
    print(f"Summary: {summary_path}")
    print(f"Config:  {config_path}")


if __name__ == "__main__":
    main()

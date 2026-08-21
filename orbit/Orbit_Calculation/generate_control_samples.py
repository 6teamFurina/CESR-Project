#!/usr/bin/env python3
"""Generate deterministic CESR corrector samples shared by Bmad and SciBmad."""

from __future__ import annotations

import argparse
import csv
import json
import math
import random
from pathlib import Path


HERE = Path(__file__).resolve().parent
ORBIT_ROOT = HERE.parent
PROJECT_ROOT = HERE.parents[1]
LATEST_CONTROL_METADATA = (
    PROJECT_ROOT
    / "Latest_Lattice"
    / "bmad_reference"
    / "control_tracking"
    / "controls.csv"
)
LEGACY_REFERENCE = ORBIT_ROOT / "reference" / "closed_orbit_response_6x119.csv"
LATEST_OUTPUT = HERE / "inputs" / "latest_cesr" / "corrector_samples.csv"
LEGACY_OUTPUT = HERE / "inputs" / "cesr_corrector_samples_1000.csv"


def control_names(reference: Path) -> list[str]:
    with reference.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.reader(stream))
    if not rows:
        raise RuntimeError(f"Empty control reference in {reference}")
    header = rows[0]
    if header and header[0] == "lord_id":
        columns = {name: index for index, name in enumerate(header)}
        required = {"lord_name", "lord_key", "variable"}
        missing = required.difference(columns)
        if missing:
            raise RuntimeError(
                f"Control metadata is missing columns {sorted(missing)}: {reference}"
            )
        names = [
            row[columns["lord_name"]]
            for row in rows[1:]
            if len(row) == len(header)
            and row[columns["lord_key"]].upper() == "OVERLAY"
            and row[columns["variable"]].upper() in {"HKICK", "VKICK"}
        ]
    elif header and header[0] in {"coordinate", "observable"}:
        # Historical response matrices remain supported when a legacy ring is
        # selected explicitly.  Their dimensions come from the header.
        names = header[1:]
    else:
        raise RuntimeError(f"Unexpected control reference header in {reference}")
    if not names or len(set(names)) != len(names):
        raise RuntimeError(f"Control names are empty or duplicated in {reference}")
    return names


def clipped_gaussian(
    generator: random.Random,
    sigma: float,
    clip_sigma: float,
) -> float:
    limit = clip_sigma * sigma
    return max(-limit, min(limit, generator.gauss(0.0, sigma)))


def make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260729)
    parser.add_argument("--sigma-rad", type=float, default=5.0e-6)
    parser.add_argument("--clip-sigma", type=float, default=3.0)
    parser.add_argument(
        "--ring",
        choices=("latest", "legacy"),
        default="latest",
        help="Select latest repaired CESR controls, or explicitly reproduce legacy inputs",
    )
    parser.add_argument("--reference", type=Path, default=None)
    parser.add_argument("--output", type=Path, default=None)
    return parser


def main() -> int:
    args = make_parser().parse_args()
    if args.samples < 2:
        raise ValueError("--samples must be at least 2")
    if not math.isfinite(args.sigma_rad) or args.sigma_rad <= 0:
        raise ValueError("--sigma-rad must be finite and positive")
    if not math.isfinite(args.clip_sigma) or args.clip_sigma <= 0:
        raise ValueError("--clip-sigma must be finite and positive")

    default_reference = (
        LATEST_CONTROL_METADATA if args.ring == "latest" else LEGACY_REFERENCE
    )
    default_output = LATEST_OUTPUT if args.ring == "latest" else LEGACY_OUTPUT
    reference = (args.reference or default_reference).expanduser().resolve()
    output = (args.output or default_output).expanduser().resolve()
    names = control_names(reference)
    generator = random.Random(args.seed)

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["sample_id", *names])
        writer.writerow([0, *([0.0] * len(names))])
        for sample_id in range(1, args.samples):
            values = [
                clipped_gaussian(generator, args.sigma_rad, args.clip_sigma)
                for _ in names
            ]
            writer.writerow([sample_id, *(f"{value:.17g}" for value in values)])

    metadata = {
        "format": "cesr-corrector-samples-v2",
        "ring": "latest_cesr" if args.ring == "latest" else "legacy",
        "samples": args.samples,
        "controls": len(names),
        "seed": args.seed,
        "distribution": "independent clipped Gaussian",
        "sigma_rad": args.sigma_rad,
        "clip_sigma": args.clip_sigma,
        "sample_0_is_zero_baseline": True,
        "control_order": names,
        "reference_response_csv": str(reference),
        "output_csv": str(output),
    }
    output.with_suffix(".metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n",
        encoding="utf-8",
    )
    print(f"Wrote {args.samples} samples x {len(names)} controls to {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

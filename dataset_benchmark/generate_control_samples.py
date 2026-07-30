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
DEFAULT_REFERENCE = (
    HERE.parent
    / "bmad_control_response_rf_on"
    / "bmad_control_response_rf_on.csv"
)


def control_names(reference: Path) -> list[str]:
    with reference.open(encoding="utf-8", newline="") as stream:
        header = next(csv.reader(stream))
    if not header or header[0] != "observable":
        raise RuntimeError(f"Unexpected response-matrix header in {reference}")
    names = header[1:]
    if len(names) != 119 or len(set(names)) != 119:
        raise RuntimeError(
            f"Expected 119 unique CESR correctors, found {len(names)}"
        )
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
    parser.add_argument("--reference", type=Path, default=DEFAULT_REFERENCE)
    parser.add_argument(
        "--output",
        type=Path,
        default=HERE / "inputs" / "cesr_corrector_samples_1000.csv",
    )
    return parser


def main() -> int:
    args = make_parser().parse_args()
    if args.samples < 2:
        raise ValueError("--samples must be at least 2")
    if not math.isfinite(args.sigma_rad) or args.sigma_rad <= 0:
        raise ValueError("--sigma-rad must be finite and positive")
    if not math.isfinite(args.clip_sigma) or args.clip_sigma <= 0:
        raise ValueError("--clip-sigma must be finite and positive")

    reference = args.reference.expanduser().resolve()
    output = args.output.expanduser().resolve()
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
        "format": "cesr-corrector-samples-v1",
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

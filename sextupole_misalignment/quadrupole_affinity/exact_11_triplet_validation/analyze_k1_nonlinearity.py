#!/usr/bin/env python3
"""Quantify the even (nonlinear) part of each symmetric K1 response scan."""

from __future__ import annotations

import argparse
import csv
import json
import tomllib
from pathlib import Path
from typing import Any

import numpy as np

HERE = Path(__file__).resolve().parent


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--scan-root", type=Path, required=True)
    result.add_argument("--output-dir", type=Path, required=True)
    result.add_argument("--bpm-noise-m", type=float, default=5.0e-6)
    result.add_argument("--material-even-to-odd", type=float, default=0.05)
    return result


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def distribution(values: list[float]) -> dict[str, float]:
    array = np.asarray(values, dtype=float)
    return {
        "minimum": float(np.min(array)),
        "median": float(np.median(array)),
        "p90": float(np.percentile(array, 90.0)),
        "maximum": float(np.max(array)),
    }


def rms(values: np.ndarray) -> float:
    return float(np.sqrt(np.mean(values**2)))


def main() -> int:
    args = parser().parse_args()
    scan_root = args.scan_root.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    output_rows: list[dict[str, Any]] = []

    for scan_dir in sorted(path for path in scan_root.iterdir() if path.is_dir()):
        metadata_path = scan_dir / "scan_metadata.toml"
        if not metadata_path.exists():
            continue
        metadata = tomllib.loads(metadata_path.read_text(encoding="utf-8"))
        labels = [
            line.strip()
            for line in (scan_dir / "scenario_labels.txt").read_text(
                encoding="utf-8"
            ).splitlines()
            if line.strip()
        ]
        truth_index = labels.index("truth")
        observations = np.load(scan_dir / "bpm_orbits.npy", mmap_mode="r")
        levels = (
            np.asarray(metadata["k2_levels"], dtype=float)
            * float(metadata["k2_step_m3"])
        )
        slope_noise = args.bpm_noise_m / np.sqrt(float(levels @ levels))
        slopes = np.einsum(
            "k,cbkmp->cbmp", levels, observations[truth_index]
        ) / float(levels @ levels)
        conditions = read_rows(scan_dir / "k1_conditions.csv")
        nominal = slopes[0]

        for candidate in metadata["candidate_quadrupoles"]:
            plus_index = next(
                index
                for index, row in enumerate(conditions)
                if row["quadrupole"] == candidate and int(row["sign"]) == 1
            )
            minus_index = next(
                index
                for index, row in enumerate(conditions)
                if row["quadrupole"] == candidate and int(row["sign"]) == -1
            )
            plus = slopes[plus_index]
            minus = slopes[minus_index]
            odd = 0.5 * (plus - minus)
            even = 0.5 * (plus + minus) - nominal
            odd_rms = rms(odd)
            even_rms = rms(even)
            even_to_odd = even_rms / odd_rms if odd_rms > 0.0 else float("inf")
            odd_snr = odd_rms / (slope_noise / np.sqrt(2.0))
            even_snr = even_rms / (slope_noise * np.sqrt(1.5))
            output_rows.append(
                {
                    "sextupole": metadata["target_sextupole"],
                    "quadrupole": candidate,
                    "quadrupole_fraction": float(metadata["quadrupole_fraction"]),
                    "odd_k1_response_rms": odd_rms,
                    "even_k1_response_rms": even_rms,
                    "even_to_odd_rms_ratio": even_to_odd,
                    "odd_response_rms_snr": odd_snr,
                    "even_response_rms_snr": even_snr,
                    "material_ratio_and_detectable": int(
                        even_to_odd >= args.material_even_to_odd and even_snr >= 1.0
                    ),
                }
            )

    if len(output_rows) != 380:
        raise RuntimeError(f"Expected 380 target-candidate rows, found {len(output_rows)}")
    output_rows.sort(key=lambda row: (str(row["sextupole"]), str(row["quadrupole"])))
    output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(output_dir / "k1_nonlinearity_by_target_candidate.csv", output_rows)
    ratios = [float(row["even_to_odd_rms_ratio"]) for row in output_rows]
    odd_snrs = [float(row["odd_response_rms_snr"]) for row in output_rows]
    even_snrs = [float(row["even_response_rms_snr"]) for row in output_rows]
    maximum_row = max(output_rows, key=lambda row: float(row["even_to_odd_rms_ratio"]))
    summary = {
        "format": "cesr-repaired-lattice-k1-symmetric-nonlinearity-v1",
        "target_count": 76,
        "target_candidate_count": 380,
        "quadrupole_fraction": float(output_rows[0]["quadrupole_fraction"]),
        "definition": {
            "odd": "0.5 * (K1_plus_response - K1_minus_response)",
            "even": "0.5 * (K1_plus_response + K1_minus_response) - nominal_response",
            "response": "five-point K2 slope of the 9-bump by 111-BPM by 2-plane closed orbit",
        },
        "even_to_odd_rms_ratio": distribution(ratios),
        "odd_response_rms_snr": distribution(odd_snrs),
        "even_response_rms_snr": distribution(even_snrs),
        "material_even_to_odd_threshold": args.material_even_to_odd,
        "detectable_even_snr_threshold": 1.0,
        "material_and_detectable_count": sum(
            int(row["material_ratio_and_detectable"]) for row in output_rows
        ),
        "maximum_ratio_pair": {
            "sextupole": maximum_row["sextupole"],
            "quadrupole": maximum_row["quadrupole"],
            "even_to_odd_rms_ratio": maximum_row["even_to_odd_rms_ratio"],
            "even_response_rms_snr": maximum_row["even_response_rms_snr"],
        },
        "limitations": [
            "The SNR calculation uses provisional independent 5 um raw BPM-plane noise.",
            "Only the symmetric three-point K1 curvature is tested; cubic K1 dependence requires additional amplitudes.",
            "Closed-orbit K2 slopes are tested here; direct launch-trajectory phase/coupling channels are not yet included.",
        ],
    }
    (output_dir / "k1_nonlinearity_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

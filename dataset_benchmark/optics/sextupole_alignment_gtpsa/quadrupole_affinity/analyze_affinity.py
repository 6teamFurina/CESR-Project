#!/usr/bin/env python3
"""Calculate nuisance-marginalized information gain and precision affinity."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import tomllib
from pathlib import Path
from typing import Any, Iterable

import numpy as np

HERE = Path(__file__).resolve().parent
LATEST_RESULTS = HERE / "results" / "scibmad_latest"

NOISE_BY_OBSERVABLE = {
    "trajectory_x": 5.0e-6,
    "trajectory_y": 5.0e-6,
    "x_from_x_probe": 5.0e-6,
    "x_from_px_probe": 5.0e-6,
    "x_from_y_probe": 5.0e-6,
    "x_from_py_probe": 5.0e-6,
    "y_from_x_probe": 5.0e-6,
    "y_from_px_probe": 5.0e-6,
    "y_from_y_probe": 5.0e-6,
    "y_from_py_probe": 5.0e-6,
}


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--response-dir", type=Path, default=LATEST_RESULTS / "responses")
    result.add_argument("--output-dir", type=Path, default=LATEST_RESULTS / "affinity")
    result.add_argument("--nuisance-rms-m", type=float, default=3.0e-4)
    result.add_argument("--k2-step-m3", type=float, default=0.01)
    result.add_argument("--k2-levels", default="-2,-1,0,1,2")
    return result


def write_csv(path: Path, rows: Iterable[dict[str, Any]]) -> None:
    rows = list(rows)
    if not rows:
        raise ValueError(f"Refusing to write empty CSV: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def observable_from_label(label: str) -> str:
    return label.split(":", 1)[1]


def slope_noise(labels: list[str], k2_step: float, levels: list[float]) -> np.ndarray:
    sum_squared_k2 = sum((level * k2_step) ** 2 for level in levels)
    if sum_squared_k2 <= 0.0:
        raise ValueError("K2 grid cannot estimate a slope")
    return np.asarray(
        [NOISE_BY_OBSERVABLE[observable_from_label(label)] / math.sqrt(sum_squared_k2) for label in labels],
        dtype=float,
    )


def positive_information(matrix: np.ndarray) -> np.ndarray:
    symmetric = 0.5 * (matrix + matrix.T)
    values, vectors = np.linalg.eigh(symmetric)
    scale = max(float(np.max(np.abs(values))), 1.0)
    floor = 1.0e-12 * scale
    return (vectors * np.maximum(values, floor)) @ vectors.T


def nuisance_inverse(
    whitened_nuisance: np.ndarray,
    block_count: int,
    nuisance_rms: float,
) -> np.ndarray:
    information = block_count * (whitened_nuisance.T @ whitened_nuisance)
    information.flat[:: information.shape[0] + 1] += 1.0 / nuisance_rms**2
    return np.linalg.inv(information)


def marginalized_information(
    target_blocks: list[np.ndarray],
    whitened_nuisance: np.ndarray,
    sigma: np.ndarray,
    inverse_nuisance_information: np.ndarray,
) -> np.ndarray:
    whitened_targets = [block / sigma[:, None] for block in target_blocks]
    f_cc = sum((target.T @ target for target in whitened_targets), start=np.zeros((2, 2)))
    f_cn = sum(
        (target.T @ whitened_nuisance for target in whitened_targets),
        start=np.zeros((2, whitened_nuisance.shape[1])),
    )
    schur = f_cc - f_cn @ inverse_nuisance_information @ f_cn.T
    return positive_information(schur)


def covariance_metrics(information: np.ndarray) -> tuple[np.ndarray, float, float]:
    covariance = np.linalg.inv(information)
    sigmas = np.sqrt(np.maximum(np.diag(covariance), 0.0))
    worst_axis = float(np.max(sigmas))
    worst_eigen = float(math.sqrt(np.max(np.linalg.eigvalsh(covariance))))
    return sigmas, worst_axis, worst_eigen


def safe_name(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", name).lower()


def text_lines(path: Path) -> list[str]:
    return [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def target_bundles(response_dir: Path) -> list[Path]:
    directory_bundles = sorted(
        path for path in (response_dir / "targets").glob("*_responses") if path.is_dir()
    )
    if directory_bundles:
        return directory_bundles
    return sorted((response_dir / "targets").glob("*_responses.npz"))


def load_target_bundle(
    path: Path, response_dir: Path
) -> tuple[str, list[str], np.ndarray, np.ndarray, list[str], np.ndarray, np.ndarray]:
    if path.is_file():
        with np.load(path) as saved:
            return (
                path.name.removesuffix("_responses.npz").upper(),
                [str(value) for value in saved["observation_labels"]],
                np.asarray(saved["target_response_nominal"], dtype=float),
                np.asarray(saved["nuisance_response_nominal"], dtype=float),
                [str(value) for value in saved["candidate_names"]],
                np.asarray(saved["target_response_candidate_plus"], dtype=float),
                np.asarray(saved["target_response_candidate_minus"], dtype=float),
            )

    target = path.name.removesuffix("_responses").upper()
    labels = text_lines(response_dir / "observation_labels.txt")
    nominal = np.load(path / "target_response_nominal.npy")
    nuisance = np.load(path / "nuisance_response_nominal.npy")
    candidates = text_lines(path / "candidate_names.txt")
    plus = np.stack(
        [np.load(path / f"candidate_{safe_name(name)}_plus.npy") for name in candidates]
    )
    minus = np.stack(
        [np.load(path / f"candidate_{safe_name(name)}_minus.npy") for name in candidates]
    )
    return target, labels, nominal, nuisance, candidates, plus, minus


def main() -> int:
    args = parser().parse_args()
    if not math.isfinite(args.nuisance_rms_m) or args.nuisance_rms_m <= 0.0:
        raise ValueError("--nuisance-rms-m must be positive and finite")
    levels = [float(value) for value in args.k2_levels.split(",")]
    response_dir = args.response_dir.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    screen_rows: dict[tuple[str, str], dict[str, str]] = {}
    with (response_dir / "quadrupole_optics_screen.csv").open(encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream):
            screen_rows[(row["sextupole"], row["quadrupole"])] = row

    rows: list[dict[str, Any]] = []
    target_files = target_bundles(response_dir)
    if not target_files:
        raise FileNotFoundError(f"No target response files under {response_dir / 'targets'}")
    nuisance_counts: set[int] = set()
    observation_counts: set[int] = set()
    for target_path in target_files:
        target, labels, nominal, nuisance, candidates, plus, minus = load_target_bundle(
            target_path, response_dir
        )
        observation_counts.add(len(labels))
        nuisance_counts.add(nuisance.shape[1])
        sigma = slope_noise(labels, args.k2_step_m3, levels)
        whitened_nuisance = nuisance / sigma[:, None]
        inverse_nuisance_one = nuisance_inverse(
            whitened_nuisance, 1, args.nuisance_rms_m
        )
        inverse_nuisance_three = nuisance_inverse(
            whitened_nuisance, 3, args.nuisance_rms_m
        )
        f0 = marginalized_information(
            [nominal], whitened_nuisance, sigma, inverse_nuisance_one
        )
        sigma0, worst_axis0, worst_eigen0 = covariance_metrics(f0)
        sign0, logdet0 = np.linalg.slogdet(f0)
        if sign0 <= 0.0:
            raise RuntimeError(f"Non-positive baseline information determinant for {target}")
        for index, quadrupole in enumerate(candidates):
            fq = marginalized_information(
                [nominal, plus[index], minus[index]],
                whitened_nuisance,
                sigma,
                inverse_nuisance_three,
            )
            sigmaq, worst_axisq, worst_eigenq = covariance_metrics(fq)
            signq, logdetq = np.linalg.slogdet(fq)
            if signq <= 0.0:
                raise RuntimeError(f"Non-positive candidate information determinant for {target}/{quadrupole}")
            source = screen_rows[(target, quadrupole)]
            rows.append(
                {
                    "sextupole": target,
                    "sextupole_s_m": source["sextupole_s_m"],
                    "quadrupole": quadrupole,
                    "quadrupole_index": source.get(
                        "quadrupole_inventory_index", source.get("quadrupole_index", "")
                    ),
                    "quadrupole_s_m": source.get("quadrupole_s_m", ""),
                    "screen_rank": source["selected_rank"],
                    "optics_leverage": source["optics_leverage"],
                    "information_gain_logdet": float(logdetq - logdet0),
                    "precision_improvement_worst_axis": float(worst_axis0 / worst_axisq),
                    "precision_improvement_worst_direction": float(worst_eigen0 / worst_eigenq),
                    "baseline_sigma_x_um": float(1.0e6 * sigma0[0]),
                    "baseline_sigma_y_um": float(1.0e6 * sigma0[1]),
                    "candidate_sigma_x_um": float(1.0e6 * sigmaq[0]),
                    "candidate_sigma_y_um": float(1.0e6 * sigmaq[1]),
                    "max_abs_tune_shift": source["max_abs_tune_shift"],
                    "max_detector_beta_beating": source["max_detector_beta_beating"],
                }
            )
    write_csv(output_dir / "quadrupole_affinity_scores.csv", rows)

    summaries: list[dict[str, Any]] = []
    for target in sorted({str(row["sextupole"]) for row in rows}):
        target_rows = [row for row in rows if row["sextupole"] == target]
        best_info = max(target_rows, key=lambda row: float(row["information_gain_logdet"]))
        best_precision = max(
            target_rows, key=lambda row: float(row["precision_improvement_worst_axis"])
        )
        summaries.append(
            {
                "sextupole": target,
                "best_information_quadrupole": best_info["quadrupole"],
                "best_information_gain_logdet": best_info["information_gain_logdet"],
                "best_precision_quadrupole": best_precision["quadrupole"],
                "best_precision_improvement": best_precision["precision_improvement_worst_axis"],
                "baseline_sigma_x_um": best_info["baseline_sigma_x_um"],
                "baseline_sigma_y_um": best_info["baseline_sigma_y_um"],
            }
        )
    write_csv(output_dir / "best_quadrupoles_by_sextupole.csv", summaries)
    if len(nuisance_counts) != 1 or len(observation_counts) != 1:
        raise RuntimeError(
            f"Inconsistent response shapes: nuisance={nuisance_counts}, observations={observation_counts}"
        )
    response_engine = "unknown"
    candidate_metadata_path = response_dir / "candidate_metadata.toml"
    if candidate_metadata_path.exists():
        response_engine = tomllib.loads(candidate_metadata_path.read_text(encoding="utf-8"))[
            "engine"
        ]
    nuisance_count = next(iter(nuisance_counts))
    metadata = {
        "format": "cesr-sextupole-quadrupole-affinity-v2",
        "target_count": len(summaries),
        "retained_pair_count": len(rows),
        "response_engine": response_engine,
        "response_directory": str(response_dir),
        "observation_count": next(iter(observation_counts)),
        "nuisance_model": f"{nuisance_count} nominal own Kn2-offset response columns from the other active sextupoles with independent Gaussian offset priors",
        "nuisance_rms_m": args.nuisance_rms_m,
        "information_definition": "logdet(F_nominal+candidate_plus+candidate_minus)-logdet(F_nominal), after Schur marginalization",
        "precision_definition": "max(sigma_x_nominal,sigma_y_nominal)/max(sigma_x_candidate,sigma_y_candidate)",
        "measurement_noise": NOISE_BY_OBSERVABLE,
        "k2_step_m3": args.k2_step_m3,
        "k2_levels": levels,
        "slope_variance": "sigma_observable^2 / sum(delta_K2_level^2)",
        "measurement_protocol": "fixed reference launch plus x=1 mm, px=0.1 mrad, y=1 mm, py=0.1 mrad probes; 5 um BPM noise per trajectory sample",
        "bump_protocol": "nominal launch pre-screen; no separate 3x3 closed-orbit bump replication included",
        "limitations": [
            "Noise levels are provisional diagonal measured-style assumptions, not an archived CESR covariance matrix.",
            "Nuisance response columns are calculated at nominal quadrupole optics and reused in candidate +/- blocks.",
            "The heatmaps rank single quadrupoles; greedy multi-knob complementarity is a subsequent calculation.",
        ],
    }
    (output_dir / "affinity_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )
    print(f"Affinity pairs: {len(rows)}")
    print(f"Targets: {len(summaries)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

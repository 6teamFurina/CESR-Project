#!/usr/bin/env python3
"""Validate the repaired-lattice SciBmad affinity response and score artifacts."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
LATEST_RESULTS = HERE / "results" / "scibmad_latest"


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--response-dir", type=Path, default=LATEST_RESULTS / "responses")
    result.add_argument("--affinity-dir", type=Path, default=LATEST_RESULTS / "affinity")
    result.add_argument("--allow-incomplete", action="store_true")
    return result


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def safe_name(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", name).lower()


def main() -> int:
    args = parser().parse_args()
    response_dir = args.response_dir.expanduser().resolve()
    affinity_dir = args.affinity_dir.expanduser().resolve()
    screen = rows(response_dir / "quadrupole_optics_screen.csv")
    selected = rows(response_dir / "selected_candidates.csv")
    sextupoles = {row["sextupole"] for row in screen}
    quadrupoles = {row["quadrupole"] for row in screen}
    assert len(sextupoles) == 76
    assert len(quadrupoles) == 113
    assert len(screen) == 76 * 113
    assert len(selected) == 76 * 15
    for sextupole in sextupoles:
        target_rows = [row for row in selected if row["sextupole"] == sextupole]
        assert sorted(int(row["selected_rank"]) for row in target_rows) == list(range(1, 16))

    labels = [
        line.strip()
        for line in (response_dir / "observation_labels.txt").read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    assert len(labels) == 1110
    assert len(set(labels)) == len(labels)
    target_directories = sorted(
        path for path in (response_dir / "targets").glob("*_responses") if path.is_dir()
    )
    if not args.allow_incomplete:
        assert len(target_directories) == 76

    complete_candidates = 0
    finite_target_norms: list[float] = []
    for directory in target_directories:
        nominal_path = directory / "target_response_nominal.npy"
        nuisance_path = directory / "nuisance_response_nominal.npy"
        if not nominal_path.exists() or not nuisance_path.exists():
            if args.allow_incomplete:
                continue
            raise FileNotFoundError(f"Incomplete nominal bundle: {directory}")
        nominal = np.load(nominal_path)
        nuisance = np.load(nuisance_path)
        assert nominal.shape == (1110, 2)
        assert nuisance.shape == (1110, 150)
        assert np.isfinite(nominal).all()
        assert np.isfinite(nuisance).all()
        finite_target_norms.append(float(np.linalg.norm(nominal)))
        candidate_names_path = directory / "candidate_names.txt"
        if not candidate_names_path.exists():
            if args.allow_incomplete:
                continue
            raise FileNotFoundError(f"Missing candidate list: {directory}")
        candidates = [
            line.strip()
            for line in candidate_names_path.read_text(encoding="utf-8").splitlines()
            if line.strip()
        ]
        assert len(candidates) == 15
        for name in candidates:
            for sign in ("plus", "minus"):
                path = directory / f"candidate_{safe_name(name)}_{sign}.npy"
                if not path.exists():
                    if args.allow_incomplete:
                        break
                    raise FileNotFoundError(path)
                response = np.load(path)
                assert response.shape == (1110, 2)
                assert np.isfinite(response).all()
            else:
                continue
            break
        else:
            complete_candidates += 1

    score_path = affinity_dir / "quadrupole_affinity_scores.csv"
    score_count = 0
    if score_path.exists():
        scores = rows(score_path)
        score_count = len(scores)
        if not args.allow_incomplete:
            assert score_count == 76 * 15
        for row in scores:
            for field in (
                "information_gain_logdet",
                "precision_improvement_worst_axis",
                "baseline_sigma_x_um",
                "baseline_sigma_y_um",
                "candidate_sigma_x_um",
                "candidate_sigma_y_um",
            ):
                assert math.isfinite(float(row[field]))

    summary = {
        "screen_pairs": len(screen),
        "selected_pairs": len(selected),
        "sextupoles": len(sextupoles),
        "independent_quadrupoles": len(quadrupoles),
        "observations": len(labels),
        "nominal_target_bundles": len(finite_target_norms),
        "complete_candidate_bundles": complete_candidates,
        "score_rows": score_count,
    }
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Compare the P=1-per-corrector SciBmad Jacobian with Bmad and wide SciBmad."""

from __future__ import annotations

import json
from pathlib import Path

from compare_corrector_jacobians import compare_pair

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"


def main() -> int:
    pairs = {
        "detector_optics": (
            "detector_optics_jacobian.csv",
            "scibmad_per_corrector_detector_optics_jacobian.csv",
        ),
        "closed_orbit": (
            "closed_orbit_jacobian_6x119.csv",
            "scibmad_per_corrector_closed_orbit_jacobian.csv",
        ),
        "ring_tunes": (
            "ring_tune_jacobian.csv",
            "scibmad_per_corrector_ring_tune_jacobian.csv",
        ),
    }
    per_dir = RESULTS / "scibmad_per_corrector"
    comparison: dict[str, object] = {
        "per_corrector_vs_bmad": {},
        "per_corrector_vs_wide_scibmad": {},
    }
    for name, (common_suffix, per_filename) in pairs.items():
        per_path = per_dir / per_filename
        comparison["per_corrector_vs_bmad"][name] = compare_pair(
            RESULTS / "bmad" / f"bmad_{common_suffix}",
            per_path,
        )
        comparison["per_corrector_vs_wide_scibmad"][name] = compare_pair(
            RESULTS / "scibmad" / f"scibmad_{common_suffix}",
            per_path,
        )

    output = per_dir / "comparison.json"
    output.write_text(
        json.dumps(comparison, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    for reference_name, matrices in comparison.items():
        print(reference_name)
        for matrix_name, details in matrices.items():
            metrics = details["overall"]
            print(
                f"  {matrix_name}: relative Frobenius="
                f"{metrics['relative_frobenius']:.6e}, correlation="
                f"{metrics['cosine_correlation']:.12f}"
            )
    print(f"Comparison: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

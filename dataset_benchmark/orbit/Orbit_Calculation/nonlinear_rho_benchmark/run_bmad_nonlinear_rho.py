#!/usr/bin/env python3
"""Run the shared nonlinear-rho orbit inputs through one persistent Tao process."""

from __future__ import annotations

import csv
import json
import platform
import resource
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Any


HERE = Path(__file__).resolve().parent
CALCULATION_DIR = HERE.parent
ORBIT_ROOT = CALCULATION_DIR.parent
PROJECT_ROOT = HERE.parents[3]
sys.path.insert(0, str(CALCULATION_DIR))
sys.path.insert(0, str(PROJECT_ROOT / "test_codes"))

import benchmark_bmad as common  # noqa: E402


def selected_ring() -> str:
    ring = "latest"
    for argument in sys.argv[1:]:
        if argument.startswith("--ring="):
            ring = argument.split("=", 1)[1].lower()
    ring in {"latest", "legacy"} or raise_value_error()
    return ring


def raise_value_error() -> None:
    raise ValueError("--ring must be latest or legacy")


RING = selected_ring()
ARTIFACT_RING = "latest_cesr" if RING == "latest" else "legacy"
INPUT_PATH = HERE / "shared_input" / ARTIFACT_RING / "nonlinear_rho_correctors.csv"
MANIFEST_PATH = HERE / "shared_input" / ARTIFACT_RING / "sample_manifest.csv"
LATTICE_PATH = (
    PROJECT_ROOT / "Latest_Lattice" / "lat.bmad"
    if RING == "latest"
    else ORBIT_ROOT / "reference" / "cesr_bmad_compatible.bmad"
)
RESULT_DIR = HERE / "results" / ARTIFACT_RING / (
    "bmad_reference" if RING == "latest" else "bmad"
)
SAMPLES_FILENAME = "bmad_rf_on_samples.csv" if RING == "latest" else "bmad_samples.csv"


def read_manifest(path: Path) -> dict[int, tuple[str, float, int]]:
    result: dict[int, tuple[str, float, int]] = {}
    with path.open(encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream):
            result[int(row["sample_id"])] = (
                row["scenario"],
                float(row["rho"]),
                int(row["trial_id"]),
            )
    return result


def write_sample_timings(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_group_timings(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    sample_ids, names, samples = common.read_samples(INPUT_PATH)
    manifest = read_manifest(MANIFEST_PATH)
    if sample_ids != list(manifest):
        raise RuntimeError("Input and manifest sample IDs differ")

    detectors, horizontal, vertical = common.response_tools.parse_cesr_layout(
        LATTICE_PATH
    )
    if names != horizontal + vertical:
        raise RuntimeError("Input control order differs from the Bmad lattice order")
    references = common.variable_references(horizontal, vertical)

    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    init_path = RESULT_DIR / "tao_nonlinear_rho.init"
    init_path.write_text(
        common.response_tools.build_tao_init(
            LATTICE_PATH,
            detectors,
            horizontal,
            vertical,
            derivative_step=1.0e-6,
        ),
        encoding="utf-8",
    )

    try:
        from pytao import Tao
    except ImportError as exc:
        raise RuntimeError("Run this script in the Ubuntu-Bmad PyTao environment") from exc

    initialize_start = time.perf_counter()
    tao = Tao(init_file=str(init_path), noplot=True)
    common.enable_rf(tao)
    common.activate_benchmark_data(tao)
    initialization_seconds = time.perf_counter() - initialize_start

    x_rows = common.response_tools.active_data(tao, "x")
    y_rows = common.response_tools.active_data(tao, "y")
    labels = [f"{row['ele_name'].upper()}:x" for row in x_rows]
    labels += [f"{row['ele_name'].upper()}:y" for row in y_rows]
    warmup_start = time.perf_counter()
    for sample in samples[1:3]:
        common.apply_sample(tao, references, sample)
        common.read_observables(tao)
    common.apply_sample(tao, references, samples[0])
    common.read_observables(tao)
    warmup_seconds = time.perf_counter() - warmup_start

    observables: list[list[float]] = []
    converged: list[bool] = []
    timing_rows: list[dict[str, Any]] = []
    group_accumulator: dict[tuple[str, float], dict[str, float]] = defaultdict(
        lambda: {
            "samples": 0.0,
            "converged": 0.0,
            "update_seconds": 0.0,
            "read_seconds": 0.0,
            "physics_seconds": 0.0,
        }
    )

    physics_start = time.perf_counter()
    for sample_index, (sample_id, sample) in enumerate(zip(sample_ids, samples)):
        update_start = time.perf_counter()
        common.apply_sample(tao, references, sample)
        update_seconds = time.perf_counter() - update_start

        read_start = time.perf_counter()
        values, good = common.read_observables(tao)
        read_seconds = time.perf_counter() - read_start
        sample_seconds = update_seconds + read_seconds

        observables.append(values)
        converged.append(good)
        scenario, rho, trial_id = manifest[sample_id]
        timing_rows.append(
            {
                "sample_id": sample_id,
                "scenario": scenario,
                "rho": f"{rho:.17g}",
                "trial_id": trial_id,
                "converged": str(good).lower(),
                "update_seconds": f"{update_seconds:.17g}",
                "read_seconds": f"{read_seconds:.17g}",
                "physics_seconds": f"{sample_seconds:.17g}",
            }
        )
        if scenario != "baseline":
            group = group_accumulator[(scenario, rho)]
            group["samples"] += 1
            group["converged"] += int(good)
            group["update_seconds"] += update_seconds
            group["read_seconds"] += read_seconds
            group["physics_seconds"] += sample_seconds
        if (sample_index + 1) % 100 == 0:
            print(f"Bmad completed {sample_index + 1}/{len(samples)} samples", flush=True)
    physics_seconds = time.perf_counter() - physics_start

    output_path = RESULT_DIR / SAMPLES_FILENAME
    timing_path = RESULT_DIR / "bmad_sample_timings.csv"
    group_path = RESULT_DIR / "bmad_group_timings.csv"
    common.write_outputs(output_path, sample_ids, labels, observables, converged)
    write_sample_timings(timing_path, timing_rows)

    group_rows: list[dict[str, Any]] = []
    for (scenario, rho), values in group_accumulator.items():
        sample_count = int(values["samples"])
        group_rows.append(
            {
                "scenario": scenario,
                "rho": f"{rho:.17g}",
                "samples": sample_count,
                "converged": int(values["converged"]),
                "update_seconds": f"{values['update_seconds']:.17g}",
                "read_seconds": f"{values['read_seconds']:.17g}",
                "physics_seconds": f"{values['physics_seconds']:.17g}",
                "samples_per_physics_second": f"{sample_count / values['physics_seconds']:.17g}",
            }
        )
    write_group_timings(group_path, group_rows)

    version = tao.cmd("show version", raises=False)
    metadata = {
        "format": "ring-nonlinear-rho-bmad-v2",
        "ring": ARTIFACT_RING,
        "engine": "Bmad/Tao/PyTao",
        "input_csv": str(INPUT_PATH),
        "lattice": str(LATTICE_PATH),
        "sample_count": len(samples),
        "converged_count": sum(converged),
        "failed_count": len(converged) - sum(converged),
        "initialization_seconds": initialization_seconds,
        "warmup_seconds": warmup_seconds,
        "physics_seconds": physics_seconds,
        "samples_per_second": len(samples) / physics_seconds,
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "tao_version_raw": version,
        "maximum_resident_set_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        "execution_model": (
            f"one persistent Tao instance; {len(names)} variable commands batched with "
            "suppress_lattice_calc=True; one model recalculation per sample"
        ),
        "timed_region": "variable update + Tao lattice recalculation + observable read",
        "convergence_indicator": (
            "PyTao orbit data good_model flags and finite model values; Tao does not "
            "expose SciBmad's explicit one-turn closure norm in this path"
        ),
    }
    (RESULT_DIR / "bmad_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )
    print(
        f"Bmad physics: {physics_seconds:.3f} s, "
        f"{len(samples) / physics_seconds:.3f} samples/s, "
        f"converged {sum(converged)}/{len(converged)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

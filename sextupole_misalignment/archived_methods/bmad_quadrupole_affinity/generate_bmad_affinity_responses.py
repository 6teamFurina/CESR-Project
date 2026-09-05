#!/usr/bin/env python3
"""Generate quadrupole-conditioned sextupole-center response dictionaries.

The calculation is deliberately staged:

1. Apply a symmetric fractional K1 change to every active quadrupole and rank
   its finite change of beta/phase at every active sextupole.
2. Retain only the top candidates for each target sextupole.
3. For every target, finite-difference the target K2 slope with respect to its
   own x/y offset at nominal and candidate-quadrupole +/- settings.
4. At nominal optics, finite-difference the same target K2 slope with respect
   to the x/y offsets of all other sextupoles.  These 150 columns are the
   explicit linear nuisance dictionary used by the affinity analysis.

The expensive response stage is resumable one target at a time.  Bmad/Tao is
used through the validated Ubuntu-Bmad PyTao environment.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import platform
import sys
import time
from pathlib import Path
from typing import Any, Iterable

import numpy as np

HERE = Path(__file__).resolve().parent
ALIGNMENT_DIR = HERE.parent.parent
PROJECT_DIR = ALIGNMENT_DIR.parent
OPTICS_DIR = PROJECT_DIR / "optics"
ORBIT_DIR = PROJECT_DIR / "orbit"
ORBIT_CALCULATION_DIR = ORBIT_DIR / "Orbit_Calculation"
CORRECTOR_JACOBIAN_DIR = OPTICS_DIR / "corrector_jacobian_benchmark"

for directory in (OPTICS_DIR, ORBIT_CALCULATION_DIR, CORRECTOR_JACOBIAN_DIR):
    sys.path.insert(0, str(directory))

from audit_bmad_extended_inventory import dictionary_float  # noqa: E402
from benchmark_bmad import (  # noqa: E402
    activate_benchmark_data,
    apply_sample,
    read_samples,
    response_tools,
    variable_references,
)
from benchmark_bmad_extended_jacobian import strength_inventory  # noqa: E402
from benchmark_bmad_optics import (  # noqa: E402
    collect_state,
    detector_inventory,
    lat_values,
    numeric_lat_values,
)

DETECTOR_COLUMNS = (
    "orbit_x",
    "orbit_y",
    "phi_1",
    "phi_2",
    "c11",
    "c12",
    "c21",
    "c22",
)
RING_COLUMNS = ("tune_1", "tune_2")


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument(
        "--inputs",
        type=Path,
        default=ORBIT_CALCULATION_DIR / "inputs" / "cesr_corrector_samples_1000.csv",
    )
    result.add_argument(
        "--lattice",
        type=Path,
        default=ORBIT_DIR / "reference" / "cesr_bmad_compatible.bmad",
    )
    result.add_argument("--output-dir", type=Path, default=HERE / "results" / "responses")
    result.add_argument("--top-k", type=int, default=15)
    result.add_argument("--quadrupole-fraction", type=float, default=0.001)
    result.add_argument("--k2-step-m3", type=float, default=0.01)
    result.add_argument("--offset-step-m", type=float, default=1.0e-4)
    result.add_argument("--max-tune-shift", type=float, default=0.01)
    result.add_argument("--max-beta-beating", type=float, default=0.20)
    result.add_argument(
        "--targets",
        default="all",
        help="Comma-separated target names or 'all'.",
    )
    result.add_argument("--overwrite", action="store_true")
    result.add_argument("--screen-only", action="store_true")
    return result


def validate_args(args: argparse.Namespace) -> None:
    if not 10 <= args.top_k <= 20:
        raise ValueError("--top-k must be between 10 and 20")
    for value, name in (
        (args.quadrupole_fraction, "--quadrupole-fraction"),
        (args.k2_step_m3, "--k2-step-m3"),
        (args.offset_step_m, "--offset-step-m"),
        (args.max_tune_shift, "--max-tune-shift"),
        (args.max_beta_beating, "--max-beta-beating"),
    ):
        if not math.isfinite(value) or value <= 0.0:
            raise ValueError(f"{name} must be positive and finite")


def write_csv(path: Path, rows: Iterable[dict[str, Any]]) -> None:
    rows = list(rows)
    if not rows:
        raise ValueError(f"Refusing to write empty CSV: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def wrap_turn_difference(high: np.ndarray, low: np.ndarray) -> np.ndarray:
    return (high - low + 0.5) % 1.0 - 0.5


def wrap_radian_difference(high: np.ndarray, low: np.ndarray) -> np.ndarray:
    return (high - low + math.pi) % (2.0 * math.pi) - math.pi


def measurement_layout(detector_names: list[str]) -> tuple[list[str], np.ndarray]:
    labels: list[str] = []
    phase_mask: list[bool] = []
    for detector in detector_names:
        for observable in DETECTOR_COLUMNS:
            labels.append(f"{detector.lower()}:{observable}")
            phase_mask.append(observable.startswith("phi_"))
    labels.extend(("ring:tune_1", "ring:tune_2"))
    phase_mask.extend((False, False))
    return labels, np.asarray(phase_mask, dtype=bool)


def flatten_measurements(state: dict[str, Any]) -> np.ndarray:
    columns = state["columns"]
    relative_phase: dict[str, list[float]] = {}
    for observable in ("phi_1", "phi_2"):
        values = np.asarray(columns[observable], dtype=float)
        relative_phase[observable] = list(wrap_turn_difference(values, values[0]))
    result: list[float] = []
    for detector in range(99):
        for observable in DETECTOR_COLUMNS:
            source = relative_phase.get(observable, columns[observable])
            result.append(float(source[detector]))
    result.extend((float(state["Q1"]), float(state["Q2"])))
    return np.asarray(result, dtype=float)


def central_state_difference(
    plus: dict[str, Any],
    minus: dict[str, Any],
    step: float,
    phase_mask: np.ndarray,
) -> np.ndarray:
    high = flatten_measurements(plus)
    low = flatten_measurements(minus)
    difference = high - low
    difference[phase_mask] = wrap_turn_difference(high[phase_mask], low[phase_mask])
    return difference / (2.0 * step)


def set_element(tao: Any, name: str, attribute: str, value: float) -> None:
    tao.cmds(
        [f"set ele {name} {attribute} = {value:.17g}"],
        suppress_lattice_calc=True,
        suppress_plotting=True,
    )


def collect(tao: Any, detector_data: dict[str, Any]) -> dict[str, Any]:
    state, _, _ = collect_state(tao, detector_data, 0.0)
    return state


def k2_slope_at_offset(
    tao: Any,
    detector_data: dict[str, Any],
    target: dict[str, Any],
    offset_element: str,
    offset_attribute: str,
    offset_value: float,
    k2_step: float,
    phase_mask: np.ndarray,
) -> np.ndarray:
    set_element(tao, offset_element, offset_attribute, offset_value)
    set_element(tao, target["name"], "K2", target["baseline"] + k2_step)
    plus = collect(tao, detector_data)
    set_element(tao, target["name"], "K2", target["baseline"] - k2_step)
    minus = collect(tao, detector_data)
    set_element(tao, target["name"], "K2", target["baseline"])
    set_element(tao, offset_element, offset_attribute, 0.0)
    return central_state_difference(plus, minus, k2_step, phase_mask)


def mixed_k2_offset_response(
    tao: Any,
    detector_data: dict[str, Any],
    target: dict[str, Any],
    offset_element: str,
    offset_attribute: str,
    k2_step: float,
    offset_step: float,
    phase_mask: np.ndarray,
) -> np.ndarray:
    plus = k2_slope_at_offset(
        tao,
        detector_data,
        target,
        offset_element,
        offset_attribute,
        offset_step,
        k2_step,
        phase_mask,
    )
    minus = k2_slope_at_offset(
        tao,
        detector_data,
        target,
        offset_element,
        offset_attribute,
        -offset_step,
        k2_step,
        phase_mask,
    )
    return (plus - minus) / (2.0 * offset_step)


def read_sextupole_optics(tao: Any, expected: int) -> dict[str, np.ndarray | list[str]]:
    names_raw, _ = lat_values(tao, "SEX_*", "ele.name", expected)
    names = [str(name).upper() for name in names_raw]
    result: dict[str, np.ndarray | list[str]] = {"names": names}
    for key, field in (
        ("beta_x", "ele.a.beta"),
        ("beta_y", "ele.b.beta"),
        ("phi_x", "ele.a.phi"),
        ("phi_y", "ele.b.phi"),
        ("s_m", "ele.s"),
    ):
        values, _ = numeric_lat_values(tao, "SEX_*", field, expected)
        result[key] = np.asarray(values, dtype=float)
    for key, field in (("phi_ref_x", "ele.a.phi"), ("phi_ref_y", "ele.b.phi")):
        values, _ = numeric_lat_values(tao, "DET_00W", field, 1)
        result[key] = np.asarray(values, dtype=float)
    return result


def normalized_relative_phases(optics: dict[str, np.ndarray | list[str]]) -> tuple[np.ndarray, np.ndarray]:
    phi_x = np.asarray(optics["phi_x"], dtype=float)
    phi_y = np.asarray(optics["phi_y"], dtype=float)
    ref_x = float(np.asarray(optics["phi_ref_x"])[0])
    ref_y = float(np.asarray(optics["phi_ref_y"])[0])
    return (
        wrap_radian_difference(phi_x, np.full_like(phi_x, ref_x)),
        wrap_radian_difference(phi_y, np.full_like(phi_y, ref_y)),
    )


def screen_quadrupoles(
    tao: Any,
    detector_data: dict[str, Any],
    quadrupoles: list[dict[str, Any]],
    sextupoles: list[dict[str, Any]],
    fraction: float,
    top_k: int,
    max_tune_shift: float,
    max_beta_beating: float,
    output_dir: Path,
) -> tuple[list[dict[str, Any]], dict[str, list[str]]]:
    baseline_state = collect(tao, detector_data)
    baseline_detector_beta_x = np.asarray(baseline_state["columns"]["beta_1"], dtype=float)
    baseline_detector_beta_y = np.asarray(baseline_state["columns"]["beta_2"], dtype=float)
    baseline_optics = read_sextupole_optics(tao, len(sextupoles))
    ring_sextupoles = list(baseline_optics["names"])
    inventory_names = [item["name"] for item in sextupoles]
    if set(ring_sextupoles) != set(inventory_names):
        raise RuntimeError("SEX_* optics inventory differs from the active sextupole inventory")
    index_by_name = {name: index for index, name in enumerate(ring_sextupoles)}

    rows: list[dict[str, Any]] = []
    for q_index, quadrupole in enumerate(quadrupoles, start=1):
        step = abs(float(quadrupole["baseline"])) * fraction
        set_element(tao, quadrupole["name"], "K1", quadrupole["baseline"] + step)
        plus_state = collect(tao, detector_data)
        plus_optics = read_sextupole_optics(tao, len(sextupoles))
        set_element(tao, quadrupole["name"], "K1", quadrupole["baseline"] - step)
        minus_state = collect(tao, detector_data)
        minus_optics = read_sextupole_optics(tao, len(sextupoles))
        set_element(tao, quadrupole["name"], "K1", quadrupole["baseline"])

        d_log_beta_x = 0.5 * (
            np.log(np.asarray(plus_optics["beta_x"]))
            - np.log(np.asarray(minus_optics["beta_x"]))
        )
        d_log_beta_y = 0.5 * (
            np.log(np.asarray(plus_optics["beta_y"]))
            - np.log(np.asarray(minus_optics["beta_y"]))
        )
        plus_phi_x, plus_phi_y = normalized_relative_phases(plus_optics)
        minus_phi_x, minus_phi_y = normalized_relative_phases(minus_optics)
        d_phi_x = 0.5 * wrap_radian_difference(plus_phi_x, minus_phi_x)
        d_phi_y = 0.5 * wrap_radian_difference(plus_phi_y, minus_phi_y)
        leverage = np.sqrt(d_log_beta_x**2 + d_log_beta_y**2 + d_phi_x**2 + d_phi_y**2)

        tune_shift = max(
            abs(float(plus_state["Q1"]) - float(baseline_state["Q1"])),
            abs(float(plus_state["Q2"]) - float(baseline_state["Q2"])),
            abs(float(minus_state["Q1"]) - float(baseline_state["Q1"])),
            abs(float(minus_state["Q2"]) - float(baseline_state["Q2"])),
        )
        beta_beating = max(
            float(np.max(np.abs(np.asarray(plus_state["columns"]["beta_1"]) / baseline_detector_beta_x - 1.0))),
            float(np.max(np.abs(np.asarray(plus_state["columns"]["beta_2"]) / baseline_detector_beta_y - 1.0))),
            float(np.max(np.abs(np.asarray(minus_state["columns"]["beta_1"]) / baseline_detector_beta_x - 1.0))),
            float(np.max(np.abs(np.asarray(minus_state["columns"]["beta_2"]) / baseline_detector_beta_y - 1.0))),
        )
        allowed = tune_shift <= max_tune_shift and beta_beating <= max_beta_beating
        for sextupole in sextupoles:
            s_index = index_by_name[sextupole["name"]]
            rows.append(
                {
                    "sextupole": sextupole["name"],
                    "sextupole_s_m": float(np.asarray(baseline_optics["s_m"])[s_index]),
                    "quadrupole": quadrupole["name"],
                    "quadrupole_index": q_index,
                    "quadrupole_k1_m2": quadrupole["baseline"],
                    "delta_k1_m2": step,
                    "optics_leverage": float(leverage[s_index]),
                    "delta_log_beta_x": float(d_log_beta_x[s_index]),
                    "delta_log_beta_y": float(d_log_beta_y[s_index]),
                    "delta_phi_x_rad": float(d_phi_x[s_index]),
                    "delta_phi_y_rad": float(d_phi_y[s_index]),
                    "max_abs_tune_shift": tune_shift,
                    "max_detector_beta_beating": beta_beating,
                    "allowed": int(allowed),
                }
            )
        if q_index % 10 == 0 or q_index == len(quadrupoles):
            print(f"Optics screen {q_index}/{len(quadrupoles)}: {quadrupole['name']}", flush=True)

    candidates: dict[str, list[str]] = {}
    for sextupole in sextupoles:
        target_rows = [row for row in rows if row["sextupole"] == sextupole["name"]]
        allowed_rows = [row for row in target_rows if bool(row["allowed"])]
        if len(allowed_rows) < top_k:
            raise RuntimeError(
                f"Only {len(allowed_rows)} allowed quadrupoles for {sextupole['name']}; need {top_k}"
            )
        allowed_rows.sort(key=lambda row: float(row["optics_leverage"]), reverse=True)
        selected = allowed_rows[:top_k]
        candidates[sextupole["name"]] = [str(row["quadrupole"]) for row in selected]
        rank = {str(row["quadrupole"]): index for index, row in enumerate(selected, start=1)}
        for row in target_rows:
            row["selected"] = int(str(row["quadrupole"]) in rank)
            row["selected_rank"] = rank.get(str(row["quadrupole"]), "")

    write_csv(output_dir / "quadrupole_optics_screen.csv", rows)
    (output_dir / "selected_candidates.json").write_text(
        json.dumps(candidates, indent=2) + "\n", encoding="utf-8"
    )
    return rows, candidates


def target_names(requested: str, sextupoles: list[dict[str, Any]]) -> list[str]:
    all_names = [item["name"] for item in sextupoles]
    if requested.strip().lower() == "all":
        return all_names
    selected = [name.strip().upper() for name in requested.split(",") if name.strip()]
    unknown = sorted(set(selected) - set(all_names))
    if unknown:
        raise ValueError(f"Unknown targets: {unknown}")
    return selected


def calculate_target_responses(
    tao: Any,
    detector_data: dict[str, Any],
    target: dict[str, Any],
    sextupoles: list[dict[str, Any]],
    quadrupole_by_name: dict[str, dict[str, Any]],
    candidates: list[str],
    k2_step: float,
    offset_step: float,
    quadrupole_fraction: float,
    labels: list[str],
    phase_mask: np.ndarray,
) -> dict[str, np.ndarray]:
    start = time.perf_counter()
    target_response = np.column_stack(
        [
            mixed_k2_offset_response(
                tao,
                detector_data,
                target,
                target["name"],
                attribute,
                k2_step,
                offset_step,
                phase_mask,
            )
            for attribute in ("X_OFFSET", "Y_OFFSET")
        ]
    )

    nuisance_labels: list[str] = []
    nuisance_columns: list[np.ndarray] = []
    other_sextupoles = [item for item in sextupoles if item["name"] != target["name"]]
    for nuisance_index, nuisance in enumerate(other_sextupoles, start=1):
        for attribute, plane in (("X_OFFSET", "x"), ("Y_OFFSET", "y")):
            nuisance_labels.append(f"{nuisance['name']}:{plane}_offset")
            nuisance_columns.append(
                mixed_k2_offset_response(
                    tao,
                    detector_data,
                    target,
                    nuisance["name"],
                    attribute,
                    k2_step,
                    offset_step,
                    phase_mask,
                )
            )
        if nuisance_index % 15 == 0 or nuisance_index == len(other_sextupoles):
            print(
                f"  {target['name']} nuisance {nuisance_index}/{len(other_sextupoles)}",
                flush=True,
            )
    nuisance_response = np.column_stack(nuisance_columns)

    candidate_plus: list[np.ndarray] = []
    candidate_minus: list[np.ndarray] = []
    delta_k1: list[float] = []
    for candidate_index, candidate_name in enumerate(candidates, start=1):
        quadrupole = quadrupole_by_name[candidate_name]
        step = abs(float(quadrupole["baseline"])) * quadrupole_fraction
        delta_k1.append(step)
        conditioned: list[np.ndarray] = []
        for sign in (1.0, -1.0):
            set_element(
                tao,
                candidate_name,
                "K1",
                float(quadrupole["baseline"]) + sign * step,
            )
            conditioned.append(
                np.column_stack(
                    [
                        mixed_k2_offset_response(
                            tao,
                            detector_data,
                            target,
                            target["name"],
                            attribute,
                            k2_step,
                            offset_step,
                            phase_mask,
                        )
                        for attribute in ("X_OFFSET", "Y_OFFSET")
                    ]
                )
            )
        set_element(tao, candidate_name, "K1", float(quadrupole["baseline"]))
        candidate_plus.append(conditioned[0])
        candidate_minus.append(conditioned[1])
        print(
            f"  {target['name']} candidate {candidate_index}/{len(candidates)}: {candidate_name}",
            flush=True,
        )

    elapsed = time.perf_counter() - start
    return {
        "observation_labels": np.asarray(labels),
        "target_response_nominal": target_response,
        "nuisance_labels": np.asarray(nuisance_labels),
        "nuisance_response_nominal": nuisance_response,
        "candidate_names": np.asarray(candidates),
        "candidate_delta_k1_m2": np.asarray(delta_k1, dtype=float),
        "target_response_candidate_plus": np.stack(candidate_plus),
        "target_response_candidate_minus": np.stack(candidate_minus),
        "calculation_seconds": np.asarray([elapsed], dtype=float),
    }


def main() -> int:
    args = parser().parse_args()
    validate_args(args)
    inputs = args.inputs.expanduser().resolve()
    lattice = args.lattice.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    _, control_names, _ = read_samples(inputs)
    detectors, horizontal, vertical = response_tools.parse_cesr_layout(lattice)
    if control_names != horizontal + vertical:
        raise RuntimeError("Input control order differs from the Bmad lattice order")
    references = variable_references(horizontal, vertical)
    init_path = output_dir / "tao_quadrupole_affinity_rf_on.init"
    init_path.write_text(
        response_tools.build_tao_init(lattice, detectors, horizontal, vertical, 1.0e-6),
        encoding="utf-8",
    )

    from pytao import Tao

    initialization_start = time.perf_counter()
    tao = Tao(init_file=str(init_path), noplot=True)
    tao.cmds(
        ["set global rf_on = T", "set particle_start pz = 0"],
        suppress_lattice_calc=True,
        suppress_plotting=True,
    )
    activate_benchmark_data(tao)
    detector_data = detector_inventory(tao, detectors)
    apply_sample(tao, references, [0.0] * len(references))
    collect(tao, detector_data)
    magnets = strength_inventory(tao, lattice)
    quadrupoles = [item for item in magnets if item["kind"] == "quadrupole"]
    sextupoles = [item for item in magnets if item["kind"] == "sextupole"]
    if (len(quadrupoles), len(sextupoles)) != (106, 76):
        raise RuntimeError(
            f"Expected 106 quadrupoles and 76 sextupoles, got {len(quadrupoles)} and {len(sextupoles)}"
        )
    labels, phase_mask = measurement_layout(detector_data["names"])
    initialization_seconds = time.perf_counter() - initialization_start

    screen_path = output_dir / "selected_candidates.json"
    if args.overwrite or not screen_path.exists():
        _, candidates = screen_quadrupoles(
            tao,
            detector_data,
            quadrupoles,
            sextupoles,
            args.quadrupole_fraction,
            args.top_k,
            args.max_tune_shift,
            args.max_beta_beating,
            output_dir,
        )
    else:
        candidates = json.loads(screen_path.read_text(encoding="utf-8"))
    if args.screen_only:
        return 0

    quadrupole_by_name = {item["name"]: item for item in quadrupoles}
    sextupole_by_name = {item["name"]: item for item in sextupoles}
    requested_targets = target_names(args.targets, sextupoles)
    response_dir = output_dir / "targets"
    response_dir.mkdir(parents=True, exist_ok=True)
    completed: list[dict[str, Any]] = []
    for target_index, name in enumerate(requested_targets, start=1):
        path = response_dir / f"{name.lower()}_responses.npz"
        if path.exists() and not args.overwrite:
            with np.load(path) as saved:
                seconds = float(saved["calculation_seconds"][0])
            completed.append({"target": name, "calculation_seconds": seconds, "reused": 1})
            print(f"Target {target_index}/{len(requested_targets)} {name}: reused", flush=True)
            continue
        print(f"Target {target_index}/{len(requested_targets)} {name}: starting", flush=True)
        arrays = calculate_target_responses(
            tao,
            detector_data,
            sextupole_by_name[name],
            sextupoles,
            quadrupole_by_name,
            list(candidates[name]),
            args.k2_step_m3,
            args.offset_step_m,
            args.quadrupole_fraction,
            labels,
            phase_mask,
        )
        temporary_path = path.with_name(f".{path.name}.tmp.npz")
        np.savez_compressed(temporary_path, **arrays)
        temporary_path.replace(path)
        seconds = float(arrays["calculation_seconds"][0])
        completed.append({"target": name, "calculation_seconds": seconds, "reused": 0})
        print(
            f"Target {target_index}/{len(requested_targets)} {name}: {seconds:.2f} s",
            flush=True,
        )

    write_csv(output_dir / "target_response_timings.csv", completed)
    metadata = {
        "format": "cesr-sextupole-quadrupole-affinity-responses-v1",
        "engine": "Bmad/Tao/PyTao",
        "method": "RF-on symmetric finite differences of target K2 slopes",
        "rf_mode": "on (6D periodic optics)",
        "sextupole_count": len(sextupoles),
        "quadrupole_count": len(quadrupoles),
        "candidate_count_per_target": args.top_k,
        "nuisance_columns_per_target": 2 * (len(sextupoles) - 1),
        "quadrupole_fraction": args.quadrupole_fraction,
        "k2_step_m3": args.k2_step_m3,
        "offset_step_m": args.offset_step_m,
        "max_tune_shift": args.max_tune_shift,
        "max_beta_beating": args.max_beta_beating,
        "detector_observables": list(DETECTOR_COLUMNS),
        "ring_observables": list(RING_COLUMNS),
        "observation_count": len(labels),
        "initialization_seconds": initialization_seconds,
        "targets_requested": requested_targets,
        "python_version": platform.python_version(),
        "tao_version_raw": tao.cmd("show version", raises=False),
        "lattice": str(lattice),
        "limitations": [
            "The 15-candidate pre-screen uses finite beta/phase leverage at the target sextupole.",
            "Other-sextupole K2-slope nuisance columns are evaluated at nominal quadrupole optics and reused for candidate +/- conditions.",
            "This is a nominal-orbit response screen; the planned 3x3 bump grid is a later exact validation stage.",
        ],
    }
    (output_dir / "response_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

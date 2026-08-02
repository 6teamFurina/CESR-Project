#!/usr/bin/env python3
"""Benchmark RF-off Bmad chromatic optics for the shared CESR samples.

One persistent Tao instance is used.  For each corrector sample, Bmad computes
periodic coasting optics at pz = 0 and at symmetric offsets pz = +/- delta.
The central differences provide the same d/ddelta columns emitted by the
SciBmad Descriptor(6, 2) benchmark.  Bmad's native ring chromaticities are
also saved and checked against the tune finite differences.
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
from typing import Any, Sequence, Union

OPTICS_DIR = Path(__file__).resolve().parent
DATASET_DIR = OPTICS_DIR.parent
ORBIT_DIR = DATASET_DIR / "orbit"
sys.path.insert(0, str(ORBIT_DIR))

from benchmark_bmad import (  # noqa: E402
    activate_benchmark_data,
    apply_sample,
    read_samples,
    response_tools,
    variable_references,
)


def make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--inputs",
        type=Path,
        default=ORBIT_DIR / "inputs" / "cesr_corrector_samples_1000.csv",
    )
    parser.add_argument(
        "--lattice",
        type=Path,
        default=ORBIT_DIR / "reference" / "cesr_bmad_compatible.bmad",
    )
    parser.add_argument("--sample-count", type=int, default=10)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=OPTICS_DIR / "results" / "chromatic_test_10",
    )
    parser.add_argument("--warmup-samples", type=int, default=1)
    parser.add_argument(
        "--delta-step",
        type=float,
        default=1.0e-5,
        help="Symmetric relative-momentum step used for d/ddelta",
    )
    return parser


def finite_float(value: Any) -> float:
    result = float(value)
    if not math.isfinite(result):
        raise RuntimeError(f"Bmad returned a non-finite value: {value!r}")
    return result


def lat_values(
    tao: Any,
    elements: str,
    candidates: Union[str, Sequence[str]],
    expected: int,
) -> tuple[list[Any], str]:
    if isinstance(candidates, str):
        candidates = (candidates,)
    errors: list[str] = []
    for who in candidates:
        try:
            values = list(
                tao.lat_list(
                    elements,
                    who,
                    which="model",
                    flags="-array_out -track_only",
                )
            )
            if len(values) != expected:
                raise RuntimeError(
                    f"returned {len(values)} values instead of {expected}"
                )
            return values, who
        except Exception as exc:  # PyTao exception classes vary by release.
            errors.append(f"{who}: {exc}")
    raise RuntimeError(
        f"Tao query failed for {elements!r}; tried {list(candidates)!r}: "
        + " | ".join(errors)
    )


def numeric_lat_values(
    tao: Any,
    elements: str,
    candidates: Union[str, Sequence[str]],
    expected: int,
) -> tuple[list[float], str]:
    values, selected = lat_values(tao, elements, candidates, expected)
    return [finite_float(value) for value in values], selected


def dictionary_float(mapping: dict[str, Any], candidates: Sequence[str]) -> float:
    lower = {str(key).lower(): value for key, value in mapping.items()}
    for candidate in candidates:
        if candidate.lower() in lower:
            return finite_float(lower[candidate.lower()])
    raise RuntimeError(
        f"None of {tuple(candidates)!r} is present in ring_general keys "
        f"{tuple(mapping)!r}"
    )


def detector_inventory(tao: Any, expected_names: Sequence[str]) -> dict[str, Any]:
    names_raw, name_field = lat_values(tao, "det_*", "ele.name", 99)
    names = [str(value).upper() for value in names_raw]
    if names != list(expected_names):
        raise RuntimeError(
            "Tao DET_* order differs from the shared detector order\n"
            f"Tao: {names}\nExpected: {list(expected_names)}"
        )
    indices, index_field = numeric_lat_values(tao, "det_*", "ele.ix_ele", 99)
    positions, s_field = numeric_lat_values(tao, "det_*", "ele.s", 99)
    all_indices, _ = numeric_lat_values(tao, "*", "ele.ix_ele", 870)
    return {
        "names": names,
        "indices": [int(value) for value in indices],
        "s": positions,
        "first_index": int(all_indices[0]),
        "last_index": int(all_indices[-1]),
        "selected_fields": {
            "name": name_field,
            "beamline_index": index_field,
            "s": s_field,
        },
    }


STATE_FIELDS: tuple[tuple[str, tuple[str, ...]], ...] = (
    ("beta_1", ("ele.a.beta", "twiss.beta_a")),
    ("alpha_1", ("ele.a.alpha", "twiss.alpha_a")),
    ("phi_1_radians", ("ele.a.phi", "twiss.phi_a")),
    ("beta_2", ("ele.b.beta", "twiss.beta_b")),
    ("alpha_2", ("ele.b.alpha", "twiss.alpha_b")),
    ("phi_2_radians", ("ele.b.phi", "twiss.phi_b")),
    ("orbit_x", ("orbit.vec.1",)),
    ("orbit_px", ("orbit.vec.2",)),
    ("orbit_y", ("orbit.vec.3",)),
    ("orbit_py", ("orbit.vec.4",)),
    ("orbit_z", ("orbit.vec.5",)),
    ("orbit_pz", ("orbit.vec.6",)),
)

# Bmad stores Twiss chromatic derivatives directly in its element structures,
# but Tao's ``lat_list`` interface does not document mode-qualified
# ``ele.a/b.dbeta_dpz`` or ``ele.a/b.dalpha_dpz`` selectors.  In Tao
# 20260726-0 those selectors silently returned the same mode array for both
# requests.  beta/alpha derivatives therefore remain the symmetric differences
# already computed from the pz=+/-step states below.  The native dispersion
# selectors are documented ``lat_list`` fields and are safe to use.
NATIVE_DERIVATIVE_FIELDS: tuple[tuple[str, tuple[str, ...]], ...] = (
    ("orbit_x", ("ele.x.eta", "twiss.eta_x")),
    ("orbit_px", ("ele.x.etap", "twiss.etap_x")),
    ("orbit_y", ("ele.y.eta", "twiss.eta_y")),
    ("orbit_py", ("ele.y.etap", "twiss.etap_y")),
    ("orbit_z", ("ele.z.eta", "twiss.eta_z")),
    ("orbit_pz", ("ele.z.etap", "twiss.etap_z")),
)

TWISS_COLUMNS = (
    "phi_1",
    "beta_1",
    "alpha_1",
    "phi_2",
    "beta_2",
    "alpha_2",
    "phi_3",
    "gamma_c",
    "c11",
    "c12",
    "c21",
    "c22",
)

ORBIT_COLUMNS = (
    "orbit_x",
    "orbit_px",
    "orbit_y",
    "orbit_py",
    "orbit_z",
    "orbit_pz",
)


def set_momentum_offset(tao: Any, pz: float) -> None:
    tao.cmds(
        [f"set particle_start pz = {pz:.17g}"],
        suppress_lattice_calc=True,
        suppress_plotting=True,
    )


def coupling_gamma(columns: dict[str, list[float]]) -> list[float]:
    result: list[float] = []
    for c11, c12, c21, c22 in zip(
        columns["c11"], columns["c12"], columns["c21"], columns["c22"]
    ):
        determinant = c11 * c22 - c12 * c21
        radicand = 1.0 - determinant
        if radicand < -1.0e-12:
            raise RuntimeError(f"Invalid Bmad coupling determinant: {determinant}")
        result.append(math.sqrt(max(0.0, radicand)))
    return result


def collect_state(
    tao: Any,
    inventory: dict[str, Any],
    pz: float,
) -> tuple[dict[str, Any], float, float]:
    """Set one energy offset and collect a complete periodic-optics state."""
    set_momentum_offset(tao, pz)
    selected_fields: dict[str, str] = {}

    # The first query triggers the one Bmad lattice recalculation at this pz.
    recalc_start = time.perf_counter()
    orbit_x, selected_fields["orbit_x"] = numeric_lat_values(
        tao, "det_*", "orbit.vec.1", 99
    )
    recalc_seconds = time.perf_counter() - recalc_start

    read_start = time.perf_counter()
    columns: dict[str, list[float]] = {"orbit_x": orbit_x}
    for output_name, candidates in STATE_FIELDS:
        if output_name == "orbit_x":
            continue
        columns[output_name], selected_fields[output_name] = numeric_lat_values(
            tao, "det_*", candidates, 99
        )

    c_values, c_field = numeric_lat_values(tao, "det_*", "ele.c_mat", 99 * 4)
    columns["c11"] = c_values[0::4]
    columns["c12"] = c_values[1::4]
    columns["c21"] = c_values[2::4]
    columns["c22"] = c_values[3::4]
    selected_fields["c11_c12_c21_c22"] = c_field
    for mode in (1, 2):
        columns[f"phi_{mode}"] = [
            value / (2.0 * math.pi)
            for value in columns.pop(f"phi_{mode}_radians")
        ]
    # A coasting closed ring has no periodic longitudinal c mode. SciBmad's
    # coasting phi_3 has zero constant term; its d/ddelta coefficient is the
    # accumulated longitudinal slip, constructed below from orbit_z.
    columns["phi_3"] = [0.0] * 99
    columns["gamma_c"] = coupling_gamma(columns)

    first_selector = str(inventory["first_index"])
    last_selector = str(inventory["last_index"])
    start_orbit: list[float] = []
    end_orbit: list[float] = []
    for coordinate in range(1, 7):
        start_value, _ = numeric_lat_values(
            tao, first_selector, f"orbit.vec.{coordinate}", 1
        )
        end_value, _ = numeric_lat_values(
            tao, last_selector, f"orbit.vec.{coordinate}", 1
        )
        start_orbit.append(start_value[0])
        end_orbit.append(end_value[0])
    transverse_closure_norm = math.sqrt(
        sum((end_orbit[index] - start_orbit[index]) ** 2 for index in range(4))
    )

    ring = tao.ring_general()
    q1 = dictionary_float(ring, ("Q_a", "q_a"))
    q2 = dictionary_float(ring, ("Q_b", "q_b"))
    chrom1 = dictionary_float(ring, ("chrom_a",))
    chrom2 = dictionary_float(ring, ("chrom_b",))
    longitudinal_advance = [
        orbit_z - start_orbit[4] for orbit_z in columns["orbit_z"]
    ]
    end_longitudinal_advance = end_orbit[4] - start_orbit[4]
    read_seconds = time.perf_counter() - read_start

    return (
        {
            "columns": columns,
            "start_orbit": start_orbit,
            "end_orbit": end_orbit,
            "transverse_closure_norm": transverse_closure_norm,
            "Q1": q1,
            "Q2": q2,
            "chrom1_native": chrom1,
            "chrom2_native": chrom2,
            "longitudinal_advance": longitudinal_advance,
            "end_longitudinal_advance": end_longitudinal_advance,
            "selected_fields": selected_fields,
        },
        recalc_seconds,
        read_seconds,
    )


def collect_native_derivatives(
    tao: Any,
) -> tuple[dict[str, list[float]], dict[str, str], dict[str, str]]:
    """Read native Bmad chromatic fields, falling back field-by-field."""
    derivatives: dict[str, list[float]] = {}
    selected: dict[str, str] = {}
    errors: dict[str, str] = {}
    for column, candidates in NATIVE_DERIVATIVE_FIELDS:
        try:
            derivatives[column], selected[column] = numeric_lat_values(
                tao, "det_*", candidates, 99
            )
        except Exception as exc:
            errors[column] = str(exc)
    return derivatives, selected, errors


def phase_central_difference(plus: float, minus: float, step: float) -> float:
    # Phase columns are in turns. Reduce only the small +/- difference to the
    # nearest turn so a 0/1 boundary does not create a spurious derivative.
    difference = (plus - minus + 0.5) % 1.0 - 0.5
    return difference / (2.0 * step)


def collect_chromatic_sample(
    tao: Any,
    inventory: dict[str, Any],
    delta_step: float,
) -> tuple[dict[str, Any], float, float]:
    base, base_recalc, base_read = collect_state(tao, inventory, 0.0)
    native_derivatives, native_fields, native_errors = collect_native_derivatives(tao)
    plus, plus_recalc, plus_read = collect_state(tao, inventory, delta_step)
    minus, minus_recalc, minus_read = collect_state(tao, inventory, -delta_step)

    derivatives: dict[str, list[float]] = {}
    for column in (*TWISS_COLUMNS, *ORBIT_COLUMNS):
        if column == "phi_3":
            plus_values = plus["longitudinal_advance"]
            minus_values = minus["longitudinal_advance"]
            derivatives[column] = [
                (high - low) / (2.0 * delta_step)
                for high, low in zip(plus_values, minus_values)
            ]
            continue
        plus_values = plus["columns"][column]
        minus_values = minus["columns"][column]
        if column in ("phi_1", "phi_2"):
            derivatives[column] = [
                phase_central_difference(high, low, delta_step)
                for high, low in zip(plus_values, minus_values)
            ]
        else:
            derivatives[column] = [
                (high - low) / (2.0 * delta_step)
                for high, low in zip(plus_values, minus_values)
            ]

    derivative_sources = {column: "symmetric finite difference" for column in derivatives}
    for column, values in native_derivatives.items():
        derivatives[column] = values
        derivative_sources[column] = f"native Bmad field {native_fields[column]}"

    q1_fd = phase_central_difference(plus["Q1"], minus["Q1"], delta_step)
    q2_fd = phase_central_difference(plus["Q2"], minus["Q2"], delta_step)
    # SciBmad reports the ring coasting slip with the opposite sign to the
    # accumulated local longitudinal phase advance.
    slip_factor = -(
        plus["end_longitudinal_advance"] - minus["end_longitudinal_advance"]
    ) / (2.0 * delta_step)
    return (
        {
            "columns": base["columns"],
            "derivatives": derivatives,
            "start_orbit": base["start_orbit"],
            "transverse_closure_norm": base["transverse_closure_norm"],
            "Q1": base["Q1"],
            "Q2": base["Q2"],
            "xi_1": base["chrom1_native"],
            "xi_2": base["chrom2_native"],
            "xi_1_finite_difference": q1_fd,
            "xi_2_finite_difference": q2_fd,
            "slip_tps_constant": 0.0,
            "slip_factor": slip_factor,
            "selected_fields": base["selected_fields"],
            "derivative_sources": derivative_sources,
            "native_derivative_query_errors": native_errors,
        },
        base_recalc + plus_recalc + minus_recalc,
        base_read + plus_read + minus_read,
    )


def detector_output_fields() -> tuple[str, ...]:
    fields = ["sample_id", "s", "beamline_index", "name"]
    for column in TWISS_COLUMNS:
        fields.extend((column, f"d{column}_ddelta"))
    for column in ORBIT_COLUMNS:
        fields.extend((column, f"d{column}_ddelta"))
    return tuple(fields)


def write_detector_output(
    path: Path,
    sample_ids: Sequence[int],
    inventory: dict[str, Any],
    results: Sequence[dict[str, Any]],
) -> None:
    fields = detector_output_fields()
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for sample_id, result in zip(sample_ids, results):
            for detector_index in range(99):
                row: dict[str, Any] = {
                    "sample_id": sample_id,
                    "s": inventory["s"][detector_index],
                    "beamline_index": inventory["indices"][detector_index],
                    "name": inventory["names"][detector_index].lower(),
                }
                for column in (*TWISS_COLUMNS, *ORBIT_COLUMNS):
                    row[column] = result["columns"][column][detector_index]
                    row[f"d{column}_ddelta"] = result["derivatives"][column][
                        detector_index
                    ]
                writer.writerow(row)


def write_ring_output(
    path: Path,
    sample_ids: Sequence[int],
    results: Sequence[dict[str, Any]],
    sample_seconds: Sequence[float],
) -> None:
    fields = (
        "sample_id",
        "Q1",
        "Q2",
        "Qx_fractional",
        "Qy_fractional",
        "slip_tps_constant",
        "xi_1",
        "xi_2",
        "slip_factor",
        "bmad_physics_seconds",
        "xi_1_finite_difference",
        "xi_2_finite_difference",
    )
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for sample_id, result, elapsed in zip(sample_ids, results, sample_seconds):
            writer.writerow(
                {
                    "sample_id": sample_id,
                    "Q1": result["Q1"],
                    "Q2": result["Q2"],
                    "Qx_fractional": result["Q1"] % 1.0,
                    "Qy_fractional": result["Q2"] % 1.0,
                    "slip_tps_constant": result["slip_tps_constant"],
                    "xi_1": result["xi_1"],
                    "xi_2": result["xi_2"],
                    "slip_factor": result["slip_factor"],
                    "bmad_physics_seconds": elapsed,
                    "xi_1_finite_difference": result["xi_1_finite_difference"],
                    "xi_2_finite_difference": result["xi_2_finite_difference"],
                }
            )


def write_start_orbits(
    path: Path,
    sample_ids: Sequence[int],
    results: Sequence[dict[str, Any]],
) -> None:
    fields = ("sample_id", "x", "px", "y", "py", "z", "pz", "closure_norm")
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(fields)
        for sample_id, result in zip(sample_ids, results):
            writer.writerow(
                [
                    sample_id,
                    *result["start_orbit"],
                    result["transverse_closure_norm"],
                ]
            )


def main() -> int:
    args = make_parser().parse_args()
    inputs = args.inputs.expanduser().resolve()
    lattice = args.lattice.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    if not math.isfinite(args.delta_step) or args.delta_step <= 0.0:
        raise ValueError("--delta-step must be a positive finite number")

    sample_ids, control_names, all_samples = read_samples(inputs)
    if not 1 <= args.sample_count <= len(all_samples):
        raise ValueError(f"--sample-count must be between 1 and {len(all_samples)}")
    if args.warmup_samples < 1:
        raise ValueError("--warmup-samples must be positive")
    sample_ids = sample_ids[: args.sample_count]
    samples = all_samples[: args.sample_count]

    expected_detectors, horizontal, vertical = response_tools.parse_cesr_layout(lattice)
    if control_names != horizontal + vertical:
        raise RuntimeError("Input control order differs from the Bmad lattice order")
    references = variable_references(horizontal, vertical)

    output_dir.mkdir(parents=True, exist_ok=True)
    init_path = output_dir / "tao_bmad_chromatic_optics_rf_off.init"
    init_path.write_text(
        response_tools.build_tao_init(
            lattice,
            expected_detectors,
            horizontal,
            vertical,
            derivative_step=args.delta_step,
        ),
        encoding="utf-8",
    )

    try:
        from pytao import Tao
    except ImportError as exc:
        raise RuntimeError(
            "Activate the Linux environment containing PyTao before running this script"
        ) from exc

    initialization_start = time.perf_counter()
    tao = Tao(init_file=str(init_path), noplot=True)
    tao.cmds(
        [
            "set global rf_on = F",
            f"set global delta_e_chrom = {args.delta_step:.17g}",
            "set particle_start pz = 0",
        ],
        suppress_lattice_calc=True,
        suppress_plotting=True,
    )
    activate_benchmark_data(tao)
    inventory = detector_inventory(tao, expected_detectors)
    initialization_seconds = time.perf_counter() - initialization_start

    warmup_count = min(args.warmup_samples, len(samples))
    warmup_start = time.perf_counter()
    for sample in samples[:warmup_count]:
        apply_sample(tao, references, sample)
        collect_chromatic_sample(tao, inventory, args.delta_step)
    warmup_seconds = time.perf_counter() - warmup_start

    results: list[dict[str, Any]] = []
    sample_seconds: list[float] = []
    variable_update_seconds = 0.0
    lattice_recalc_seconds = 0.0
    output_query_seconds = 0.0
    physics_start = time.perf_counter()
    for row, (sample_id, sample) in enumerate(zip(sample_ids, samples), start=1):
        sample_start = time.perf_counter()
        update_start = time.perf_counter()
        apply_sample(tao, references, sample)
        variable_update_seconds += time.perf_counter() - update_start
        result, recalc_seconds, read_seconds = collect_chromatic_sample(
            tao, inventory, args.delta_step
        )
        lattice_recalc_seconds += recalc_seconds
        output_query_seconds += read_seconds
        elapsed = time.perf_counter() - sample_start
        results.append(result)
        sample_seconds.append(elapsed)
        print(
            f"Bmad chromatic optics {row}/{len(samples)} "
            f"(sample_id={sample_id}): {elapsed:.3f} s"
        )
    physics_seconds = time.perf_counter() - physics_start

    detector_path = output_dir / "bmad_detector_chromatic_twiss.csv"
    ring_path = output_dir / "bmad_ring_chromatic_twiss.csv"
    orbit_path = output_dir / "bmad_start_closed_orbits.csv"
    metadata_path = output_dir / "bmad_chromatic_optics_metadata.json"
    write_start = time.perf_counter()
    write_detector_output(detector_path, sample_ids, inventory, results)
    write_ring_output(ring_path, sample_ids, results, sample_seconds)
    write_start_orbits(orbit_path, sample_ids, results)
    write_seconds = time.perf_counter() - write_start

    xi1_disagreement = [
        abs(result["xi_1"] - result["xi_1_finite_difference"])
        for result in results
    ]
    xi2_disagreement = [
        abs(result["xi_2"] - result["xi_2_finite_difference"])
        for result in results
    ]
    metadata = {
        "format": "cesr-chromatic-optics-benchmark-v1",
        "engine": "Bmad/Tao/PyTao",
        "method": (
            "persistent RF-off Tao model; periodic optics at pz=0,+step,-step; "
            "symmetric finite differences for local d/ddelta; native Bmad "
            "ring chromaticities"
        ),
        "rf_mode": "off (coasting)",
        "delta_step": args.delta_step,
        "derivative_method": "second-order accurate symmetric finite difference",
        "ring_chromaticity_source": "Bmad ring_general chrom_a/chrom_b",
        "local_derivative_sources": results[0]["derivative_sources"],
        "native_derivative_query_errors": results[0]["native_derivative_query_errors"],
        "input_csv": str(inputs),
        "lattice": str(lattice),
        "output_directory": str(output_dir),
        "sample_count": len(samples),
        "control_count": len(control_names),
        "detector_count": len(inventory["names"]),
        "detector_row_count": len(samples) * len(inventory["names"]),
        "initialization_seconds": initialization_seconds,
        "warmup_sample_count": warmup_count,
        "warmup_seconds": warmup_seconds,
        "physics_seconds": physics_seconds,
        "variable_update_seconds": variable_update_seconds,
        "lattice_recalculation_seconds": lattice_recalc_seconds,
        "output_query_seconds": output_query_seconds,
        "seconds_per_sample": sample_seconds,
        "seconds_per_sample_mean": sum(sample_seconds) / len(sample_seconds),
        "samples_per_second": len(samples) / physics_seconds,
        "write_seconds": write_seconds,
        "valid_count": sum(
            math.isfinite(result["transverse_closure_norm"])
            for result in results
        ),
        "maximum_transverse_closure_norm": max(
            result["transverse_closure_norm"] for result in results
        ),
        "maximum_abs_xi1_native_minus_fd": max(xi1_disagreement),
        "maximum_abs_xi2_native_minus_fd": max(xi2_disagreement),
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "tao_version_raw": tao.cmd("show version", raises=False),
        "selected_tao_fields": results[0]["selected_fields"],
        "detector_location": (
            "exit of zero-length DET_* marker; physically identical to marker beginning"
        ),
        "persistent_tao_instance": True,
        "model_load_count": 1,
        "lattice_recalculations_per_sample": 3,
        "timed_region": (
            "119-variable update + Bmad optics at pz=0,+step,-step + all output queries"
        ),
    }
    metadata_path.write_text(
        json.dumps(metadata, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )

    print(
        f"Bmad chromatic-optics physics: {physics_seconds:.3f} s, "
        f"{len(samples) / physics_seconds:.3f} samples/s"
    )
    print(f"Detector chromatic optics: {detector_path}")
    print(f"Ring chromatic optics:     {ring_path}")
    print(f"Start orbits:              {orbit_path}")
    print(f"Metadata:                  {metadata_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

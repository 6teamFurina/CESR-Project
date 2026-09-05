#!/usr/bin/env python3
"""Time maintained kernels on saved latest-lattice SciBmad observations.

No physical-machine acquisition rate is assumed. This measures one machine,
one noise realization, and all 76 targets, with reusable operators separated
from online processing. Production analysis and result files are unchanged.
"""

import json
import platform
import sys
import time
import tomllib
from pathlib import Path

import numpy as np
import scipy

import analyze_state_space_bpm_gtpsa_inverse as study


OUT = Path(__file__).resolve().parent / "results" / "runtime_breakdown_20260905"


def timed(function):
    start = time.perf_counter()
    value = function()
    return value, time.perf_counter() - start


def repeated(function, count=11):
    function()
    samples = [timed(function)[1] for _ in range(count)]
    return {
        "median_seconds": float(np.median(samples)),
        "min_seconds": min(samples),
        "max_seconds": max(samples),
        "samples_seconds": samples,
    }


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    source = study.SCAN_ROOT / study.DEFAULT_CASE
    metadata = tomllib.loads((source / "scan_metadata.toml").read_text())
    targets = (source / "target_names.txt").read_text().splitlines()
    bpms = (source / "bpm_names.txt").read_text().splitlines()
    commands = np.array([
        [float(row["bump_x_command_m"]), float(row["bump_y_command_m"])]
        for row in study.read_rows(source / "bump_points.csv")
    ])
    dk = np.array(metadata["k2_delta_m3"])
    zb = int(np.flatnonzero(np.all(commands == 0, axis=1))[0])
    zk = int(np.flatnonzero(dk == 0)[0])
    setup = {}
    model, setup["load_and_build_local_orbit_model_seconds"] = timed(
        lambda: study.subset_orbit_model(
            study.DEFAULT_MODEL, study.DEFAULT_KNOBS, targets, bpms, commands
        )
    )
    templates, setup["source_templates_from_saved_gtpsa_maps_seconds"] = timed(
        lambda: study.derivative.source_templates(study.DEFAULT_MODEL, 0.272)
    )
    protocol, setup["build_protocol_seconds"] = timed(
        lambda: study.build_protocol(commands, dk, 3072, 256)
    )
    variance = 1e-10 / (3072 * 8 - 1)
    brownian, setup["simulation_only_brownian_covariance_seconds"] = timed(
        lambda: study.covariance_sqrt(
            study.brownian_functional_covariance(protocol, variance)
        )
    )
    static = np.load(source / "observable_bpm_readbacks.npy", mmap_mode="r")[:1]
    drift = np.load(source / "observable_drift_bpm_readbacks.npy", mmap_mode="r")[:1]
    response = study.recover_forward_drift_response(
        static, drift, metadata["drift_halfwidth_m"]
    )
    static = static.reshape(1, len(targets), 5, 3, -1)
    operators = []
    observations = []
    operator_times = []
    simulation_times = []
    for target in range(len(targets)):
        operator, elapsed = timed(lambda: study.build_state_space_operator(
            study.nominal_drift_matrix(model, target, commands),
            protocol, 5e-6, 32, variance,
        ))
        operators.append(operator)
        operator_times.append(elapsed)
        observed, elapsed = timed(lambda: study.simulate_forward_target_observables(
            static[:, target], response[:, target], operator, protocol, zk,
            brownian, 5e-6, 32, 1, 20261230, target,
        ))
        observations.append(observed)
        simulation_times.append(elapsed)
    setup["state_space_operators_all_targets_seconds"] = sum(operator_times)
    setup["simulation_only_observable_generation_seconds"] = sum(simulation_times)

    def local_orbits():
        return [study.reconstruct_target_local_orbits(
            observations[t].calibration_readbacks, t, zb, model
        ) for t in range(len(targets))]

    local = local_orbits()

    def filter_observations():
        result = []
        for t, op in enumerate(operators):
            posterior = study.hidden_state_filtered_averages(
                observations[t].projected_reference_observations, op
            )
            result.append(observations[t].unfiltered_core_readbacks
                          - posterior @ op.nominal_bpm_drift.T)
        return result

    filtered = filter_observations()
    # Exercise the full 156 x 222 reference projection, which production's
    # sufficient-statistic simulator supplies already projected. These inputs
    # preserve its local reference observations in the nominal response space.
    reference_residuals = [
        observations[t].projected_reference_observations
        @ operators[t].nominal_bpm_drift.T for t in range(len(targets))
    ]

    def project_references():
        return [reference_residuals[t] @ operators[t].bpm_to_local_drift.T
                for t in range(len(targets))]

    def inverse():
        results = []
        for t in range(len(targets)):
            slopes, bumps = study.k2_slopes_from_core(filtered[t], protocol, dk)
            centers = study.fixed_template_centers(slopes, local[t][0][:, bumps],
                                                   templates[t])
            results.append(centers + local[t][1])
        return results

    kernels = {
        "bpm_to_sextupole_local_and_reference_orbit": repeated(local_orbits),
        "reference_projection": repeated(project_references),
        "state_smoother_and_bpm_correction": repeated(filter_observations),
        "slopes_fixed_template_inverse_and_absolute_offset": repeated(inverse),
    }
    np.save(OUT / "one_machine_fixed_template_offsets.npy", np.array(inverse()))

    # Time the current complete comparative analysis at the same one-machine,
    # one-noise-realization scale; this includes static multi-start fitting and
    # both filtered/unfiltered profiled and fixed-template alternatives.
    def comparison():
        return [study.run_machine_facing_target(
            static[:, t], observations[t], operators[t], t, model, dk, zb, zk,
            protocol, 5e-6, 2.5e-3, templates[t]
        ) for t in range(len(targets))]

    compared, comparative_seconds = timed(comparison)
    expected = np.array([x.filtered_fixed_template_absolute_offsets[0]
                         for x in compared])
    actual = np.array(inverse())
    np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-16)
    output = {
        "date": "2026-09-05",
        "python": sys.executable,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "source": str(source),
        "lattice": metadata["lattice"],
        "machines": 1,
        "targets": len(targets),
        "stochastic_realizations": 1,
        "online_timing_semantics": "warm in-memory single-machine 76-target passes; excludes hardware, raw streaming accumulation, file IO and model/operator setup",
        "setup": setup,
        "kernels": kernels,
        "all_comparison_estimators_seconds": comparative_seconds,
        "fixed_path_vs_production_max_abs_m": float(np.max(np.abs(actual - expected))),
        "acquisition": {
            "signal_reads_per_target": 3072 * 8,
            "reference_reads_per_target": len(protocol.reference_times),
            "calibration_reads_per_target": len(commands) * 32,
            "total_reads_per_target": protocol.total_acquisitions + len(commands) * 32,
            "total_reads_all_targets": (protocol.total_acquisitions + len(commands) * 32) * len(targets),
            "full_bpm_vectors_per_acquisition": 1,
            "real_hardware_effective_rate_hz": None,
            "settling_time_seconds": None,
        },
        "limitations": [
            "One latent machine and one stochastic realization; timings are workstation-specific.",
            "Reference-projection timing uses nominal-space residuals reconstructed from saved simulation statistics.",
            "Latest-lattice straight multipoles in a curved reference remain an approximation; no girder pitch is varied.",
        ],
    }
    (OUT / "python_timing.json").write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()

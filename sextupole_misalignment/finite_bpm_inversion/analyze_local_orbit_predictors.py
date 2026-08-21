#!/usr/bin/env python3
"""Compare three leakage-safe predictors of relative sextupole local orbit.

Predictions use only nominal latest-lattice model artifacts, known corrector
commands, and nominal-K2 BPM orbit differences. Exact target-local orbits are
loaded only after all predictions and hyperparameter selection are complete.
"""

from __future__ import annotations

import argparse
import csv
import tomllib
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
DEFAULT_SOURCE = (
    HERE.parent
    / "direct_observable_nuisance_ablation"
    / "results"
    / "all_76_orbit_protocol"
)
DEFAULT_MODEL = HERE / "results" / "local_orbit_model"
DEFAULT_KNOBS = (
    HERE.parent
    / "quadrupole_affinity"
    / "exact_11_triplet_validation"
    / "results"
    / "bump_knobs"
    / "local_bump_knobs.csv"
)


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write empty CSV: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def summarize(errors_m: np.ndarray) -> dict[str, float]:
    radial_um = np.linalg.norm(errors_m, axis=-1) * 1e6
    return {
        "x_rmse_um": float(np.sqrt(np.mean(errors_m[..., 0] ** 2)) * 1e6),
        "y_rmse_um": float(np.sqrt(np.mean(errors_m[..., 1] ** 2)) * 1e6),
        "rmse_2d_um": float(np.sqrt(np.mean(radial_um**2))),
        "median_2d_um": float(np.median(radial_um)),
        "p90_2d_um": float(np.percentile(radial_um, 90)),
        "p99_2d_um": float(np.percentile(radial_um, 99)),
        "max_2d_um": float(np.max(radial_um)),
    }


def ridge_operator(matrix: np.ndarray, ratio: float) -> np.ndarray:
    """Return the Tikhonov operator mapping observations to coefficients."""
    u, singular, vt = np.linalg.svd(matrix, full_matrices=False)
    if singular.size == 0 or singular[0] <= 0:
        raise ValueError("Response matrix has no positive singular value")
    if ratio == 0.0:
        cutoff = np.finfo(float).eps * max(matrix.shape) * singular[0]
        factors = np.divide(
            1.0,
            singular,
            out=np.zeros_like(singular),
            where=singular > cutoff,
        )
    else:
        ridge = ratio * singular[0]
        factors = singular / (singular * singular + ridge * ridge)
    return (vt.T * factors) @ u.T


def forward_map(
    cumulative_from: np.ndarray,
    cumulative_to: np.ndarray,
    wraps: bool,
    one_turn: np.ndarray,
) -> np.ndarray:
    middle = one_turn if wraps else np.eye(6)
    return cumulative_to @ middle @ np.linalg.inv(cumulative_from)


def build_two_sided_maps(
    bpm_maps: np.ndarray,
    target_maps: np.ndarray,
    one_turn: np.ndarray,
    bpm_line_indices: np.ndarray,
    target_line_indices: np.ndarray,
) -> tuple[np.ndarray, list[dict[str, object]]]:
    """Map upstream/downstream BPM position residuals to each target position."""
    target_count = len(target_line_indices)
    maps = np.zeros((target_count, 2, 4))
    rows: list[dict[str, object]] = []
    transverse = np.array([0, 1, 2, 3])
    position = np.array([0, 2])
    momentum = np.array([1, 3])
    for target in range(target_count):
        line_index = target_line_indices[target]
        before = np.flatnonzero(bpm_line_indices < line_index)
        after = np.flatnonzero(bpm_line_indices > line_index)
        upstream = int(before[np.argmax(bpm_line_indices[before])]) if before.size else int(np.argmax(bpm_line_indices))
        downstream = int(after[np.argmin(bpm_line_indices[after])]) if after.size else int(np.argmin(bpm_line_indices))

        map_up_down = forward_map(
            bpm_maps[upstream],
            bpm_maps[downstream],
            bool(bpm_line_indices[upstream] > bpm_line_indices[downstream]),
            one_turn,
        )[np.ix_(transverse, transverse)]
        map_up_target = forward_map(
            bpm_maps[upstream],
            target_maps[target],
            bool(bpm_line_indices[upstream] > line_index),
            one_turn,
        )[np.ix_(transverse, transverse)]

        down_from_position = map_up_down[np.ix_(position, position)]
        down_from_momentum = map_up_down[np.ix_(position, momentum)]
        target_from_position = map_up_target[np.ix_(position, position)]
        target_from_momentum = map_up_target[np.ix_(position, momentum)]
        momentum_inverse = np.linalg.pinv(down_from_momentum, rcond=1e-12)
        maps[target, :, :2] = (
            target_from_position
            - target_from_momentum @ momentum_inverse @ down_from_position
        )
        maps[target, :, 2:] = target_from_momentum @ momentum_inverse
        singular = np.linalg.svd(down_from_momentum, compute_uv=False)
        rows.append(
            {
                "target_index": target + 1,
                "upstream_bpm_index": upstream + 1,
                "downstream_bpm_index": downstream + 1,
                "momentum_block_condition": float(singular[0] / singular[-1]),
                "momentum_block_sigma_min": float(singular[-1]),
            }
        )
    return maps, rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--model-dir", type=Path, default=DEFAULT_MODEL)
    parser.add_argument("--knobs", type=Path, default=DEFAULT_KNOBS)
    parser.add_argument(
        "--ridge-ratios",
        default="0,1e-8,3e-8,1e-7,3e-7,1e-6,3e-6,1e-5,3e-5,1e-4,3e-4,1e-3,3e-3,1e-2,3e-2,1e-1,3e-1,1",
    )
    parser.add_argument("--cv-folds", type=int, default=5)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=HERE / "results" / "local_orbit_predictors",
    )
    args = parser.parse_args()
    source = args.input_dir.resolve()
    model_dir = args.model_dir.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)

    with (source / "scan_metadata.toml").open("rb") as stream:
        scan_metadata = tomllib.load(stream)
    bpm_names = (source / "bpm_names.txt").read_text(encoding="utf-8").splitlines()
    target_names = (source / "target_names.txt").read_text(encoding="utf-8").splitlines()
    bpm_rows = read_rows(model_dir / "bpm_locations.csv")
    target_rows = read_rows(model_dir / "target_locations.csv")
    control_rows = read_rows(model_dir / "control_inventory.csv")
    if bpm_names != [row["bpm"] for row in bpm_rows]:
        raise ValueError("Scan and model BPM inventories do not match")
    if target_names != [row["target"] for row in target_rows]:
        raise ValueError("Scan and model target inventories do not match")

    bpm_response = np.load(model_dir / "bpm_control_response.npy")
    target_response_flat = np.load(model_dir / "target_control_response.npy")
    bpm_maps = np.load(model_dir / "bpm_cumulative_maps.npy")
    target_maps = np.load(model_dir / "target_cumulative_maps.npy")
    one_turn = np.load(model_dir / "one_turn_map.npy")
    nt = len(target_names)
    nd = len(bpm_names)
    nc = len(control_rows)
    if bpm_response.shape != (2 * nd, nc):
        raise ValueError(f"Unexpected BPM response shape: {bpm_response.shape}")
    if target_response_flat.shape != (2 * nt, nc):
        raise ValueError(f"Unexpected target response shape: {target_response_flat.shape}")
    target_response = target_response_flat.reshape(nt, 2, nc)

    bump_rows = read_rows(source / "bump_points.csv")
    bump_commands = np.asarray(
        [
            (float(row["bump_x_command_m"]), float(row["bump_y_command_m"]))
            for row in bump_rows
        ]
    )
    zero_candidates = np.flatnonzero(np.all(bump_commands == 0.0, axis=1))
    if zero_candidates.size != 1:
        raise ValueError("Expected exactly one zero bump")
    zero_bump = int(zero_candidates[0])
    levels = np.asarray(scan_metadata["k2_levels"], dtype=float)
    nominal_candidates = np.flatnonzero(levels == 0.0)
    if nominal_candidates.size != 1:
        raise ValueError("Expected exactly one nominal K2 level")
    nominal_k2 = int(nominal_candidates[0])

    control_lookup = {
        (row["corrector"], row["field"]): index
        for index, row in enumerate(control_rows)
    }
    knob_x = np.zeros((nt, nc))
    knob_y = np.zeros((nt, nc))
    target_lookup = {name: index for index, name in enumerate(target_names)}
    for row in read_rows(args.knobs.resolve()):
        target = target_lookup[row["target_sextupole"]]
        control = control_lookup[(row["corrector"], row["field"])]
        knob_x[target, control] = float(row["field_per_x_bump_m"])
        knob_y[target, control] = float(row["field_per_y_bump_m"])
    command_vectors = (
        knob_x[:, None, :] * bump_commands[None, :, 0, None]
        + knob_y[:, None, :] * bump_commands[None, :, 1, None]
    )
    model_bpm = np.einsum("oc,tbc->tbo", bpm_response, command_vectors)
    model_target = np.einsum("toc,tbc->tbo", target_response, command_vectors)
    command_consistency_um = float(
        np.max(np.abs(model_target - bump_commands[None, :, :])) * 1e6
    )

    # Predictor inputs.  The target-local truth is deliberately not loaded
    # until command-only, two-sided, and global-MAP predictions are complete.
    bpm_orbits = np.load(source / "bpm_orbits.npy", mmap_mode="r")
    if bpm_orbits.shape[0] != nt or bpm_orbits.shape[-2:] != (nd, 2):
        raise ValueError(f"Unexpected BPM tensor shape: {bpm_orbits.shape}")
    observed = np.array(
        bpm_orbits[:, :, :, nominal_k2, :, :], dtype=float, copy=True
    )
    observed -= observed[:, :, zero_bump : zero_bump + 1, :, :]
    observed_flat = observed.reshape(*observed.shape[:3], 2 * nd)
    nr, nb = observed.shape[1:3]
    command_prediction = np.broadcast_to(
        model_target[:, None, :, :], (nt, nr, nb, 2)
    ).copy()

    two_sided_maps, neighbor_rows = build_two_sided_maps(
        bpm_maps,
        target_maps,
        one_turn,
        np.asarray([int(row["line_index"]) for row in bpm_rows]),
        np.asarray([int(row["line_index"]) for row in target_rows]),
    )
    for row, target_name in zip(neighbor_rows, target_names):
        row["target"] = target_name
        row["upstream_bpm"] = bpm_names[int(row["upstream_bpm_index"]) - 1]
        row["downstream_bpm"] = bpm_names[int(row["downstream_bpm_index"]) - 1]
    residual = observed_flat - model_bpm[:, None, :, :]
    two_sided_prediction = command_prediction.copy()
    for target in range(nt):
        upstream = int(neighbor_rows[target]["upstream_bpm_index"]) - 1
        downstream = int(neighbor_rows[target]["downstream_bpm_index"]) - 1
        channels = np.array(
            [2 * upstream, 2 * upstream + 1, 2 * downstream, 2 * downstream + 1]
        )
        nearby_residual = np.take(residual[target], channels, axis=-1)
        correction = nearby_residual @ two_sided_maps[target].T
        two_sided_prediction[target] += correction

    # Global command-conditioned MAP.  Effective control corrections are
    # standardized by the RMS size of the existing 0.5-mm local-bump commands.
    nonzero_bumps = np.arange(nb) != zero_bump
    control_scale = np.sqrt(
        np.mean(command_vectors[:, nonzero_bumps, :] ** 2, axis=(0, 1))
    )
    positive_scale = control_scale[control_scale > 0]
    scale_floor = np.median(positive_scale) * 1e-8
    control_scale = np.maximum(control_scale, scale_floor)
    standardized_response = bpm_response * control_scale[None, :]

    residual_samples = residual[:, :, nonzero_bumps, :].reshape(-1, 2 * nd)
    ratios = [float(value) for value in args.ridge_ratios.split(",")]
    if any(value < 0 or not np.isfinite(value) for value in ratios):
        raise ValueError("Ridge ratios must be finite and nonnegative")
    if args.cv_folds < 2:
        raise ValueError("At least two BPM cross-validation folds are required")
    cv_sse = np.zeros(len(ratios))
    cv_count = np.zeros(len(ratios), dtype=int)
    bpm_indices = np.arange(nd)
    for fold in range(args.cv_folds):
        validation_bpms = bpm_indices % args.cv_folds == fold
        validation_channels = np.repeat(validation_bpms, 2)
        training_channels = ~validation_channels
        train_matrix = standardized_response[training_channels]
        validation_matrix = standardized_response[validation_channels]
        train_values = residual_samples[:, training_channels]
        validation_values = residual_samples[:, validation_channels]
        for index, ratio in enumerate(ratios):
            operator = ridge_operator(train_matrix, ratio)
            coefficients = train_values @ operator.T
            prediction = coefficients @ validation_matrix.T
            cv_sse[index] += float(np.sum((prediction - validation_values) ** 2))
            cv_count[index] += prediction.size
    cv_rmse_um = np.sqrt(cv_sse / cv_count) * 1e6
    chosen_index = int(np.argmin(cv_rmse_um))
    chosen_ratio = ratios[chosen_index]
    full_operator = ridge_operator(standardized_response, chosen_ratio)
    all_coefficients = residual.reshape(-1, 2 * nd) @ full_operator.T
    all_coefficients = all_coefficients.reshape(nt, nr, nb, nc)
    global_prediction = command_prediction.copy()
    scaled_target_response = target_response * control_scale[None, None, :]
    for target in range(nt):
        global_prediction[target] += (
            all_coefficients[target] @ scaled_target_response[target].T
        )

    predictions = {
        "command_only": command_prediction,
        "two_sided_transport": two_sided_prediction,
        "global_map": global_prediction,
    }
    for name, prediction in predictions.items():
        if not np.all(np.isfinite(prediction)):
            raise ValueError(f"Non-finite {name} prediction")
        # Persist the leakage-safe predictor output before evaluation truth is
        # loaded. Downstream center-inversion studies consume these arrays and
        # do not need access to target-local truth.
        np.save(output / f"{name}_local_orbits.npy", prediction)

    # Evaluation-only truth begins here.
    target_orbits = np.load(source / "target_orbits.npy", mmap_mode="r")
    truth = np.array(
        target_orbits[:, :, :, nominal_k2, :], dtype=float, copy=True
    )
    truth -= truth[:, :, zero_bump : zero_bump + 1, :]

    summary_rows: list[dict[str, object]] = []
    per_target_rows: list[dict[str, object]] = []
    prediction_rows: list[dict[str, object]] = []
    for method, prediction in predictions.items():
        errors = prediction[:, :, nonzero_bumps, :] - truth[:, :, nonzero_bumps, :]
        summary_rows.append({"method": method, **summarize(errors)})
        for target, target_name in enumerate(target_names):
            per_target_rows.append(
                {
                    "method": method,
                    "target": target_name,
                    "target_index": target + 1,
                    **summarize(errors[target]),
                }
            )
        for target, target_name in enumerate(target_names):
            for realization in range(nr):
                for bump in np.flatnonzero(nonzero_bumps):
                    error = prediction[target, realization, bump] - truth[target, realization, bump]
                    prediction_rows.append(
                        {
                            "method": method,
                            "target": target_name,
                            "target_index": target + 1,
                            "realization": realization + 1,
                            "bump_index": bump + 1,
                            "bump_x_command_um": bump_commands[bump, 0] * 1e6,
                            "bump_y_command_um": bump_commands[bump, 1] * 1e6,
                            "truth_dx_um": truth[target, realization, bump, 0] * 1e6,
                            "truth_dy_um": truth[target, realization, bump, 1] * 1e6,
                            "prediction_dx_um": prediction[target, realization, bump, 0] * 1e6,
                            "prediction_dy_um": prediction[target, realization, bump, 1] * 1e6,
                            "error_x_um": error[0] * 1e6,
                            "error_y_um": error[1] * 1e6,
                            "error_2d_um": np.linalg.norm(error) * 1e6,
                        }
                    )

    write_rows(output / "summary.csv", summary_rows)
    write_rows(output / "per_target_summary.csv", per_target_rows)
    write_rows(output / "per_prediction_errors.csv", prediction_rows)
    write_rows(output / "two_sided_neighbors.csv", neighbor_rows)
    write_rows(
        output / "global_map_cv.csv",
        [
            {
                "ridge_ratio": ratio,
                "heldout_bpm_channel_rmse_um": cv_rmse_um[index],
                "selected": index == chosen_index,
            }
            for index, ratio in enumerate(ratios)
        ],
    )

    method_colors = {
        "command_only": "#7f7f7f",
        "two_sided_transport": "#ed7d31",
        "global_map": "#4472c4",
    }
    fig, ax = plt.subplots(figsize=(12, 5.5))
    x = np.arange(nt) + 1
    for method in predictions:
        rows = [row for row in per_target_rows if row["method"] == method]
        ax.plot(
            x,
            [float(row["rmse_2d_um"]) for row in rows],
            marker="o",
            ms=2.5,
            lw=1,
            color=method_colors[method],
            label=method.replace("_", " "),
        )
    ax.set_xlabel("Sextupole inventory index")
    ax.set_ylabel("Relative local-orbit 2D RMSE [micrometers]")
    ax.set_title("Finite-BPM prediction at all 76 sextupoles")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output / "per_target_local_orbit_rmse.png", dpi=180)
    plt.close(fig)

    summary_lines = "\n".join(
        f"| {row['method']} | {row['x_rmse_um']:.3f} | {row['y_rmse_um']:.3f} | "
        f"{row['rmse_2d_um']:.3f} | {row['median_2d_um']:.3f} | "
        f"{row['p90_2d_um']:.3f} | {row['max_2d_um']:.3f} |"
        for row in summary_rows
    )
    worst_condition = max(float(row["momentum_block_condition"]) for row in neighbor_rows)
    report = f"""# Relative local-orbit predictor comparison

All three predictors use nominal-K2 BPM orbit differences and known corrector
commands. Exact target-local orbit is loaded only after prediction and BPM-only
ridge selection, and is used solely for scoring.

- targets / latent realizations: {nt} / {nr} per target
- evaluated nonzero bumps: {int(np.sum(nonzero_bumps))} per realization
- hidden machine errors: target offset, all-other-sextupole offsets, and
  independent quadrupole strength errors from the frozen all-76 tensor
- BPM noise/offset/gain errors: none
- maximum nominal target-command consistency error: {command_consistency_um:.6g} micrometers
- global-MAP ridge ratio selected by {args.cv_folds}-fold held-out-BPM residual:
  {chosen_ratio:g}
- largest two-sided transverse momentum-block condition number:
  {worst_condition:.6g}

| method | x RMSE [um] | y RMSE [um] | 2D RMSE [um] | median [um] | P90 [um] | max [um] |
|---|---:|---:|---:|---:|---:|---:|
{summary_lines}

`command_only` uses the nominal SciBmad control-to-target response.
`two_sided_transport` adds a correction inferred from the residual x/y orbit
at the nearest upstream and downstream BPMs using nominal local transport.
`global_map` fits a regularized effective-corrector correction to all BPM
residuals, with the prior centered on the known commanded bump.

## Two-sided BPM method: required quantities

The local-orbit predictor itself runs only at nominal K2. For each target and
each bump, the machine-facing inputs are:

| quantity | source | directly available on a machine? |
|---|---|---|
| x/y orbit at the nearest upstream BPM | zero-bump and current-bump BPM readbacks | yes |
| x/y orbit at the nearest downstream BPM | zero-bump and current-bump BPM readbacks | yes |
| corrector settings defining the bump | setpoints/readbacks, converted to model fields | yes, subject to calibration |
| target K2 state | sextupole setpoint/readback, used to select nominal K2 | yes, subject to calibration |
| corrector-to-two-BPM response | nominal latest-lattice SciBmad model | model-derived; BPM part can also be measured |
| corrector-to-target response | nominal latest-lattice SciBmad model | not directly measurable at the sextupole |
| upstream-to-target and upstream-to-downstream 4D transport | nominal latest-lattice SciBmad model | model-derived, not a direct readback |

Only four BPM channels per state are consumed: upstream x/y and downstream
x/y. The full K2 scan is needed later by the magnetic-center inverse, but it is
not needed to predict these nominal-K2 local bump coordinates.

## Two-sided BPM method: implementation

For each nonzero bump, the code first subtracts the zero-bump state:

`delta y_b = y_BPM(b, nominal K2) - y_BPM(0, nominal K2)`.

Known corrector commands are propagated through nominal SciBmad responses to
obtain `delta y_model` at the two BPMs and `delta z_command` at the target.
The measured-minus-model residuals at the upstream and downstream BPMs are
called `r_u` and `r_d`.

The nominal transverse map from upstream to downstream is partitioned so that

`r_d = A_ud r_u + B_ud p_u`,

where `p_u = (delta px, delta py)`. The two unmeasured momenta are inferred by

`p_u = pinv(B_ud) (r_d - A_ud r_u)`.

The residual is then transported to the target:

`delta z_residual = A_us r_u + B_us p_u`,

and the final prediction is

`delta z_prediction = delta z_command + delta z_residual`.

The implementation builds this as one target-specific `2 x 4` matrix acting
on `[r_ux, r_uy, r_dx, r_dy]`. Neighbor selection is circular around the ring,
and the one-turn map is used when the upstream/downstream interval crosses the
lattice boundary.

## Code and artifacts

The method lives in the `finite_bpm_inversion/` directory:

- [`analyze_local_orbit_predictors.py`](../../analyze_local_orbit_predictors.py):
  constructs the two-sided matrices, predicts all relative local orbits, and
  scores them only after prediction;
- [`generate_local_orbit_models.jl`](../../generate_local_orbit_models.jl):
  generates the latest-lattice SciBmad corrector responses, cumulative maps,
  one-turn map, and element inventories;
- [`validate_local_orbit_predictor_results.py`](../../validate_local_orbit_predictor_results.py):
  independently recomputes the summary statistics and checks result counts;
- [`two_sided_neighbors.csv`](two_sided_neighbors.csv): selected neighboring
  BPMs and transport conditioning for every target;
- [`per_prediction_errors.csv`](per_prediction_errors.csv): simulation-only
  truth comparison for every method, target, realization, and bump.

## Does inference require an unmeasurable machine quantity?

No exact target-local orbit, true sextupole offset, true quadrupole errors, or
other-sextupole misalignments enter the predictor. Those hidden values are not
required by the machine-facing inference. `target_orbits.npy` is used only in
simulation after prediction to calculate the error reported above; a real
machine cannot provide this direct scoring truth.

The method does depend on quantities at the target that are model-derived
rather than directly measured: the corrector-to-target response and local
transport matrices. It also presently assumes that corrector readbacks are
converted to physical kicks with the nominal calibration, and that BPM gains,
rolls, noise, and missing channels are absent. Stable BPM offsets cancel in the
zero-bump difference, but gain/roll/calibration errors do not. Therefore the
present result is implementable from machine readbacks plus a calibrated
SciBmad model, but it is not a model-free measurement of the internal orbit.

These standalone local-orbit estimates have now been propagated through the
K2-slope center fit. See the maintained end-to-end result in
[`../two_sided_center_inversion/SUMMARY.md`](../two_sided_center_inversion/SUMMARY.md).
"""
    (output / "SUMMARY.md").write_text(report, encoding="utf-8")
    print(report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

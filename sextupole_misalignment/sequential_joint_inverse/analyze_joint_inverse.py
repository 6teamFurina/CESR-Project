#!/usr/bin/env python3
"""Benchmark sequential excitation with a shared all-target inverse.

The exact SciBmad artifact is indexed by latent machine first.  Every machine
contains all target scans, so train/validation/test splitting by machine is
leakage-safe and a joint model can consume context from all 76 scans.

The required physics baseline is the fixed-template full-BPM GLS inverse.  Two
small residual learners are then compared: a shared target-local ridge model
and a shared joint ridge model whose compact context is formed from all target
scan residuals.  A deterministic one-hidden-layer random-feature residual
model tests whether a modest nonlinearity helps without adding a heavy ML
dependency to the reproducibility path.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
import time
import tomllib
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
STUDY_ROOT = HERE.parent
GTPSA_INVERSE = STUDY_ROOT / "gtpsa_derivative_stochastic_inverse"
sys.path.insert(0, str(GTPSA_INVERSE))
import analyze_stochastic_inverse as physics  # noqa: E402


WITHOUT = "without_quadrupole_misalignment"
WITH = "with_quadrupole_misalignment"
CASES = (WITHOUT, WITH)
CASE_LABELS = {
    WITHOUT: "No quadrupole alignment drift",
    WITH: "50 um/plane RMS quadrupole alignment drift",
}


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write empty CSV: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0])
    for row in rows[1:]:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def summarize(errors_m: np.ndarray) -> dict[str, float]:
    return physics.summarize(np.asarray(errors_m))


def case_abbreviation(case: str) -> str:
    return "no_quad_align" if case == WITHOUT else "quad_align_50um"


def load_metadata(directory: Path) -> dict[str, object]:
    with (directory / "scan_metadata.toml").open("rb") as stream:
        return tomllib.load(stream)


@dataclass
class CaseData:
    name: str
    metadata: dict[str, object]
    target_names: list[str]
    bpm_names: list[str]
    bump_commands: np.ndarray
    delta_k2: np.ndarray
    measured_bpm: np.ndarray
    drift_measured_bpm: np.ndarray
    reference_bpm: np.ndarray
    reference_target: np.ndarray
    truth: np.ndarray
    drift_response: np.ndarray


def load_case(root: Path, case: str, model_target_names: list[str]) -> CaseData:
    source = root / case
    metadata = load_metadata(source)
    target_names = (source / "target_names.txt").read_text(encoding="utf-8").splitlines()
    bpm_names = (source / "bpm_names.txt").read_text(encoding="utf-8").splitlines()
    if target_names != model_target_names[: len(target_names)]:
        raise ValueError(f"Target order mismatch in {case}")
    bump_rows = read_rows(source / "bump_points.csv")
    bump_commands = np.asarray(
        [
            (float(row["bump_x_command_m"]), float(row["bump_y_command_m"]))
            for row in bump_rows
        ]
    )
    delta_k2 = np.asarray(metadata["k2_delta_m3"], dtype=float)
    bpm = np.asarray(np.load(source / "bpm_orbits.npy", mmap_mode="r"), dtype=float)
    drift_bpm = np.asarray(
        np.load(source / "drift_bpm_orbits.npy", mmap_mode="r"), dtype=float
    )
    reference_bpm = np.asarray(
        np.load(source / "reference_bpm_orbits.npy", mmap_mode="r"), dtype=float
    )
    reference_target = np.asarray(
        np.load(source / "reference_target_orbits.npy", mmap_mode="r"), dtype=float
    )
    latent_root = root / "paired_latents"
    gains = np.asarray(np.load(latent_root / "bpm_gain_errors.npy"), dtype=float)
    sextupole_offsets = np.asarray(
        np.load(latent_root / "sextupole_offsets.npy"), dtype=float
    )[:, : len(target_names)]
    expected = (
        int(metadata["machine_count"]),
        len(target_names),
        int(metadata["bump_count"]),
        int(metadata["k2_count"]),
        len(bpm_names),
        2,
    )
    if bpm.shape != expected or drift_bpm.shape != expected:
        raise ValueError(f"Unexpected {case} BPM shapes: {bpm.shape}, {drift_bpm.shape}")
    if gains.shape != (expected[0], expected[4], 2):
        raise ValueError(f"Unexpected BPM gain shape: {gains.shape}")
    gain_factor = 1.0 + gains[:, None, None, None, :, :]
    measured = bpm * gain_factor
    drift_measured = drift_bpm * gain_factor
    measured_reference = reference_bpm * (1.0 + gains)
    truth = sextupole_offsets - reference_target
    drift_response = physics.recover_drift_response(
        measured,
        drift_measured,
        float(metadata["drift_halfwidth_m"]),
    )
    arrays = (measured, drift_measured, measured_reference, truth, drift_response)
    if not all(np.all(np.isfinite(array)) for array in arrays):
        raise ValueError(f"Non-finite values in {case}")
    return CaseData(
        name=case,
        metadata=metadata,
        target_names=target_names,
        bpm_names=bpm_names,
        bump_commands=bump_commands,
        delta_k2=delta_k2,
        measured_bpm=measured,
        drift_measured_bpm=drift_measured,
        reference_bpm=measured_reference.reshape(expected[0], -1),
        reference_target=reference_target,
        truth=truth,
        drift_response=drift_response,
    )


def machine_split(machine_count: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if machine_count < 4:
        raise ValueError("At least four exact machines are required for leakage-safe splits")
    train_count = max(2, int(round(0.625 * machine_count)))
    val_count = max(1, int(round(0.1875 * machine_count)))
    if train_count + val_count >= machine_count:
        train_count = machine_count - 2
        val_count = 1
    indices = np.arange(machine_count)
    return (
        indices[:train_count],
        indices[train_count : train_count + val_count],
        indices[train_count + val_count :],
    )


def state_mean_random_walk_covariance(
    repeats: int,
    state_count: int,
    endpoint_rms_m: float,
) -> np.ndarray:
    """Covariance of per-state means for an interleaved scalar random walk."""
    read_count = repeats * state_count
    step_variance = endpoint_rms_m**2 / max(read_count - 1, 1)
    covariance = np.zeros((state_count, state_count))
    for p in range(state_count):
        for q in range(state_count):
            total = 0.0
            for repeat in range(repeats):
                n = state_count * repeat + p
                last_before = min(repeats - 1, math.floor((n - q) / state_count))
                count_before = max(last_before + 1, 0)
                arithmetic = (
                    state_count * last_before * (last_before + 1) / 2.0
                    + count_before * (q + 1)
                    if count_before
                    else 0.0
                )
                total += arithmetic + (repeats - count_before) * (n + 1)
            covariance[p, q] = step_variance * total / repeats**2
    covariance = 0.5 * (covariance + covariance.T)
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    return (eigenvectors * np.maximum(eigenvalues, 0.0)) @ eigenvectors.T


def covariance_square_root(covariance: np.ndarray) -> np.ndarray:
    values, vectors = np.linalg.eigh(0.5 * (covariance + covariance.T))
    return vectors @ np.diag(np.sqrt(np.maximum(values, 0.0)))


def augmented_gradients(
    data: CaseData,
    machine_indices: np.ndarray,
    augmentations: int,
    repeats: int,
    bpm_noise_rms_m: float,
    drift_endpoint_rms_m: float,
    seed: int,
) -> np.ndarray:
    measured = data.measured_bpm[machine_indices]
    deterministic, _ = physics.parity_gradients(
        measured, data.delta_k2, data.bump_commands
    )
    deterministic = np.asarray(deterministic)
    machine_count, target_count, _, channel_count = deterministic.shape
    schedule = physics.signed_state_indices(data.bump_commands, data.delta_k2)
    normalization = float(np.ptp(data.delta_k2)) * 2.0 * float(
        np.max(np.abs(data.bump_commands))
    )
    white_std = 2.0 * bpm_noise_rms_m / (np.sqrt(repeats) * normalization)
    walk_covariance = state_mean_random_walk_covariance(
        repeats, len(schedule), drift_endpoint_rms_m
    )
    walk_sqrt = covariance_square_root(walk_covariance)
    response = data.drift_response[machine_indices]
    rng = np.random.default_rng(seed)
    result = np.empty(
        (augmentations, machine_count, target_count, 2, channel_count), dtype=float
    )
    for augmentation in range(augmentations):
        white = white_std * rng.standard_normal(deterministic.shape)
        state_means = rng.standard_normal(
            (machine_count, target_count, len(schedule))
        ) @ walk_sqrt.T
        drift = np.zeros_like(deterministic)
        for position, (block, sign, bump, k2) in enumerate(schedule):
            drift[:, :, block, :] += (
                sign
                * state_means[:, :, position, None]
                * response[:, :, bump, k2, :]
                / normalization
            )
        result[augmentation] = deterministic + white + drift
    return result


def physics_estimates(
    gradients: np.ndarray,
    design: np.ndarray,
    left_inverses: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    right = np.concatenate((gradients[..., 0, :], gradients[..., 1, :]), axis=-1)
    estimates = np.einsum("tij,amtj->amti", left_inverses, right)
    predicted = np.einsum("tji,amti->amtj", design, estimates)
    return estimates, right - predicted


@dataclass
class PCAState:
    mean: np.ndarray
    scale: np.ndarray
    components: np.ndarray


def fit_pca(values: np.ndarray, modes: int) -> PCAState:
    values = np.asarray(values, dtype=float)
    mean = np.mean(values, axis=0)
    scale = np.std(values, axis=0)
    scale = np.where(scale > 1.0e-14 * max(float(np.max(scale)), 1.0), scale, 1.0)
    normalized = (values - mean) / scale
    _, _, vt = np.linalg.svd(normalized, full_matrices=False)
    retained = min(modes, vt.shape[0], vt.shape[1])
    return PCAState(mean=mean, scale=scale, components=vt[:retained])


def transform_pca(values: np.ndarray, state: PCAState) -> np.ndarray:
    return ((values - state.mean) / state.scale) @ state.components.T


@dataclass
class FeaturePipeline:
    design: np.ndarray
    left_inverses: np.ndarray
    residual_pca: PCAState
    reference_pca: PCAState
    global_pca: PCAState
    target_count: int

    @classmethod
    def fit(
        cls,
        gradients: np.ndarray,
        reference_bpm: np.ndarray,
        design: np.ndarray,
        left_inverses: np.ndarray,
        residual_modes: int,
        reference_modes: int,
        global_modes: int,
    ) -> "FeaturePipeline":
        estimates, residuals = physics_estimates(gradients, design, left_inverses)
        a, m, t, channels = residuals.shape
        residual_pca = fit_pca(residuals.reshape(-1, channels), residual_modes)
        residual_scores = transform_pca(
            residuals.reshape(-1, channels), residual_pca
        ).reshape(a, m, t, -1)
        reference_pca = fit_pca(reference_bpm, reference_modes)
        global_raw = np.concatenate((estimates / 3.0e-4, residual_scores), axis=-1)
        global_pca = fit_pca(global_raw.reshape(a * m, -1), global_modes)
        return cls(
            design=design,
            left_inverses=left_inverses,
            residual_pca=residual_pca,
            reference_pca=reference_pca,
            global_pca=global_pca,
            target_count=t,
        )

    def transform(
        self,
        gradients: np.ndarray,
        reference_bpm: np.ndarray,
        truth: np.ndarray,
    ) -> dict[str, np.ndarray]:
        estimates, residuals = physics_estimates(
            gradients, self.design, self.left_inverses
        )
        a, m, t, channels = residuals.shape
        if t != self.target_count:
            raise ValueError("Target count changed across feature transforms")
        residual_scores = transform_pca(
            residuals.reshape(-1, channels), self.residual_pca
        ).reshape(a, m, t, -1)
        reference_scores = transform_pca(reference_bpm, self.reference_pca)
        reference_scores = np.broadcast_to(
            reference_scores[None, :, None, :],
            (a, m, t, reference_scores.shape[-1]),
        )
        target_identity = np.broadcast_to(
            np.eye(t)[None, None, :, :], (a, m, t, t)
        )
        local = np.concatenate(
            (estimates / 3.0e-4, residual_scores, reference_scores, target_identity),
            axis=-1,
        )
        global_raw = np.concatenate((estimates / 3.0e-4, residual_scores), axis=-1)
        global_scores = transform_pca(
            global_raw.reshape(a * m, -1), self.global_pca
        ).reshape(a, m, -1)
        global_scores = np.broadcast_to(
            global_scores[:, :, None, :], (a, m, t, global_scores.shape[-1])
        )
        joint = np.concatenate((local, global_scores), axis=-1)
        repeated_truth = np.broadcast_to(truth[None], estimates.shape)
        correction_um = (repeated_truth - estimates) * 1.0e6
        return {
            "physics": estimates,
            "local_features": local.reshape(-1, local.shape[-1]),
            "joint_features": joint.reshape(-1, joint.shape[-1]),
            "correction_um": correction_um.reshape(-1, 2),
            "truth": repeated_truth,
        }

    def save(self, path: Path) -> None:
        np.savez_compressed(
            path,
            design=self.design,
            left_inverses=self.left_inverses,
            residual_mean=self.residual_pca.mean,
            residual_scale=self.residual_pca.scale,
            residual_components=self.residual_pca.components,
            reference_mean=self.reference_pca.mean,
            reference_scale=self.reference_pca.scale,
            reference_components=self.reference_pca.components,
            global_mean=self.global_pca.mean,
            global_scale=self.global_pca.scale,
            global_components=self.global_pca.components,
        )


@dataclass
class RidgeModel:
    mean: np.ndarray
    scale: np.ndarray
    weights: np.ndarray
    alpha: float

    def standardized(self, values: np.ndarray) -> np.ndarray:
        return (values - self.mean) / self.scale

    def predict(self, values: np.ndarray) -> np.ndarray:
        normalized = self.standardized(values)
        augmented = np.concatenate((normalized, np.ones((len(normalized), 1))), axis=1)
        return augmented @ self.weights


def fit_ridge(
    train_x: np.ndarray,
    train_y: np.ndarray,
    val_x: np.ndarray,
    val_y: np.ndarray,
    alphas: np.ndarray,
) -> tuple[RidgeModel, list[dict[str, float]]]:
    mean = np.mean(train_x, axis=0)
    scale = np.std(train_x, axis=0)
    scale = np.where(scale > 1.0e-12, scale, 1.0)
    x = (train_x - mean) / scale
    xv = (val_x - mean) / scale
    x = np.concatenate((x, np.ones((len(x), 1))), axis=1)
    xv = np.concatenate((xv, np.ones((len(xv), 1))), axis=1)
    gram = x.T @ x
    rhs = x.T @ train_y
    penalty = np.eye(gram.shape[0])
    penalty[-1, -1] = 0.0
    rows: list[dict[str, float]] = []
    best: tuple[float, np.ndarray, float] | None = None
    for alpha in alphas:
        weights = np.linalg.solve(gram + float(alpha) * penalty, rhs)
        residual = xv @ weights - val_y
        rmse = float(np.sqrt(np.mean(np.sum(residual * residual, axis=-1))))
        rows.append({"alpha": float(alpha), "validation_rmse_2d_um": rmse})
        if best is None or rmse < best[0]:
            best = (rmse, weights, float(alpha))
    assert best is not None
    return RidgeModel(mean=mean, scale=scale, weights=best[1], alpha=best[2]), rows


@dataclass
class RandomFeatureModel:
    base: RidgeModel
    projection: np.ndarray
    bias: np.ndarray
    hidden_scale: float

    def features(self, values: np.ndarray) -> np.ndarray:
        normalized = (values - self.base.mean[: values.shape[1]]) / self.base.scale[
            : values.shape[1]
        ]
        hidden = np.tanh(
            self.hidden_scale * (normalized @ self.projection) + self.bias
        )
        return np.concatenate((values, hidden), axis=1)

    def predict(self, values: np.ndarray) -> np.ndarray:
        return self.base.predict(self.features(values))


def fit_random_feature_model(
    train_x: np.ndarray,
    train_y: np.ndarray,
    val_x: np.ndarray,
    val_y: np.ndarray,
    alphas: np.ndarray,
    hidden_units: int,
    seed: int,
) -> tuple[RandomFeatureModel, list[dict[str, float]]]:
    input_mean = np.mean(train_x, axis=0)
    input_scale = np.std(train_x, axis=0)
    input_scale = np.where(input_scale > 1.0e-12, input_scale, 1.0)
    normalized_train = (train_x - input_mean) / input_scale
    normalized_val = (val_x - input_mean) / input_scale
    rng = np.random.default_rng(seed)
    projection = rng.standard_normal((train_x.shape[1], hidden_units)) / np.sqrt(
        train_x.shape[1]
    )
    bias = rng.uniform(-0.5, 0.5, hidden_units)
    selection_rows: list[dict[str, float]] = []
    best: tuple[float, RandomFeatureModel] | None = None
    for hidden_scale in (0.5, 1.0, 2.0):
        train_hidden = np.tanh(hidden_scale * (normalized_train @ projection) + bias)
        val_hidden = np.tanh(hidden_scale * (normalized_val @ projection) + bias)
        extended_train = np.concatenate((train_x, train_hidden), axis=1)
        extended_val = np.concatenate((val_x, val_hidden), axis=1)
        ridge, rows = fit_ridge(
            extended_train, train_y, extended_val, val_y, alphas
        )
        for row in rows:
            selection_rows.append(
                {
                    "hidden_scale": hidden_scale,
                    "alpha": row["alpha"],
                    "validation_rmse_2d_um": row["validation_rmse_2d_um"],
                }
            )
        score = min(row["validation_rmse_2d_um"] for row in rows)
        candidate = RandomFeatureModel(
            base=ridge,
            projection=projection,
            bias=bias,
            hidden_scale=hidden_scale,
        )
        if best is None or score < best[0]:
            best = (score, candidate)
    assert best is not None
    return best[1], selection_rows


@dataclass
class TrainedBundle:
    case: str
    pipeline: FeaturePipeline
    local_ridge: RidgeModel
    joint_ridge: RidgeModel
    joint_random: RandomFeatureModel


def train_bundle(
    case: CaseData,
    train_indices: np.ndarray,
    val_indices: np.ndarray,
    design: np.ndarray,
    left_inverses: np.ndarray,
    args: argparse.Namespace,
) -> tuple[TrainedBundle, list[dict[str, object]]]:
    train_gradients = augmented_gradients(
        case,
        train_indices,
        args.train_augmentations,
        args.repeats,
        args.bpm_noise_rms_m,
        args.drift_endpoint_rms_m,
        args.measurement_seed,
    )
    val_gradients = augmented_gradients(
        case,
        val_indices,
        args.validation_augmentations,
        args.repeats,
        args.bpm_noise_rms_m,
        args.drift_endpoint_rms_m,
        args.measurement_seed + 10_000,
    )
    pipeline = FeaturePipeline.fit(
        train_gradients,
        case.reference_bpm[train_indices],
        design,
        left_inverses,
        args.residual_modes,
        args.reference_modes,
        args.global_modes,
    )
    train_view = pipeline.transform(
        train_gradients, case.reference_bpm[train_indices], case.truth[train_indices]
    )
    val_view = pipeline.transform(
        val_gradients, case.reference_bpm[val_indices], case.truth[val_indices]
    )
    alphas = np.logspace(-4, 4, 17)
    local_ridge, local_rows = fit_ridge(
        train_view["local_features"],
        train_view["correction_um"],
        val_view["local_features"],
        val_view["correction_um"],
        alphas,
    )
    joint_ridge, joint_rows = fit_ridge(
        train_view["joint_features"],
        train_view["correction_um"],
        val_view["joint_features"],
        val_view["correction_um"],
        alphas,
    )
    joint_random, random_rows = fit_random_feature_model(
        train_view["joint_features"],
        train_view["correction_um"],
        val_view["joint_features"],
        val_view["correction_um"],
        alphas,
        args.hidden_units,
        args.model_seed,
    )
    selection: list[dict[str, object]] = []
    for model, rows in (
        ("shared_target_local_ridge", local_rows),
        ("shared_joint_ridge", joint_rows),
        ("shared_joint_random_feature", random_rows),
    ):
        for row in rows:
            selection.append({"training_case": case.name, "model": model, **row})
    return (
        TrainedBundle(
            case=case.name,
            pipeline=pipeline,
            local_ridge=local_ridge,
            joint_ridge=joint_ridge,
            joint_random=joint_random,
        ),
        selection,
    )


def evaluation_rows(
    bundle: TrainedBundle | None,
    evaluation: CaseData,
    test_indices: np.ndarray,
    design: np.ndarray,
    left_inverses: np.ndarray,
    args: argparse.Namespace,
) -> tuple[list[dict[str, object]], dict[str, np.ndarray], list[dict[str, object]]]:
    gradients = augmented_gradients(
        evaluation,
        test_indices,
        args.test_augmentations,
        args.repeats,
        args.bpm_noise_rms_m,
        args.drift_endpoint_rms_m,
        args.measurement_seed + 20_000,
    )
    if bundle is None:
        estimates, _ = physics_estimates(gradients, design, left_inverses)
        truth = np.broadcast_to(evaluation.truth[test_indices][None], estimates.shape)
        predictions = {"physics_gls": estimates}
        training_case = "fixed_nominal_physics"
    else:
        view = bundle.pipeline.transform(
            gradients,
            evaluation.reference_bpm[test_indices],
            evaluation.truth[test_indices],
        )
        estimates = view["physics"]
        truth = view["truth"]
        shape = estimates.shape
        predictions = {
            "shared_target_local_ridge": estimates
            + bundle.local_ridge.predict(view["local_features"]).reshape(shape) * 1.0e-6,
            "shared_joint_ridge": estimates
            + bundle.joint_ridge.predict(view["joint_features"]).reshape(shape) * 1.0e-6,
            "shared_joint_random_feature": estimates
            + bundle.joint_random.predict(view["joint_features"]).reshape(shape) * 1.0e-6,
        }
        training_case = bundle.case

    summary_rows: list[dict[str, object]] = []
    target_rows: list[dict[str, object]] = []
    saved: dict[str, np.ndarray] = {"truth": truth}
    for model, estimate in predictions.items():
        errors = estimate - truth
        metrics = summarize(errors)
        per_machine = np.sqrt(np.mean(np.sum(errors * errors, axis=-1), axis=-1)) * 1.0e6
        per_target_rmse = np.sqrt(np.mean(np.sum(errors * errors, axis=-1), axis=(0, 1))) * 1.0e6
        summary_rows.append(
            {
                "training_case": training_case,
                "evaluation_case": evaluation.name,
                "model": model,
                "held_out_machine_count": len(test_indices),
                "measurement_augmentations": args.test_augmentations,
                "fit_count": int(np.prod(errors.shape[:-1])),
                **metrics,
                "fraction_below_50um": float(
                    np.mean(np.linalg.norm(errors, axis=-1) * 1.0e6 < 50.0)
                ),
                "median_machine_rmse_2d_um": float(np.median(per_machine)),
                "worst_machine_rmse_2d_um": float(np.max(per_machine)),
                "worst_target_rmse_2d_um": float(np.max(per_target_rmse)),
                "preferred_30um_aggregate_gate": bool(metrics["rmse_2d_um"] < 30.0),
                "hard_50um_rmse_p99_all_target_gate": bool(
                    metrics["rmse_2d_um"] < 50.0
                    and metrics["p99_2d_um"] < 50.0
                    and float(np.max(per_target_rmse)) < 50.0
                ),
            }
        )
        for target, name in enumerate(evaluation.target_names):
            target_rows.append(
                {
                    "training_case": training_case,
                    "evaluation_case": evaluation.name,
                    "model": model,
                    "target": name,
                    "target_index": target + 1,
                    **summarize(errors[:, :, target]),
                }
            )
        saved[f"{model}_estimate"] = estimate
        saved[f"{model}_error"] = errors
    return summary_rows, saved, target_rows


def save_ridge(path: Path, model: RidgeModel) -> None:
    np.savez_compressed(
        path,
        mean=model.mean,
        scale=model.scale,
        weights=model.weights,
        alpha=np.asarray(model.alpha),
    )


def main() -> int:
    started = time.perf_counter()
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--scan-root", type=Path, default=HERE / "results" / "exact_joint_machines"
    )
    parser.add_argument(
        "--model-dir",
        type=Path,
        default=STUDY_ROOT / "finite_bpm_inversion" / "results" / "local_orbit_model",
    )
    parser.add_argument("--sextupole-length-m", type=float, default=0.272)
    parser.add_argument("--repeats", type=int, default=3072)
    parser.add_argument("--bpm-noise-rms-m", type=float, default=5.0e-6)
    parser.add_argument("--drift-endpoint-rms-m", type=float, default=1.0e-5)
    parser.add_argument("--train-augmentations", type=int, default=8)
    parser.add_argument("--validation-augmentations", type=int, default=8)
    parser.add_argument("--test-augmentations", type=int, default=32)
    parser.add_argument("--residual-modes", type=int, default=24)
    parser.add_argument("--reference-modes", type=int, default=8)
    parser.add_argument("--global-modes", type=int, default=16)
    parser.add_argument("--hidden-units", type=int, default=96)
    parser.add_argument("--measurement-seed", type=int, default=20260923)
    parser.add_argument("--model-seed", type=int, default=20261023)
    parser.add_argument(
        "--output-dir", type=Path, default=HERE / "results" / "joint_inverse_analysis"
    )
    args = parser.parse_args()
    scan_root = args.scan_root.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)

    model_rows = read_rows(args.model_dir.resolve() / "target_locations.csv")
    model_target_names = [row["target"] for row in model_rows]
    data = {case: load_case(scan_root, case, model_target_names) for case in CASES}
    machine_count = int(data[WITHOUT].metadata["machine_count"])
    if int(data[WITH].metadata["machine_count"]) != machine_count:
        raise ValueError("Paired cases have different machine counts")
    train_indices, val_indices, test_indices = machine_split(machine_count)
    if data[WITHOUT].target_names != data[WITH].target_names:
        raise ValueError("Paired cases have different target inventories")

    templates = physics.source_templates(args.model_dir.resolve(), args.sextupole_length_m)
    target_count = len(data[WITHOUT].target_names)
    design = physics.center_design(templates[:target_count])
    left_inverses = np.asarray(
        [np.linalg.pinv(matrix, rcond=1.0e-12) for matrix in design]
    )

    bundles: dict[str, TrainedBundle] = {}
    selection_rows: list[dict[str, object]] = []
    for case in CASES:
        bundle, rows = train_bundle(
            data[case], train_indices, val_indices, design, left_inverses, args
        )
        bundles[case] = bundle
        selection_rows.extend(rows)
        prefix = output / f"model_{case_abbreviation(case)}"
        bundle.pipeline.save(prefix.with_name(prefix.name + "_pipeline.npz"))
        save_ridge(prefix.with_name(prefix.name + "_local_ridge.npz"), bundle.local_ridge)
        save_ridge(prefix.with_name(prefix.name + "_joint_ridge.npz"), bundle.joint_ridge)
        np.savez_compressed(
            prefix.with_name(prefix.name + "_joint_random_feature.npz"),
            projection=bundle.joint_random.projection,
            bias=bundle.joint_random.bias,
            hidden_scale=np.asarray(bundle.joint_random.hidden_scale),
            base_mean=bundle.joint_random.base.mean,
            base_scale=bundle.joint_random.base.scale,
            base_weights=bundle.joint_random.base.weights,
            base_alpha=np.asarray(bundle.joint_random.base.alpha),
        )

    summary_rows: list[dict[str, object]] = []
    target_rows: list[dict[str, object]] = []
    prediction_payload: dict[str, np.ndarray] = {}
    for evaluation_case in CASES:
        rows, saved, per_target = evaluation_rows(
            None,
            data[evaluation_case],
            test_indices,
            design,
            left_inverses,
            args,
        )
        summary_rows.extend(rows)
        target_rows.extend(per_target)
        for key, value in saved.items():
            prediction_payload[
                f"physics__{case_abbreviation(evaluation_case)}__{key}"
            ] = value
        for training_case in CASES:
            rows, saved, per_target = evaluation_rows(
                bundles[training_case],
                data[evaluation_case],
                test_indices,
                design,
                left_inverses,
                args,
            )
            summary_rows.extend(rows)
            target_rows.extend(per_target)
            prefix = (
                f"train_{case_abbreviation(training_case)}"
                f"__eval_{case_abbreviation(evaluation_case)}"
            )
            for key, value in saved.items():
                prediction_payload[f"{prefix}__{key}"] = value

    write_rows(output / "summary.csv", summary_rows)
    write_rows(output / "per_target_summary.csv", target_rows)
    write_rows(output / "model_selection.csv", selection_rows)
    np.savez_compressed(output / "held_out_predictions.npz", **prediction_payload)

    latent_root = scan_root / "paired_latents"
    alignment_normals = np.load(
        latent_root / "quadrupole_alignment_standard_normals.npy"
    )
    alignment_rms = float(data[WITH].metadata["quadrupole_alignment_rms_m_per_plane"])
    actual_offsets = alignment_rms * alignment_normals
    paired_reference_delta = data[WITH].reference_bpm - data[WITHOUT].reference_bpm
    paired_truth_delta = data[WITH].truth - data[WITHOUT].truth
    diagnostics = {
        "machine_count": machine_count,
        "target_count": target_count,
        "train_machine_indices_zero_based": train_indices.tolist(),
        "validation_machine_indices_zero_based": val_indices.tolist(),
        "test_machine_indices_zero_based": test_indices.tolist(),
        "quadrupole_alignment_requested_rms_um_per_plane": alignment_rms * 1.0e6,
        "quadrupole_alignment_realized_x_rms_um": float(
            np.sqrt(np.mean(actual_offsets[..., 0] ** 2)) * 1.0e6
        ),
        "quadrupole_alignment_realized_y_rms_um": float(
            np.sqrt(np.mean(actual_offsets[..., 1] ** 2)) * 1.0e6
        ),
        "paired_reference_bpm_change_rms_um": float(
            np.sqrt(np.mean(paired_reference_delta**2)) * 1.0e6
        ),
        "paired_reference_bpm_change_p99_um": float(
            np.percentile(np.abs(paired_reference_delta) * 1.0e6, 99)
        ),
        "paired_beam_relative_truth_change_rms_2d_um": float(
            np.sqrt(np.mean(np.sum(paired_truth_delta**2, axis=-1))) * 1.0e6
        ),
        "without_truth_fraction_outside_bump_radius": float(
            np.mean(
                np.linalg.norm(data[WITHOUT].truth, axis=-1)
                > float(data[WITHOUT].metadata["bump_amplitude_m"])
            )
        ),
        "with_truth_fraction_outside_bump_radius": float(
            np.mean(
                np.linalg.norm(data[WITH].truth, axis=-1)
                > float(data[WITH].metadata["bump_amplitude_m"])
            )
        ),
    }
    (output / "diagnostics.json").write_text(
        json.dumps(diagnostics, indent=2) + "\n", encoding="utf-8"
    )

    same_case_rows = [
        row
        for row in summary_rows
        if (
            row["training_case"] == "fixed_nominal_physics"
            or row["training_case"] == row["evaluation_case"]
        )
    ]
    models = [
        "physics_gls",
        "shared_target_local_ridge",
        "shared_joint_ridge",
        "shared_joint_random_feature",
    ]
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 5.0), sharey=True)
    for axis, case in zip(axes, CASES):
        values = []
        for model in models:
            candidates = [
                row
                for row in same_case_rows
                if row["evaluation_case"] == case and row["model"] == model
            ]
            if len(candidates) != 1:
                raise ValueError(f"Missing same-case result for {case}/{model}")
            values.append(float(candidates[0]["rmse_2d_um"]))
        axis.bar(np.arange(len(models)), values, color=["#777777", "#4c78a8", "#59a14f", "#e15759"])
        axis.axhline(50.0, color="#b22222", ls="--", lw=1.2, label="50 um reference gate")
        axis.set_xticks(
            np.arange(len(models)),
            ["Physics\nGLS", "Local\nridge", "Joint\nridge", "Joint random-\nfeature"],
        )
        axis.set_title(CASE_LABELS[case])
        axis.grid(axis="y", alpha=0.25)
    axes[0].set_ylabel("Held-out beam-relative center 2D RMSE [um]")
    axes[1].legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(output / "held_out_model_comparison.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.8))
    for case, color in ((WITHOUT, "#4c78a8"), (WITH, "#e15759")):
        radial = np.linalg.norm(data[case].truth, axis=-1).ravel() * 1.0e6
        axes[0].hist(
            radial,
            bins=36,
            histtype="step",
            lw=1.7,
            color=color,
            label=CASE_LABELS[case],
        )
    axes[0].axvline(
        float(data[WITH].metadata["bump_amplitude_m"]) * 1.0e6,
        color="#222222",
        ls="--",
        lw=1.2,
        label="1.5 mm bump radius",
    )
    axes[0].set_xlabel("Beam-relative center radius [um]")
    axes[0].set_ylabel("Machine-target count")
    axes[0].set_title("Excitation-domain exposure")
    axes[0].legend(fontsize=8)
    axes[0].grid(alpha=0.2)
    reference_rms_by_machine = np.sqrt(
        np.mean(paired_reference_delta.reshape(machine_count, -1) ** 2, axis=1)
    ) * 1.0e6
    truth_rms_by_machine = np.sqrt(
        np.mean(np.sum(paired_truth_delta**2, axis=-1), axis=1)
    ) * 1.0e6
    axes[1].scatter(reference_rms_by_machine, truth_rms_by_machine, color="#e15759")
    for machine, (x_value, y_value) in enumerate(
        zip(reference_rms_by_machine, truth_rms_by_machine), start=1
    ):
        axes[1].annotate(str(machine), (x_value, y_value), fontsize=7, xytext=(3, 2), textcoords="offset points")
    axes[1].set_xlabel("Paired reference BPM change RMS [um]")
    axes[1].set_ylabel("Beam-relative truth change 2D RMS [um]")
    axes[1].set_title("Uncorrected quadrupole-drift consequence")
    axes[1].grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(output / "quadrupole_drift_domain_diagnostic.png", dpi=180)
    plt.close(fig)

    primary_table_rows = []
    for case in CASES:
        for model in models:
            candidates = [
                row
                for row in same_case_rows
                if row["evaluation_case"] == case and row["model"] == model
            ]
            row = candidates[0]
            primary_table_rows.append(
                f"| {CASE_LABELS[case]} | {model} | "
                f"{float(row['rmse_2d_um']):.3f} | {float(row['p90_2d_um']):.3f} | "
                f"{float(row['p99_2d_um']):.3f} | {float(row['worst_target_rmse_2d_um']):.3f} | "
                f"{100.0*float(row['fraction_below_50um']):.2f}% |"
            )
    ood_candidates = [
        row
        for row in summary_rows
        if row["training_case"] == WITHOUT
        and row["evaluation_case"] == WITH
        and row["model"] == "shared_joint_ridge"
    ]
    ood_row = ood_candidates[0]
    best_by_case = {}
    for case in CASES:
        candidates = [
            row
            for row in same_case_rows
            if row["evaluation_case"] == case and row["model"] != "physics_gls"
        ]
        best_by_case[case] = min(candidates, key=lambda row: float(row["rmse_2d_um"]))
    joint_comparison = {}
    for case in CASES:
        local = next(
            row
            for row in same_case_rows
            if row["evaluation_case"] == case
            and row["model"] == "shared_target_local_ridge"
        )
        joint = next(
            row
            for row in same_case_rows
            if row["evaluation_case"] == case and row["model"] == "shared_joint_ridge"
        )
        joint_comparison[case] = {
            "relative_change_percent": 100.0
            * (float(joint["rmse_2d_um"]) / float(local["rmse_2d_um"]) - 1.0),
        }

    report = f"""# Sequential excitation and all-target joint-inverse pilot

## Protocol

The exact dataset contains {machine_count} latest-lattice SciBmad latent
machines.  In each machine, all {target_count} sextupole offsets, BPM gains,
corrector gains, K2 gains, quadrupole strength errors, and quadrupole rolls are
fixed while the {target_count} sextupoles are excited one at a time.  The
paired quadrupole-alignment case adds independent Gaussian x/y displacement to
each physical quadrupole with `50 um` RMS per plane, coherent across its slices
and fixed throughout all target scans.

Each physical case contains
`{int(data[WITHOUT].metadata['total_exact_states']):,}` exact SciBmad state
lanes: 15 scan states and 15 paired drift-secant states per target and machine.

Machines are split, without target-level leakage, into {len(train_indices)}
training, {len(val_indices)} validation, and {len(test_indices)} held-out test
machines.  Each test machine has {args.test_augmentations} independent
measurement realizations.  Measurement augmentation uses `5 um` RMS white
noise per BPM plane/read, {args.repeats} repeated balanced eight-state cycles,
and a scalar random walk with `10 um` RMS endpoint change propagated through
the paired exact SciBmad drift secant.

The learned target is the two-component beam-relative sextupole magnetic
center.  All models receive full-ring BPM parity contrasts and nominal command
values; latent errors and exact target-local orbit are evaluation-only.

## Held-out result

| quadrupole-alignment input | inverse | 2D RMSE [um] | P90 [um] | P99 [um] | worst-target RMSE [um] | below 50 um |
|---|---|---:|---:|---:|---:|---:|
{chr(10).join(primary_table_rows)}

The best learned same-distribution result without quadrupole alignment is
`{best_by_case[WITHOUT]['model']}` at
`{float(best_by_case[WITHOUT]['rmse_2d_um']):.3f} um` RMSE.  With the paired
50-um/plane quadrupole drift enabled, the best learned result is
`{best_by_case[WITH]['model']}` at
`{float(best_by_case[WITH]['rmse_2d_um']):.3f} um` RMSE.  Training the joint
ridge only on the no-alignment distribution and evaluating it on the aligned
case gives `{float(ood_row['rmse_2d_um']):.3f} um` RMSE; this is the explicit
distribution-shift check and must not be replaced by the in-distribution row.

Adding all-target context changes ridge RMSE relative to the matched local
shared model by `{joint_comparison[WITHOUT]['relative_change_percent']:+.3f}%`
without quadrupole alignment and
`{joint_comparison[WITH]['relative_change_percent']:+.3f}%` with the 50-um
drift.  A negative value is the predeclared evidence that joint context helped;
a positive value means this ensemble does not support the added joint-model
complexity.  The strict reference gate requires aggregate RMSE, P99, and every
target-level RMSE to remain below 50 um.
The best learned no-alignment row reports
`{'PASS' if best_by_case[WITHOUT]['hard_50um_rmse_p99_all_target_gate'] else 'FAIL'}`
for that gate; the best learned 50-um/plane row reports
`{'PASS' if best_by_case[WITH]['hard_50um_rmse_p99_all_target_gate'] else 'FAIL'}`.

## Quadrupole-drift interpretation

The requested 50-um setting is interpreted as a per-plane RMS.  The finite
draw realizes `{diagnostics['quadrupole_alignment_realized_x_rms_um']:.3f} um`
in x and `{diagnostics['quadrupole_alignment_realized_y_rms_um']:.3f} um` in y.
If the facility value instead denotes an isotropic two-dimensional radial RMS,
the corresponding per-plane RMS would be about `35.4 um`; correlated girder
motion is also a distinct prior.  Neither alternative is silently folded into
the primary 50-um/plane row.
Before any orbit correction, it changes the paired full-ring reference BPM
orbit by `{diagnostics['paired_reference_bpm_change_rms_um']:.3f} um` RMS and
the beam-relative center truth by
`{diagnostics['paired_beam_relative_truth_change_rms_2d_um']:.3f} um` 2D RMS.
The fraction of truths outside the maintained 1.5-mm bump radius changes from
`{100.0*diagnostics['without_truth_fraction_outside_bump_radius']:.3f}%` to
`{100.0*diagnostics['with_truth_fraction_outside_bump_radius']:.3f}%`.

This is intentionally the uncorrected residual-drift input requested for the
paired model test.  If the orbit excursion, tail error, or out-of-range
fraction dominates the aligned case, the next physical protocol must perform
and record a BPM-only orbit correction relative to the yearly nominal orbit
before the sextupole scans; a neural model must not be credited with replacing
that machine operation.

## Model definitions

- `physics_gls`: fixed latest-lattice covariance-uniform full-BPM source
  inverse, applied independently to each target block.
- `shared_target_local_ridge`: one parameter-sharing residual model for all
  targets, using only the target scan, target identity, and common baseline
  orbit modes.
- `shared_joint_ridge`: the same residual model plus a compact context derived
  from all {target_count} target scans in the machine; one inference call
  returns all `{2*target_count}` center coordinates.
- `shared_joint_random_feature`: the joint inputs plus a fixed nonlinear tanh
  feature layer, with only the output ridge weights fitted.

The response residual, baseline-orbit, and all-target context projections are
fit on training machines only.  Model and PCA files are saved with the result
so validation/test information cannot silently enter feature construction.

## Scope and limitations

This is a synthetic SciBmad pilot, not a CESR position-precision claim.  The
sample contains only {machine_count} independent static machines, the error
priors other than the user-provided quadrupole drift are maintained sensitivity
settings rather than measured CESR distributions, BPM white noise is
independent, and the drift is one target-local scalar mode with an independent
balanced-cycle realization for each target scan rather than one continuous
trajectory spanning all 76 scans.  There is no
actuator hysteresis, K2 polarity asymmetry, missing/outlier BPM process, or
sim-to-real validation.  The learned models return point estimates rather than
a calibrated posterior covariance or OOD probability.  Exact target orbit is
used only to form evaluation truth.

The latest lattice emits the straight-multipole-in-curved-reference warning.
No girder pitch is varied in this experiment, so the documented curved-DQX
girder-pitch discrepancy is not an excitation here, but remains part of the
lattice provenance.
"""
    (output / "SUMMARY.md").write_text(report, encoding="utf-8")
    metadata = {
        "format": "cesr-sequential-joint-inverse-analysis-v1",
        "generated_seconds": time.perf_counter() - started,
        "scan_root": str(scan_root),
        "machine_count": machine_count,
        "target_count": target_count,
        "train_indices_zero_based": train_indices.tolist(),
        "validation_indices_zero_based": val_indices.tolist(),
        "test_indices_zero_based": test_indices.tolist(),
        "repeats_per_signed_state": args.repeats,
        "bpm_noise_rms_m_per_read": args.bpm_noise_rms_m,
        "drift_endpoint_rms_m": args.drift_endpoint_rms_m,
        "train_augmentations": args.train_augmentations,
        "validation_augmentations": args.validation_augmentations,
        "test_augmentations": args.test_augmentations,
        "residual_modes": args.residual_modes,
        "reference_modes": args.reference_modes,
        "global_modes": args.global_modes,
        "hidden_units": args.hidden_units,
        "local_feature_count": int(bundles[WITHOUT].local_ridge.mean.size),
        "joint_feature_count": int(bundles[WITHOUT].joint_ridge.mean.size),
        "measurement_seed": args.measurement_seed,
        "model_seed": args.model_seed,
    }
    (output / "analysis_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )
    print(report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

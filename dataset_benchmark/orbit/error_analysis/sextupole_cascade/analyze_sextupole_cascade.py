#!/usr/bin/env python3

"""Analyze vector closure of the CESR global sextupole-strength scan."""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path


def norm(values: list[float]) -> float:
    return math.sqrt(sum(value * value for value in values))


def dot(left: list[float], right: list[float]) -> float:
    return sum(a * b for a, b in zip(left, right))


def subtract(left: list[float], right: list[float]) -> list[float]:
    return [a - b for a, b in zip(left, right)]


def solve_three(matrix: list[list[float]], rhs: list[float]) -> list[float]:
    augmented = [row[:] + [value] for row, value in zip(matrix, rhs)]
    for column in range(3):
        pivot = max(range(column, 3), key=lambda row: abs(augmented[row][column]))
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        scale = augmented[column][column]
        if abs(scale) < 1e-15:
            raise ValueError("Singular quadratic strength-scan design")
        augmented[column] = [value / scale for value in augmented[column]]
        for row in range(3):
            if row == column:
                continue
            factor = augmented[row][column]
            augmented[row] = [
                value - factor * pivot_value
                for value, pivot_value in zip(augmented[row], augmented[column])
            ]
    return [augmented[row][3] for row in range(3)]


def fit_quadratic(lambdas: list[float], values: list[float]) -> list[float]:
    powers = [[1.0, value, value * value] for value in lambdas]
    normal = [
        [sum(row[i] * row[j] for row in powers) for j in range(3)]
        for i in range(3)
    ]
    rhs = [sum(row[i] * value for row, value in zip(powers, values)) for i in range(3)]
    return solve_three(normal, rhs)


def percentile(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    position = probability * (len(ordered) - 1)
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    fraction = position - lower
    return ordered[lower] * (1 - fraction) + ordered[upper] * fraction


def load_vectors(path: Path) -> tuple[list[float], list[int], dict[tuple[float, int], list[float]]]:
    by_key: dict[tuple[float, int], list[tuple[str, float]]] = defaultdict(list)
    with path.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            by_key[(float(row["lambda2"]), int(row["trial"]))].append(
                (row["detector"], float(row["c3_y_m"]))
            )
    lambdas = sorted({key[0] for key in by_key})
    trials = sorted({key[1] for key in by_key})
    vectors: dict[tuple[float, int], list[float]] = {}
    detector_order: list[str] | None = None
    for key, entries in by_key.items():
        names = [entry[0] for entry in entries]
        if detector_order is None:
            detector_order = names
        elif names != detector_order:
            raise ValueError(f"Detector order differs at {key}")
        vectors[key] = [entry[1] for entry in entries]
    for lambda2 in lambdas:
        for trial in trials:
            if (lambda2, trial) not in vectors:
                raise ValueError(f"Missing vector for lambda={lambda2}, trial={trial}")
    if lambdas[0] != 0.0 or lambdas[-1] != 1.0:
        raise ValueError("Scan must include lambda2=0 and lambda2=1")
    return lambdas, trials, vectors


def analyze(lambdas: list[float], trials: list[int], vectors):
    direction_rows: list[dict[str, float]] = []
    lambda_rows: list[dict[str, float]] = []
    zero_fractions: list[float] = []
    zero_projections: list[float] = []
    fitted_components: dict[int, tuple[list[float], list[float], list[float]]] = {}

    for trial in trials:
        by_lambda = [vectors[(lambda2, trial)] for lambda2 in lambdas]
        components = [[], [], []]
        for detector in range(len(by_lambda[0])):
            fitted = fit_quadratic(
                lambdas, [values[detector] for values in by_lambda]
            )
            for order in range(3):
                components[order].append(fitted[order])
        fitted_components[trial] = (components[0], components[1], components[2])

    for trial in trials:
        zero = vectors[(0.0, trial)]
        baseline = vectors[(1.0, trial)]
        sext = subtract(baseline, zero)
        denominator = norm(baseline)
        if denominator == 0:
            raise ValueError(f"Zero nominal C3 vector for trial {trial}")
        zero_fraction = norm(zero) / denominator
        zero_projection = dot(zero, baseline) / denominator**2
        sext_fraction = norm(sext) / denominator
        sext_projection = dot(sext, baseline) / denominator**2
        constant, linear, quadratic = fitted_components[trial]
        zero_fractions.append(zero_fraction)
        zero_projections.append(zero_projection)
        direction_rows.append(
            {
                "trial": trial,
                "zero_norm_fraction": zero_fraction,
                "zero_signed_projection": zero_projection,
                "sextupole_difference_norm_fraction": sext_fraction,
                "sextupole_difference_signed_projection": sext_projection,
                "lambda_constant_norm_fraction": norm(constant) / denominator,
                "lambda_constant_signed_projection": dot(constant, baseline) / denominator**2,
                "lambda_linear_norm_fraction": norm(linear) / denominator,
                "lambda_linear_signed_projection": dot(linear, baseline) / denominator**2,
                "lambda_quadratic_norm_fraction": norm(quadratic) / denominator,
                "lambda_quadratic_signed_projection": dot(quadratic, baseline) / denominator**2,
            }
        )

    for lambda2 in lambdas:
        measured_all: list[float] = []
        baseline_all: list[float] = []
        pure_residual_all: list[float] = []
        anchored_residual_all: list[float] = []
        full_residual_all: list[float] = []
        pure_direction: list[float] = []
        anchored_direction: list[float] = []
        full_direction: list[float] = []
        for trial in trials:
            measured = vectors[(lambda2, trial)]
            zero = vectors[(0.0, trial)]
            baseline = vectors[(1.0, trial)]
            scale = lambda2**2
            pure = [scale * value for value in baseline]
            anchored = [z + scale * (b - z) for z, b in zip(zero, baseline)]
            constant, linear, quadratic = fitted_components[trial]
            full = [
                c0 + lambda2 * c1 + scale * c2
                for c0, c1, c2 in zip(constant, linear, quadratic)
            ]
            baseline_norm = norm(baseline)
            pure_residual = subtract(measured, pure)
            anchored_residual = subtract(measured, anchored)
            full_residual = subtract(measured, full)
            pure_direction.append(norm(pure_residual) / baseline_norm)
            anchored_direction.append(norm(anchored_residual) / baseline_norm)
            full_direction.append(norm(full_residual) / baseline_norm)
            measured_all.extend(measured)
            baseline_all.extend(baseline)
            pure_residual_all.extend(pure_residual)
            anchored_residual_all.extend(anchored_residual)
            full_residual_all.extend(full_residual)
        baseline_norm_all = norm(baseline_all)
        lambda_rows.append(
            {
                "lambda2": lambda2,
                "lambda2_squared": lambda2**2,
                "global_c3_norm_ratio": norm(measured_all) / baseline_norm_all,
                "pure_lambda2_squared_relative_residual": norm(pure_residual_all) / baseline_norm_all,
                "anchored_constant_plus_lambda2_squared_relative_residual": norm(anchored_residual_all) / baseline_norm_all,
                "full_quadratic_lambda_relative_residual": norm(full_residual_all) / baseline_norm_all,
                "anchored_direction_p10": percentile(anchored_direction, 0.10),
                "anchored_direction_median": percentile(anchored_direction, 0.50),
                "anchored_direction_p90": percentile(anchored_direction, 0.90),
                "full_quadratic_direction_p10": percentile(full_direction, 0.10),
                "full_quadratic_direction_median": percentile(full_direction, 0.50),
                "full_quadratic_direction_p90": percentile(full_direction, 0.90),
            }
        )
    baseline_all = [value for trial in trials for value in vectors[(1.0, trial)]]
    baseline_norm = norm(baseline_all)
    global_components = {}
    for order, label in enumerate(("constant", "linear", "quadratic")):
        component = [
            value
            for trial in trials
            for value in fitted_components[trial][order]
        ]
        global_components[label] = {
            "norm_fraction": norm(component) / baseline_norm,
            "signed_projection": dot(component, baseline_all) / baseline_norm**2,
        }
    return direction_rows, lambda_rows, zero_fractions, zero_projections, global_components


def write_csv(path: Path, rows: list[dict[str, float]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def render_svg(path: Path, rows: list[dict[str, float]], trials: int) -> None:
    width, height = 900, 560
    left, right, top, bottom = 90, 830, 90, 470
    plot_width, plot_height = right - left, bottom - top
    y_max = max(1.05, 1.1 * max(row["global_c3_norm_ratio"] for row in rows))

    def xp(value: float) -> float:
        return left + plot_width * value

    def yp(value: float) -> float:
        return bottom - plot_height * value / y_max

    measured = " ".join(
        f'{xp(row["lambda2"]):.2f},{yp(row["global_c3_norm_ratio"]):.2f}' for row in rows
    )
    pure = " ".join(
        f'{xp(row["lambda2"]):.2f},{yp(row["lambda2_squared"]):.2f}' for row in rows
    )
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<g font-family="Arial, Helvetica, sans-serif" fill="#202020">',
        '<text x="450" y="31" text-anchor="middle" font-size="21" font-weight="600">Vertical cubic response under global sextupole scaling</text>',
        f'<text x="450" y="55" text-anchor="middle" font-size="13" fill="#555">Concatenated detector vectors across {trials} fixed directions</text>',
    ]
    for index in range(6):
        value = y_max * index / 5
        y = yp(value)
        parts.append(f'<line x1="{left}" y1="{y:.2f}" x2="{right}" y2="{y:.2f}" stroke="#E2E5E7"/>')
        parts.append(f'<text x="{left-10}" y="{y+4:.2f}" text-anchor="end" font-size="11">{value:.2f}</text>')
    for index in range(5):
        value = index / 4
        x = xp(value)
        parts.append(f'<line x1="{x:.2f}" y1="{top}" x2="{x:.2f}" y2="{bottom}" stroke="#ECEEEF"/>')
        parts.append(f'<text x="{x:.2f}" y="{bottom+21}" text-anchor="middle" font-size="11">{value:g}</text>')
    parts.extend(
        [
            f'<polyline points="{pure}" fill="none" stroke="#777" stroke-width="2" stroke-dasharray="7 5"/>',
            f'<polyline points="{measured}" fill="none" stroke="#0072B2" stroke-width="2.8"/>',
        ]
    )
    for row in rows:
        parts.append(
            f'<circle cx="{xp(row["lambda2"]):.2f}" cy="{yp(row["global_c3_norm_ratio"]):.2f}" r="4" fill="white" stroke="#0072B2" stroke-width="2"/>'
        )
    parts.extend(
        [
            f'<rect x="{left}" y="{top}" width="{plot_width}" height="{plot_height}" fill="none" stroke="#555"/>',
            '<line x1="250" y1="76" x2="286" y2="76" stroke="#0072B2" stroke-width="2.8"/><text x="296" y="81" font-size="12">Measured C3 norm</text>',
            '<line x1="490" y1="76" x2="526" y2="76" stroke="#777" stroke-width="2" stroke-dasharray="7 5"/><text x="536" y="81" font-size="12">Pure lambda2 squared</text>',
            f'<text x="{(left+right)/2}" y="{bottom+52}" text-anchor="middle" font-size="13">Global sextupole multiplier, lambda2</text>',
            '<text x="24" y="280" text-anchor="middle" font-size="13" transform="rotate(-90 24 280)">C3 vector norm / nominal norm</text>',
            '</g></svg>',
        ]
    )
    path.write_text("\n".join(parts), encoding="utf-8")


def write_report(path: Path, rows, zero_fractions, zero_projections, components, trials: int) -> None:
    interior = [row for row in rows if 0.0 < row["lambda2"] < 1.0]
    worst_anchored = max(
        (row["anchored_constant_plus_lambda2_squared_relative_residual"] for row in interior),
        default=0.0,
    )
    worst_full = max(row["full_quadratic_lambda_relative_residual"] for row in rows)
    lines = [
        "# Sextupole-cascade strength-scan result",
        "",
        f"The scan uses {trials} fixed vertical-corrector directions. Every lattice variant has its own nominal closed orbit and linear detector response.",
        "",
        "## Global result",
        "",
        f"- With all order-2 multipoles removed, the C3 norm fraction has median {percentile(zero_fractions, 0.5):.4f} and P10--P90 [{percentile(zero_fractions, 0.1):.4f}, {percentile(zero_fractions, 0.9):.4f}] across directions.",
        f"- Its signed projection onto the nominal C3 vector has median {percentile(zero_projections, 0.5):.4f} and P10--P90 [{percentile(zero_projections, 0.1):.4f}, {percentile(zero_projections, 0.9):.4f}].",
        f"- The largest global vector residual of the anchored model C3(lambda)=C3(0)+lambda^2[C3(1)-C3(0)] over interior scan points is {worst_anchored:.4%} of the nominal norm.",
        f"- A full vector fit C3(lambda)=A0+A1*lambda+A2*lambda^2 has maximum scan-point residual {worst_full:.4%} of the nominal norm.",
        f"- The fitted component norm fractions (signed projections onto nominal in parentheses) are: A0 {components['constant']['norm_fraction']:.4f} ({components['constant']['signed_projection']:.4f}), A1 {components['linear']['norm_fraction']:.4f} ({components['linear']['signed_projection']:.4f}), and A2 {components['quadratic']['norm_fraction']:.4f} ({components['quadratic']['signed_projection']:.4f}).",
        "",
        "A material A2 component supports a two-sextupole cascade. A material A1 component instead indicates one sextupole interaction combined with a fixed nonlinear source or strength-dependent feed-down. A nonzero A0 proves that sextupoles are not the only source. The next discriminating control is the nonlinear wiggler model.",
        "",
        "## Lambda scan",
        "",
        "| lambda2 | measured norm / nominal | pure lambda2^2 residual | constant + lambda2^2 residual | full quadratic residual | full-fit direction median [P10, P90] |",
        "|---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            f'| {row["lambda2"]:.2f} | {row["global_c3_norm_ratio"]:.6f} | '
            f'{row["pure_lambda2_squared_relative_residual"]:.6f} | '
            f'{row["anchored_constant_plus_lambda2_squared_relative_residual"]:.6f} | '
            f'{row["full_quadratic_lambda_relative_residual"]:.6f} | '
            f'{row["full_quadratic_direction_median"]:.6f} [{row["full_quadratic_direction_p10"]:.6f}, {row["full_quadratic_direction_p90"]:.6f}] |'
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("vectors", type=Path)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    output_dir = args.output_dir or args.vectors.parent
    output_dir.mkdir(parents=True, exist_ok=True)
    lambdas, trials, vectors = load_vectors(args.vectors)
    direction_rows, lambda_rows, zero_fractions, zero_projections, components = analyze(
        lambdas, trials, vectors
    )
    write_csv(output_dir / "sextupole_cascade_direction_attribution.csv", direction_rows)
    write_csv(output_dir / "sextupole_cascade_lambda_closure.csv", lambda_rows)
    render_svg(output_dir / "sextupole_cascade_strength_scan.svg", lambda_rows, len(trials))
    write_report(
        output_dir / "SEXTUPOLE_CASCADE_RESULTS.md",
        lambda_rows,
        zero_fractions,
        zero_projections,
        components,
        len(trials),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

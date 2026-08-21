#!/usr/bin/env python3
"""Inspect Bmad's standard and field-tracking maps for ID_S1A."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from pytao import Tao


HERE = Path(__file__).resolve().parent
LATTICE_DIR = HERE.parent
LATTICE = LATTICE_DIR / "lat.bmad"
SLAVES = (997, 999, 1001)


def map_data(tao: Tao, index: int) -> tuple[np.ndarray, np.ndarray]:
    result = tao.ele(
        f"0>>{index}", which="model", defaults=False, mat6=True, warn=False
    ).mat6
    return np.asarray(result.mat6, dtype=float), np.asarray(result.vec0, dtype=float)


def block_map(tao: Tao) -> tuple[np.ndarray, np.ndarray]:
    matrix = np.eye(6)
    vector = np.zeros(6)
    for index in range(997, 1002):
        local_matrix, local_vector = map_data(tao, index)
        vector = local_matrix @ vector + local_vector
        matrix = local_matrix @ matrix
    return matrix, vector


def show_commands(tao: Tao, commands: tuple[str, ...]) -> None:
    for command in commands:
        print(f"command={command}")
        try:
            output = tao.cmd(command)
            for line in output:
                print(line)
        except Exception as exc:  # diagnostic: report every attempted method
            print(f"error={type(exc).__name__}: {exc}")


def main() -> None:
    tao = Tao(
        lattice_file=str(LATTICE),
        noplot=True,
        noinit=True,
        nostartup=True,
    )
    show_commands(
        tao,
        (
            "show element 0>>997",
            "show element ID_S1A",
        ),
    )
    standard_matrix, standard_vector = block_map(tao)
    print(f"standard_R12={standard_matrix[0, 1]:.16e}")
    print(f"standard_vec0={standard_vector.tolist()}")

    for tracking_method in ("runge_kutta", "time_runge_kutta", "symp_lie_ptc"):
        candidate = Tao(
            lattice_file=str(LATTICE),
            noplot=True,
            noinit=True,
            nostartup=True,
        )
        commands = []
        for index in SLAVES:
            commands.extend(
                (
                    f"set element 0>>{index} tracking_method = {tracking_method}",
                    f"set element 0>>{index} mat6_calc_method = tracking",
                )
            )
        show_commands(candidate, tuple(commands))
        try:
            matrix, vector = block_map(candidate)
            print(f"method={tracking_method}")
            print(f"R12={matrix[0, 1]:.16e}")
            print(f"max_abs_dR_from_standard={np.max(np.abs(matrix - standard_matrix)):.12e}")
            print(f"vec0={vector.tolist()}")
        except Exception as exc:
            print(f"method_error={tracking_method}: {type(exc).__name__}: {exc}")


if __name__ == "__main__":
    main()

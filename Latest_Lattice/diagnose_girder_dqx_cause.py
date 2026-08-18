#!/usr/bin/env python3
"""Isolate the DQX quadrupole contribution to one girder-pitch response."""

from __future__ import annotations

import json
from pathlib import Path

from pytao import Tao


HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "bmad_girder_reference" / "dqx4b_pitch_isolation.json"
STEP = 1.0e-6
DQX = "0>>22"
GIRDER = "0>>1209"


def orbit_at(tao: Tao, index: int) -> list[float]:
    return [
        float(
            tao.lat_list(
                "*", f"orbit.vec.{coordinate}", ix_branch=0,
                which="model", flags="-array_out -track_only",
            )[index]
        )
        for coordinate in range(1, 7)
    ]


def response(tao: Tao) -> list[float]:
    tao.cmd(f"set element {GIRDER} y_pitch = {STEP:.17g}")
    plus = orbit_at(tao, 22)
    tao.cmd(f"set element {GIRDER} y_pitch = {-STEP:.17g}")
    minus = orbit_at(tao, 22)
    tao.cmd(f"set element {GIRDER} y_pitch = 0")
    return [(p - m) / (2 * STEP) for p, m in zip(plus, minus)]


def main() -> None:
    tao = Tao(
        lattice_file=str(HERE / "lat.bmad"),
        noplot=True,
        noinit=True,
        nostartup=True,
    )
    tao.cmd("set branch 0 geometry = open")
    # DQX's quadrupole component is stored in the bend multipole coefficient
    # B1, not in the generic K1 attribute. A fresh Tao process makes restoring
    # the original B1 unnecessary after this diagnostic.
    nominal = response(tao)
    tao.cmd(f"set element {DQX} tracking_method = runge_kutta")
    runge_kutta = response(tao)
    tao.cmd(f"set element {DQX} tracking_method = bmad_standard")
    tao.cmd(f"set element {DQX} b1 = 0")
    no_k1 = response(tao)
    OUTPUT.write_text(
        json.dumps(
            {
                "step": STEP,
                "nominal": nominal,
                "runge_kutta": runge_kutta,
                "no_k1": no_k1,
            },
            indent=2,
        ) + "\n",
        encoding="utf-8",
    )
    print(OUTPUT)


if __name__ == "__main__":
    main()

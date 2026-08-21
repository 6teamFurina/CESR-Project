#!/usr/bin/env python3
"""Export one-pass Bmad tracking derivatives for all girder parameters."""

from __future__ import annotations

import csv
import math
from collections import defaultdict
from pathlib import Path

from pytao import Tao


HERE = Path(__file__).resolve().parent
LATTICE_DIR = HERE.parent
REFERENCE_DIR = LATTICE_DIR / "bmad_reference"
INVENTORY_DIR = REFERENCE_DIR / "inventory"
LATTICE = LATTICE_DIR / "lat.bmad"
LORDS_CSV = INVENTORY_DIR / "bmad_control_lords.csv"
MEMBERS_CSV = REFERENCE_DIR / "girder" / "members.csv"
OUTPUT = REFERENCE_DIR / "girder" / "tracking_response.csv"
PARTICLE_START = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
PARAMETERS = ("X_OFFSET", "Y_OFFSET", "Z_OFFSET", "X_PITCH", "Y_PITCH", "TILT")


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def tracking_arrays(tao: Tao) -> list[list[float]]:
    result = [
        [
            float(value)
            for value in tao.lat_list(
                "*", f"orbit.vec.{coordinate}", ix_branch=0,
                which="model", flags="-array_out -track_only",
            )
        ]
        for coordinate in range(1, 7)
    ]
    if {len(values) for values in result} != {1178}:
        raise RuntimeError("Unexpected branch-0 tracking-array length")
    return result


def main() -> None:
    girders = {
        row["name"]: row["element_id"]
        for row in read_csv(LORDS_CSV)
        if row["key"] == "Girder"
    }
    members: dict[str, list[int]] = defaultdict(list)
    for row in read_csv(MEMBERS_CSV):
        members[row["girder"]].append(int(row["member_index"]))

    tao = Tao(
        lattice_file=str(LATTICE), noplot=True, noinit=True, nostartup=True,
    )
    tao.cmd("set branch 0 geometry = open")
    names = ("x", "px", "y", "py", "z", "pz")
    tao.cmds(
        [
            f"set particle_start {name} = {value:.17g}"
            for name, value in zip(names, PARTICLE_START)
        ],
        suppress_lattice_calc=True,
        suppress_plotting=True,
    )

    rows: list[dict[str, object]] = []
    calculation = 0
    for girder, element_id in girders.items():
        observations = sorted(set(members[girder]) | {1177})
        for parameter in PARAMETERS:
            calculation += 1
            step = 1.0e-6 if parameter in {"X_PITCH", "Y_PITCH", "TILT"} else 1.0e-5
            tao.cmd(f"set element {element_id} {parameter} = {step:.17g}")
            plus = tracking_arrays(tao)
            tao.cmd(f"set element {element_id} {parameter} = {-step:.17g}")
            minus = tracking_arrays(tao)
            tao.cmd(f"set element {element_id} {parameter} = 0")
            for index in observations:
                row: dict[str, object] = {
                    "girder": girder,
                    "parameter": parameter.lower(),
                    "step": step,
                    "observation_index": index,
                }
                for coordinate, name in enumerate(names):
                    derivative = (plus[coordinate][index] - minus[coordinate][index]) / (2 * step)
                    if not math.isfinite(derivative):
                        raise RuntimeError(f"Non-finite {girder} {parameter} response")
                    row[f"d{name}"] = derivative
                rows.append(row)
            print(f"Bmad girder tracking {calculation}/72 {girder} {parameter}")

    fields = [
        "girder", "parameter", "step", "observation_index",
        "dx", "dpx", "dy", "dpy", "dz", "dpz",
    ]
    with OUTPUT.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {OUTPUT}")
    print(f"Response observations: {len(rows)}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Finite-difference every exported Bmad control at its affected slave exits."""

from __future__ import annotations

import csv
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Any

from pytao import Tao


HERE = Path(__file__).resolve().parent
LATTICE = HERE / "lat.bmad"
LORDS_CSV = HERE / "bmad_control_lords.csv"
RELATIONS_CSV = HERE / "bmad_control_relations.csv"
INVENTORY_CSV = HERE / "bmad_element_inventory.csv"
OUTPUT_DIR = HERE / "bmad_control_tracking_reference"

PARTICLE_START = (1.0e-3, 2.0e-4, -8.0e-4, 1.5e-4, 5.0e-4, 2.0e-4)
STEP_BY_VARIABLE = {
    "HKICK": 1.0e-6,
    "VKICK": 1.0e-6,
    "K1": 1.0e-5,
    "A1": 1.0e-5,
    "COMMAND": 1.0e-3,
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def tracking_arrays(tao: Tao) -> list[list[float]]:
    arrays: list[list[float]] = []
    for coordinate in range(1, 7):
        values = tao.lat_list(
            "*",
            f"orbit.vec.{coordinate}",
            ix_branch=0,
            which="model",
            flags="-array_out -track_only",
        )
        arrays.append([float(value) for value in values])
    lengths = {len(values) for values in arrays}
    if lengths != {1178}:
        raise RuntimeError(f"Expected six 1178-position arrays, found {lengths}")
    return arrays


def set_control(tao: Tao, element_id: str, variable: str, value: float) -> None:
    tao.cmd(f"set element {element_id} {variable} = {value:.17g}")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    lords = [
        row for row in read_csv(LORDS_CSV) if row["key"] in {"Overlay", "Group"}
    ]
    relations = read_csv(RELATIONS_CSV)
    inventory = read_csv(INVENTORY_CSV)

    tracking_by_base_name: dict[str, list[int]] = defaultdict(list)
    for element in inventory:
        if element["branch"] != "0" or element["is_tracking"] != "True":
            continue
        name = element["name"].upper()
        base_name = name.split("#", 1)[0]
        tracking_by_base_name[base_name].append(int(element["index"]))

    relations_by_lord: dict[str, list[dict[str, str]]] = defaultdict(list)
    for relation in relations:
        relations_by_lord[relation["lord_id"]].append(relation)

    tao = Tao(
        lattice_file=str(LATTICE),
        noplot=True,
        noinit=True,
        nostartup=True,
    )
    # Compare a deterministic one-pass map rather than the periodic closed
    # orbit. Tao otherwise rejects non-pz particle_start coordinates on a ring.
    tao.cmd("set branch 0 geometry = open")
    coordinate_names = ("x", "px", "y", "py", "z", "pz")
    tao.cmds(
        [
            f"set particle_start {name} = {value:.17g}"
            for name, value in zip(coordinate_names, PARTICLE_START)
        ],
        suppress_lattice_calc=True,
        suppress_plotting=True,
    )
    baseline = tracking_arrays(tao)

    output_rows: list[dict[str, Any]] = []
    control_rows: list[dict[str, Any]] = []
    for control_index, lord in enumerate(lords, 1):
        variables = lord["control_variables"].split(";")
        if len(variables) != 1 or "=" not in variables[0]:
            raise RuntimeError(f"Unexpected control variables for {lord}: {variables}")
        variable, baseline_text = variables[0].split("=", 1)
        baseline_value = float(baseline_text)
        if baseline_value != 0.0:
            raise RuntimeError(f"Expected zero control baseline for {lord['name']}")
        step = STEP_BY_VARIABLE[variable]

        observation_indices: set[int] = {1177}
        lord_relations = relations_by_lord[lord["element_id"]]
        for relation in lord_relations:
            slave_index = int(relation["slave_id"].split(">>", 1)[1])
            if slave_index <= 1177:
                observation_indices.add(slave_index)
            else:
                slices = tracking_by_base_name[relation["slave_name"].upper()]
                if not slices:
                    raise RuntimeError(f"No tracking slices for {relation['slave_name']}")
                observation_indices.update(slices)

        set_control(tao, lord["element_id"], variable, step)
        plus = tracking_arrays(tao)
        set_control(tao, lord["element_id"], variable, -step)
        minus = tracking_arrays(tao)
        set_control(tao, lord["element_id"], variable, 0.0)

        for index in sorted(observation_indices):
            row: dict[str, Any] = {
                "lord_id": lord["element_id"],
                "lord_name": lord["name"],
                "lord_key": lord["key"],
                "variable": variable,
                "step": step,
                "observation_index": index,
            }
            for coordinate, name in enumerate(coordinate_names):
                derivative = (plus[coordinate][index] - minus[coordinate][index]) / (2 * step)
                if not math.isfinite(derivative):
                    raise RuntimeError(f"Non-finite response for {lord['name']} at {index}")
                row[f"d{name}"] = derivative
            output_rows.append(row)

        control_rows.append(
            {
                "lord_id": lord["element_id"],
                "lord_name": lord["name"],
                "lord_key": lord["key"],
                "variable": variable,
                "step": step,
                "relation_count": len(lord_relations),
                "observation_count": len(observation_indices),
            }
        )
        print(
            f"Bmad control tracking {control_index}/{len(lords)} "
            f"{lord['name']} ({len(lord_relations)} relations, "
            f"{len(observation_indices)} observations)"
        )

    response_path = OUTPUT_DIR / "control_tracking_response.csv"
    fields = [
        "lord_id", "lord_name", "lord_key", "variable", "step",
        "observation_index", "dx", "dpx", "dy", "dpy", "dz", "dpz",
    ]
    with response_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(output_rows)

    control_path = OUTPUT_DIR / "controls.csv"
    with control_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(control_rows[0]))
        writer.writeheader()
        writer.writerows(control_rows)

    metadata = {
        "format": "latest-cesr-bmad-control-tracking-response-v1",
        "tao_version": str(tao.version()),
        "lattice": str(LATTICE),
        "branch": 0,
        "tracking_element_count": 1177,
        "control_count": len(lords),
        "relationship_count": len(relations),
        "response_row_count": len(output_rows),
        "particle_start": dict(zip(coordinate_names, PARTICLE_START)),
        "steps": STEP_BY_VARIABLE,
        "baseline_end_orbit": [values[-1] for values in baseline],
    }
    metadata_path = OUTPUT_DIR / "metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {response_path}")
    print(f"Wrote {control_path}")
    print(f"Wrote {metadata_path}")


if __name__ == "__main__":
    main()

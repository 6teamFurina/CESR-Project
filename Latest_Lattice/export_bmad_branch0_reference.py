#!/usr/bin/env python3
"""Export an unambiguous Tao/Bmad reference for CESR branch 0 only."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any

from pytao import Tao


HERE = Path(__file__).resolve().parent
LATTICE = HERE / "lat.bmad"
OUTPUT_DIR = HERE / "bmad_reference_branch0"


def jsonable(value: Any) -> Any:
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError(f"Non-finite value in Bmad reference: {value}")
        return value
    if hasattr(value, "tolist"):
        return jsonable(value.tolist())
    if hasattr(value, "item"):
        return jsonable(value.item())
    if isinstance(value, dict):
        return {str(key): jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(item) for item in value]
    if hasattr(value, "model_dump"):
        return jsonable(value.model_dump())
    return str(value)


def query(tao: Tao, field: str, count: int) -> list[Any]:
    result = list(
        jsonable(
            tao.lat_list(
                "*",
                field,
                ix_branch=0,
                which="model",
                flags="-array_out -track_only",
            )
        )
    )
    if len(result) != count:
        raise RuntimeError(f"{field} returned {len(result)} values, expected {count}")
    return result


def normalize_map(result: Any) -> dict[str, Any]:
    if result is None:
        raise RuntimeError("Tao returned no map")
    if isinstance(result, dict):
        matrix = result.get("mat6", result.get("matrix"))
        vector = result.get("vec0", result.get("vector"))
        error = result.get("symplectic_error")
    else:
        matrix = getattr(result, "mat6", None)
        vector = getattr(result, "vec0", None)
        error = getattr(result, "symplectic_error", None)
    if matrix is None or vector is None:
        raise RuntimeError(f"Incomplete Tao map result: {result!r}")
    return {
        "mat6": jsonable(matrix),
        "vec0": jsonable(vector),
        "symplectic_error": jsonable(error),
    }


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    tao = Tao(
        lattice_file=str(LATTICE),
        noplot=True,
        noinit=True,
        nostartup=True,
    )

    indices = [
        int(value)
        for value in query(
            tao,
            "ele.ix_ele",
            len(
                tao.lat_list(
                    "*",
                    "ele.ix_ele",
                    ix_branch=0,
                    which="model",
                    flags="-array_out -track_only",
                )
            ),
        )
    ]
    if indices != list(range(indices[-1] + 1)):
        raise RuntimeError("Branch-0 tracking indices are not contiguous from BEGINNING")
    count = len(indices)

    fields = {
        "name": "ele.name",
        "key": "ele.key",
        "s": "ele.s",
        "length": "ele.l",
        "p0c": "ele.p0c",
        "x": "orbit.vec.1",
        "px": "orbit.vec.2",
        "y": "orbit.vec.3",
        "py": "orbit.vec.4",
        "z": "orbit.vec.5",
        "pz": "orbit.vec.6",
        "beta_a": "ele.a.beta",
        "alpha_a": "ele.a.alpha",
        "phi_a": "ele.a.phi",
        "beta_b": "ele.b.beta",
        "alpha_b": "ele.b.alpha",
        "phi_b": "ele.b.phi",
        "eta_x": "ele.x.eta",
        "etap_x": "ele.x.etap",
        "eta_y": "ele.y.eta",
        "etap_y": "ele.y.etap",
    }
    columns = {name: query(tao, field, count) for name, field in fields.items()}

    elements: list[dict[str, Any]] = []
    identity_map = {
        "mat6": [[1.0 if row == column else 0.0 for column in range(6)] for row in range(6)],
        "vec0": [0.0] * 6,
        "symplectic_error": 0.0,
    }
    for position, index in enumerate(indices):
        element_id = f"0>>{index}"
        element = {"index": index, "element_id": element_id}
        for name, values in columns.items():
            element[name] = values[position]

        if index == 0:
            element["local_map"] = identity_map
            element["cumulative_map"] = identity_map
        else:
            local = tao.ele(
                element_id,
                which="model",
                defaults=False,
                mat6=True,
                warn=False,
            )
            element["local_map"] = normalize_map(local.mat6)
            element["cumulative_map"] = normalize_map(
                tao.matrix("0>>0", element_id)
            )
        elements.append(element)
        if index == 0 or index % 100 == 0 or index == indices[-1]:
            print(f"Exported branch 0 element {index}/{indices[-1]}")

    one_turn = normalize_map(tao.matrix("0>>0", f"0>>{indices[-1]}"))
    payload = {
        "format": "cesr-bmad-branch0-reference-v1",
        "lattice": str(LATTICE),
        "tao_version": jsonable(tao.version()),
        "branch": jsonable(tao.lat_branch_list()[0]),
        "tracking_position_count_including_beginning": count,
        "tracking_element_count": count - 1,
        "coordinate_order": ["x", "px", "y", "py", "z", "pz"],
        "ring_general": jsonable(tao.ring_general()),
        "one_turn_map": one_turn,
        "elements": elements,
        "errors": [],
    }

    json_path = OUTPUT_DIR / "bmad_reference_branch0.json"
    json_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    csv_fields = ["index", "element_id", *fields]
    with (OUTPUT_DIR / "element_index.csv").open(
        "w", encoding="utf-8", newline=""
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=csv_fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(elements)

    map_fields = [
        "index",
        "name",
        "key",
        "length",
        "x",
        "px",
        "y",
        "py",
        "z",
        "pz",
        *[f"r{row}{column}" for row in range(1, 7) for column in range(1, 7)],
        *[f"v{row}" for row in range(1, 7)],
    ]
    with (OUTPUT_DIR / "local_maps.csv").open(
        "w", encoding="utf-8", newline=""
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=map_fields)
        writer.writeheader()
        for element in elements:
            row = {field: element[field] for field in map_fields[:10]}
            matrix = element["local_map"]["mat6"]
            vector = element["local_map"]["vec0"]
            for i in range(6):
                for j in range(6):
                    row[f"r{i + 1}{j + 1}"] = matrix[i][j]
                row[f"v{i + 1}"] = vector[i]
            writer.writerow(row)

    print(f"Wrote {json_path}")
    print(f"Tao version: {payload['tao_version']}")
    print(f"Tracking elements: {payload['tracking_element_count']}")
    print(f"Errors: {len(payload['errors'])}")


if __name__ == "__main__":
    main()

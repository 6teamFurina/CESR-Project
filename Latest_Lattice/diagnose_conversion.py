#!/usr/bin/env python3
"""Inventory the latest CESR Bmad lattice and its raw SciBmad export.

Run this script with the PyTao interpreter from the Ubuntu-Bmad environment.
It intentionally performs no checksum or other cryptographic validation.
"""

from __future__ import annotations

import csv
import json
import re
from collections import Counter
from pathlib import Path
from typing import Any

from pytao import Tao


HERE = Path(__file__).resolve().parent
LATTICE = HERE / "lat.bmad"
RAW_EXPORT = HERE / "latest_cesr_scibmad_bmad_20260814.jl"
WRITE_LOG = HERE / "tao_write_scibmad_bmad_20260814.log"


def scalar(value: Any) -> Any:
    if hasattr(value, "item"):
        return value.item()
    return value


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def make_tao() -> Tao:
    try:
        return Tao(
            lattice_file=str(LATTICE),
            noplot=True,
            noinit=True,
            nostartup=True,
        )
    except TypeError:
        return Tao(f'-noinit -noplot -nostartup -lat "{LATTICE}"')


def collect_bmad_inventory(tao: Tao) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    elements: list[dict[str, Any]] = []
    branches: list[dict[str, Any]] = []

    for branch_info in tao.lat_branch_list():
        info = dict(branch_info)
        branch = int(info["index"])
        branches.append(
            {
                "branch": branch,
                "branch_name": info["branch_name"],
                "n_ele_track": int(info["n_ele_track"]),
                "n_ele_max": int(info["n_ele_max"]),
            }
        )

        fields = {
            "index": "ele.ix_ele",
            "name": "ele.name",
            "key": "ele.key",
            "s": "ele.s",
            "length": "ele.l",
        }
        values = {
            field: list(
                tao.lat_list(
                    "*",
                    query,
                    ix_branch=branch,
                    flags="-array_out",
                )
            )
            for field, query in fields.items()
        }
        counts = {field: len(items) for field, items in values.items()}
        if len(set(counts.values())) != 1:
            raise RuntimeError(f"Inconsistent branch {branch} list lengths: {counts}")

        tracking_indices = {
            int(scalar(value))
            for value in tao.lat_list(
                "*",
                "ele.ix_ele",
                ix_branch=branch,
                flags="-array_out -track_only",
            )
        }

        for position in range(next(iter(counts.values()))):
            index = int(scalar(values["index"][position]))
            elements.append(
                {
                    "branch": branch,
                    "branch_name": info["branch_name"],
                    "index": index,
                    "element_id": f"{branch}>>{index}",
                    "name": str(scalar(values["name"][position])),
                    "key": str(scalar(values["key"][position])),
                    "s": float(scalar(values["s"][position])),
                    "length": float(scalar(values["length"][position])),
                    "is_tracking": index in tracking_indices,
                }
            )

    return branches, elements


SLAVE_RE = re.compile(
    r"^\s*(\S+)\s+\((\d+>>\d+)\)\s+(\S+)\s+"
    r"([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+(.*?)\s*$"
)


def collect_controls(
    tao: Tao, elements: list[dict[str, Any]]
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    lords: list[dict[str, Any]] = []
    relations: list[dict[str, Any]] = []

    for element in elements:
        if element["key"] not in {"Overlay", "Group", "Girder"}:
            continue
        element_id = element["element_id"]
        try:
            variables = tao.ele_control_var(element_id)
        except Exception:
            variables = {}
        lords.append(
            {
                "element_id": element_id,
                "name": element["name"],
                "key": element["key"],
                "control_variables": ";".join(
                    f"{name}={value:.17g}" for name, value in variables.items()
                ),
            }
        )

        for line in tao.show(f"element -all {element_id}"):
            match = SLAVE_RE.match(line)
            if not match:
                continue
            slave, slave_id, attribute, value, expression_value, expression = match.groups()
            relations.append(
                {
                    "lord_id": element_id,
                    "lord_name": element["name"],
                    "lord_key": element["key"],
                    "slave_id": slave_id,
                    "slave_name": slave,
                    "attribute": attribute,
                    "attribute_value": value,
                    "expression_value": expression_value,
                    "expression": expression,
                }
            )

    return lords, relations


def analyze_raw_export() -> dict[str, Any]:
    text = RAW_EXPORT.read_text(encoding="utf-8")
    lines = text.splitlines()
    log_lines = WRITE_LOG.read_text(encoding="utf-8", errors="replace").splitlines()

    constructors: Counter[str] = Counter()
    constructor_rows: list[dict[str, Any]] = []
    assignment_re = re.compile(r"^\s+([A-Za-z0-9_!]+)\s*=\s*([A-Za-z0-9_]+)\(")
    for line_number, line in enumerate(lines, 1):
        match = assignment_re.match(line)
        if match:
            name, constructor = match.groups()
            constructors[constructor] += 1
            constructor_rows.append(
                {"line": line_number, "name": name, "constructor": constructor}
            )

    beamline_lines: dict[str, int] = {}
    for line_number, line in enumerate(lines, 1):
        match = re.match(r"^([A-Za-z0-9_!]+)\s*=\s*Beamline\(", line)
        if match:
            beamline_lines[match.group(1)] = line_number

    forward_forks: list[dict[str, Any]] = []
    for row in constructor_rows:
        if row["constructor"] != "Fork":
            continue
        line = lines[row["line"] - 1]
        match = re.search(r"to_line\s*=\s*([A-Za-z0-9_!]+)", line)
        target = match.group(1) if match else ""
        target_line = beamline_lines.get(target)
        forward_forks.append(
            {
                "fork": row["name"],
                "fork_line": row["line"],
                "target": target,
                "target_line": target_line,
                "is_forward_reference": target_line is None or target_line > row["line"],
            }
        )

    warning_counts = {
        "wiggler_not_translated": sum(
            "Wiggler:" in line and "cannot yet be translated" in line for line in log_lines
        ),
        "girder_not_translated": sum(
            "GIRDER ELEMENTS CANNOT YET BE TRANSLATED" in line for line in log_lines
        ),
    }

    line_elements = [
        row for row in constructor_rows if row["constructor"] == "LineElement"
    ]
    return {
        "constructors": dict(sorted(constructors.items())),
        "forward_forks": forward_forks,
        "line_element_placeholders": line_elements,
        "warning_counts": warning_counts,
        "defexpr_assignment_count": sum("DefExpr(" in line for line in lines),
    }


def main() -> None:
    tao = make_tao()
    branches, elements = collect_bmad_inventory(tao)
    lords, relations = collect_controls(tao, elements)
    raw_export = analyze_raw_export()

    write_csv(
        HERE / "bmad_branches.csv",
        branches,
        ["branch", "branch_name", "n_ele_track", "n_ele_max"],
    )
    write_csv(
        HERE / "bmad_element_inventory.csv",
        elements,
        [
            "branch",
            "branch_name",
            "index",
            "element_id",
            "name",
            "key",
            "s",
            "length",
            "is_tracking",
        ],
    )
    write_csv(
        HERE / "bmad_control_lords.csv",
        lords,
        ["element_id", "name", "key", "control_variables"],
    )
    write_csv(
        HERE / "bmad_control_relations.csv",
        relations,
        [
            "lord_id",
            "lord_name",
            "lord_key",
            "slave_id",
            "slave_name",
            "attribute",
            "attribute_value",
            "expression_value",
            "expression",
        ],
    )

    main_branch = [element for element in elements if element["branch"] == 0]
    summary = {
        "bmad_version": scalar(tao.version()),
        "lattice": str(LATTICE),
        "branch_count": len(branches),
        "all_element_count": len(elements),
        "main_branch_element_count": len(main_branch),
        "main_branch_key_counts": dict(
            sorted(Counter(element["key"] for element in main_branch).items())
        ),
        "control_lord_counts": dict(sorted(Counter(row["key"] for row in lords).items())),
        "control_relation_count": len(relations),
        "raw_export": raw_export,
    }
    (HERE / "conversion_inventory.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

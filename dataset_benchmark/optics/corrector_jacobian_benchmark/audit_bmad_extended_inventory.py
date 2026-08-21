#!/usr/bin/env python3
"""Audit active Bmad quadrupole/sextupole strengths used by the benchmark."""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
OPTICS_DIR = HERE.parent
DATASET_DIR = OPTICS_DIR.parent
ORBIT_DIR = DATASET_DIR / "orbit"
ORBIT_CALCULATION_DIR = ORBIT_DIR / "Orbit_Calculation"
sys.path.insert(0, str(ORBIT_CALCULATION_DIR))

from benchmark_bmad import activate_benchmark_data, response_tools  # noqa: E402

ELEMENT_RE = re.compile(
    r"(?im)^\s*([A-Za-z][A-Za-z0-9_]*)\s*:\s*(Quadrupole|Sextupole)\b"
)


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument(
        "--lattice",
        type=Path,
        default=ORBIT_DIR / "reference" / "cesr_bmad_compatible.bmad",
    )
    result.add_argument("--output-dir", type=Path, default=HERE / "results" / "inventory")
    return result


def dictionary_float(mapping: dict[str, Any], key: str) -> float:
    lower = {str(name).lower(): value for name, value in mapping.items()}
    return float(lower[key.lower()])


def main() -> int:
    args = parser().parse_args()
    lattice = args.lattice.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    text = response_tools.strip_bmad_comments(lattice.read_text(encoding="utf-8"))
    candidates = [(match.group(1).upper(), match.group(2).lower()) for match in ELEMENT_RE.finditer(text)]

    detectors, horizontal, vertical = response_tools.parse_cesr_layout(lattice)
    output_dir.mkdir(parents=True, exist_ok=True)
    init_path = output_dir / "tao_extended_inventory_rf_off.init"
    init_path.write_text(
        response_tools.build_tao_init(lattice, detectors, horizontal, vertical, 1.0e-6),
        encoding="utf-8",
    )
    from pytao import Tao

    tao = Tao(init_file=str(init_path), noplot=True)
    tao.cmds(
        ["set global rf_on = F", "set particle_start pz = 0"],
        suppress_lattice_calc=True,
        suppress_plotting=True,
    )
    activate_benchmark_data(tao)

    inventory: list[dict[str, Any]] = []
    for name, kind in candidates:
        attribute = "K1" if kind == "quadrupole" else "K2"
        attributes = tao.ele_gen_attribs(name)
        strength = dictionary_float(attributes, attribute)
        if strength == 0.0:
            continue
        inventory.append(
            {"name": name, "kind": kind, "attribute": attribute, "strength": strength}
        )

    inventory.sort(key=lambda item: (item["kind"], item["name"]))
    path = output_dir / "bmad_extended_strength_inventory.json"
    path.write_text(json.dumps(inventory, indent=2) + "\n", encoding="utf-8")
    quadrupoles = [item for item in inventory if item["kind"] == "quadrupole"]
    sextupoles = [item for item in inventory if item["kind"] == "sextupole"]
    print(f"Active quadrupole strengths: {len(quadrupoles)}")
    print(f"Active sextupole strengths:  {len(sextupoles)}")
    print(f"Inventory: {path}")
    for item in (quadrupoles[0], quadrupoles[-1], sextupoles[0], sextupoles[-1]):
        print(item)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

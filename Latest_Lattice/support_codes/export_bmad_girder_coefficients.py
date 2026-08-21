#!/usr/bin/env python3
"""Linearize all Bmad girder rigid transforms into SciBmad alignment fields."""

from __future__ import annotations

import csv
import re
from pathlib import Path

from pytao import Tao


HERE = Path(__file__).resolve().parent
LATTICE_DIR = HERE.parent
REFERENCE_DIR = LATTICE_DIR / "bmad_reference"
INVENTORY_DIR = REFERENCE_DIR / "inventory"
LATTICE = LATTICE_DIR / "lat.bmad"
LORDS_CSV = INVENTORY_DIR / "bmad_control_lords.csv"
INVENTORY_CSV = INVENTORY_DIR / "bmad_element_inventory.csv"
OUTPUT_DIR = REFERENCE_DIR / "girder"

GIRDER_PARAMETERS = (
    "X_OFFSET", "Y_OFFSET", "Z_OFFSET", "X_PITCH", "Y_PITCH", "TILT",
)
SCI_PROPERTIES = ("x_offset", "y_offset", "z_offset", "x_rot", "y_rot", "tilt")
OFFSET_STEP = 1.0e-5
ANGLE_STEP = 1.0e-6
SLAVE_RE = re.compile(r"^\s*(.*?)\s+\((0>>\d+)\)\s*([A-Za-z_]+)\s+[-+0-9.Ee]+\s*$")


def read_girders() -> list[dict[str, str]]:
    with LORDS_CSV.open("r", encoding="utf-8", newline="") as stream:
        return [row for row in csv.DictReader(stream) if row["key"] == "Girder"]


def tracking_slices_by_base() -> dict[str, list[dict[str, object]]]:
    result: dict[str, list[dict[str, object]]] = {}
    with INVENTORY_CSV.open("r", encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream):
            if row["branch"] != "0" or row["is_tracking"] != "True":
                continue
            name = row["name"].upper()
            base = name.split("#", 1)[0]
            result.setdefault(base, []).append(
                {
                    "member_id": row["element_id"],
                    "member_index": int(row["index"]),
                    "member_name": row["name"],
                    "member_key": row["key"],
                }
            )
    return result


def girder_slaves(tao: Tao, element_id: str) -> list[dict[str, object]]:
    result: list[dict[str, object]] = []
    in_slaves = False
    for line in tao.show(f"element -all {element_id}"):
        if line.strip() == "Slaves:":
            in_slaves = True
            continue
        if not in_slaves or line.lstrip().startswith("Name"):
            continue
        match = SLAVE_RE.match(line)
        if match:
            name, member_id, key = match.groups()
            result.append(
                {
                    "member_id": member_id,
                    "member_index": int(member_id.split(">>", 1)[1]),
                    "member_name": name.strip(),
                    "member_key": key,
                }
            )
        elif result and not line.strip():
            break
    if not result:
        raise RuntimeError(f"No girder slaves parsed for {element_id}")
    return result


def sci_alignment(attributes: dict[str, object]) -> dict[str, float]:
    # Bmad x/y_pitch describe slopes of the local z axis. Beamlines names the
    # same rotations by their right-handed rotation axes.
    roll = attributes.get("ROLL_TOT", attributes.get("TILT_TOT", 0.0))
    return {
        "x_offset": float(attributes.get("X_OFFSET_TOT", 0.0)),
        "y_offset": float(attributes.get("Y_OFFSET_TOT", 0.0)),
        "z_offset": float(attributes.get("Z_OFFSET_TOT", 0.0)),
        "x_rot": -float(attributes.get("Y_PITCH_TOT", 0.0)),
        "y_rot": float(attributes.get("X_PITCH_TOT", 0.0)),
        "tilt": float(roll),
    }


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    girders = read_girders()
    tao = Tao(
        lattice_file=str(LATTICE),
        noplot=True,
        noinit=True,
        nostartup=True,
    )

    tracking_slices = tracking_slices_by_base()
    members_by_girder: dict[str, list[dict[str, object]]] = {}
    for girder in girders:
        expanded: list[dict[str, object]] = []
        for member in girder_slaves(tao, girder["element_id"]):
            if int(member["member_index"]) <= 1177:
                expanded.append(member)
                continue
            slices = tracking_slices.get(str(member["member_name"]).upper(), [])
            if not slices:
                raise RuntimeError(f"No tracking slices for girder member {member}")
            expanded.extend(slices)
        unique = {int(member["member_index"]): member for member in expanded}
        members_by_girder[girder["name"]] = [unique[index] for index in sorted(unique)]
    coefficient_rows: list[dict[str, object]] = []
    member_rows: list[dict[str, object]] = []
    for girder_index, girder in enumerate(girders, 1):
        name = girder["name"]
        members = members_by_girder[name]
        for member in members:
            member_rows.append({"girder": name, **member})

        for parameter in GIRDER_PARAMETERS:
            step = ANGLE_STEP if parameter in {"X_PITCH", "Y_PITCH", "TILT"} else OFFSET_STEP
            tao.cmd(f"set element {girder['element_id']} {parameter} = {step:.17g}")
            plus = {
                str(member["member_id"]): sci_alignment(
                    tao.ele_gen_attribs(str(member["member_id"]))
                )
                for member in members
            }
            tao.cmd(f"set element {girder['element_id']} {parameter} = {-step:.17g}")
            minus = {
                str(member["member_id"]): sci_alignment(
                    tao.ele_gen_attribs(str(member["member_id"]))
                )
                for member in members
            }
            tao.cmd(f"set element {girder['element_id']} {parameter} = 0")

            for member in members:
                member_id = str(member["member_id"])
                for prop in SCI_PROPERTIES:
                    coefficient_rows.append(
                        {
                            "girder": name,
                            "member_index": member["member_index"],
                            "member_name": member["member_name"],
                            "member_key": member["member_key"],
                            "girder_parameter": parameter.lower(),
                            "scibmad_property": prop,
                            "coefficient": (plus[member_id][prop] - minus[member_id][prop]) / (2 * step),
                        }
                    )
        print(f"Exported girder {girder_index}/{len(girders)} {name} ({len(members)} members)")

    member_path = OUTPUT_DIR / "members.csv"
    with member_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(member_rows[0]))
        writer.writeheader()
        writer.writerows(member_rows)

    coefficient_path = OUTPUT_DIR / "alignment_coefficients.csv"
    with coefficient_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(coefficient_rows[0]))
        writer.writeheader()
        writer.writerows(coefficient_rows)

    print(f"Wrote {member_path}")
    print(f"Wrote {coefficient_path}")
    print(f"Girders: {len(girders)}")
    print(f"Members: {len(member_rows)}")
    print(f"Coefficients: {len(coefficient_rows)}")


if __name__ == "__main__":
    main()

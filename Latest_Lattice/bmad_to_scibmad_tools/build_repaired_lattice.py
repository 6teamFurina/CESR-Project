#!/usr/bin/env python3
"""Build a loadable charged-particle SciBmad lattice from Tao's raw export."""

from __future__ import annotations

import re
from pathlib import Path


HERE = Path(__file__).resolve().parent
LATTICE_DIR = HERE.parent
SOURCE = LATTICE_DIR / "bmad_reference" / "raw_exports" / "latest_cesr_scibmad_bmad_20260814.jl"
OUTPUT = LATTICE_DIR / "latest_cesr_scibmad_repaired.jl"

PHOTON_REGISTRY = r'''

# Photon branch registry. SciBmad/Beamlines currently has no Fork or Mirror
# tracking constructor, so branches are explicit and independently runnable.
const latest_photon_branches = Dict(
    "S4B_LINE" => s4b_line,
    "S7A_LINE" => s7a_line,
    "S7B1_S7B2_LINE" => s7b1_s7b2_line,
    "S1A1_S1A2_S1A3_LINE" => s1a1_s1a2_s1a3_line,
    "S1A1_LINE" => s1a1_line,
    "S1A2_S1A3_LINE" => s1a2_s1a3_line,
    "S2A_LINE" => s2a_line,
    "S2B_LINE" => s2b_line,
    "S3A_LINE" => s3a_line,
    "S3B_LINE" => s3b_line,
    "S4A_LINE" => s4a_line,
)

const latest_photon_fork_targets = Dict(
    "LS_4B" => "S4B_LINE",
    "LS_7A" => "S7A_LINE",
    "LS_7B1_7B2" => "S7B1_S7B2_LINE",
    "LS_1A1_1A2_1A3" => "S1A1_S1A2_S1A3_LINE",
    "LS_2A" => "S2A_LINE",
    "LS_2B" => "S2B_LINE",
    "LS_3A" => "S3A_LINE",
    "LS_3B" => "S3B_LINE",
    "LS_4A" => "S4A_LINE",
    "FORK_S1A1" => "S1A1_LINE",
    "FORK_S1A2_S1A3" => "S1A2_S1A3_LINE",
)

const latest_photon_mirror_placeholders = ("MIRROR_S7", "MIRROR_S1A")

# Both archived Bmad mirrors are flat reference mirrors at 10 keV with no
# reflectivity table or curvature. The reflected reference geometry is already
# encoded in the corresponding Bmad photon branch, so a reference ray is zero
# in each branch's local coordinates on both sides of the mirror.
const latest_photon_mirrors = Dict(
    "MIRROR_S7" => (ref_tilt=-pi / 2, graze_angle=0.004, photon_energy=1.0e4),
    "MIRROR_S1A" => (ref_tilt=-pi / 2, graze_angle=0.004, photon_energy=1.0e4),
)

"""A paraxial photon ray expressed in a Bmad photon branch's local frame."""
struct LatestPhotonRay
    x::Float64
    xprime::Float64
    y::Float64
    yprime::Float64
    path_length::Float64
end

LatestPhotonRay(; x=0.0, xprime=0.0, y=0.0, yprime=0.0, path_length=0.0) =
    LatestPhotonRay(x, xprime, y, yprime, path_length)

function latest_photon_branch(name::AbstractString)
    key = uppercase(String(name))
    haskey(latest_photon_branches, key) || throw(ArgumentError(
        "Unknown photon branch $name; expected one of $(sort!(collect(keys(latest_photon_branches))))",
    ))
    return latest_photon_branches[key]
end

function latest_photon_branch_for_fork(name::AbstractString)
    key = uppercase(String(name))
    haskey(latest_photon_fork_targets, key) || throw(ArgumentError(
        "Unknown photon fork $name; expected one of $(sort!(collect(keys(latest_photon_fork_targets))))",
    ))
    return latest_photon_branch(latest_photon_fork_targets[key])
end

"""
    track_latest_photon_branch(name; ray=LatestPhotonRay())

Propagate a paraxial ray through one of the eleven drift/marker photon lines.
Coordinates are local to the already-reflected Bmad branch reference frame;
therefore the flat reference mirrors are identity operations here. This helper
does not model mirror reflectivity, finite apertures, curvature, or off-axis
specular scattering.
"""
function track_latest_photon_branch(
    name::AbstractString;
    ray::LatestPhotonRay=LatestPhotonRay(),
)
    x, xp = ray.x, ray.xprime
    y, yp = ray.y, ray.yprime
    path_length = ray.path_length
    for ele in latest_photon_branch(name).line
        length = Float64(ele.L)
        x += length * xp
        y += length * yp
        path_length += length
    end
    return LatestPhotonRay(x, xp, y, yp, path_length)
end
'''


def replace_exactly(text: str, old: str, new: str, count: int) -> str:
    actual = text.count(old)
    if actual != count:
        raise RuntimeError(f"Expected {count} occurrences of {old!r}, found {actual}")
    return text.replace(old, new)


def main() -> None:
    text = SOURCE.read_text(encoding="utf-8")
    text = replace_exactly(
        text,
        "# Translated from Bmad lattice file: lat.bmad\n\nusing Beamlines\n",
        "# Repaired from the Tao 20260814-0 raw SciBmad export of lat.bmad.\n"
        "# See bmad_to_scibmad_tools/CONVERSION_REPORT.md and\n"
        "# bmad_to_scibmad_tools/build_repaired_lattice.py.\n\n"
        "using Beamlines\n"
        "include(joinpath(@__DIR__, \"essential_supports\", \"latest_lattice_support.jl\"))\n",
        1,
    )

    fork_pattern = re.compile(r"Fork\(to_line\s*=\s*[A-Za-z0-9_!]+\)")
    text, fork_count = fork_pattern.subn("Marker()", text)
    if fork_count != 11:
        raise RuntimeError(f"Expected 11 Fork replacements, found {fork_count}")

    text = replace_exactly(text, "Mirror()", "Marker()", 2)

    replacements = {
        "id_s1a!s1 = LineElement(L =  8.84599999998708864E-001)":
            "id_s1a!s1 = LatestWigglerSegment(L =  8.84599999998708864E-001, "
            "s_offset = 0.0, n_steps = 90)",
        "id_s1a!s2 = LineElement(L =  2.92999999999892680E-001)":
            "id_s1a!s2 = LatestWigglerSegment(L =  2.92999999999892680E-001, "
            "s_offset = 8.84599999998708864E-001, n_steps = 30)",
        "id_s1a!s3 = LineElement(L =  1.17740000000139844E+000)":
            "id_s1a!s3 = LatestWigglerSegment(L =  1.17740000000139844E+000, "
            "s_offset = 1.17759999999860154E+000, n_steps = 120)",
    }
    for old, new in replacements.items():
        text = replace_exactly(text, old, new, 1)

    # Apply Bmad's reference-time convention on the existing zero-length fork
    # marker. BeamTracking does not allow PatchParams and a four-potential on
    # the same LineElement, and a magnetostatic downstream segment is
    # insensitive to where this constant reference-time shift is placed.
    text = replace_exactly(
        text,
        "ls_1a1_1a2_1a3 = Marker()",
        "ls_1a1_1a2_1a3 = Marker(dt = LATEST_WIGGLER_REFERENCE_DT)",
        1,
    )

    # Bmad tracks each of the twelve combined-function DQX bends with
    # NUM_STEPS=100. Tao's SciBmad writer currently omits this attribute and
    # SciBmad's one-step default is far too coarse for their strong Kn1 field.
    dqx_pattern = re.compile(r"(?m)^(  dqx[1-6][bd] = SBend\()")
    text, dqx_count = dqx_pattern.subn(
        r"\1tracking_method = Yoshida(order = 6, n_steps = 100), ",
        text,
    )
    if dqx_count != 12:
        raise RuntimeError(f"Expected 12 DQX tracking-method repairs, found {dqx_count}")

    # Tao writes four Overlay expressions against super-lord element names
    # that are not emitted into the SciBmad file after superposition splits
    # the physical elements. The controlled quantity is the local field
    # coefficient, so every longitudinal slice receives the same expression.
    split_superlords = {
        "hs3b": ("hs3b!s1", "hs3b!s2"),
        "hs4b": ("hs4b!s1", "hs4b!s2"),
        "b48w": ("b48w!s1", "b48w!s2"),
        "b48e": ("b48e!s1", "b48e!s2"),
    }
    for superlord, slices in split_superlords.items():
        pattern = re.compile(
            rf"(?m)^{re.escape(superlord)}\.([A-Za-z0-9_]+)\s*=\s*(DefExpr\(.*\))$"
        )
        match = pattern.search(text)
        if match is None:
            raise RuntimeError(f"Missing control expression for split super-lord {superlord}")
        attribute, expression = match.groups()
        expanded = "\n".join(
            f"{slice_name}.{attribute} = {expression}" for slice_name in slices
        )
        text, count = pattern.subn(expanded, text, count=1)
        if count != 1:
            raise RuntimeError(f"Could not expand control expression for {superlord}")

    # The writer emits only the main-ring end marker even though every photon
    # branch Beamline references its own branch-end element.
    branch_end_markers = "\n".join(
        f"  end_b{branch} = Marker()" for branch in range(1, 12)
    )
    text = replace_exactly(
        text,
        "end    # @elements",
        f"{branch_end_markers}\nend    # @elements",
        1,
    )

    if re.search(r"\b(Fork|Mirror)\(", text):
        raise RuntimeError("Unsupported Fork or Mirror constructor remains")
    if "id_s1a!s1 = LineElement" in text or "id_s1a!s2 = LineElement" in text or "id_s1a!s3 = LineElement" in text:
        raise RuntimeError("A field-free ID_S1A placeholder remains")

    if "latest_photon_branches" in text:
        raise RuntimeError("Photon registry unexpectedly exists in the raw export")
    text = text.rstrip() + PHOTON_REGISTRY + "\n"

    OUTPUT.write_text(text, encoding="utf-8")
    print(f"Wrote {OUTPUT}")
    print(
        "Replacements: 11 Fork->Marker, 2 Mirror->Marker, "
        "3 wiggler segments plus reference-time patch, 12 DQX tracking methods, "
        "4 split-superlord controls, 11 branch-end markers, photon registry"
    )


if __name__ == "__main__":
    main()

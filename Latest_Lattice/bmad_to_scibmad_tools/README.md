# Bmad-to-SciBmad conversion and validation tools

This directory documents the current CHESS-U 6 GeV CESR lattice and its first
reproducible Bmad-to-SciBmad conversion audit.

The concise acceptance table is in
`../scibmad_validation/FINAL_VALIDATION_SUMMARY.md`; the detailed
writer-defect and repair analysis is in `CONVERSION_REPORT.md`.

The parent directory keeps the current SciBmad lattice, its runtime support,
and the two-file Bmad lattice. This directory holds generators, exporters,
diagnostics, validation programs, and their explanatory documentation. None
of these tools is required for an ordinary load of the repaired SciBmad
lattice; runtime dependencies live in `../essential_supports/`.

## Inputs

- `../lat.bmad`: top-level `CHESS-U_6000MEV_20251020_S1` lattice file.
- `../chess-u_02.bmad`: base CESR/CHESS-U layout called by `lat.bmad`.

The two input files were copied from the former `CESR Project/Archive`
directory. Their duplicate Archive copies were compared byte-for-byte before
that directory was removed. No cryptographic hashes are calculated in this
project.

## Raw exports

- `../bmad_reference/raw_exports/latest_cesr_scibmad_bmad_20260814.jl`:
  authoritative raw export produced
  by Tao `20260814-0`, the latest conda-forge Bmad build found on 2026-08-16.
- `../bmad_reference/raw_exports/tao_write_scibmad_bmad_20260814.log`: complete
  latest-export log.
- `../bmad_reference/raw_exports/latest_cesr_scibmad.jl` and
  `../bmad_reference/raw_exports/tao_write_scibmad.log`: retained comparison
  output from the pre-existing local Tao `20260801-1` environment.

The authoritative raw export was generated with:

```text
tao -noinit -noplot -nostartup -lattice_file lat.bmad
Tao> write scibmad bmad_reference/raw_exports/latest_cesr_scibmad_bmad_20260814.jl
Tao> quit
```

The raw file is deliberately unedited. It does not yet load in the current
SciBmad project; see `CONVERSION_REPORT.md`.

## Repaired executable main ring

`../latest_cesr_scibmad_repaired.jl` is the loadable charged-particle branch-0
model. It is generated deterministically from the authoritative raw export by
`build_repaired_lattice.py`; do not edit the generated file by hand.

The repair preserves all 1,177 branch-0 tracking elements and the exact
`768.4378690000005 m` circumference. It restores the three phase-continuous
`ID_S1A` wiggler slices, Bmad's wiggler reference-time convention, the omitted
100-step integration settings on all 12 strong DQX combined-function bends,
the four control expressions that target split super-lords, and the missing
branch-end markers. Unsupported photon forks and mirrors are represented by
zero-length markers in the generated Beamlines objects. An explicit registry
and branch-local paraxial reference-ray helper make all eleven archived photon
lines queryable and runnable; mirror reflectivity and off-axis scattering are
not modeled because current SciBmad/Beamlines has no photon `Fork`/`Mirror`
tracker.

Load it from Julia with:

```julia
include("CESR Project/Latest_Lattice/latest_cesr_scibmad_repaired.jl")
# The resulting Beamline is named `cesr`.
```

To regenerate and smoke-test it from PowerShell:

```powershell
& "C:/Users/JoeyN/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe" `
  "CESR Project/Latest_Lattice/bmad_to_scibmad_tools/build_repaired_lattice.py"
julia --startup-file=no --project="CESR Project" `
  "CESR Project/Latest_Lattice/bmad_to_scibmad_tools/load_repaired_lattice.jl"
```

## Numerical validation

- `../scibmad_validation/LOCAL_MAP_COMPARISON.md` and
  `../scibmad_validation/scibmad_local_map_comparison.csv` contain the
  1,177-element Bmad/SciBmad local-map scan. Outside the split wiggler block,
  the largest local matrix discrepancy is `4.29e-13`.
- `compare_wiggler_block.jl` shows that the complete wiggler block has an
  affine/exit-orbit mismatch of `7.12e-15`. Its remaining `5.89e-6` `R12`
  difference is also produced when Bmad is switched from its standard matrix
  approximation to Runge-Kutta field tracking.
- `RING_OPTICS_COMPARISON.md` records RF-on closed-orbit agreement to
  `5.25e-14`; the three one-turn eigenphase tune differences are approximately
  `5.9e-14`, `2.4e-9`, and `5.5e-8`.
- `../scibmad_validation/FULL_CONTROL_VALIDATION.md` compares one-pass tracking derivatives for all
  124 Overlay/Group lords, all 347 lord-to-slave relationships, and 475
  observation points. The median per-control maximum relative discrepancy is
  `3.74e-4`; the worst informative response is `1.84%`, with cosine
  `0.999834` and small absolute norm.
- `../scibmad_validation/GIRDER_VALIDATION.md` covers all 12 girders, all six alignment parameters,
  150 member tracking elements/slices, and 972 response observations. The
  supplied `set_latest_girder!` interface reproduces Bmad's linearized member
  geometry. Offset responses agree closely, while pitch through strong DQX
  combined-function bends differs by as much as `21.7%` because the current
  SciBmad tracker uses straight multipoles in a curved reference system.
- `../scibmad_validation/PHOTON_BRANCH_VALIDATION.md` verifies exact element counts and lengths for
  all 11 photon lines, fork lookup, mirror metadata, and branch-local paraxial
  reference-ray propagation.

Nominal charged-particle branch-0 tracking, optics, and the complete exported
control graph are therefore validated. The remaining quantified limitations
are exact curved-coordinate multipole response under girder pitch and general
photon mirror optics; see `CONVERSION_REPORT.md`.

## Reproduce the inventory

From Windows PowerShell:

```powershell
wsl.exe -d Ubuntu-Bmad -- `
  /home/joeyfurina/.cache/codex-bmad-20260814/bin/python `
  "/mnt/d/Ring_Design_Development/CESR Project/Latest_Lattice/bmad_to_scibmad_tools/diagnose_conversion.py"
```

This writes:

- `../bmad_reference/inventory/bmad_branches.csv`;
- `../bmad_reference/inventory/bmad_element_inventory.csv`;
- `../bmad_reference/inventory/bmad_control_lords.csv`;
- `../bmad_reference/inventory/bmad_control_relations.csv`;
- `../bmad_reference/inventory/conversion_inventory.json`.

The raw SciBmad load smoke test is:

```powershell
julia --startup-file=no --project="CESR Project" `
  "CESR Project/Latest_Lattice/bmad_to_scibmad_tools/load_raw_export.jl"
```

Its captured failure is in `../scibmad_validation/scibmad_raw_load.log`.

# CESR Digital Twin in SciBmad

This project creates a digital twin of the Cornell Electron Storage Ring
(CESR) in SciBmad. The model preserves the physical structure and numerical
settings of the real ring represented by the Bmad CESR lattice. It also
reproduces the Bmad control structure with SciBmad deferred expressions
(`DefExpr`), allowing control knobs to modify dependent element attributes
without changing the static lattice definition.

The project also implements a planar wiggler element for SciBmad. Its tracking
model is based on the magnetic four-potential and supports GTPSA transport-map
calculation.

## Core Model Files

| Path | Purpose |
|---|---|
| `cesr.bmad` | Original Bmad CESR lattice and the source of the element settings, Overlay definitions, and Group control coefficients. |
| `cesr.jl` | Static SciBmad lattice. It contains the physical ring elements, their baseline numerical settings, the consolidated solenoid model, the CESR wigglers, kicker strengths, and RF-cavity helpers. It does not contain mutable control settings. |
| `cesr_controls.jl` | Deferred-expression control layer generated from `cesr.bmad`. It implements the Bmad-style Overlays and Groups and propagates their contributions to physical element attributes. |
| `cesr_model.jl` | Main model entry point. `load_cesr_model()` creates an independent copy of the static ring and attaches the complete control layer. It can also explicitly enable or disable the RF cavities. |
| `Project.toml` | Julia project definition and direct package dependencies. |
| `.gitignore` | Excludes the local Julia `Manifest.toml` and Jupyter checkpoint directories from version control. |
| `.gitattributes` | Repository text-file settings. |

The initialized model returned by `load_cesr_model()` contains two fields:
`model.ring`, the independent CESR beamline, and `model.controls`, the mutable
collection of Overlay and Group knobs.

## Directories

### `wigglers/`

This directory contains the CESR wiggler implementation and its derivation:

- `wiggler.jl` defines the planar wiggler field, magnetic four-potential,
  SciBmad element construction, and GTPSA transport-map utilities.
- `cesr_wiggler_transport_map.ipynb` explains the field and four-potential
  formulas, derives the transport-map calculation, and demonstrates the
  corresponding implementation.

### `bmad_comparison/`

This directory contains the Bmad/Tao reference outputs and the numerical
comparison between the Bmad and SciBmad models:

- `bmad_reference_output.tar.gz` and
  `bmad_reference_rf_on_output.tar.gz` contain the exported Bmad reference
  results.
- `bmad_scibmad_rf_on_comparison.csv` contains the element-by-element RF-on
  comparison.
- `bmad_scibmad_rf_on_comparison_summary.md` summarizes the closed-orbit,
  tune, local-map, cumulative-map, and element-orbit differences.
- `run_bmad_reference.sh` and `export_bmad_reference.py` generate the Bmad
  reference data on a Linux system with Bmad/Tao installed.
- `compare_rf_on_optics.jl` compares the RF-on closed orbit and tunes.
- `BMAD_REFERENCE_EXPORT.md` documents the Bmad export workflow.

### `test_codes/`

This directory contains development and regression utilities used to verify
that the ring and its control structure run correctly:

- `generate_cesr_controls.py` reads the Overlay and Group definitions from
  `cesr.bmad` and regenerates the coefficient tables in `cesr_controls.jl`.
- `test_cesr_controls.jl` verifies zero-control consistency, direct Overlays,
  Group-to-Overlay propagation, sextupole controls, and controls on
  superimposed element slices.
- `test_rf_on_twiss.jl` verifies the six-dimensional RF-on closed orbit and
  Twiss calculation.
- `test_bmad_scibmad.jl` performs the GTPSA element-by-element comparison with
  the exported Bmad transfer maps.

## Current Agreement with Bmad

The comparison must be described by several metrics rather than a single
accuracy number. With the RF cavities enabled, the current results are:

- maximum local element transfer-matrix difference: `2.476e-4`;
- CESR wiggler local matrix differences: approximately `3.8e-5` to `4.0e-5`;
- solenoid-quadrupole overlap local matrix differences: approximately
  `9.1e-5`;
- maximum isolated-element exit-orbit difference: `2.526e-6`;
- maximum closed-orbit coordinate difference: `6.404e-7`;
- maximum cumulative transfer-matrix difference around the ring: `1.267e-2`.

Therefore, the characteristic local agreement is in the `10^-5` to `10^-4`
range, while the worst local discrepancy is currently at the `10^-4` level.
The cumulative one-turn difference is larger because the local discrepancies
accumulate around the full ring. Detailed results and the reproduction commands
are available in `bmad_comparison/bmad_scibmad_rf_on_comparison_summary.md`.

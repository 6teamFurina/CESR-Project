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

- `bmad_optics_outputs/` contains the exported RF-off and RF-on Bmad
  optics/map references: `bmad_reference_rf_off_output.tar.gz` and
  `bmad_reference_rf_on_output.tar.gz`, respectively. 
- `optic_maps_rf_off/` and `optic_maps_rf_on/` contain the corresponding
  SciBmad comparison results and mode-specific comparison programs.
- `optic_maps_rf_off/bmad_scibmad_rf_off_comparison.csv` and
  `optic_maps_rf_on/bmad_scibmad_rf_on_comparison.csv` contain the
  element-by-element RF-off and RF-on comparisons.
- `bmad_scibmad_rf_on_comparison_summary.md` summarizes the closed-orbit,
  tune, local-map, cumulative-map, and element-orbit differences.
- `test_control_response_scibmad.jl` computes the complete SciBmad
  control-to-closed-orbit response matrix with first-order GTPSA parameters
  and compares it with the labeled Tao/Bmad matrix.
- `bmad_scibmad_control_response_summary.md` summarizes the RF-on and RF-off
  control-response comparisons. Detailed matrices and reports are stored in
  `bmad_control_response_rf_on/` and `bmad_control_response_rf_off/`.
- `bmad_optics_outputs/run_bmad_reference.sh` and
  `bmad_optics_outputs/export_bmad_reference.py` generate the Bmad reference
  data on a Linux system with Bmad/Tao installed.
- `optic_maps_rf_on/compare_rf_on_optics.jl` compares the RF-on closed orbit
  and tunes.
- `bmad_optics_outputs/BMAD_REFERENCE_EXPORT.md` documents the Bmad export
  workflow.

### `test_codes/`

This directory contains development and regression utilities used to verify
that the ring and its control structure run correctly:

- `generate_cesr_controls.py` reads the Overlay and Group definitions from
  `cesr.bmad` and regenerates the coefficient tables in `cesr_controls.jl`.
- `test_cesr_controls.jl` verifies zero-control consistency, direct Overlays,
  Group-to-Overlay propagation, sextupole controls, controls on superimposed
  element slices, and the conversion of laboratory-frame Bmad kicks into
  tilted SciBmad element coordinates. The current suite contains 24 passing
  tests.
- `test_rf_on_twiss.jl` verifies the six-dimensional RF-on closed orbit and
  Twiss calculation.
- `test_rf_off_closed_orbit.jl` verifies the patched four-dimensional RF-off
  coasting closed orbit against an independent Float64 central-difference
  Newton solve.
- `test_bmad_scibmad.jl` performs the GTPSA element-by-element comparison with
  the exported Bmad transfer maps.
- `test_control_response_tao.py` runs in the Linux Bmad/PyTao environment and
  exports the labeled CESR control-to-closed-orbit response matrix for either
  the RF-off or RF-on lattice.

## RF-Off Coasting Closed-Orbit Patch

`scibmad_coasting_forwarddiff_patch.jl` fixes the CESR RF-off closed-orbit
failure without modifying the installed SciBmad or BeamTracking packages. When
`coasting_beam == true`, the physical closed-orbit problem contains only the
four transverse unknowns `(x, px, y, py)`, while `z` and `pz` remain fixed.

SciBmad already selects this 4D formulation for a coasting beam. The failure
occurred because the current BeamTracking implicit-integrator ForwardDiff path
assumes that each `Dual` value carries six partial derivatives, whereas the
native coasting path seeds only four. The external adapter therefore evaluates
the one-turn residual with six ForwardDiff directions and extracts the
transverse `4 x 4` Jacobian. Only this Jacobian construction is patched; the
closed-orbit iteration, convergence checks, singular-matrix handling, and
return code use SciBmad's native `BatchSolve.newton!` implementation.

The patched RF-off orbit can be calculated with:

```julia
include("scibmad_coasting_forwarddiff_patch.jl")
using .SciBmadCoastingForwardDiffPatch

solution = find_closed_orbit_coasting_forwarddiff(
    ring;
    coasting_beam=true,
)
```

The regression comparison is run from the project root with:

```console
julia --project=. test_codes/test_rf_off_closed_orbit.jl
```

For the current CESR lattice, both the patched ForwardDiff/BatchSolve method
and the Float64 central-difference reference converge in three Newton
iterations. Their closed-orbit vectors differ by at most `8.37e-16`; the
patched orbit has a one-turn transverse closure residual of `3.03e-15`, and
the maximum row-sum difference between the two `4 x 4` residual Jacobians is
`1.79e-8`.

## Control-Response Validation

The control-response test differentiates the horizontal and vertical closed
orbit at 99 `DET_*` markers with respect to 119 Bmad-compatible CESR Overlay
controls (58 horizontal and 61 vertical). The resulting matrix has shape
`198 x 119` and units of `m/rad`.

The Bmad reference is generated with
`test_codes/test_control_response_tao.py`. The SciBmad comparison is run from
the project root with:

```console
julia --project=. bmad_comparison/test_control_response_scibmad.jl --mode=both
```

All 119 SciBmad control derivatives are calculated simultaneously using
first-order GTPSA parameters. If the one-turn map derivatives are `A = dF/dz`
and `B = dF/dk`, the closed-orbit response is obtained from
`d z_closed/dk = (I - A)^-1 B`. RF-on uses the full 6D closed orbit. RF-off
uses the corresponding 4D coasting-beam equation with fixed `z = pz = 0`.

The current comparison results are:

| Mode | Relative Frobenius difference | Maximum absolute difference (m/rad) | Full-matrix correlation | GTPSA closure residual |
|---|---:|---:|---:|---:|
| RF-on | `0.203229%` | `6.7472e-2` | `0.999997991471` | `7.105e-15` |
| RF-off | `0.201201%` | `6.7263e-2` | `0.999998031423` | `7.105e-15` |

Agreement by response block is:

| Mode | `x <- H` | `x <- V` | `y <- H` | `y <- V` |
|---|---:|---:|---:|---:|
| RF-on | `0.04864%` | `0.49898%` | `0.18615%` | `0.24958%` |
| RF-off | `0.04779%` | `0.49911%` | `0.18526%` | `0.24716%` |

The block values are relative Frobenius differences. The overall agreement is
therefore approximately `0.2%`, with correlations above `0.9999979` in both
RF configurations.

Bmad `HKICK` and `VKICK` controls specify laboratory-frame orbit kicks, while
SciBmad normal and skew dipole multipoles rotate with the element alignment.
For alignment tilt `t`, the control layer applies
`HKICK -> (Kn0L, Ks0L) = (-cos(t), -sin(t)) HKICK` and
`VKICK -> (Kn0L, Ks0L) = (-sin(t), cos(t)) VKICK`. This conversion is required
for the correct sign of horizontal responses and for vertical correctors
implemented by 45-degree tilted elements.

The stored control-response comparison retains its original Float64
finite-difference baseline solve for reproducibility. The external coasting
patch described above now provides a ForwardDiff alternative using SciBmad's
native `BatchSolve.newton!`. In both workflows, the reported 119 control
derivatives are calculated with GTPSA and do not use control finite
differences.

Each RF-mode output directory contains:

- `scibmad_control_response_<mode>.csv`: the labeled SciBmad `198 x 119`
  response matrix;
- `bmad_scibmad_control_response_entries_<mode>.csv`: all 23,562 individual
  Bmad/SciBmad entry comparisons and differences;
- `bmad_scibmad_control_response_columns_<mode>.csv`: one relative 2-norm
  error for each of the 119 control columns;
- `bmad_scibmad_control_response_summary_<mode>.md`: detailed mode-specific
  metrics, response blocks, worst entries, worst control, and singular values.

The root `bmad_scibmad_control_response_summary.md` combines both modes and
records the coordinate convention and RF-off solver details. It can be
regenerated from existing matrices without rerunning the GTPSA tracking:

```console
julia --project=. bmad_comparison/test_control_response_scibmad.jl --mode=summary
```

## Current Agreement with Bmad

The comparison must be described by several metrics rather than a single
accuracy number. Matrix percentages are normalized by the corresponding Bmad
matrix scale, using `max(abs(R_SciBmad - R_Bmad)) / max(abs(R_Bmad))`; this
avoids unstable percentages from individual Bmad matrix entries that are zero
or nearly zero. With the RF cavities enabled, the current results are:

- maximum local element transfer-matrix difference: `2.476e-4`, corresponding
  to a maximum normalized relative difference of `0.007645%`;
- CESR wiggler local matrix differences: approximately `3.8e-5` to `4.0e-5`,
  corresponding to `0.00160%` to `0.00168%` relative difference;
- solenoid-quadrupole overlap local matrix differences: approximately
  `9.1e-5`, or approximately `0.00584%` relative difference;
- maximum isolated-element exit-orbit difference: `2.526e-6`, or `0.02697%`
  relative to the largest Bmad orbit coordinate at that element;
- maximum closed-orbit coordinate difference: `6.404e-7`; the complete orbit
  vectors differ by `0.0390%` in relative 2-norm;
- maximum cumulative transfer-matrix difference around the ring: `1.267e-2`,
  equal to `0.05221%` of the Bmad cumulative-map scale at that location; the
  maximum normalized cumulative difference over all locations is `0.16772%`.

Therefore, the characteristic local relative agreement is at the `10^-5`
level, with a worst local relative discrepancy of `7.645e-5` (`0.007645%`).
The cumulative normalized discrepancy reaches the `10^-3` level because the
local differences accumulate around the full ring. Detailed absolute and
relative results, including percentage differences for the closed orbit and
tunes, are available in
`bmad_comparison/bmad_scibmad_rf_on_comparison_summary.md`.

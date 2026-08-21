# CESR Digital Twin in SciBmad

This project develops a SciBmad digital twin of the Cornell Electron Storage
Ring (CESR). The default model is the validated latest CHESS-U 6 GeV lattice in
`Latest_Lattice/latest_cesr_scibmad_repaired.jl`. The earlier maintained CESR
ring, its controls, tests, and comparison artifacts are retained under
`older_ring_version/` for historical reproduction.

The project also implements a planar wiggler element for SciBmad. Its tracking
model is based on the magnetic four-potential and supports GTPSA transport-map
calculation.

## Current lattice files

| Path | Purpose |
|---|---|
| `Latest_Lattice/latest_cesr_scibmad_repaired.jl` | Validated latest SciBmad lattice and the default for all new CESR calculations. |
| `Latest_Lattice/lat.bmad` | Latest Bmad entry lattice used only for conversion/reference validation. |
| `Latest_Lattice/chess-u_02.bmad` | Base Bmad layout loaded by `lat.bmad`; it is part of the same two-file Bmad lattice. |

Support programs and their documentation are in
`Latest_Lattice/support_codes/`, SciBmad validation outputs are in
`Latest_Lattice/scibmad_validation/`, and Bmad reference/load outputs are in
`Latest_Lattice/bmad_reference/`.

## Historical ring files

| Path | Purpose |
|---|---|
| `older_ring_version/cesr.bmad` | Original Bmad CESR lattice and the source of the element settings, Overlay definitions, and Group control coefficients. |
| `older_ring_version/cesr.jl` | Historical static SciBmad lattice. It contains the physical ring elements, baseline numerical settings, consolidated solenoid model, CESR wigglers, kicker strengths, and RF-cavity helpers. |
| `older_ring_version/cesr_controls.jl` | Deferred-expression control layer generated from the historical `cesr.bmad`. |
| `older_ring_version/cesr_model.jl` | Historical model entry point. `load_cesr_model()` creates an independent copy of that ring and attaches its control layer. |

The initialized model returned by `load_cesr_model()` contains two fields:
`model.ring`, the independent CESR beamline, and `model.controls`, the mutable
collection of Overlay and Group knobs.

## Dataset Vision

The validated CESR model is intended to support a reproducible
parameter-to-observable dataset and later uncertainty-aware inverse problems.
This section records the current dataset design; it is a research plan, not a
claim that the complete generator or every proposed observable has already
been implemented and validated.

### Dataset semantics

The first dataset uses the following simplified parameter and observation
model:

```text
command = nominal
physical = apply_element_errors(command, error)
physics_output = SciBmad(physical)
observable_readback = measurement_model(physics_output, noise)
```

- `nominal` and `command` are equal in the first version.
- `command` contains the known settings supplied to the accelerator model.
- `physical` contains the latent values that actually enter SciBmad. It can
  differ from `command` because of errors attached to physical accelerator
  elements, including strength, calibration, alignment, and field errors.
- `physics_output` contains the noise-free simulator results.
- `observable_readback` contains the simulated beam-diagnostic observations,
  optionally with measurement noise, missing channels, or other diagnostic
  effects.

Control-system setpoint and indicated-setting errors are deliberately out of
scope for the first dataset. The qualified name `observable_readback` is used
because accelerator control systems also commonly use "readback" for a
device's indicated setting.

Every parameter and observation must retain a stable device or element
identifier, physical units, coordinate convention, availability mask,
uncertainty information, random seed, lattice/configuration version, and
generation provenance. The stored physical values remain available as
simulation truth even when selected values are hidden from a particular
training task.

### Candidate generated data

The parameter side can include quadrupole and sextupole strengths, corrector
kicks, RF settings, element alignment, calibration or field errors, and other
explicitly supported mutable lattice quantities. Error realizations may be
independent, family-correlated, spatially correlated, slowly varying, or
discrete faults, provided their distributions and seeds are recorded.

The initial observable and label set can include:

- closed orbit around the ring and horizontal/vertical orbit at the 99
  `DET_*` markers;
- transverse and longitudinal tunes where defined;
- Twiss functions, dispersion, chromaticity, and RF-off momentum derivatives
  of selected Twiss quantities;
- local, cumulative, one-turn, orbit-response, and selected optics-response
  maps;
- solver convergence, closure residual, linear stability, and physical or
  control-limit flags;
- noise-free physics outputs together with noisy or partially missing
  observable readbacks.

Expensive nonlinear outputs are planned as a separate, linked high-fidelity
data product rather than as mandatory fields in every base sample. Candidate
outputs include dynamic aperture, momentum aperture, tune footprint, particle
survival, survival turns, and loss location. A dynamic-aperture record should
retain the boundary or survival map as a function of launch angle, momentum
offset, and tracking turns, together with the tracking and aperture
configuration; a single scalar aperture area is not a sufficient primary
label.

### Planned learning tasks

The same generated machine states can expose different task-specific views:

1. **Forward observable prediction:** use all `command` values to predict
   `observable_readback`. With zero element error this is a deterministic
   nominal forward model. If latent errors vary while only commands are
   provided, the target is a conditional distribution rather than a unique
   deterministic value.
2. **Masked physical-parameter reconstruction:** provide observable readbacks,
   an explicit parameter mask, and any unmasked physical-parameter context,
   then infer the masked physical values. Masks may cover individual devices,
   element families, contiguous lattice regions, parameter types, or
   combinations of faults.
3. **Physical-error localization and estimation:** provide `command` and
   `observable_readback`, then infer the locations, types, magnitudes, and
   uncertainty of the latent physical errors that separate `physical` from
   `command`.
4. **Later intervention and diagnosis tasks:** combine multiple observations
   from the same latent machine state under known corrector, quadrupole, RF, or
   diagnostic excitations to improve identifiability and test
   counterfactual predictions.

Evaluation should include physical consistency: inferred parameters should be
placed back into SciBmad and judged by their ability to reproduce the observed
machine state and predict held-out interventions, not only by componentwise
parameter error.

### Open research-design questions

Two issues require explicit study before the corresponding benchmarks are
finalized:

1. **Inverse degeneracy and identifiability.** Different element errors may
   produce indistinguishable or nearly indistinguishable orbit, tune, or Twiss
   observations. Candidate treatments include sensitivity-Jacobian
   rank/singular-value analysis, identifiable parameter groups, multiple
   diagnostic excitations, physically motivated priors, posterior or
   multi-hypothesis predictions, and evaluation in the observable or
   identifiable subspace rather than requiring an unjustified unique inverse.
2. **Dynamic-aperture and nonlinear-data generation.** The launch-coordinate
   convention, momentum slices, tracking turns, physical apertures, radiation
   settings, boundary-search method, fidelity hierarchy, storage
   representation, and Bmad comparison protocol must be fixed. The nonlinear
   labels must then be validated before they are described as high-confidence
   CESR training targets.

## Directories

### `wigglers/` and `older_ring_version/wigglers/`

The shared `wigglers/wiggler.jl` file contains the wiggler implementation used
by both lattice generations. Historical derivation and experiment materials
are grouped in `older_ring_version/wigglers/`:

- `wiggler.jl` defines the planar wiggler field, magnetic four-potential,
  SciBmad element construction, and GTPSA transport-map utilities.
- `cesr_wiggler_transport_map.ipynb` explains the field and four-potential
  formulas, derives the transport-map calculation, and demonstrates the
  corresponding implementation.

### `older_ring_version/bmad_comparison/`

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
- `nonlinear_twiss/compute_nonlinear_twiss.jl` computes the RF-off second-order
  GTPSA Twiss table, exports chromatic and alpha/beta momentum derivatives to
  CSV, and plots the three quantities around the ring. With RF cavities on,
  `pz` (`delta`) is a dynamical variable rather than a controllable parameter,
  so it cannot serve as the independent variable for these derivatives. The
  nonlinear Twiss study therefore uses RF cavities off.
- `bmad_optics_outputs/BMAD_REFERENCE_EXPORT.md` documents the Bmad export
  workflow.

### `dataset_benchmark/`

This independent directory contains the matched RF-on dataset-throughput
benchmark. Bmad and SciBmad consume the same deterministic 1000-sample,
119-control input and produce the same 198 labeled detector coordinates. It
contains the runners, shared inputs, archived source package, Bmad-compatible
reference lattice, preliminary runs, formal 1000-sample results, and numerical
comparison reports.

### `older_ring_version/test_codes/`

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

`older_ring_version/scibmad_coasting_forwarddiff_patch.jl` fixes the historical
CESR RF-off closed-orbit
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
include("older_ring_version/scibmad_coasting_forwarddiff_patch.jl")
using .SciBmadCoastingForwardDiffPatch

solution = find_closed_orbit_coasting_forwarddiff(
    ring;
    coasting_beam=true,
)
```

The regression comparison is run from the project root with:

```console
julia --project=. older_ring_version/test_codes/test_rf_off_closed_orbit.jl
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
`older_ring_version/test_codes/test_control_response_tao.py`. The SciBmad comparison is run from
the project root with:

```console
julia --project=. older_ring_version/bmad_comparison/test_control_response_scibmad.jl --mode=both
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

The RF-off control-response workflow now uses the external coasting patch for
its baseline orbit: six-direction ForwardDiff supplies the transverse `4 x 4`
Jacobian and SciBmad's native `BatchSolve.newton!` solves the 4D closed-orbit
equation. The reported 119 control derivatives are calculated simultaneously
with GTPSA and do not use control finite differences.

For second- and higher-order corrector-response coefficients, the project
default is likewise GTPSA parameterization, combined with implicit
differentiation of the closed-orbit fixed-point equation when the orbit itself
depends on the controls. Direction ensembles are reported as
`median [P10, P90]`. Signed finite differences are retained as independent
validation and higher-order-contamination checks, while direct nonlinear
solutions define the amplitude range in which a truncated response expansion
is valid. The adopted quadratic mixed-block calculation and its four-sign
validation are documented in
`dataset_benchmark/orbit/error_analysis/mixed_terms/`.

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
julia --project=. older_ring_version/bmad_comparison/test_control_response_scibmad.jl --mode=summary
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
`older_ring_version/bmad_comparison/bmad_scibmad_rf_on_comparison_summary.md`.

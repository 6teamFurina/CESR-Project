# Relative local-orbit predictor comparison

All three predictors use nominal-K2 BPM orbit differences and known corrector
commands. Exact target-local orbit is loaded only after prediction and BPM-only
ridge selection, and is used solely for scoring.

- targets / latent realizations: 76 / 8 per target
- evaluated nonzero bumps: 4 per realization
- hidden machine errors: target offset, all-other-sextupole offsets, and
  independent quadrupole strength errors from the frozen all-76 tensor
- BPM noise/offset/gain errors: none
- maximum nominal target-command consistency error: 4.01155e-12 micrometers
- global-MAP ridge ratio selected by 5-fold held-out-BPM residual:
  0.01
- largest two-sided transverse momentum-block condition number:
  4.28986

| method | x RMSE [um] | y RMSE [um] | 2D RMSE [um] | median [um] | P90 [um] | max [um] |
|---|---:|---:|---:|---:|---:|---:|
| command_only | 10.302 | 20.451 | 22.899 | 13.026 | 37.964 | 91.689 |
| two_sided_transport | 0.237 | 0.185 | 0.301 | 0.049 | 0.250 | 4.385 |
| global_map | 7.048 | 13.665 | 15.375 | 8.469 | 24.095 | 96.678 |

`command_only` uses the nominal SciBmad control-to-target response.
`two_sided_transport` adds a correction inferred from the residual x/y orbit
at the nearest upstream and downstream BPMs using nominal local transport.
`global_map` fits a regularized effective-corrector correction to all BPM
residuals, with the prior centered on the known commanded bump.

## Two-sided BPM method: required quantities

The local-orbit predictor itself runs only at nominal K2. For each target and
each bump, the machine-facing inputs are:

| quantity | source | directly available on a machine? |
|---|---|---|
| x/y orbit at the nearest upstream BPM | zero-bump and current-bump BPM readbacks | yes |
| x/y orbit at the nearest downstream BPM | zero-bump and current-bump BPM readbacks | yes |
| corrector settings defining the bump | setpoints/readbacks, converted to model fields | yes, subject to calibration |
| target K2 state | sextupole setpoint/readback, used to select nominal K2 | yes, subject to calibration |
| corrector-to-two-BPM response | nominal latest-lattice SciBmad model | model-derived; BPM part can also be measured |
| corrector-to-target response | nominal latest-lattice SciBmad model | not directly measurable at the sextupole |
| upstream-to-target and upstream-to-downstream 4D transport | nominal latest-lattice SciBmad model | model-derived, not a direct readback |

Only four BPM channels per state are consumed: upstream x/y and downstream
x/y. The full K2 scan is needed later by the magnetic-center inverse, but it is
not needed to predict these nominal-K2 local bump coordinates.

## Two-sided BPM method: implementation

For each nonzero bump, the code first subtracts the zero-bump state:

`delta y_b = y_BPM(b, nominal K2) - y_BPM(0, nominal K2)`.

Known corrector commands are propagated through nominal SciBmad responses to
obtain `delta y_model` at the two BPMs and `delta z_command` at the target.
The measured-minus-model residuals at the upstream and downstream BPMs are
called `r_u` and `r_d`.

The nominal transverse map from upstream to downstream is partitioned so that

`r_d = A_ud r_u + B_ud p_u`,

where `p_u = (delta px, delta py)`. The two unmeasured momenta are inferred by

`p_u = pinv(B_ud) (r_d - A_ud r_u)`.

The residual is then transported to the target:

`delta z_residual = A_us r_u + B_us p_u`,

and the final prediction is

`delta z_prediction = delta z_command + delta z_residual`.

The implementation builds this as one target-specific `2 x 4` matrix acting
on `[r_ux, r_uy, r_dx, r_dy]`. Neighbor selection is circular around the ring,
and the one-turn map is used when the upstream/downstream interval crosses the
lattice boundary.

## Code and artifacts

The method lives in the `finite_bpm_inversion/` directory:

- [`analyze_local_orbit_predictors.py`](../../analyze_local_orbit_predictors.py):
  constructs the two-sided matrices, predicts all relative local orbits, and
  scores them only after prediction;
- [`generate_local_orbit_models.jl`](../../generate_local_orbit_models.jl):
  generates the latest-lattice SciBmad corrector responses, cumulative maps,
  one-turn map, and element inventories;
- [`validate_local_orbit_predictor_results.py`](../../validate_local_orbit_predictor_results.py):
  independently recomputes the summary statistics and checks result counts;
- [`two_sided_neighbors.csv`](two_sided_neighbors.csv): selected neighboring
  BPMs and transport conditioning for every target;
- [`per_prediction_errors.csv`](per_prediction_errors.csv): simulation-only
  truth comparison for every method, target, realization, and bump.

## Does inference require an unmeasurable machine quantity?

No exact target-local orbit, true sextupole offset, true quadrupole errors, or
other-sextupole misalignments enter the predictor. Those hidden values are not
required by the machine-facing inference. `target_orbits.npy` is used only in
simulation after prediction to calculate the error reported above; a real
machine cannot provide this direct scoring truth.

The method does depend on quantities at the target that are model-derived
rather than directly measured: the corrector-to-target response and local
transport matrices. It also presently assumes that corrector readbacks are
converted to physical kicks with the nominal calibration, and that BPM gains,
rolls, noise, and missing channels are absent. Stable BPM offsets cancel in the
zero-bump difference, but gain/roll/calibration errors do not. Therefore the
present result is implementable from machine readbacks plus a calibrated
SciBmad model, but it is not a model-free measurement of the internal orbit.

These standalone local-orbit estimates have now been propagated through the
K2-slope center fit. See the maintained end-to-end result in
[`../two_sided_center_inversion/SUMMARY.md`](../two_sided_center_inversion/SUMMARY.md).

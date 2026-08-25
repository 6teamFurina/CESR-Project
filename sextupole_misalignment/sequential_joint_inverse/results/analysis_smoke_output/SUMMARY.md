# Sequential excitation and all-target joint-inverse pilot

## Protocol

The exact dataset contains 4 latest-lattice SciBmad latent
machines.  In each machine, all 4 sextupole offsets, BPM gains,
corrector gains, K2 gains, quadrupole strength errors, and quadrupole rolls are
fixed while the 4 sextupoles are excited one at a time.  The
paired quadrupole-alignment case adds independent Gaussian x/y displacement to
each physical quadrupole with `50 um` RMS per plane, coherent across its slices
and fixed throughout all target scans.

Machines are split, without target-level leakage, into 2
training, 1 validation, and 1 held-out test
machines.  Each test machine has 2 independent
measurement realizations.  Measurement augmentation uses `5 um` RMS white
noise per BPM plane/read, 3072 repeated balanced eight-state cycles,
and a scalar random walk with `10 um` RMS endpoint change propagated through
the paired exact SciBmad drift secant.

The learned target is the two-component beam-relative sextupole magnetic
center.  All models receive full-ring BPM parity contrasts and nominal command
values; latent errors and exact target-local orbit are evaluation-only.

## Held-out result

| quadrupole-alignment input | inverse | 2D RMSE [um] | P90 [um] | P99 [um] | worst-target RMSE [um] | below 50 um |
|---|---|---:|---:|---:|---:|---:|
| No quadrupole alignment drift | physics_gls | 25.336 | 36.429 | 41.519 | 31.467 | 100.00% |
| No quadrupole alignment drift | shared_target_local_ridge | 23.526 | 35.586 | 37.546 | 28.179 | 100.00% |
| No quadrupole alignment drift | shared_joint_ridge | 26.010 | 36.579 | 42.878 | 33.224 | 100.00% |
| No quadrupole alignment drift | shared_joint_random_feature | 25.278 | 37.164 | 41.308 | 30.706 | 100.00% |
| 50 um/plane RMS quadrupole alignment drift | physics_gls | 52.155 | 80.687 | 90.379 | 84.116 | 75.00% |
| 50 um/plane RMS quadrupole alignment drift | shared_target_local_ridge | 60.540 | 79.289 | 88.601 | 82.576 | 37.50% |
| 50 um/plane RMS quadrupole alignment drift | shared_joint_ridge | 60.729 | 76.747 | 85.759 | 79.775 | 37.50% |
| 50 um/plane RMS quadrupole alignment drift | shared_joint_random_feature | 60.727 | 77.633 | 86.904 | 80.911 | 37.50% |

The best learned same-distribution result without quadrupole alignment is
`shared_target_local_ridge` at
`23.526 um` RMSE.  With the paired
50-um/plane quadrupole drift enabled, the best learned result is
`shared_target_local_ridge` at
`60.540 um` RMSE.  Training the joint
ridge only on the no-alignment distribution and evaluating it on the aligned
case gives `51.184 um` RMSE; this is the explicit
distribution-shift check and must not be replaced by the in-distribution row.

## Quadrupole-drift interpretation

The requested 50-um setting is interpreted as a per-plane RMS.  The finite
draw realizes `50.355 um`
in x and `49.585 um` in y.
Before any orbit correction, it changes the paired full-ring reference BPM
orbit by `836.545 um` RMS and
the beam-relative center truth by
`1264.510 um` 2D RMS.
The fraction of truths outside the maintained 1.5-mm bump radius changes from
`0.000%` to
`31.250%`.

This is intentionally the uncorrected residual-drift input requested for the
paired model test.  If the orbit excursion, tail error, or out-of-range
fraction dominates the aligned case, the next physical protocol must perform
and record a BPM-only orbit correction relative to the yearly nominal orbit
before the sextupole scans; a neural model must not be credited with replacing
that machine operation.

## Model definitions

- `physics_gls`: fixed latest-lattice covariance-uniform full-BPM source
  inverse, applied independently to each target block.
- `shared_target_local_ridge`: one parameter-sharing residual model for all
  targets, using only the target scan, target identity, and common baseline
  orbit modes.
- `shared_joint_ridge`: the same residual model plus a compact context derived
  from all 4 target scans in the machine; one inference call
  returns all `8` center coordinates.
- `shared_joint_random_feature`: the joint inputs plus a fixed nonlinear tanh
  feature layer, with only the output ridge weights fitted.

The response residual, baseline-orbit, and all-target context projections are
fit on training machines only.  Model and PCA files are saved with the result
so validation/test information cannot silently enter feature construction.

## Scope and limitations

This is a synthetic SciBmad pilot, not a CESR position-precision claim.  The
sample contains only 4 independent static machines, the error
priors other than the user-provided quadrupole drift are maintained sensitivity
settings rather than measured CESR distributions, BPM white noise is
independent, and the drift is one target-local scalar mode.  There is no
actuator hysteresis, K2 polarity asymmetry, missing/outlier BPM process, or
sim-to-real validation.  Exact target orbit is used only to form evaluation
truth.

The latest lattice emits the straight-multipole-in-curved-reference warning.
No girder pitch is varied in this experiment, so the documented curved-DQX
girder-pitch discrepancy is not an excitation here, but remains part of the
lattice provenance.

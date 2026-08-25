# Sequential excitation and all-target joint-inverse pilot

## Protocol

The exact dataset contains 16 latest-lattice SciBmad latent
machines.  In each machine, all 76 sextupole offsets, BPM gains,
corrector gains, K2 gains, quadrupole strength errors, and quadrupole rolls are
fixed while the 76 sextupoles are excited one at a time.  The
paired quadrupole-alignment case adds independent Gaussian x/y displacement to
each physical quadrupole with `50 um` RMS per plane, coherent across its slices
and fixed throughout all target scans.

Each physical case contains
`36,480` exact SciBmad state
lanes: 15 scan states and 15 paired drift-secant states per target and machine.

Machines are split, without target-level leakage, into 10
training, 3 validation, and 3 held-out test
machines.  Each test machine has 32 independent
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
| No quadrupole alignment drift | physics_gls | 33.923 | 54.069 | 86.559 | 60.070 | 87.25% |
| No quadrupole alignment drift | shared_target_local_ridge | 33.078 | 52.491 | 84.159 | 60.873 | 87.88% |
| No quadrupole alignment drift | shared_joint_ridge | 33.094 | 52.196 | 84.217 | 62.313 | 88.20% |
| No quadrupole alignment drift | shared_joint_random_feature | 33.139 | 52.370 | 84.479 | 62.631 | 88.28% |
| 50 um/plane RMS quadrupole alignment drift | physics_gls | 100.970 | 146.142 | 376.013 | 395.166 | 45.94% |
| 50 um/plane RMS quadrupole alignment drift | shared_target_local_ridge | 99.796 | 149.080 | 379.330 | 361.084 | 46.57% |
| 50 um/plane RMS quadrupole alignment drift | shared_joint_ridge | 99.618 | 148.595 | 377.844 | 362.513 | 46.82% |
| 50 um/plane RMS quadrupole alignment drift | shared_joint_random_feature | 99.119 | 147.288 | 371.918 | 354.570 | 46.08% |

The best learned same-distribution result without quadrupole alignment is
`shared_target_local_ridge` at
`33.078 um` RMSE.  With the paired
50-um/plane quadrupole drift enabled, the best learned result is
`shared_joint_random_feature` at
`99.119 um` RMSE.  Training the joint
ridge only on the no-alignment distribution and evaluating it on the aligned
case gives `765.859 um` RMSE; this is the explicit
distribution-shift check and must not be replaced by the in-distribution row.

Adding all-target context changes ridge RMSE relative to the matched local
shared model by `+0.046%`
without quadrupole alignment and
`-0.178%` with the 50-um
drift.  A negative value is the predeclared evidence that joint context helped;
a positive value means this ensemble does not support the added joint-model
complexity.  The strict reference gate requires aggregate RMSE, P99, and every
target-level RMSE to remain below 50 um.
The best learned no-alignment row reports
`FAIL`
for that gate; the best learned 50-um/plane row reports
`FAIL`.

## Quadrupole-drift interpretation

The requested 50-um setting is interpreted as a per-plane RMS.  The finite
draw realizes `51.086 um`
in x and `49.988 um` in y.
If the facility value instead denotes an isotropic two-dimensional radial RMS,
the corresponding per-plane RMS would be about `35.4 um`; correlated girder
motion is also a distinct prior.  Neither alternative is silently folded into
the primary 50-um/plane row.
Before any orbit correction, it changes the paired full-ring reference BPM
orbit by `860.882 um` RMS and
the beam-relative center truth by
`1335.740 um` 2D RMS.
The fraction of truths outside the maintained 1.5-mm bump radius changes from
`0.000%` to
`28.618%`.

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
  from all 76 target scans in the machine; one inference call
  returns all `152` center coordinates.
- `shared_joint_random_feature`: the joint inputs plus a fixed nonlinear tanh
  feature layer, with only the output ridge weights fitted.

The response residual, baseline-orbit, and all-target context projections are
fit on training machines only.  Model and PCA files are saved with the result
so validation/test information cannot silently enter feature construction.

## Scope and limitations

This is a synthetic SciBmad pilot, not a CESR position-precision claim.  The
sample contains only 16 independent static machines, the error
priors other than the user-provided quadrupole drift are maintained sensitivity
settings rather than measured CESR distributions, BPM white noise is
independent, and the drift is one target-local scalar mode with an independent
balanced-cycle realization for each target scan rather than one continuous
trajectory spanning all 76 scans.  There is no
actuator hysteresis, K2 polarity asymmetry, missing/outlier BPM process, or
sim-to-real validation.  The learned models return point estimates rather than
a calibrated posterior covariance or OOD probability.  Exact target orbit is
used only to form evaluation truth.

The latest lattice emits the straight-multipole-in-curved-reference warning.
No girder pitch is varied in this experiment, so the documented curved-DQX
girder-pitch discrepancy is not an excitation here, but remains part of the
lattice provenance.

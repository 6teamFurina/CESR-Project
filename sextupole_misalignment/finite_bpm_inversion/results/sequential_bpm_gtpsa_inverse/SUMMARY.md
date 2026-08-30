# Sequential BPM/GTPSA local-orbit and sextupole-center inverse

## Result

This run replaces exact target-local orbit coordinates with a machine-facing
two-sided reconstruction.  For every one of 76 sextupoles in each
of 16 paired latest-lattice machines, the estimator uses only the
nearest upstream/downstream BPM observable readbacks, known local-bump and K2
commands, and nominal order-one SciBmad/GTPSA cumulative and one-turn maps.
All target-local SciBmad orbit arrays and latent sextupole offsets are loaded
only after the BPM/GTPSA estimates have been saved, and are evaluation-only.

| acquisition | local-orbit quantity | 2D RMSE [um] | P90 [um] | P99 [um] | maximum [um] |
|---|---|---:|---:|---:|---:|
| deterministic_static_readback | relative_local_orbit_nonzero_bumps | 14.306 | 23.529 | 36.333 | 47.449 |
| deterministic_static_readback | absolute_reference_orbit | 7.866 | 6.937 | 35.456 | 88.474 |
| stochastic_15_state_means | relative_local_orbit_nonzero_bumps | 13.851 | 21.895 | 38.210 | 43.664 |
| stochastic_15_state_means | absolute_reference_orbit | 7.555 | 11.714 | 23.117 | 38.378 |

| acquisition | center method | beam-relative RMSE [um] | relative P99 [um] | absolute-offset RMSE [um] | absolute P99 [um] | bound hits |
|---|---|---:|---:|---:|---:|---:|
| deterministic_static_readback | bpm_gtpsa_two_sided | 21.658 | 65.609 | 23.346 | 72.376 | 0 |
| deterministic_static_readback | exact_local_orbit_oracle | 20.912 | 63.617 | 20.912 | 63.617 | 0 |
| stochastic_15_state_means | bpm_gtpsa_two_sided | 60.667 | 165.486 | 60.654 | 165.333 | 0 |

On deterministic static readbacks, replacing the exact local orbit changes the
beam-relative center RMSE from
`20.912 um` for the evaluation-only
oracle to `21.658 um` for the
BPM/GTPSA reconstruction.  Adding the separately reconstructed absolute
reference orbit gives an absolute sextupole-offset RMSE of
`23.346 um`.

The stochastic result uses the last 3 machines, 32 independent measurement realizations, 5.0 micrometers RMS white noise per BPM plane/read averaged over 3,072 reads, and a 10.0-micrometer endpoint-RMS scalar random walk over a repeated 15-state acquisition.  On the same three-machine subset, the deterministic beam-relative and absolute center RMSE values are 20.695 and 21.243 micrometers; the full 16-machine deterministic rows are therefore not the paired stochastic baseline.

## Method boundary

The local predictor first subtracts the zero-bump, nominal-K2 BPM readback.
Nominal control-to-BPM and control-to-target responses supply the commanded
bump prediction.  The measured-minus-model x/y residual at the nearest BPM on
each side supplies four position coordinates; the GTPSA transport maps infer
the missing upstream transverse momenta and transport the residual to the
target.  The same two-sided operator acts on the absolute BPM-minus-nominal
reference orbit to recover the target reference orbit.  The K2 slopes of all
111 BPMs are then fit against these reconstructed local coordinates while two
unknown propagation vectors are profiled out.

The exact nonlinear RF-on closed orbits remain the SciBmad forward simulation
used to synthesize BPM readbacks.  They are not inputs to the machine-facing
inverse.  The direct high-order K2/offset GTPSA-map limitation is also not
invoked: this route uses the stable order-one transport maps.

## Scope

The static physical ensemble includes all-sextupole offsets, 1% RMS local
corrector gains, 1% RMS K2 gains, 1% RMS BPM gains, independent quadrupole
strength errors within +/-1%, 1 mrad RMS quadrupole rolls, and 50-micrometer
per-plane quadrupole offsets after the recorded noisy-BPM GTPSA-ORM baseline
correction.  BPM offsets, rolls, missing/outlier channels, actuator hysteresis,
unknown model-to-machine calibration mismatch, and sim-to-real validation are
not included.  The 15-state random-walk schedule is a transparent sensitivity
model, not an optimized CESR acquisition order.

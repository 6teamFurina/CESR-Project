# CESR sextupole magnetic-center studies

New or extended calculations in this study use the repaired latest SciBmad
lattice at `Latest_Lattice/latest_cesr_scibmad_repaired.jl`.

## Maintained studies

- `direct_observable_nuisance_ablation/` contains the scan-conditioned
  bump-by-`K2` magnetic-center estimator and its nuisance/protocol ablations.
- `finite_bpm_inversion/` studies the replacement of exact target-local orbit
  coordinates by finite BPM information.
- `real_machine_nuisance_ablation/` adds BPM/corrector gain, K2 calibration,
  quadrupole strength/roll/misalignment, scan-time drift, and BPM noise one at
  a time to the finite-BPM center inverse.
- `interleaved_measurement_protocol/` compares blocked and interleaved
  `0,+,0,-,0` acquisition, repeated per-point averaging, correlated random-walk
  drift, and BPM white noise while all sextupole offsets remain hidden and
  quadrupole misalignment is disabled.
- `quadrupole_affinity/` studies quadrupole selection and nuisance leverage on
  the latest repaired lattice.
- `sextupole_misalignment_only_bpm_taylor_map/` is the isolated 76-sextupole
  misalignment benchmark comparing the maintained finite-difference source
  inverse, direct observable-derivative inverses, scan-fitted high-order
  observation Taylor maps, and a separately qualified direct GTPSA subset.
- `sextupole_excitation_validity_envelope/` expands the target-sextupole
  entrance/exit orbit in the two model-based corrector bump knobs and target
  `delta K2`, then uses 49,476 exact latest-lattice SciBmad states to report a
  signed, per-target order-two/order-four last-pass and first-fail envelope.
  The limits quantify Taylor/model validity and deliberately do not claim
  corrector, power-supply, aperture, lifetime, or operator safety limits.
- `gtpsa_derivative_stochastic_inverse/` fixes the two local `dO/dK2` source
  templates with latest-lattice SciBmad/GTPSA transport, then treats BPM white
  noise and random-walk drift through parity contrasts and analytic
  covariance.  Its all-76 noise-plus-drift-only benchmark passes the
  50-micrometer RMSE/P99 acceptance gate without invoking the failing direct
  high-order map.  A newer paired compound test that also enables BPM,
  corrector, K2, quadrupole-strength, and quadrupole-roll errors (but excludes
  quadrupole misalignment) finds no nonlinear error explosion, yet fails the
  full tail gate because quadrupole-strength sensitivity raises P99 and two
  target-level RMSEs above 50 micrometers.

## `archived_methods/`

This folder preserves the two historical studies built around the older
`cesr_model.jl` lattice and nominal/conditioned response dictionaries:

- `response_map/`: the 76-sextupole `Kn2`/offset GTPSA response map and local
  SVD baseline;
- `targeted_bump_k2_inversion/`: the P0--P3 and P1/P2 source-reconstruction
  experiments that consume the historical response map.

They remain useful for provenance and method comparison, but are not the
default starting point for new sextupole-alignment calculations.

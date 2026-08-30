# CESR sextupole magnetic-center studies

New or extended calculations in this study use the repaired latest SciBmad
lattice at `Latest_Lattice/latest_cesr_scibmad_repaired.jl`.

For a poster-style overview of the measurement protocol, finite-BPM local-orbit
inference, center inverse, stochastic-noise mitigation, results, and the
unresolved quadrupole-misalignment problem, see [`EXPLAINER.md`](EXPLAINER.md).

## Maintained studies

- `direct_observable_nuisance_ablation/` contains the scan-conditioned
  bump-by-`K2` magnetic-center estimator and its nuisance/protocol ablations.
- `finite_bpm_inversion/` studies the replacement of exact target-local orbit
  coordinates by finite BPM information.  Its full-error extension uses one
  nominal theoretical GTPSA correction ORM, two-sided order-one GTPSA local-
  orbit transport, balanced signed K2 states, and a periodic-reference random-
  walk state-space filter across 16 machines and all 76 targets.  Exact target
  orbit and every realized offset/gain/magnet error are evaluation-only.  The
  filtered reconstructed-orbit fixed-template result has 27.081-micrometer
  beam-relative center RMSE and 28.783-micrometer absolute-offset RMSE.  The
  filter reduces BPM time-state error from 2.632 to 0.320 micrometers, while
  the balanced parity contrast makes its final center-RMSE change negligible.
  The deterministic static absolute-offset RMSE is 23.232 micrometers.
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
- `sextupole_cross_response/` uses compact order-one SciBmad/GTPSA transport,
  one 62-corrector first-derivative calculation, and analytic local-sextupole
  source composition to form all 76-by-76 bump and K2--bump--center orbit
  propagation matrices without a global Hessian.  The result rejects a
  target-only observation assumption and identifies a compact shared spatial
  basis worth testing while showing that the target axis must be retained for
  a multiple-center inverse; selected paired exact SciBmad scans validate its
  sign and scale.
- `sequential_joint_inverse/` makes one complete latent machine the dataset
  unit, holds every static nuisance fixed while all 76 sextupoles are scanned
  one at a time, and compares the fixed physics inverse with shared local and
  all-target joint residual models.  Its paired cases omit or enable an
  independent 50-micrometer/plane RMS quadrupole alignment drift relative to
  the yearly nominal geometry; full-ring BPM readbacks, not oracle sextupole
  orbit, are the model inputs.  The 16-machine production pilot finds 33.078
  micrometers held-out 2D RMSE without that drift and no material all-target
  context gain; the uncorrected drift case reaches 99.119 micrometers while
  28.618% of truths leave the 1.5-millimeter excitation radius.
- `quadrupole_orbit_correction/` performs that recorded BPM-reference step on
  the same 16 paired latest-lattice machines.  It uses the 103 normal steering
  Overlay controls and only BPM readbacks plus a measured response matrix; the
  latent quadrupole offsets and sextupole-local orbit remain evaluation-only.
  A response matrix stored in the zero-offset reference state reduces the
  aggregate BPM-coordinate RMS difference from 860.882 to 126.034 micrometers
  and the target-sextupole 2D orbit difference from 1,335.740 to 70.709
  micrometers.  The fraction outside the 1.5-millimeter excitation radius
  returns from 28.618% to zero.  Holding those baseline commands fixed during
  all 76 scans reduces the best held-out center RMSE from 99.119 micrometers in
  the uncorrected-offset protocol to 33.444 micrometers, only 1.104% above the
  33.078-micrometer zero-offset benchmark.  P99 remains 84.486 micrometers, so
  the correction restores the workflow but does not pass the strict tail gate.
  A matched follow-up replaces the finite-difference ORM with the first-order
  SciBmad/GTPSA periodic closed-orbit Jacobian and adds independent noisy means
  to the stored reference and current correction orbits.  With 5-micrometer
  per-read BPM noise averaged over 3,072 reads, the best RMSE is 33.416
  micrometers; the corrected BPM orbit differs from the finite-difference,
  noiseless-correction case by only 0.077 micrometers RMS.  This establishes
  consistency at that high-repeat noise level, not low-repeat robustness.  A
  stricter production case now uses one nominal 222-by-103 GTPSA ORM with no
  realized gain/error scaling and no finite-difference ORM.  Six independently
  loaded Julia worker models generate the fixed-error one-at-a-time 76-target
  scans; the recorded thread-versus-serial orbit difference is exactly zero.
  That correction reduces BPM RMS from 860.882 to 126.037 micrometers and
  target-orbit 2D RMS from 1,335.740 to 70.945 micrometers.
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
`older_ring_version/cesr_model.jl` lattice and nominal/conditioned response dictionaries:

- `response_map/`: the 76-sextupole `Kn2`/offset GTPSA response map and local
  SVD baseline;
- `targeted_bump_k2_inversion/`: the P0--P3 and P1/P2 source-reconstruction
  experiments that consume the historical response map.

They remain useful for provenance and method comparison, but are not the
default starting point for new sextupole-alignment calculations.

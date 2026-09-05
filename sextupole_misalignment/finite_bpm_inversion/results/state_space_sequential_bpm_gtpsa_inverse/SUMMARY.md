# Full-error state-space BPM/GTPSA sextupole inverse

This latest-lattice SciBmad experiment uses 16 fixed latent
machines and all 76 sequential sextupole scans.  The baseline orbit
is corrected from observable BPM readbacks with one nominal theoretical GTPSA
ORM.  The correction response is not remeasured by central difference and is
not scaled by any realized BPM/corrector gain or other latent machine error.

The inverse process opens only pre-materialized observable BPM readings,
commands, the nominal
order-one GTPSA response/transport, and the declared stochastic priors.  Eight
balanced K2-sign/bump-sign signal states are repeated
3,072 times.  Every 256
cycles and at the endpoint, same-bump K2=0 references observe a hidden
two-plane local-orbit random walk.  Their finite 32-read
calibration errors are marginalized rather than treated as exact.  The
profiled comparison estimator supplies the optimizer with an exact analytic
variable-projection Jacobian; it does not use SciPy's numerical-difference
default.

| acquisition | inverse | beam-relative RMSE [um] | relative P99 [um] | absolute-offset RMSE [um] | absolute P99 [um] |
|---|---|---:|---:|---:|---:|
| deterministic static | noise-floor profiled BPM/GTPSA | 19.470 | 58.138 | 21.343 | 64.823 |
| balanced time series | unfiltered profiled BPM/GTPSA | 53.185 | 142.699 | 53.875 | 142.932 |
| periodic-reference time series | filtered profiled BPM/GTPSA | 53.179 | 142.685 | 53.869 | 142.845 |
| deterministic static | reconstructed-orbit fixed GTPSA template | 21.126 | 63.269 | 23.232 | 66.133 |
| balanced time series | unfiltered fixed GTPSA template | 27.081 | 69.184 | 28.783 | 75.207 |
| periodic-reference time series | filtered fixed GTPSA template | 27.081 | 69.184 | 28.783 | 75.207 |

The aggregate BPM-state deviation from the no-time-error observable states is
2.632 um before and
0.320 um after hidden-state
correction.  The finite-calibration BPM/GTPSA local-orbit RMSE is
14.407 um and the absolute reference-orbit
RMSE is 7.974 um.

The state correction is therefore active, but the fixed-template center RMSE
is unchanged at 0.001-um reporting precision because the balanced signed-state
contrast already rejects first-order drift.  White noise and static
source/model mismatch dominate the remaining center error.  The filtered
absolute aggregate RMSE passes the maintained 30-um gate, but its P99 fails
the strict 50-um tail gate.

All sextupole offsets and all realized measurement/magnet errors remain unknown
to the correction response and inverse.  Exact target-local orbits and latent
offsets enter only below the persisted-estimate boundary for these metrics.
This is a synthetic full-error SciBmad experiment, not demonstrated CESR
machine precision.  The latest lattice emits its documented straight-
multipole-in-curved-reference warning; this study does not vary girder pitch.

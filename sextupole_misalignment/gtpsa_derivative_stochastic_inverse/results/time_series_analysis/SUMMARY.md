# Eight-state time-series inverse result

The deterministic signal comes from the latest-lattice SciBmad scan for all
76 active normal sextupoles and 4 hidden
all-sextupole-offset realizations per target.  Every signal-state read has
5.0 um RMS BPM white noise.  A scalar physical
orbit drift evolves continuously as a random walk; its endpoint RMS over the
core-only eight-state scan is 10.0
um.  Adding references lengthens the equal-cadence scan, so its endpoint RMS
is conservatively increased to
10.032 um.

The eight signal states are unchanged.  Every
256 cycles (and at the final endpoint), each fixed
bump uses an interleaved `K2=0,+,0,-,0` block.  Other cycles contain only the
eight signal states.  The four K2=0 baseline means are calibrated with
32 reads each and retained as
finite-uncertainty nuisance states rather than treated as exact.

| case | 2D RMSE [um] | median [um] | P90 [um] | P99 [um] | maximum [um] |
|---|---:|---:|---:|---:|---:|
| clean | 12.761 | 8.323 | 20.627 | 38.946 | 40.779 |
| bpm_white_noise_matched_filter | 20.385 | 15.846 | 31.464 | 49.146 | 90.159 |
| balanced_8state_combined | 21.108 | 16.565 | 32.503 | 50.502 | 90.304 |
| reference_filtered_drift | 12.764 | 8.270 | 20.729 | 39.448 | 41.512 |
| reference_filtered_combined | 20.297 | 15.817 | 31.289 | 48.681 | 94.707 |

- selected signal reads/state: 3072
- balanced eight-state acquisitions/target:
  24576
- interleaved time-series acquisitions/target:
  24732
- interleaved reference cycles/target:
  13
- separate reference-calibration acquisitions/target:
  128
- drift stochastic component:
  5.757 ->
  0.319 um
- filtered worst target-level RMSE: 44.561 um
- hard gate (aggregate RMSE, P99, and every target RMSE <
  50.0 um): PASS
- preferred aggregate RMSE < 30.0 um:
  PASS
- RMSE relative to the requested 300 um scale:
  6.766%
- stochastic analysis wall time: 9.588 s

White-noise suppression is the fixed-template, covariance-matched full-BPM
GLS estimator.  With the present equal independent BPM noise this is
numerically the same matched filter as OLS, but the implementation accepts a
nonuniform diagonal or full measured BPM covariance.  Drift suppression uses
the actual acquisition order and elapsed-step random walk.  Later references
update the already accumulated center-error functional, which is equivalent
to smoothing the time series for the final center without storing the full
BPM tensor.

These are synthetic SciBmad sensitivity results, not demonstrated CESR
precision.  The drift has one calibrated spatial mode per latent machine,
the acquisition cadence is uniform, and the BPM/drift magnitudes are assumed
rather than measured.  Machine deployment still requires measured BPM
covariance and cadence, calibrated drift modes, K2/corrector readbacks,
settling masks, missing/outlier BPM handling, and tests with multidirectional
drift.

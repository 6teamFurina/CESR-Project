# Compound nuisance revalidation with the eight-state inverse

All deterministic states use the current latest-lattice SciBmad protocol:
signed +/-1.500 mm local bumps and delta-K2 extrema
-0.100/0.100 m^-3.  The paired benchmark covers
all 76 targets and 4 hidden machines per target.
Quadrupole misalignment is deliberately excluded.  The fully combined case
simultaneously activates 1% BPM gain, 1% corrector gain, 1% K2 gain,
independent +/-1% quadrupole strength, 1 mrad RMS quadrupole roll, 5 um/read
BPM white noise, and the maintained 10 um endpoint random-walk drift.

| case | 2D RMSE [um] | increment vs clean reference RMS [um] | P90 [um] | P99 [um] |
|---|---:|---:|---:|---:|
| Reference: sextupole offsets only | 12.761 | 0.000 | 20.627 | 38.946 |
| BPM gain | 12.752 | 0.369 | 20.631 | 38.919 |
| Corrector gain | 14.476 | 6.421 | 23.939 | 38.265 |
| K2 calibration gain | 13.394 | 3.016 | 22.196 | 39.134 |
| Quadrupole strength | 23.825 | 19.384 | 36.008 | 69.750 |
| Quadrupole roll | 13.953 | 5.172 | 22.613 | 37.016 |
| Linear sum of paired increments | 25.873 | 21.443 | 40.176 | 68.293 |
| All static nuisances combined | 25.897 | 21.458 | 41.059 | 65.552 |
| Reference + white noise + drift inverse | 20.332 | 15.803 | 31.315 | 48.714 |
| All nuisances except quadrupole misalignment | 30.334 | 26.640 | 47.316 | 73.738 |

The component-matched deterministic decomposition gives:

- quadrature of five component increments:
  21.283 um
- vector sum of paired component increments:
  21.443 um
- actual compound increment:
  21.458 um
- nonlinear compound interaction RMS:
  1.565 um
- nonlinear interaction P99:
  4.834 um
- nonlinear interaction maximum:
  7.174 um at
  SEX_17W, realization
  1
- nonlinear interaction / actual increment:
  7.293%

For the complete stochastic case:

- aggregate 2D RMSE: 30.334 um
- P99: 73.738 um
- worst target-level RMSE: 66.560 um
- targets at or above 50 um RMSE:
  2 (SEX_09AW, SEX_38E)
- compound white-noise component: 15.814 um
- compound filtered-drift component: 0.324 um
- compound total stochastic component: 15.817 um
- clean-reference total stochastic component: 15.817 um
- proxy relative RMSE on a 300 um scale:
  10.111%
- hard aggregate/P99/all-target 50 um gate: FAIL
- preferred aggregate 30 um gate: FAIL

BPM gain is applied to both signal and reference readbacks but is not supplied
to the nominal center template.  Corrector and K2 gains likewise alter the
physical SciBmad scan while the inverse uses commanded spans.  The compound
drift response is independently recovered from a paired exact SciBmad secant
with identical static nuisance realizations.  White noise and random-walk
histories are then propagated through the repeated acquisition sequence.

These remain synthetic sensitivities, not measured CESR error priors.  The
drift basis is one calibrated scalar mode, BPM white noise is temporally
independent, actuator hysteresis and polarity asymmetry are absent, and the
nuisance ensemble has only four hidden machines per target.

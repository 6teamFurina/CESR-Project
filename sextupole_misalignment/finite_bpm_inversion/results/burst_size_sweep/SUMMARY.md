# Full-error burst-size sweep

This paired latest-lattice SciBmad study changes only the number of
consecutive BPM turns acquired during one fixed magnet-state visit.
All 16 latent machines, 76 sextupoles,
32 stochastic realizations, 3,072 turns per signed state,
13 reference cycles, finite calibration, nominal GTPSA transport,
state-space filtering, and the fixed-template inverse are retained.

| burst | state visits / target | visit reduction | absolute RMSE [um] | RMSE delta vs B=1 [um] | absolute P99 [um] | worst-target RMSE [um] | BPM-state RMSE [um] |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 24732 | 1.00x | 28.783 | +0.000 | 75.207 | 72.529 | 0.320 |
| 2 | 12444 | 1.99x | 28.782 | -0.001 | 75.046 | 72.530 | 0.320 |
| 4 | 6300 | 3.93x | 28.796 | +0.013 | 75.153 | 72.540 | 0.320 |
| 8 | 3228 | 7.66x | 28.809 | +0.026 | 74.976 | 72.528 | 0.319 |
| 16 | 1692 | 14.62x | 28.867 | +0.084 | 75.199 | 71.902 | 0.317 |
| 32 | 924 | 26.77x | 29.205 | +0.422 | 75.871 | 72.549 | 0.315 |
| 64 | 540 | 45.80x | 30.289 | +1.506 | 77.608 | 72.598 | 0.311 |
| 128 | 348 | 71.07x | 33.773 | +4.990 | 87.981 | 72.812 | 0.307 |

All burst rows contain the same 24,576 signal turns and 156 periodic
reference observations per target. Any wall-time reduction therefore
comes from fewer physical state visits. The current random-walk model
advances per BPM acquisition and contains no corrector or sextupole
settling interval; a facility-time conclusion requires measured settling
and drift spectra.

`B=16` reduces visits by 14.62x while changing filtered absolute RMSE by only +0.084 micrometers and P99 by -0.008 micrometers. It is the conservative modeled operating candidate. `B=32` reduces visits by 26.77x, with RMSE and P99 penalties of +0.422 and +0.664 micrometers; it is an aggressive candidate requiring target-wise and machine-data checks. `B=64` reaches 30.289 micrometers and crosses the 30-micrometer aggregate gate, so it is not the current default.

The scalar filtered BPM-state RMSE decreases slightly with burst size,
but final center error rises from `B=32` onward. Consecutive-state
clustering changes how residual drift projects onto the signed parity
contrasts used by the inverse, so BPM-state RMS alone is not a valid
selection metric. The baseline P99 and worst-target RMSE already exceed
the stricter 50-micrometer tail gate; burst acquisition does not resolve
that pre-existing model/template limitation.

The latest lattice emits the documented straight-multipole-in-curved-
reference warning. This study does not vary girder pitch.

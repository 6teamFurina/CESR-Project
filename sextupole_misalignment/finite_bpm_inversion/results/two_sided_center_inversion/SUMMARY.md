# Two-sided-BPM sextupole-center inversion

This end-to-end experiment replaces exact target-local coordinates with the
relative local orbit predicted from the nearest upstream/downstream BPM pair.
The maintained K2-slope center fit then consumes all 111 BPM x/y response
channels. Exact target-local orbit and target offset are evaluation-only.

- targets / latent realizations: 76 / 8 per target
- completed center fits: 608
- bump protocol: five-point axial cross at `0.5 mm`
- K2 protocol: `(-0.02, 0, +0.02) m^-3`
- target offset: independent x/y uniform within `+/-350 micrometers`
- other 75 sextupoles: independent x/y `300 micrometer` RMS offsets
- all 113 quadrupoles: independent physical strength errors within `+/-1%`
- BPM noise/offset/gain errors and missing channels: none
- internal target orbit used by the two-sided inverse: **no**
- aggregate x/y/2D center RMSE: **3.826 / 4.444 / 5.864 micrometers**
- aggregate median / P90 / P99 / maximum: **3.645 / 9.418 / 17.897 / 25.270 micrometers**
- per-target RMSE median / P90 / maximum: **4.325 / 7.968 / 17.569 micrometers**
- correlation between per-case local-orbit prediction RMSE and center error:
  **-0.014288**
- two-sided versus oracle center-error-vector RMS / median / P90 / maximum:
  **0.192 / 0.025 / 0.181 / 1.818 micrometers**

For context, the same frozen tensor gave `13.913 micrometers` beam-relative
2D RMSE when commanded bump coordinates were used directly, while the frozen
oracle-local-orbit fit gave `5.870 micrometers`. The oracle value is the same
under an exact common coordinate translation because the source fit depends
only on local orbit minus center.

## Ten largest per-target center RMSE values

| rank | target | RMSE [micrometers] | P90 [micrometers] | max [micrometers] |
|---:|---|---:|---:|---:|
| 1 | SEX_13E | 17.569 | 24.757 | 25.270 |
| 2 | SEX_24W | 13.385 | 16.770 | 20.125 |
| 3 | SEX_28E | 12.249 | 19.403 | 20.098 |
| 4 | SEX_17W | 12.032 | 15.254 | 15.966 |
| 5 | SEX_21W | 9.954 | 12.630 | 13.691 |
| 6 | SEX_44W | 9.045 | 12.981 | 13.020 |
| 7 | SEX_17E | 8.656 | 11.884 | 13.205 |
| 8 | SEX_39E | 8.052 | 11.494 | 11.898 |
| 9 | SEX_12E | 7.885 | 10.052 | 10.138 |
| 10 | SEX_18E | 7.818 | 9.914 | 10.678 |

The result remains a noise-free SciBmad study. It tests propagation of the
two-sided local-orbit estimate through the physical center inverse, not real-
machine BPM or corrector calibration accuracy.

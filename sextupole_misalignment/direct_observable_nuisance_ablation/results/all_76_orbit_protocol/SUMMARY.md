# All-76 economical orbit protocol

The latest repaired SciBmad CESR lattice was evaluated for all 76 active
normal sextupoles. Each target has 8 independent latent realizations and
uses five axial-cross bumps, three K2 levels `(-2K,0,+2K)`, and nominal K1
commands. Every tensor contains unknown offsets on the other 75 sextupoles and
independent physical quadrupole strength errors bounded by ±1%. Target-local
orbit coordinates are exact and no BPM noise is added.

- exact SciBmad states: 9120
- aggregate fits: 608
- aggregate 2D RMSE: **5.870 µm**
- aggregate median / P90 / P99: **3.664 / 9.395 / 17.897 µm**
- maximum realization error: **25.274 µm**
- per-target RMSE median / P90 / maximum: **4.300 / 7.933 / 17.574 µm**
- realized other-sextupole x/y RMS: **301.142 / 300.103 µm**
- realized quadrupole-error range: **-1.0000% to 1.0000%**
- SciBmad generation wall time: **553.6 s**

## Ten largest per-target RMSE values

| rank | target | RMSE [µm] | P90 [µm] | max [µm] |
|---:|---|---:|---:|---:|
| 1 | SEX_13E | 17.574 | 24.762 | 25.274 |
| 2 | SEX_24W | 13.436 | 16.873 | 20.086 |
| 3 | SEX_28E | 12.242 | 19.390 | 20.087 |
| 4 | SEX_17W | 12.018 | 15.195 | 15.877 |
| 5 | SEX_21W | 9.990 | 12.706 | 13.821 |
| 6 | SEX_44W | 9.131 | 13.121 | 13.121 |
| 7 | SEX_17E | 8.662 | 11.913 | 13.231 |
| 8 | SEX_39E | 8.044 | 11.479 | 11.900 |
| 9 | SEX_18E | 7.823 | 9.922 | 10.679 |
| 10 | SEX_12E | 7.674 | 9.484 | 9.540 |

These values measure the present noise-free shared thin-sextupole source fit.
They do not include target-local-orbit uncertainty, BPM noise, missing BPMs, or
measured covariance, and therefore are not predicted machine accuracy.

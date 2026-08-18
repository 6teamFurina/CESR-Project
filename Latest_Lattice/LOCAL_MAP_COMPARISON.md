# Latest CESR local-map comparison

- Bmad/Tao reference: `20260814-0`, branch 0, RF on.
- Compared elements: `1177`.
- Name mismatches after `#` to `!s` normalization: `1`.
- Maximum element-length mismatch: `0.000000e+00 m`.
- Maximum local-matrix entry mismatch: `1.139998e-04`.
- Maximum affine-vector mismatch: `1.139996e-04`.

## Largest local matrix discrepancies

| Index | Bmad name | Bmad key | SciBmad kind | max abs dR | relative Frobenius | entry |
|---:|---|---|---|---:|---:|---|
| 997 | `ID_S1A#1` | `Wiggler` | `Wiggler` | 1.139998e-04 | 5.864555e-05 | R16 |
| 999 | `ID_S1A#2` | `Wiggler` | `Wiggler` | 1.139995e-04 | 6.490548e-05 | R16 |
| 1001 | `ID_S1A#3` | `Wiggler` | `Wiggler` | 2.944379e-06 | 9.985692e-07 | R12 |
| 22 | `DQX4B` | `SBend` | `SBend` | 4.289902e-13 | 9.841753e-14 | R12 |
| 42 | `DQX4D` | `SBend` | `SBend` | 4.289902e-13 | 9.841753e-14 | R12 |
| 72 | `DQX5B` | `SBend` | `SBend` | 4.289902e-13 | 9.841753e-14 | R12 |
| 92 | `DQX5D` | `SBend` | `SBend` | 4.289902e-13 | 9.841753e-14 | R12 |
| 118 | `DQX6B` | `SBend` | `SBend` | 4.289902e-13 | 9.841753e-14 | R12 |
| 138 | `DQX6D` | `SBend` | `SBend` | 4.289902e-13 | 9.841753e-14 | R12 |
| 1014 | `DQX1B` | `SBend` | `SBend` | 4.289902e-13 | 9.841753e-14 | R12 |
| 1034 | `DQX1D` | `SBend` | `SBend` | 4.289902e-13 | 9.841753e-14 | R12 |
| 1074 | `DQX2B` | `SBend` | `SBend` | 4.289902e-13 | 9.841753e-14 | R12 |
| 1094 | `DQX2D` | `SBend` | `SBend` | 4.289902e-13 | 9.841753e-14 | R12 |
| 1136 | `DQX3B` | `SBend` | `SBend` | 4.289902e-13 | 9.841753e-14 | R12 |
| 1156 | `DQX3D` | `SBend` | `SBend` | 4.289902e-13 | 9.841753e-14 | R12 |
| 1000 | `LS_1A1_1A2_1A3` | `Photon_Fork` | `Marker` | 1.423769e-14 | 5.812512e-15 | R56 |
| 553 | `B47AW` | `SBend` | `SBend` | 5.773160e-15 | 1.158787e-15 | R12 |
| 617 | `B48E#2` | `SBend` | `SBend` | 5.773160e-15 | 1.680129e-15 | R34 |
| 623 | `B47AE` | `SBend` | `SBend` | 5.773160e-15 | 1.158787e-15 | R12 |
| 545 | `B47W` | `SBend` | `SBend` | 5.329071e-15 | 1.044170e-15 | R12 |
| 633 | `B47E` | `SBend` | `SBend` | 5.329071e-15 | 1.044170e-15 | R12 |
| 331 | `Q23W` | `Quadrupole` | `Quadrupole` | 5.107026e-15 | 2.790610e-15 | R11 |
| 909 | `Q15E` | `Quadrupole` | `Quadrupole` | 5.107026e-15 | 2.793897e-15 | R11 |
| 104 | `D043#1` | `Drift` | `Drift` | 4.884981e-15 | 1.750082e-15 | R12 |
| 106 | `D043#2` | `Drift` | `Drift` | 4.884981e-15 | 1.761441e-15 | R12 |

## Largest affine-vector discrepancies

| Index | Bmad name | Bmad key | SciBmad kind | max abs dvec0 | exit-orbit mismatch |
|---:|---|---|---|---:|---:|
| 997 | `ID_S1A#1` | `Wiggler` | `Wiggler` | 1.139996e-04 | 1.139996e-04 |
| 999 | `ID_S1A#2` | `Wiggler` | `Wiggler` | 1.139993e-04 | 1.139993e-04 |
| 1000 | `LS_1A1_1A2_1A3` | `Photon_Fork` | `Marker` | 1.962917e-06 | 1.962917e-06 |
| 1001 | `ID_S1A#3` | `Wiggler` | `Wiggler` | 9.814586e-07 | 9.814586e-07 |
| 234 | `B13W` | `SBend` | `SBend` | 3.552714e-15 | 3.552714e-15 |
| 256 | `B15W` | `SBend` | `SBend` | 3.552714e-15 | 3.552714e-15 |
| 523 | `B44W` | `SBend` | `SBend` | 3.552714e-15 | 3.552714e-15 |
| 539 | `B46W` | `SBend` | `SBend` | 3.552714e-15 | 3.552714e-15 |
| 639 | `B46E` | `SBend` | `SBend` | 3.552714e-15 | 3.552714e-15 |
| 655 | `B44E` | `SBend` | `SBend` | 3.552714e-15 | 3.552714e-15 |

## Split-wiggler interpretation

The three `ID_S1A` slices must not be interpreted independently: the continuous four-potential uses canonical boundary terms that cancel across the complete five-element block. Run `compare_wiggler_block.jl` for the physical block comparison. With the Bmad reference-time patch, the block affine/exit mismatch is about `7.12e-15`; the remaining `5.89e-6` `R12` difference is reproduced by Bmad Runge-Kutta field tracking and is absent only from Bmad's standard wiggler approximation.

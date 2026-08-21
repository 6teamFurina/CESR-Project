# Bmad–SciBmad CESR control-response comparison (RF-ON)

The response matrix maps 119 Bmad-compatible CESR control Overlay knobs to horizontal and vertical closed orbit at 99 `DET_*` markers. SciBmad computes all control derivatives simultaneously with a first-order GTPSA parameter map; no control finite differences are used.

- Matrix shape: `198 x 119`
- Units: `m/rad`
- Closed-orbit model: `6D RF-confined beam`
- SciBmad baseline closed orbit: `[-1.668519860e-05, 2.390112503e-03, 1.054311223e-06, 1.796930342e-06, -3.993611670e-04, -7.803721196e-06]`
- GTPSA closure-equation residual: `7.105427e-15`
- Runtime: `44.444 s`

## Overall agreement

- Relative Frobenius difference: `2.032285710e-03` (`0.203229%`)
- Maximum absolute entry difference: `6.747159290e-02 m/rad`
- Maximum difference normalized by the Bmad matrix maximum: `2.869110121e-03` (`0.286911%`)
- Full-matrix cosine correlation: `0.999997991471`
- Singular-value relative 2-norm difference: `1.356418312e-03`
- Worst entry: `DET_37W:y` versus `V37W`
- Worst column relative 2-norm: `3.585490324e-03` for `V05E`

## Plane blocks

| Block | Relative Frobenius | Max absolute difference (m/rad) | Max-normalized difference | Correlation |
|---|---:|---:|---:|---:|
| xH | 4.863691725e-04 | 9.701908092e-03 | 4.682599538e-04 | 0.999999883348 |
| xV | 4.989785932e-03 | 1.632960501e-03 | 2.754215569e-03 | 0.999987636600 |
| yH | 1.861519971e-03 | 1.061189099e-03 | 9.573234053e-04 | 0.999998467461 |
| yV | 2.495793330e-03 | 6.747159290e-02 | 2.869110121e-03 | 0.999997004037 |

## Leading singular values

| Index | Bmad | SciBmad | Difference |
|---:|---:|---:|---:|
| 1 | 3.463845586e+02 | 3.472920405e+02 | 9.074818841e-01 |
| 2 | 2.826703719e+02 | 2.821896717e+02 | -4.807002240e-01 |
| 3 | 2.791338868e+02 | 2.791331102e+02 | -7.765866646e-04 |
| 4 | 2.734030868e+02 | 2.734354674e+02 | 3.238060377e-02 |
| 5 | 2.360411106e+02 | 2.360528709e+02 | 1.176034254e-02 |
| 6 | 1.895053067e+02 | 1.895113176e+02 | 6.010954684e-03 |
| 7 | 1.591911379e+02 | 1.592121421e+02 | 2.100423326e-02 |
| 8 | 1.315240600e+02 | 1.315026155e+02 | -2.144455441e-02 |
| 9 | 1.015229092e+02 | 1.015332557e+02 | 1.034655249e-02 |
| 10 | 9.684016579e+01 | 9.684129072e+01 | 1.124927242e-03 |
| 11 | 8.022142241e+01 | 8.022074269e+01 | -6.797278143e-04 |
| 12 | 7.340932246e+01 | 7.340838182e+01 | -9.406435151e-04 |

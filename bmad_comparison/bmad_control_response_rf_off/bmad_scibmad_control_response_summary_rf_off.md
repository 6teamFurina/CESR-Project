# Bmad–SciBmad CESR control-response comparison (RF-OFF)

The response matrix maps 119 Bmad-compatible CESR control Overlay knobs to horizontal and vertical closed orbit at 99 `DET_*` markers. SciBmad computes all control derivatives simultaneously with a first-order GTPSA parameter map; no control finite differences are used.

- Matrix shape: `198 x 119`
- Units: `m/rad`
- Closed-orbit model: `4D coasting beam with fixed z and pz`
- Baseline solver: `six-direction ForwardDiff 4 x 4 Jacobian adapter + SciBmad BatchSolve.newton!`
- SciBmad baseline closed orbit: `[-1.663730483e-05, 2.389857926e-03, 1.065879859e-06, 1.782430008e-06, 0.000000000e+00, 0.000000000e+00]`
- GTPSA closure-equation residual: `7.105427e-15`
- Runtime: `41.676 s`

## Overall agreement

- Relative Frobenius difference: `2.012009559e-03` (`0.201201%`)
- Maximum absolute entry difference: `6.726326268e-02 m/rad`
- Maximum difference normalized by the Bmad matrix maximum: `2.860669026e-03` (`0.286067%`)
- Full-matrix cosine correlation: `0.999998031423`
- Singular-value relative 2-norm difference: `1.344825645e-03`
- Worst entry: `DET_37W:y` versus `V37W`
- Worst column relative 2-norm: `3.542456655e-03` for `V05E`

## Plane blocks

| Block | Relative Frobenius | Max absolute difference (m/rad) | Max-normalized difference | Correlation |
|---|---:|---:|---:|---:|
| xH | 4.779468679e-04 | 9.651253882e-03 | 4.669368153e-04 | 0.999999887334 |
| xV | 4.991117850e-03 | 1.595391846e-03 | 2.691981276e-03 | 0.999987625217 |
| yH | 1.852555466e-03 | 1.050569292e-03 | 9.481541443e-04 | 0.999998480541 |
| yV | 2.471638363e-03 | 6.726326268e-02 | 2.860669026e-03 | 0.999997062105 |

## Leading singular values

| Index | Bmad | SciBmad | Difference |
|---:|---:|---:|---:|
| 1 | 3.462975926e+02 | 3.471963919e+02 | 8.987993517e-01 |
| 2 | 2.826903148e+02 | 2.822117818e+02 | -4.785329428e-01 |
| 3 | 2.791944789e+02 | 2.791987105e+02 | 4.231575444e-03 |
| 4 | 2.733388024e+02 | 2.733700476e+02 | 3.124518290e-02 |
| 5 | 2.359996623e+02 | 2.360118959e+02 | 1.223358848e-02 |
| 6 | 1.894837390e+02 | 1.894879703e+02 | 4.231363973e-03 |
| 7 | 1.590537358e+02 | 1.590748782e+02 | 2.114242570e-02 |
| 8 | 1.314223077e+02 | 1.314003459e+02 | -2.196183702e-02 |
| 9 | 1.015198946e+02 | 1.015299989e+02 | 1.010429436e-02 |
| 10 | 9.683809713e+01 | 9.683905395e+01 | 9.568134896e-04 |
| 11 | 8.022104485e+01 | 8.022032952e+01 | -7.153264273e-04 |
| 12 | 7.341191115e+01 | 7.341118643e+01 | -7.247176451e-04 |

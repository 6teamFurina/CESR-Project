# SciBmad `6 x 119` response-initial-guess benchmark

Date: 2026-07-30  
Machine: local Windows, AMD Ryzen 9 5900HX, Julia 1.12.6, one Julia thread  
Workload: the same 1000 corrector samples and `1000 x 198` detector output  
Tolerances: `reltol=1e-8`, `abstol=1e-10`

## Method

At the zero-control nominal closed orbit `z0`, GTPSA calculates the one-turn
map derivatives

```text
A = dF/dz
B = dF/dk
R = dz_closed/dk = (I - A)^(-1) B
```

where `R` has shape `6 x 119`. Sample `i` starts modified Newton from its own
linear prediction

```text
z_initial[i] = z0 + R * delta_k[i].
```

The nominal `6 x 6` residual Jacobian is LU-factorized once and reused for all
samples and iterations. Every final one-turn closure is checked; any failed
lane is automatically rerun with the full-AD Newton solver.

## Results

| Initial guess | Jacobian | Physics (s) | Solve (s) | Track (s) | Samples/s | Iterations min / median / mean / max | Closure max | Fallbacks |
|---|---|---:|---:|---:|---:|---|---:|---:|
| fixed nominal `z0` | frozen nominal + full-AD fallback | 8.163 | 7.601 | 0.560 | 122.506 | 0 / 3 / 2.994 / 3 | `9.802e-11` | 0 |
| per-sample `z0 + R*delta_k` | frozen nominal + full-AD fallback | 6.855 | 6.266 | 0.587 | 145.885 | 0 / 2 / 1.995 / 2 | `8.104e-11` | 0 |

The response predictor reduced the recurring physics time by `16.0%`
(`1.191x` throughput speedup) in these runs and removed approximately one
modified-Newton update per nonzero sample.

## Response construction and cold start

- parameterized GTPSA model setup: `0.156 s`
- one-turn map and `6 x 119` response solve: `2.232 s`
- total post-compilation response construction: `2.389 s`
- first-process warmup/compilation: `100.764 s`
- response-equation residual maximum: `7.105e-15`

Consequently, the result supports a recurring-batch digital-twin claim when
the nominal response is cached and reused. It does not support a cold-start
speedup claim for a process that constructs the matrix and generates only one
1000-sample batch.

## Numerical checks

- converged: `1000 / 1000`
- full-AD fallbacks: `0`
- detector RMSE versus fixed-`z0` frozen result: `9.411e-12 m`
- detector maximum difference versus fixed-`z0` frozen result:
  `5.563e-10 m`
- detector RMSE versus Bmad: `2.268158e-6 m`
- detector maximum difference versus Bmad: `8.138494e-6 m`

The complete response matrix is stored in
`scibmad_response_initial_frozen_fallback_bmad_tolerance/closed_orbit_response_6x119.csv`.

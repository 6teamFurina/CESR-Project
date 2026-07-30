# SciBmad frozen-Jacobian experiment

Date: 2026-07-30  
Machine: local Windows PC, AMD Ryzen 9 5900HX  
Julia: 1.12.6, one Julia thread  
Dataset: the same 1,000 CESR corrector samples, 119 controls, RF on  
Tolerances: `reltol = 1e-8`, `abstol = 1e-10`  
Initial guess: the nominal closed orbit `z0` for every sample

## Method

The full-Newton reference recomputes the closed-orbit residual Jacobian with
automatic differentiation at each iteration. The experimental solver:

1. solves the nominal closed orbit;
2. extracts its `6 x 6` residual Jacobian;
3. computes one pivoted LU factorization;
4. reuses that factorization for every sample and modified-Newton iteration;
5. recomputes the full nonlinear one-turn residual at every iteration;
6. independently checks final one-turn closure while collecting detector data.

This freezes only the derivative. It does not replace nonlinear CESR tracking
with a linear model.

## Controlled result

| Solver | Converged | Closed-orbit solve | Solve + detector tracking | Throughput | Iterations (min / median / mean / max) |
|---|---:|---:|---:|---:|---:|
| Full AD Jacobian each iteration | 1000/1000 | 21.655 s | 22.247 s | 44.950 samples/s | 0 / 2 / 1.998 / 2 |
| Frozen nominal Jacobian | 1000/1000 | 6.986 s | 7.535 s | 132.709 samples/s | 0 / 3 / 2.994 / 3 |

With the same initial guess and tolerances, frozen Jacobian is `2.952x` faster
in the timed physics region. It needs one extra modified-Newton iteration, but
ordinary residual tracking is much cheaper than rebuilding the full batched AD
Jacobian.

The one-time LU factorization took `6.6e-6 s`. Its cost is negligible. The
speedup comes from eliminating repeated full-ring derivative tracking, not from
making the `6 x 6` linear solve itself faster.

The nominal Jacobian condition number was approximately `2247`, which was
manageable for this dataset.

## Closure and numerical agreement

The final six-coordinate one-turn residual norms were:

- median: `4.107e-12`
- maximum: `9.802e-11`

The maximum is below the configured SciBmad absolute residual tolerance of
`1e-10`.

Relative to the full-AD, nominal-`z0`, low-tolerance result:

- compared detector coordinates: 198,000
- convergence flags: identical
- maximum absolute detector-orbit difference: `5.566e-10 m`
- detector-orbit RMSE: `7.639e-12 m`

Relative to the existing Bmad output, the frozen-Jacobian dataset retains the
same practical agreement:

- global RMSE: `2.268e-6 m`
- maximum absolute difference: `8.138e-6 m`
- correlation: `0.999999966415`
- median per-sample relative 2-norm difference: `0.0338%`

## Interpretation and limitation

Corrector changes primarily alter the affine kick term of the one-turn map.
For a linear lattice, the orbit Jacobian is exactly shared by all corrector
settings. CESR sextupoles and RF make the map nonlinear, but these sampled
corrector changes are small enough that one nominal Jacobian remains accurate.

The current implementation is experimental and calls private SciBmad residual
helpers. It now checks every final closure and reruns only failed lanes in a
full-AD Newton sub-batch. The fallback-enabled 1,000-sample repeat used
`8.163 s`, converged 1000/1000, and required zero fallbacks. A forced
single-lane test triggered the duplicate-lane path required by BatchParam,
recovered 1/1 lanes, and produced a final closure norm of `7.439e-12`.
A production implementation should still expose a supported residual API.

The existing remote Bmad result used 67.370 s while this local frozen-Jacobian
run used 7.535 s. The resulting `8.941x` throughput ratio is cross-machine and
must not be presented as a controlled Bmad/SciBmad speedup. The same run should
be repeated on `lnx201`.

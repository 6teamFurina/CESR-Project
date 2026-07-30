# SciBmad closed-orbit tolerance comparison

Date: 2026-07-30  
Machine: local Windows PC, AMD Ryzen 9 5900HX  
Julia: 1.12.6, one Julia thread  
Dataset: the same 1,000 CESR corrector samples, 119 controls, RF on  
Initial guess: six-dimensional zero vector

## Controlled SciBmad result

| Tolerances `(reltol, abstol)` | Converged | Closed-orbit solve | Solve + detector tracking | Throughput | Newton iterations (min / median / mean / max) |
|---|---:|---:|---:|---:|---:|
| `(1e-13, 1e-13)` | 1000/1000 | 63.786 s | 64.356 s | 15.539 samples/s | 3 / 3 / 4.034 / 12 |
| `(1e-8, 1e-10)` | 1000/1000 | 25.894 s | 26.457 s | 37.798 samples/s | 3 / 3 / 3.000 / 3 |

Matching Bmad's default numerical tolerance values reduced the SciBmad timed
physics region by 58.9%, a `2.432x` speedup. The closed-orbit solve alone was
`2.463x` faster. All samples now stop after three Newton iterations rather than
letting the full batch continue until its former maximum of 12.

The convergence criteria are not mathematically identical: Bmad checks
component-wise one-turn closure, whereas SciBmad checks residual or step norms.
This experiment therefore aligns the tolerance values, not the exact stopping
rule.

## Effect on the calculated dataset

Comparing the lower-tolerance result with the original `1e-13` SciBmad result:

- compared observable entries: 198,000
- convergence flags: identical
- maximum absolute orbit difference: `1.608e-13 m`
- orbit RMSE: `1.533e-14 m`

The relaxed tolerance produced no practically meaningful change in this
dataset. The additional iterations required by `1e-13` were refining the
fixed-point residual near floating-point limits without changing the detector
orbits at a relevant scale.

## Context against the existing Bmad run

The existing Bmad result used 67.370 s (`14.843 samples/s`) on
`lnx201.classe.cornell.edu`. The lower-tolerance local SciBmad result used
26.457 s (`37.798 samples/s`), giving an observed cross-machine throughput
ratio of `2.546x` in SciBmad's favor. This is not a controlled speedup claim:
the engines ran on different machines and Bmad's component-wise convergence
rule differs from SciBmad's norm-based rule.

Numerical agreement with Bmad is unchanged to the displayed precision:

- global RMSE: `2.268e-6 m`
- maximum absolute difference: `8.138e-6 m`
- correlation: `0.999999966415`
- median per-sample relative 2-norm difference: `0.0338%`

The next defensible comparison is to run both engines repeatedly on `lnx201`
with the same tolerance values, one thread, CPU accounting, and controlled
host load.

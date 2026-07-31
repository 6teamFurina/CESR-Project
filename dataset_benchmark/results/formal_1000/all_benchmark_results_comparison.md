# CESR 1,000-sample benchmark: unified results

Date: 2026-07-30  
Workload: 119 correctors to horizontal and vertical closed orbit at 99
detectors (`1000 x 198` output table), RF on.

## All requested results in one table

| Result | Machine | Tolerances `(rel, abs)` | Initial guess and Jacobian strategy | Converged | Warmup (s) | Physics (s) | Component total (s) | Physics samples/s | Observed rate / Bmad | RMSE vs Bmad (m) |
|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|
| Bmad/Tao/PyTao | `lnx201` Linux | defaults `(1e-8, 1e-10)`* | previous orbit; one-turn matrix reused and conditionally refreshed | 1000/1000 | 0.295 | 67.370 | 69.238 | 14.843 | 1.000 | 0 |
| SciBmad, VM high precision | `lnx201` Linux | `(1e-13, 1e-13)` | zero; full batched AD Jacobian each Newton iteration | 1000/1000 | 439.887 | 280.486 | 722.564 | 3.565 | 0.240 | `2.268158e-6` |
| SciBmad, VM response initial guess + frozen Jacobian + fallback | `lnx201` Linux | `(1e-8, 1e-10)` | per-sample `z0 + (dz/dk) delta-k`; cached response; nominal Jacobian reused; failed lanes use full AD | 1000/1000 | 386.686 | 33.178 | 423.495**** | 30.140 | 2.031 | `2.268158e-6` |
| SciBmad, local high precision | Windows, Ryzen 9 5900HX | `(1e-13, 1e-13)` | zero; full batched AD Jacobian each Newton iteration | 1000/1000 | 74.790 | 64.356 | 139.962 | 15.539 | 1.047 | `2.268158e-6` |
| SciBmad, local normal precision | Windows, Ryzen 9 5900HX | `(1e-8, 1e-10)` | zero; full batched AD Jacobian each Newton iteration | 1000/1000 | 75.016 | 26.457 | 102.270 | 37.798 | 2.546 | `2.268158e-6` |
| SciBmad, local frozen Jacobian + fallback | Windows, Ryzen 9 5900HX | `(1e-8, 1e-10)` | nominal `z0`; one nominal Jacobian reused; failed lanes use full AD | 1000/1000 | 59.235 | 8.163 | 68.123** | 122.506 | 8.253 | `2.268158e-6` |
| SciBmad, response initial guess + frozen Jacobian + fallback | Windows, Ryzen 9 5900HX | `(1e-8, 1e-10)` | per-sample `z0 + (dz/dk) delta-k`; nominal Jacobian reused; failed lanes use full AD | 1000/1000 | 100.764 | 6.855 | 110.602*** | 145.885 | 9.829 | `2.268158e-6` |

\* Bmad and SciBmad use different mathematical stopping rules even when the
two tolerance values match. Bmad checks component-wise one-turn closure;
SciBmad checks residual or step norms.

\** The frozen-Jacobian component total includes nominal model setup and
nominal closed-orbit/Jacobian calculation (`0.216 s`), regular model setup,
warmup, physics, and writing. Other totals are the sums of the components
available in their metadata. They are not external wall-clock measurements.

\*** The response-initialized component total additionally includes construction
of the parameterized GTPSA model and its `6 x 119` closed-orbit response
matrix (`2.389 s` after compilation). Its `100.764 s` warmup includes
first-process compilation of this GTPSA path. The matrix can be reused by a
long-lived digital twin; the `6.855 s` physics time is the recurring-batch
quantity, not the cold-start total.

\**** The optimized `lnx201` component total includes warmup, nominal model and
closed orbit, cached response loading (`0.000730 s`), model setup, physics, and
writing. GNU `time` measured `8:56.57` total wall time because Julia startup
and compilation outside the recurring physics region remained substantial.

## Interpretation

Only the four local SciBmad rows are same-machine comparisons. The
`Observed rate / Bmad` column is included to place all requested results in one
table, but ratios involving local Windows SciBmad and Linux Bmad are
cross-machine and are not controlled speedup claims.

The Bmad and VM high-precision SciBmad rows share `lnx201`, but the runs were
made at different times on a shared host. GNU `time` reported only 49% average
CPU utilization for the VM SciBmad run, so its `0.240` rate ratio should also
be repeated under controlled load.

The optimized SciBmad `lnx201` run used `33.178 s` for physics versus Bmad's
`67.370 s`, a measured `2.031x` same-host throughput advantage. It converged
1000/1000 lanes with zero full-AD fallbacks, maximum closure norm
`8.104e-11`, maximum output difference from Bmad `8.138494e-6 m`, and
correlation `0.999999966415499`. These runs were not back-to-back, and the
optimized SciBmad process again averaged only 49% CPU, so shared-host load
remains a qualification.

For the frozen-Jacobian implementation:

- final closure residual norm median: `4.107e-12`
- final closure residual norm maximum: `9.802e-11`
- automatic full-AD fallbacks in the formal 1,000-sample run: `0`
- maximum detector-orbit difference from the matching full-AD result:
  `5.566e-10 m`
- forced single-lane fallback test: `1/1` recovered successfully

A separate same-initial-guess test isolates the Jacobian strategy: local
full-AD Newton with nominal `z0` used `22.247 s`; the fallback-enabled frozen
repeat used `8.163 s`, an observed `2.725x` speedup.

The new `6 x 119` response predictor reduced the frozen solver's iteration
count from median/mean `3 / 2.994` to `2 / 1.995`. Physics time fell from
`8.163 s` to `6.855 s`, an observed `1.191x` recurring-batch speedup
(`16.0%` less time). All 1000 lanes passed the final closure check, with no
full-AD fallback; the maximum closure norm was `8.104e-11`. Relative to the
fixed-`z0` frozen output, detector-orbit RMSE was `9.411e-12 m` and the maximum
difference was `5.563e-10 m`.

## Numerical agreement

All SciBmad variants preserve essentially the same agreement with Bmad:

- global correlation: `0.999999966415`
- global maximum absolute difference: approximately `8.138e-6 m`
- median per-sample relative 2-norm difference: approximately `0.0338%`

The differences between the high-precision, normal-precision, and
frozen-Jacobian SciBmad detector tables are far smaller than the existing
SciBmad/Bmad model difference.

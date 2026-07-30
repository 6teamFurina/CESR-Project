# CESR 1,000-sample benchmark: unified results

Date: 2026-07-30  
Workload: 119 correctors to horizontal and vertical closed orbit at 99
detectors (`1000 x 198` output table), RF on.

## All requested results in one table

| Result | Machine | Tolerances `(rel, abs)` | Initial guess and Jacobian strategy | Converged | Warmup (s) | Physics (s) | Component total (s) | Physics samples/s | Observed rate / Bmad | RMSE vs Bmad (m) |
|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|
| Bmad/Tao/PyTao | `lnx201` Linux | defaults `(1e-8, 1e-10)`* | previous orbit; one-turn matrix reused and conditionally refreshed | 1000/1000 | 0.295 | 67.370 | 69.238 | 14.843 | 1.000 | 0 |
| SciBmad, VM high precision | `lnx201` Linux | `(1e-13, 1e-13)` | zero; full batched AD Jacobian each Newton iteration | 1000/1000 | 439.887 | 280.486 | 722.564 | 3.565 | 0.240 | `2.268158e-6` |
| SciBmad, local high precision | Windows, Ryzen 9 5900HX | `(1e-13, 1e-13)` | zero; full batched AD Jacobian each Newton iteration | 1000/1000 | 74.790 | 64.356 | 139.962 | 15.539 | 1.047 | `2.268158e-6` |
| SciBmad, local normal precision | Windows, Ryzen 9 5900HX | `(1e-8, 1e-10)` | zero; full batched AD Jacobian each Newton iteration | 1000/1000 | 75.016 | 26.457 | 102.270 | 37.798 | 2.546 | `2.268158e-6` |
| SciBmad, local frozen Jacobian + fallback | Windows, Ryzen 9 5900HX | `(1e-8, 1e-10)` | nominal `z0`; one nominal Jacobian reused; failed lanes use full AD | 1000/1000 | 59.235 | 8.163 | 68.123** | 122.506 | 8.253 | `2.268158e-6` |

\* Bmad and SciBmad use different mathematical stopping rules even when the
two tolerance values match. Bmad checks component-wise one-turn closure;
SciBmad checks residual or step norms.

\** The frozen-Jacobian component total includes nominal model setup and
nominal closed-orbit/Jacobian calculation (`0.216 s`), regular model setup,
warmup, physics, and writing. Other totals are the sums of the components
available in their metadata. They are not external wall-clock measurements.

## Interpretation

Only the three local SciBmad rows are same-machine comparisons. The
`Observed rate / Bmad` column is included to place all requested results in one
table, but ratios involving local Windows SciBmad and Linux Bmad are
cross-machine and are not controlled speedup claims.

The Bmad and VM high-precision SciBmad rows share `lnx201`, but the runs were
made at different times on a shared host. GNU `time` reported only 49% average
CPU utilization for the VM SciBmad run, so its `0.240` rate ratio should also
be repeated under controlled load.

For the frozen-Jacobian implementation:

- final closure residual norm median: `4.107e-12`
- final closure residual norm maximum: `9.802e-11`
- automatic full-AD fallbacks in the formal 1,000-sample run: `0`
- maximum detector-orbit difference from the matching full-AD result:
  `5.566e-10 m`
- forced single-lane fallback test: `1/1` recovered successfully

A separate same-initial-guess test isolates the Jacobian strategy: local
full-AD Newton with nominal `z0` used `22.247 s`; the fallback-enabled frozen
repeat used `8.163 s`, an observed `2.725x` speedup. The earlier frozen run
without fallback used `7.535 s`, indicating ordinary run-to-run timing
variation of several percent.

## Numerical agreement

All SciBmad variants preserve essentially the same agreement with Bmad:

- global correlation: `0.999999966415`
- global maximum absolute difference: approximately `8.138e-6 m`
- median per-sample relative 2-norm difference: approximately `0.0338%`

The differences between the high-precision, normal-precision, and
frozen-Jacobian SciBmad detector tables are far smaller than the existing
SciBmad/Bmad model difference.

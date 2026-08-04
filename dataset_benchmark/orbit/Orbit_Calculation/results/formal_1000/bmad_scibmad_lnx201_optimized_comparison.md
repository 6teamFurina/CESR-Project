# Optimized SciBmad versus Bmad on lnx201

Date: 2026-07-31  
Host: `lnx201.classe.cornell.edu`  
Workload: 1000 samples, 119 correctors, 99 detectors, 198 outputs per sample

| Result | Physics time | Throughput | Converged | Maximum output residual vs Bmad | Correlation vs Bmad | Maximum closure residual |
|---|---:|---:|---:|---:|---:|---:|
| Bmad/Tao/PyTao | 67.370 s | 14.843 samples/s | 1000/1000 | 0 | 1 | not recorded |
| SciBmad response + frozen Jacobian + fallback | **33.178 s** | **30.140 samples/s** | 1000/1000 | `8.138494e-6 m` | `0.999999966415499` | `8.104e-11` |

The optimized SciBmad recurring physics region was `2.031x` faster than the
existing Bmad result on the same host.

SciBmad details:

- closed-orbit solve: `30.022 s`
- detector tracking: `3.104 s`
- cached `6 x 119` response load: `0.000730 s`
- Newton iterations, min / median / mean / max: `0 / 2 / 1.995 / 2`
- full-AD fallbacks: `0`
- detector-output RMSE versus Bmad: `2.268158e-6 m`
- Julia warmup/compilation batch: `386.686 s`
- total process wall time: `8:56.57`
- average CPU reported by GNU time: `49%`
- maximum resident set: `1,317,652 KiB`
- exit status: `0`

The physics comparison is same-host but not back-to-back. `lnx201` is shared,
and the SciBmad process averaged only 49% CPU, so a paired repetition remains
recommended before presenting the ratio as a stable hardware-normalized
speedup. Startup/warmup is reported separately and is not included in the
recurring physics throughput claim.

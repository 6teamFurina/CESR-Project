# CESR same-host Bmad-SciBmad benchmark on lnx201

## Reproducibility

- Matched samples: `1000`
- Jointly converged samples: `1000`
- Controls per sample: `119`
- Observables per sample: `198`
- Bmad timed region: `variable update + lattice recalculation + observable read`
- SciBmad timed region: `closed-orbit solve + detector tracking`

## Numerical agreement

| Metric | Value |
|---|---:|
| Global RMSE (m) | `2.268157923e-06` |
| Global maximum absolute difference (m) | `8.138493963e-06` |
| Global correlation | `0.999999966415` |
| Median per-sample relative 2-norm difference | `3.377253707e-04` |
| Maximum per-sample relative 2-norm difference | `4.363242817e-04` |

## Timing

| Engine | Init/setup (s) | Warmup (s) | Physics (s) | Write (s) | Recorded total (s) | Physics samples/s |
|---|---:|---:|---:|---:|---:|---:|
| Bmad/Tao/PyTao | `1.030968` | `0.295234` | `67.370023` | `0.541796` | `69.238021` | `14.843397` |
| SciBmad | `0.607271` | `439.886528` | `280.485605` | `1.584338` | `722.563742` | `3.565245` |

Steady-state SciBmad/Bmad physics-throughput ratio: `0.240191`.
The recorded one-shot totals include initialization or model setup, warmup or
compilation, the timed physics region, and file writing.

The throughput comparison uses the same host, `lnx201.classe.cornell.edu`, one
host thread, the same input table, RF state, output schema, and convergence
tolerances. The runs occurred at different times on a shared host: GNU `time`
reports that the SciBmad process received only `49%` average CPU during
`13:21.11` elapsed wall time. This is therefore a same-host observation, but it
should be repeated under controlled host load before being presented as a
stable machine-level speed ratio.

The Linux and Windows SciBmad `1000 x 198` output tables are exactly equal
(`0` differing numerical entries), independently confirming platform
reproducibility for this dataset.
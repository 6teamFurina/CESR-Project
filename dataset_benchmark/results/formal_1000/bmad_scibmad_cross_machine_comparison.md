# CESR matched-dataset Bmad-SciBmad benchmark

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
| SciBmad | `0.169961` | `75.547633` | `65.154823` | `0.298656` | `141.171074` | `15.348058` |

Steady-state SciBmad/Bmad physics-throughput ratio: `1.033999`.
The recorded one-shot totals include initialization or model setup, warmup or
compilation, the timed physics region, and file writing.

The throughput ratio is meaningful only when both runs use the same hardware,
CPU thread count, convergence tolerances, input file, and output schema. The
present Bmad and SciBmad results were measured on different machines and must
not be presented as a controlled same-hardware speedup.

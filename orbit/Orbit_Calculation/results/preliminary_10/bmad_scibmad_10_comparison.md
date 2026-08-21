# CESR matched-dataset Bmad-SciBmad benchmark

## Reproducibility

- Matched samples: `10`
- Jointly converged samples: `10`
- Controls per sample: `119`
- Observables per sample: `198`
- Bmad timed region: `variable update + lattice recalculation + observable read`
- SciBmad timed region: `closed-orbit solve + detector tracking`

## Numerical agreement

| Metric | Value |
|---|---:|
| Global RMSE (m) | `2.203905650e-06` |
| Global maximum absolute difference (m) | `7.667968080e-06` |
| Global correlation | `0.999999970202` |
| Median per-sample relative 2-norm difference | `3.257344136e-04` |
| Maximum per-sample relative 2-norm difference | `3.497296292e-04` |

## Timing

| Engine | Init/setup (s) | Warmup (s) | Physics (s) | Write (s) | Recorded total (s) | Physics samples/s |
|---|---:|---:|---:|---:|---:|---:|
| Bmad/Tao/PyTao | `1.300531` | `0.304282` | `0.568615` | `0.003415` | `2.176843` | `17.586590` |
| SciBmad | `0.159022` | `78.965123` | `47.298132` | `0.228355` | `126.650632` | `0.211425` |

Steady-state SciBmad/Bmad physics-throughput ratio: `0.012022`.
The recorded one-shot totals include initialization or model setup, warmup or
compilation, the timed physics region, and file writing.

The throughput ratio is meaningful only when both runs use the same hardware,
CPU thread count, convergence tolerances, input file, and output schema. The
present Bmad and SciBmad results were measured on different machines and must
not be presented as a controlled same-hardware speedup.

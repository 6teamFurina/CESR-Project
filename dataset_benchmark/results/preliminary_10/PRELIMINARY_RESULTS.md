# Preliminary 10-sample matched benchmark

These runs validate the matched-dataset pipeline on a small RF-on CESR batch.
They are engineering measurements rather than the formal comparison. Every
run uses the same 10 samples, 119 corrector controls, and 198 detector outputs
(horizontal and vertical closed orbit at 99 `DET_*` markers).

## Measured results

| Engine/execution | Threads | Converged | Init/setup (s) | Warmup (s) | Timed physics (s) | Samples/s |
|---|---:|---:|---:|---:|---:|---:|
| Bmad/Tao/PyTao | n/a | 10/10 | 1.301 | 0.304 | 0.569 | 17.587 |
| SciBmad CPU | 1 | 10/10 | 0.159 | 78.965 | 47.298 | 0.211 |
| SciBmad CPU | 16 | 10/10 | 0.186 | 75.844 | 46.329 | 0.216 |
| SciBmad CUDA exploration | 1 host thread | 10/10 | 0.166 | 246.577 | 29.885 | 0.335 |

The matched Bmad/SciBmad CPU tables have global correlation
`0.999999970202`, global RMSE `2.204e-6 m`, and median per-sample relative
2-norm difference `0.0326%`.

The 1-thread and 16-thread CPU output tables are exactly equal. The maximum
absolute coordinate difference between the 1-thread CPU and exploratory CUDA
tables is `9.01e-14 m`.

## Interpretation

- Increasing Julia from 1 to 16 threads did not materially improve this
  current batched solver path.
- The CUDA experiment had higher timed throughput for this small run, but its
  compilation/warmup cost was much larger. Because the GPU path is not
  currently reliable enough to reproduce cleanly, it is archived here and is
  not part of the maintained benchmark.
- A 10-sample result is too small for a production-throughput claim because
  fixed compilation costs dominate SciBmad. Use `../formal_1000/` for the
  formal dataset result.
- The engines were measured on different machines, so even the formal
  throughput ratio is not a controlled same-hardware speedup measurement.

The JSON and TOML metadata files are the machine-readable timing sources.
Their absolute input/output paths retain the original run locations for
provenance.

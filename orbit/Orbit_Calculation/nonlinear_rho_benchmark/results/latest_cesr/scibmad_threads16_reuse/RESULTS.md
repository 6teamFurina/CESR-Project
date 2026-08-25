# Reused-Batch-Model SciBmad Benchmark

The latest-CESR nonlinear-rho benchmark was rerun with 16 Julia threads, one
BLAS thread, and one reusable 600-lane batch model. The lattice model was built
once for the first cell; the same model was retained while the 103 ordered
corrector `BatchParam` values were rebound for the remaining 14 cells.

## Timing

| Interval | Time [s] | Nonzero states/s |
|---|---:|---:|
| Physics | 17.989 | 500.318 |
| One batch-model construction | 4.658 | n/a |
| Fourteen control-batch updates | 0.195 | n/a |
| Physics including reusable-model setup | 22.842 | 394.018 |
| Shared cached-response initial-guess preparation | 5.804 | n/a |
| All runtime setup plus physics | 28.646 | 314.179 |

Compilation warmup (`67.276 s`) and CSV output are excluded. The primary
setup-inclusive interval follows the original benchmark convention: it includes
the batch-model preparation needed by the timed calculation but excludes the
separately measured response-initial preparation. The response matrix was
loaded from the existing GTPSA cache and was not recomputed.

Against the explicitly configured 16-OpenMP-thread persistent native-Tao
physics interval of `134.859 s`, the SciBmad ratios are `7.497x` for physics,
`5.904x` including reusable-model setup, and `4.708x` including the separately
measured cached-response initial preparation.

## Numerical validation

All `9,001/9,001` states converged. The maximum frozen-Jacobian iteration count
was seven, the maximum closure was `9.9921e-11`, and no full-AD fallback was
used. Comparison with the earlier 16-thread run covered all `2,592,288` saved
observables and found zero convergence or iteration mismatches, zero maximum
observable difference, zero relative-L2 observable difference, and zero maximum
closure difference.

The machine and runtime qualification remains unchanged: SciBmad ran in Windows
Julia, whereas the native Tao reference ran in the validated Ubuntu-Bmad WSL
environment.

# Nonlinear-rho benchmark: latest CESR ring

Status: `complete` for the matched 9,001-state latest-ring execution on
2026-08-23.

## Experiment definition

The shared input contains one zero-kick baseline and 9,000 nonzero states:
600 reusable random directions for each combination of `all`, `horizontal`,
and `vertical` controls with
`rho = 1.13, 3.2, 4.53, 6.4, 9.05`. The 103-control registry contains 58
horizontal and 45 vertical controls. The base kick is 5 microradians, the
direction seed is `20260803` (`20260804` and `20260805` for the horizontal and
vertical subsets), and the same directions are reused at every rho within a
scenario. The ordered input and manifest are under
[`shared_input/latest_cesr/`](shared_input/latest_cesr/).

Both calculations use the RF-on branch-0 latest CESR charged-particle model.
SciBmad loads
[`Latest_Lattice/latest_cesr_scibmad_repaired.jl`](../../../Latest_Lattice/latest_cesr_scibmad_repaired.jl),
and the explicitly labeled Bmad reference loads
[`Latest_Lattice/lat.bmad`](../../../Latest_Lattice/lat.bmad). Outputs contain
144 detector locations and matched x/y blocks, for 288 observables per state.

SciBmad uses the cached 6-by-103 GTPSA closed-orbit response to form the
nominal first-order initial estimate. It then evaluates the exact nonlinear
one-turn residual with a frozen nominal 6-by-6 phase-space Jacobian and one
shared LU factorization per 600-lane group. The run used `reltol = 1e-8`,
`abstol = 1e-10`, `maxiter = 100`, explicit closure checks, and the maintained
full-AD fallback for exceptional lanes. Bmad uses one persistent Tao/PyTao
process, its native closed-orbit calculation, one model recalculation per
sample, and `good_model` plus finite detector values as the available
convergence indicator.

## Completed results

- SciBmad converged 9,001/9,001 states. The maximum explicit one-turn closure
  norm was `9.9974057e-11`; no lane used the full-AD fallback. The maximum
  frozen-Jacobian iteration count was seven.
- Bmad/Tao returned 9,001/9,001 usable states. This PyTao path does not expose
  a SciBmad-equivalent one-turn closure norm.
- All 9,000 nonzero states converged in both engines and were included in the
  paired comparison.
- Summed nonzero physics time was 105.703 s for SciBmad and 166.608 s for
  Bmad/Tao, a 1.576x SciBmad physics-only advantage. Adding the 15 SciBmad
  batch-model constructions gives 273.372 s; adding the shared response-initial
  preparation gives 281.879 s. On this workflow, the setup-inclusive ratios
  are therefore 0.609x and 0.591x, respectively, relative to Bmad.
- For the `all` scenario, the largest baseline-subtracted relative RMSE over
  the five radii is 0.0290% in x and 0.0262% in y. The vertical-only
  principal-plane y relative RMSE reaches 0.0198%; its cross-plane x relative
  RMSE is approximately 0.553--0.558% with correlation at least 0.9999858.
  Horizontal-only y is essentially zero in both models, so its nominal 100%
  relative RMSE is not a meaningful agreement statistic.
- The zero-input baseline differences are at numerical scale: x RMSE
  `2.2308e-15 m` and y RMSE `3.1306e-16 m`.

The detailed result table is
[`results/latest_cesr/comparison/RESULTS.md`](results/latest_cesr/comparison/RESULTS.md).
Raw and diagnostic outputs are ring-scoped under:

```text
results/latest_cesr/scibmad/
results/latest_cesr/bmad_reference/
results/latest_cesr/comparison/
```

The timing comparison is application-level: SciBmad ran in Windows Julia and
Bmad/Tao ran inside Ubuntu-Bmad under WSL on the same physical machine.
Compilation warmup and file I/O are excluded from the summarized timing. The
latest lattice emits the documented straight-multipole/curved-reference
warning. This corrector-kick experiment does not vary girder pitch or photon
optics, so the curved-DQX pitch and photon-branch limitations are not exercised.

## Explicit 16-thread tracking follow-up

A separate follow-up enables BeamTracking CPU multithreading inside every
frozen-Jacobian residual evaluation and the final detector tracking. It uses
16 Julia threads and one BLAS thread; the original single-thread result above
is retained unchanged for provenance. On a warmed 600-state `all`,
`rho = 9.05` cell, the same process measured 4.954 s without threaded tracking
and 1.319 s with threaded tracking, a 3.754x speedup. The observables, final
six-dimensional orbits, closures, convergence flags, and iteration counts were
identical in this controlled cell test.

The complete threaded run converged 9,001/9,001 states with maximum closure
`9.9921e-11`, maximum seven iterations, and zero fallbacks. Summed nonzero
physics time was 47.746 s, or 188.497 samples/s. This is 2.214x faster than the
recorded single-thread physics time and 3.489x faster in wall clock than the
persistent native-Tao physics interval. Batch-model setup took another
120.676 s, giving 168.423 s for setup plus physics versus 166.608 s for Tao;
including shared response-initial preparation gives 174.171 s. The two
setup-inclusive ratios are 0.989x and 0.957x, respectively.

Across all 2,592,288 saved observables, the threaded and retained single-thread
runs differ by `1.56e-9` relative L2 with a maximum absolute difference of
`4.19e-10 m`. Two samples crossed the stopping test one iteration apart after
changing the BLAS thread setting, but both runs converged and no convergence
flag changed. The 3.489x number is therefore a 16-thread SciBmad versus native
Tao application-level wall-clock comparison, not a per-core speedup.

A matched reusable-model follow-up removed the 15 independent lattice
constructions from that orchestration path. One 600-lane batch model was built
once and reused by rebinding the 103 control batches between cells. The run
again converged 9,001/9,001 states with zero fallbacks and reproduced every
saved observable, iteration count, closure, and convergence flag exactly.
Physics took 17.989 s. One model construction took 4.658 s and the 14 control
updates took 0.195 s, giving 22.842 s for physics including reusable-model
setup. The cached-response initial-guess preparation remained separate at
5.804 s; including it gives 28.646 s.

The Tao reference was then rerun with `OMP_NUM_THREADS=16` and Tao
`global n_threads=16` recorded explicitly. One persistent process evaluated the
states sequentially while native OpenMP regions could use up to 16 threads. The
9,000 nonzero states took 136.797 s in the retained repeat, and all 9,002 output
CSV lines were text-identical to the first explicit-thread run. Against this explicit-thread
reference, the reusable-model SciBmad physics and setup-inclusive speedups are
7.605x and 5.989x.

Threaded outputs and comparisons are under:

```text
results/latest_cesr/scibmad_threads16/
results/latest_cesr/scibmad_threads16_reuse/
results/latest_cesr/bmad_reference_threads16/
results/latest_cesr/bmad_reference_threads16_run1/
results/latest_cesr/thread_scaling/
results/latest_cesr/comparison_threads16/
```

The historical root-level
[`results/comparison/RESULTS.md`](results/comparison/RESULTS.md) belongs to the
older CESR export and must remain labeled legacy. The former fixed-dimension
workflow is preserved in [`README_archived.md`](README_archived.md).

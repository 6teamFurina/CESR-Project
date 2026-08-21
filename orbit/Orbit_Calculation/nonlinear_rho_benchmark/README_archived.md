# Nonlinear-rho SciBmad/Bmad orbit benchmark

This benchmark compares the maintained SciBmad and Bmad/Tao closed-orbit
paths on identical corrector inputs in the moderate nonlinear regime.

The committed configuration uses the same direction-generation convention as
the response-rho sweep:

- scenarios: `all`, `horizontal`, and `vertical`;
- radii: `1.13`, `3.2`, `4.53`, `6.4`, and `9.05`;
- 600 fixed Gaussian unit-RMS directions per scenario and radius;
- base kick: `5e-6 rad`;
- seed: `20260803`, incremented by scenario as in the original sweep;
- one shared zero-control baseline.

This gives 9,001 states. SciBmad uses the first-order response initial guess,
one frozen nominal 6x6 Jacobian/LU factorization, exact nonlinear one-turn
residuals, closure checks, and full-AD fallback. Bmad uses one persistent Tao
process and its standard RF-on closed-orbit calculation, warm-starting from
the previous successful orbit as Tao normally does.

From `CESR Project` on Windows:

```powershell
julia --project=. orbit/Orbit_Calculation/nonlinear_rho_benchmark/generate_inputs.jl
julia --threads=auto --project=. orbit/Orbit_Calculation/nonlinear_rho_benchmark/run_scibmad_nonlinear_rho.jl
```

Run Bmad in the local `Ubuntu-Bmad` WSL distribution and then compare:

```powershell
wsl.exe -d Ubuntu-Bmad -- /home/joeyfurina/miniforge3/envs/bmad/bin/python <mounted-path>/run_bmad_nonlinear_rho.py
python orbit/Orbit_Calculation/nonlinear_rho_benchmark/compare_nonlinear_rho.py
```

Generated inputs and results are kept under `shared_input/` and `results/`.

## Recorded result (2026-08-06)

Both engines converged for all 9,001 states. Across the 9,000 nonzero states,
SciBmad required 28.204 s of steady-state physics time versus 102.725 s for
Bmad/Tao, a 3.642x speedup. Including all SciBmad runtime preparation (the
shared first-order initial guesses and all 15 batch-model constructions) gives
30.907 s and a 3.324x speedup. Compilation warmup and file I/O are excluded
from both steady-state comparisons.

The frozen-Jacobian solver used at most five iterations, every explicit
one-turn closure norm was below approximately `1e-10`, and no sample required
the full-AD fallback. After subtracting each engine's zero-input orbit, the
main-plane response RMSE is approximately 0.05% in x for all/horizontal inputs
and 0.24--0.28% in y for all/vertical inputs. The zero-input x orbits themselves
differ by 3.105 micrometers RMS, so absolute-orbit and baseline-subtracted
comparisons are both retained.

See [`results/comparison/RESULTS.md`](results/comparison/RESULTS.md) for the
per-cell table and timing interpretation, and
[`results/comparison/comparison_summary.csv`](results/comparison/comparison_summary.csv)
for the full machine-readable metrics.

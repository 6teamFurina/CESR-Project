# Latest CESR rho-sweep FD-validation smoke result

Status: `smoke/validation artifact`, generated 2026-08-20. This is an
integration check using the checked-in central-difference response cache, not
a GTPSA-default production run or a statistically meaningful response-validity
boundary. The maintained runner default is now `--response-method=gtpsa`;
recompute when method metadata does not match, or explicitly select central
difference.

## Configuration

- Lattice: [`latest_cesr_scibmad_repaired.jl`](../../../../../../Latest_Lattice/latest_cesr_scibmad_repaired.jl)
- Engine: SciBmad 0.4.1; Julia 1.12.6; RF on; branch 0
- Registry: 103 steering controls (58 H, 45 V), 144 detectors, 288 observables
- Directions: Gaussian, normalized to unit RMS over each active group; seed
  20260803 with scenario seeds 20260803/20260804/20260805
- Grid: `rho={0,0.1}`, base kick `5e-6 rad`, two trials per nonzero scenario
- Solver: frozen nominal Jacobian plus full-AD fallback,
  `reltol=1e-12`, `abstol=1e-13`, `maxiter=100`
- First-order response: SciBmad central difference, `h=1e-7 rad`, eight
  controls/chunk

All 7 nonlinear states converged, no fallback was used, and the maximum final
closure norm was `8.6186e-14`. Warmup took 53.240 s, the exact nonlinear
solve-and-track run took 50.338 s, and matrix response evaluation took 0.350 s.

At `rho=0.1`, the two-trial mean detector RMSEs were:

| Active controls | Mean x RMSE | Mean y RMSE |
|---|---:|---:|
| all | `9.25e-9 m` | `1.03e-8 m` |
| horizontal | `4.50e-9 m` | `1.49e-20 m` |
| vertical | `4.96e-9 m` | `3.41e-12 m` |

With only two directions, these values verify plumbing and symmetry behavior;
they do not estimate population statistics or define a validity threshold.

Artifacts:

- [`rho_sweep_summary.csv`](rho_sweep_summary.csv)
- [`rho_sweep_trial_errors.csv`](rho_sweep_trial_errors.csv)
- [`rho_sweep_metadata.toml`](rho_sweep_metadata.toml)

The run emitted the known straight-multipole/curved-reference warning for the
latest DQX representation. No girder-pitch perturbation was applied here, but
that limitation must accompany any result that exercises such errors.

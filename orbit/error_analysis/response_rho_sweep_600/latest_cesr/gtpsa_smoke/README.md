# Latest CESR rho-sweep GTPSA smoke result

Status: `smoke`, generated 2026-08-20. This is an end-to-end integration check
of the maintained GTPSA response path, not the 600-trial production study and
not a statistically meaningful response-validity boundary.

## Configuration

- Lattice: [`latest_cesr_scibmad_repaired.jl`](../../../../../Latest_Lattice/latest_cesr_scibmad_repaired.jl)
- Engine: SciBmad 0.4.1; Julia 1.12.6; RF on; branch 0
- Registry: 103 steering controls (58 H, 45 V), 144 detectors, 288 observables
- Directions: Gaussian, normalized to unit RMS over each active group; seed
  20260803 with scenario seeds 20260803/20260804/20260805
- Grid: `rho={0,0.1}`, base kick `5e-6 rad`, two trials per nonzero scenario
- Solver: frozen nominal Jacobian plus full-AD fallback,
  `reltol=1e-12`, `abstol=1e-13`, `maxiter=100`
- First-order response: the paired SciBmad GTPSA cache documented in
  [`reference/latest_cesr/gtpsa/`](../../../../reference/latest_cesr/gtpsa/README.md)

All 7 nonlinear states converged, no fallback was used, and the maximum final
closure norm was `8.8058e-14`. Warmup took 58.280 s, the exact nonlinear
solve-and-track run took 21.730 s, and matrix response evaluation took 0.276 s.

At `rho=0.1`, the two-trial mean detector RMSEs were:

| Active controls | Mean x RMSE | Mean y RMSE |
|---|---:|---:|
| all | `9.25e-9 m` | `1.03e-8 m` |
| horizontal | `4.50e-9 m` | `2.21e-20 m` |
| vertical | `4.96e-9 m` | `3.41e-12 m` |

With only two directions, these values verify plumbing and expected symmetry
behavior; they do not estimate population statistics or define a validity
threshold.

Artifacts:

- [`rho_sweep_summary.csv`](rho_sweep_summary.csv)
- [`rho_sweep_trial_errors.csv`](rho_sweep_trial_errors.csv)
- [`rho_sweep_metadata.toml`](rho_sweep_metadata.toml)

The run emitted the known straight-multipole/curved-reference warning for the
latest DQX representation. No girder-pitch perturbation was applied here, but
that qualification must accompany any result that exercises such errors.

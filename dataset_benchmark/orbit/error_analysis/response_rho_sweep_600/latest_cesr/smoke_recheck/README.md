# Latest CESR rho-sweep post-audit smoke

Status: `smoke`, generated 2026-08-20 after the response-radius runner and
chunk/renderer provenance audit. This is a 7-state integration check, not the
600-trial production study or a statistical validity boundary.

- Lattice: [`latest_cesr_scibmad_repaired.jl`](../../../../../../Latest_Lattice/latest_cesr_scibmad_repaired.jl)
- Engine: SciBmad 0.4.1; Julia 1.12.6; RF on; branch 0
- Runtime registry: 103 steering controls (58 H, 45 V), 144 detectors, 288
  ordered `x/y` observables
- Grid: `rho={0,0.1}`, two Gaussian directions per nonzero scenario, seed
  `20260803`, base kick `5e-6 rad`
- Solver: frozen nominal Jacobian with full-AD fallback,
  `reltol=1e-12`, `abstol=1e-13`, `maxiter=100`
- First-order response: GTPSA closed-orbit/detector pair; the pair id and
  lattice/control provenance are recorded in `rho_sweep_metadata.toml`

All 7 states converged with no fallback. The maximum final closure norm was
`8.8058e-14`; warmup took `67.833 s`, exact nonlinear solve-and-track took
`18.760 s`, and response evaluation took `0.2967 s`. These times include
compilation effects and are not throughput claims.

Artifacts:

- [`rho_sweep_summary.csv`](rho_sweep_summary.csv)
- [`rho_sweep_trial_errors.csv`](rho_sweep_trial_errors.csv)
- [`rho_sweep_metadata.toml`](rho_sweep_metadata.toml)
- [`figures/scibmad_orbit_response_error.svg`](figures/scibmad_orbit_response_error.svg)

The run emitted the known straight-multipole/curved-reference warning for the
latest DQX representation. No girder-pitch perturbation was applied; the
curved-DQX girder-pitch limitation remains a qualification for studies that
exercise that feature.

# Unified-runner rho smoke

Status: `smoke`, generated 2026-08-20 through
[`run_latest_cesr_experiments.ps1`](../../../../high_throughput_nonlinear_closed_orbit_calculation_and_response_error_analysis/run_latest_cesr_experiments.ps1).
This is an orchestration/integration check, not the 600-direction paper result.

- SciBmad latest repaired CESR lattice, RF on, branch 0
- 103 dynamically discovered steering controls (58 horizontal, 45 vertical)
- 144 detectors and 288 ordered x/y observables
- `rho={0,0.1}`, two directions for each nonzero all/H/V scenario
- GTPSA response pair `cdaa3f56-a96d-4198-b687-37f9120ffdd1`
- 7/7 exact nonlinear states converged, zero fallback
- maximum closed-orbit closure norm `8.8058e-14`

Warmup (`62.493 s`) and exact-reference (`18.804 s`) timings include a fresh
Julia process and are not throughput claims. The result-local figure is
[`figures/scibmad_orbit_response_error.svg`](figures/scibmad_orbit_response_error.svg).

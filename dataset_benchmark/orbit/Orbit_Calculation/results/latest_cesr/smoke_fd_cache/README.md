# Cache-loaded latest CESR FD-validation smoke

Status: `smoke/validation artifact`, generated 2026-08-20. This records an
end-to-end check of the production-style closed-orbit path using the historical
central-difference cache; the maintained runner default is now
`--response-method=gtpsa`.

- Lattice: [`latest_cesr_scibmad_repaired.jl`](../../../../../../Latest_Lattice/latest_cesr_scibmad_repaired.jl)
- Engine: SciBmad 0.4.1; Julia 1.12.6; RF on; branch 0
- Input/output: 2 samples, 103 controls, 144 detectors, 288 observables
- Initial guess: cached `6 x 103` SciBmad central-difference response
- Solver: frozen nominal `6 x 6` Jacobian with full-AD fallback enabled,
  `reltol=1e-8`, `abstol=1e-10`
- Convergence: 2/2; fallback 0; iterations min/median/mean/max = 0/1/1/2
- Timings: warmup 99.490 s; nominal setup/solve 5.409/3.095 s; model setup
  5.840 s; solve 3.113 s; tracking 0.024 s; timed physics 3.137 s
- Maximum final closure norm: `8.1375e-11`
- Cache provenance loaded from the sidecar: central difference,
  `h=1e-7 rad`, 8 controls/chunk, 13 chunks

Artifacts:

- [`scibmad_rf_on_samples.csv`](scibmad_rf_on_samples.csv)
- [`scibmad_rf_on_metadata.toml`](scibmad_rf_on_metadata.toml)
- [`closed_orbit_response_6x103.csv`](closed_orbit_response_6x103.csv)

This is a two-sample integration/FD-validation test, not a sustained-throughput
result or a GTPSA-default production result. Recompute if response metadata
does not match the requested method, or explicitly request central difference.
The curved-reference straight-multipole warning is the documented DQX
limitation.

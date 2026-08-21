# Orbit studies: latest CESR ring

Status: `implementation_smoke_complete` for the shared closed-orbit runner and
the current calculation paper's three retained latest-ring experiments as of
2026-08-20. Production ensembles and the manuscript rebuild remain `not_run`.
Parity/cubic, response-block, and correction-paper studies are outside this
completion gate. No numerical result under an unscoped legacy result directory
should be read as a latest-ring result.

The default lattice for new work is
[`Latest_Lattice/latest_cesr_scibmad_repaired.jl`](../../Latest_Lattice/latest_cesr_scibmad_repaired.jl).
The repository ring id for generated artifacts is `latest_cesr`.

## Result and provenance convention

New artifacts should use ring-scoped paths so that a rerun cannot silently
mix lattice versions:

```text
Orbit_Calculation/inputs/latest_cesr/
reference/latest_cesr/
Orbit_Calculation/results/latest_cesr/
<error-study>/results/latest_cesr/
<manuscript>/figures/latest_cesr/
```

Every latest-ring result README should record the lattice path and ring id,
branch, dynamically discovered control and observable registries, RF state,
seed, solver tolerances, SciBmad version, output paths, convergence/closure
statistics, and a status of `not_run`, `smoke`, or `production`.

The charged-particle branch of the latest lattice currently has 1,177
tracking elements and 144 `DET_*` marker definitions. These are inventory
facts for the current export, not interface constants: runners must discover
the actual controls, detector subset, ordering, and output layout at runtime.

The latest-lattice validation reports a known limitation for girder pitch
through curved DQX multipoles. It must be mentioned whenever a study exercises
that feature. Photon-branch limitations are outside the closed-orbit studies
unless a study explicitly uses those branches.

## Verified smoke scope

Runtime discovery found 124 writable controls: 119 Overlay and 5 Group
controls. The default orbit-steering subset contains 103 normal correctors
(58 horizontal and 45 vertical), and the 144 detectors produce 288 ordered
`x/y` observables.

- The cache-loaded two-sample benchmark converged 2/2 states with 3.137 s in
  the timed solve-and-track region and a maximum closure norm of `8.14e-11`
  at `reltol=1e-8`, `abstol=1e-10`. See
  [`Orbit_Calculation/results/latest_cesr/README.md`](Orbit_Calculation/results/latest_cesr/README.md).
- The checked-in first-order SciBmad validation cache has shapes `6 x 103` and
  `288 x 103`. It was generated in 13 bounded BatchParam central-difference
  chunks with `h=1e-7 rad`; the maximum response-lane closure norm was
  `9.95e-14`. The maintained runner default is now
  `--response-method=gtpsa`; use this FD cache only as a labeled smoke/
  validation artifact, or explicitly select central difference. See
  [`reference/latest_cesr/README.md`](reference/latest_cesr/README.md).
- The minimal response-radius run used `rho={0,0.1}` and two directions per
  nonzero scenario. All 7 states converged with no fallback and maximum
  closure norm `8.81e-14`. See
  [`error_analysis/response_rho_sweep_600/README.md`](error_analysis/response_rho_sweep_600/README.md).
- Complete-element attribution covered all 1,177 elements. The one-direction
  x/y integration smokes close their summed nonlinear targets to `4.23e-15`
  and `1.28e-14`; family projections partition the all-element result. See
  [`error_analysis/thick_element_sextupole_sourcing/README.md`](error_analysis/thick_element_sextupole_sourcing/README.md).
- The matched two-direction sextupole pipeline retained 76 dynamically
  discovered active normal sextupoles. Its all-element x/y reconstruction
  closes to `1.38e-14`, and the associated nominal/direction optics and
  beta/phase predictor tables are finite and cardinality-consistent. See
  [`error_analysis/sextupole_detector_contributions/README.md`](error_analysis/sextupole_detector_contributions/README.md)
  and
  [`error_analysis/sextupole_beta_phase_correlation/README.md`](error_analysis/sextupole_beta_phase_correlation/README.md).

These are integration tests, not production timing or statistical claims. The
checked-in first-order cache is a bounded SciBmad BatchParam central-difference
smoke artifact, but this is not a conclusion that GTPSA must be replaced. An
earlier many-parameter GTPSA diagnostic stopped at `SEX_14W` because the ring
adapter represented unselected controls as TPS quantities, triggering the
lower-level `sqrt(0)` domain error. With `zero_value=0.0` and only the selected
GTPSA controls represented as parameters, the full 1,177-element map succeeds.
Against the central-difference cache, the reported relative-L2 differences
were `2.78e-8` for the `6 x 103` response and `1.68e-8` for the `288 x 103`
detector response. The checked-in central-difference cache is retained as a
labeled smoke/validation artifact; recompute it when its method metadata does
not match the requested method, or explicitly select central difference for a
reproducibility comparison.

## Study map

- [`Orbit_Calculation/README.md`](Orbit_Calculation/README.md): generic
  closed-orbit calculation and throughput workflow.
- [`Orbit_Calculation/nonlinear_rho_benchmark/README.md`](Orbit_Calculation/nonlinear_rho_benchmark/README.md):
  matched nonlinear-rho benchmark.
- [`error_analysis/README.md`](error_analysis/README.md): response-error and
  nonlinear-order studies.
- [`error_analysis/response_rho_sweep_600/README.md`](error_analysis/response_rho_sweep_600/README.md):
  radius-sweep result location.
- [`high_throughput_nonlinear_closed_orbit_calculation_and_response_error_analysis/README.md`](high_throughput_nonlinear_closed_orbit_calculation_and_response_error_analysis/README.md):
  calculation manuscript migration status.
- [`hessian_svd_nonlinear_closed_orbit_correction/README.md`](hessian_svd_nonlinear_closed_orbit_correction/README.md):
  correction-paper migration status.
- [`reference/README.md`](reference/README.md): ring-scoped response-cache
  convention.

## Historical artifacts

The previous study descriptions remain at
[`README_archived.md`](README_archived.md) in each study directory. Existing
CSV, TOML, JSON, SVG, PDF, archive packages, and custom `RESULTS.md` reports
have deliberately not been renamed in this documentation pass. They were
generated from the older `cesr.jl`/Bmad-compatible export and are historical
until a result README explicitly identifies them as a latest-ring artifact.

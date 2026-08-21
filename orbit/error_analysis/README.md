# Orbit-response error analysis: latest CESR ring

Status: `implementation_smoke_complete` for the three error-analysis threads
retained by the current paper. Production ensembles have not been run.

The frozen paper scope is: nonlinear response error versus `rho`, attribution
of the summed leading nonlinear error to complete element types/families, and
normal-sextupole element plus source--beta--phase attribution. Detailed
`hh`/`hv`/`vv` block shares, their attribution, and cubic/third-order parity
analysis are follow-up studies and are not completion gates.

All new error-analysis calculations use
[`Latest_Lattice/latest_cesr_scibmad_repaired.jl`](../../Latest_Lattice/latest_cesr_scibmad_repaired.jl)
and the repository ring id `latest_cesr`. New outputs should live below a
ring-scoped result directory, for example:

```text
response_rho_sweep_600/latest_cesr/
thick_element_sextupole_sourcing/<plane>_results/latest_cesr/
sextupole_detector_contributions/results/latest_cesr/
sextupole_beta_phase_correlation/results/latest_cesr/
```

The response approximation compares a dynamically discovered first-order
observable response with converged nonlinear SciBmad closed orbits. Control
groups, detector planes, labels, and dimensions must come from the ring
configuration; the old prefix/count fallback is only for archived data.

The bounded SciBmad BatchParam central-difference cache (`h=1e-7 rad`, eight
controls per chunk) is a labeled response validation artifact, not the current
method default. The maintained rho runner defaults to
`--response-method=gtpsa` and uses the
method-scoped pair under
`reference/latest_cesr/gtpsa/{closed_orbit_response,detector_response}.csv`.
Use `--response-method=central-difference` explicitly to reuse the archived
central-difference pair, or `--recompute-response=true` to atomically rebuild
the selected pair without deleting the other method's cache. This is not a
claim that GTPSA must be replaced. An earlier many-parameter GTPSA diagnostic
stopped at `SEX_14W` because the adapter TPS-ified unselected controls, causing
the lower-level `sqrt(0)` domain error. With `zero_value=0.0` and only selected
GTPSA controls represented as parameters, the full 1,177-element map succeeds.
The reported relative-L2 differences against the cache are `2.78e-8` for
`6 x 103` and `1.68e-8` for `288 x 103`. The latest GTPSA rho smoke converged
7/7 exact states with no fallback. Complete-element x/y reconstruction and the
matched two-direction sextupole/beta-phase pipeline also pass their closure and
finite-value audits. These small samples validate implementation, not paper
statistics.

Each completed study README must report the lattice and branch, control and
observable registries, RF state, direction distribution, seed, amplitudes and
units, solver tolerances, convergence/closure/fallback statistics, and whether
the result is `smoke` or `production`. A result is not a CESR machine-error
budget unless it is separately calibrated and validated against measurements.

The existing response sweeps, parity tables, mixed-term tables, sextupole
attributions, figures, and custom result reports were generated before the
latest-lattice migration. They remain available through their original paths
and are described as legacy in the archived study descriptions.

Required paper-study latest-ring pages:

- [`response_rho_sweep_600/README.md`](response_rho_sweep_600/README.md)
- [`thick_element_sextupole_sourcing/README.md`](thick_element_sextupole_sourcing/README.md)
- [`sextupole_detector_contributions/README.md`](sextupole_detector_contributions/README.md)
- [`sextupole_beta_phase_correlation/README.md`](sextupole_beta_phase_correlation/README.md)

Retained follow-up pages (not current-paper gates):

- [`mixed_terms/README.md`](mixed_terms/README.md)
- [`vertical_parity/README.md`](vertical_parity/README.md)
- [`quadratic_x_attribution/README.md`](quadratic_x_attribution/README.md)
- [`sextupole_cascade/README.md`](sextupole_cascade/README.md)

The previous fixed-dimension study description is preserved in
[`README_archived.md`](README_archived.md).

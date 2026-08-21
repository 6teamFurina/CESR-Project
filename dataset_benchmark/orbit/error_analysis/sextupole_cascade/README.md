# Nonlinear multipole-cascade study: latest CESR ring

Status: `not_run`.

The latest-ring cascade and optional wiggler-corner experiments will use
[`Latest_Lattice/latest_cesr_scibmad_repaired.jl`](../../../../Latest_Lattice/latest_cesr_scibmad_repaired.jl)
and recompute the nominal orbit, response, and detector observables for every
modified lattice. Active order-2 fields, detector labels, planes, and control
groups must be discovered rather than assumed.

New outputs should use `results/latest_cesr/` and
`wiggler_corner_results/latest_cesr/`. Each result page must state the strength
grid, direction/radius grid, seed, RF state, units, closure and convergence,
and the vector-level residual for every decomposition. Any SciBmad capability
limitation exercised by a wiggler or curved element must be visible in the
result.

The existing
[`results/SEXTUPOLE_CASCADE_RESULTS.md`](results/SEXTUPOLE_CASCADE_RESULTS.md)
and
[`wiggler_corner_results/WIGGLER_CORNER_RESULTS.md`](wiggler_corner_results/WIGGLER_CORNER_RESULTS.md)
are legacy reports. The former study description is preserved in
[`README_archived.md`](README_archived.md).

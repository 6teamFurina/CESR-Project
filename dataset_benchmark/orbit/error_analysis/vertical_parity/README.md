# Signed-parity response study: latest CESR ring

Status: `not_run`.

This study will repeat the positive/negative excitation experiment with the
latest SciBmad ring
[`Latest_Lattice/latest_cesr_scibmad_repaired.jl`](../../../../Latest_Lattice/latest_cesr_scibmad_repaired.jl).

The new run must obtain the selected control group, detector planes, labels,
and response matrix from runtime configuration. It must record the direction
count and seed, radius grid, base kick and units, RF state, descriptor and
solver tolerances, convergence/closure/fallback counts, and the definitions of
the even and odd residuals. New artifacts belong under `results/latest_cesr/`.

The first-order cache defaults to GTPSA with a method-scoped pair under
`reference/latest_cesr/gtpsa/`. Use
`--response-method=central-difference` for explicit finite-difference
validation, or `--recompute-response=true` to atomically regenerate the
selected pair without touching the other method's cache.

The existing `results/` CSV, figures, and `VERTICAL_PARITY_RESULTS.md` are
legacy old-ring outputs. They are retained without renaming. The old narrative
is preserved in [`README_archived.md`](README_archived.md).

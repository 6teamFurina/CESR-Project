# Nonlinear-rho benchmark: latest CESR ring

Status: `not_run` for the 9,001-state production execution; the latest-ring
input manifest and response-cache path are prepared and syntax-checked.

This study will compare exact nonlinear closed-orbit calculations with the
local first-order response model over reusable random directions and selected
input radii. The latest SciBmad lattice is
[`Latest_Lattice/latest_cesr_scibmad_repaired.jl`](../../../../Latest_Lattice/latest_cesr_scibmad_repaired.jl).

For `--ring=latest`, the SciBmad runner reads the 103-control latest input at
[`shared_input/latest_cesr/nonlinear_rho_correctors.csv`](shared_input/latest_cesr/nonlinear_rho_correctors.csv),
its paired manifest, and the GTPSA closed-orbit response under
`reference/latest_cesr/gtpsa/`. The response method is passed explicitly as
GTPSA; the root-level central-difference cache is not used by the latest run.

Use ring-scoped paths for the new run:

```text
shared_input/latest_cesr/
results/latest_cesr/scibmad/
results/latest_cesr/bmad/       # optional labeled reference
results/latest_cesr/comparison/
```

The direction generator and solver must obtain control names, horizontal or
vertical groups, detector labels, and output planes from the runtime ring
configuration. The README for a completed run must state the number of
directions, radii, samples, RF mode, base kick and units, seed, convergence
counts, closure threshold, fallback count, and SciBmad/Bmad provenance.

No latest-ring numerical result is present yet. The existing
[`results/comparison/RESULTS.md`](results/comparison/RESULTS.md) and raw files
are historical outputs from the older CESR export; they are retained without
renaming and must remain labeled legacy in downstream comparisons.

The former fixed-dimension workflow is preserved in
[`README_archived.md`](README_archived.md).

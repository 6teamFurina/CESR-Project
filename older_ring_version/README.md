# Historical CESR ring version

This directory preserves the earlier maintained CESR SciBmad/Bmad model for
historical reproduction. It is not the default lattice for new calculations;
use `../Latest_Lattice/latest_cesr_scibmad_repaired.jl` for current work.

- `cesr_model.jl` is the historical model entry point.
- `cesr.jl`, `cesr_controls.jl`, and `cesr.bmad` are the associated lattice and
  control sources.
- `test_codes/`, `bmad_comparison/`, and `wigglers/` contain version-specific
  tests, comparisons, and experiments.
- `packages/` contains frozen historical bundles.
- The shared wiggler implementation remains at `../wigglers/wiggler.jl` because
  it is also required by the current latest lattice.

Run historical Julia workflows from the project root with `--project=.` and an
`older_ring_version/...` script path.

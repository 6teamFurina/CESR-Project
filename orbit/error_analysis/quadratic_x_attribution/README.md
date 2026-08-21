# Quadratic-response attribution: latest CESR ring

Status: `not_run`.

The latest-ring version will attribute the selected quadratic detector response
using runtime-discovered normal multipoles, complete elements, controls, and
observable planes. Its default lattice is
[`Latest_Lattice/latest_cesr_scibmad_repaired.jl`](../../../Latest_Lattice/latest_cesr_scibmad_repaired.jl).

Use `results/latest_cesr/` and `element_results/latest_cesr/` for new data and
figures. Preserve the source-to-detector vector convention, signed versus
unsigned normalization, direction seed, RF state, descriptor, element
inventory, and vector closure residual in the result README. A source ranking
must not be presented as a complete attribution unless the all-element vector
closure is reported.

The existing
[`results/RESULTS.md`](results/RESULTS.md) and
[`element_results/ELEMENT_EXPOSURE_RESULTS.md`](element_results/ELEMENT_EXPOSURE_RESULTS.md)
are legacy reports from the older ring. The old study description is preserved
in [`README_archived.md`](README_archived.md).

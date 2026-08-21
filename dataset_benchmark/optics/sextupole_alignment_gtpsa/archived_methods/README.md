# Archived sextupole-alignment methods

These directories are frozen historical methods based on the older
`cesr_model.jl` lattice:

- `response_map/` generated the 76-sextupole mixed `Kn2`/offset GTPSA response
  map, finite-difference validation, and local SVD summaries.
- `targeted_bump_k2_inversion/` used that map for the P0--P3 and P1/P2
  bump-by-`K2` inversion benchmarks.

The move into this folder changes organization, not the provenance or meaning
of the saved numerical results. Absolute paths embedded in frozen result
metadata still record their original generation locations. New work must use
`Latest_Lattice/latest_cesr_scibmad_repaired.jl` and a maintained study folder.

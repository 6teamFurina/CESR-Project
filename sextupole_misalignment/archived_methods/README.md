# Archived sextupole-alignment methods

These directories retain historical methods and completed protocol comparisons.
Archival status does not imply that every study used the same lattice or engine.

- [`response_map/`](response_map/README.md) generated the 76-sextupole mixed
  `Kn2`/offset GTPSA response map and local SVD summaries with the older
  `older_ring_version/cesr_model.jl` lattice.
- [`targeted_bump_k2_inversion/`](targeted_bump_k2_inversion/README.md) used
  that older-lattice map for the P0--P3 and P1/P2 inversion benchmarks.
- [`bmad_quadrupole_affinity/`](bmad_quadrupole_affinity/README.md) contains
  the old Bmad/Tao affinity screen and explicit saved-response comparison with
  the archived SciBmad/GTPSA map.
- [`direct_observable_k1_pilots/`](direct_observable_k1_pilots/README.md)
  contains the early observable-selection and seven-condition K1 pilots,
  generated with the repaired latest SciBmad lattice.
- [`interleaved_measurement_protocol/`](interleaved_measurement_protocol/README.md)
  preserves the blocked/interleaved and repeated-read protocol comparison,
  based on latest-lattice SciBmad physical scans.

The move into this folder changes organization, not the provenance or meaning
of the saved numerical results. Absolute paths embedded in frozen result
metadata still record their original generation locations. New work must use
`Latest_Lattice/latest_cesr_scibmad_repaired.jl` and a maintained study folder.

## Relocation completed on 2026-09-05

The last three studies above were archived with user authorization. The exact
11 source/destination paths are recorded in
[`archive_moves_2026-09-05.tsv`](archive_moves_2026-09-05.tsv). All 139 moved
files retained their identity, size, and modification time during the move;
only source paths, imports, and explanatory documentation were then updated.
Saved numerical outputs and generation metadata were preserved.

The all-target oracle tensor, shared `analyze_protocol_subsampling.py` fitting
functions, latest SciBmad affinity calculation, and bump knobs remain in their
maintained locations. The archived pilots and interleaved analysis still depend
on those shared studies; their relative paths have been updated accordingly.
The isolated Taylor-map benchmark, K1 triplet/exact-11 study, and older
sequential/ORM comparison implementations were not relocated in this step.

Migration verification is recorded in
[`ARCHIVE_VALIDATION_2026-09-05.json`](ARCHIVE_VALIDATION_2026-09-05.json).
The archived interleaved validator passes; the saved Bmad-versus-SciBmad
comparison reproduces its original tables; seven Python CLI/import checks and
both Julia generator include checks pass. The Julia checks loaded the shared
latest-lattice model and verified bump-knob paths without generating new scans.
All pre-existing numerical result files across the study were preserved.

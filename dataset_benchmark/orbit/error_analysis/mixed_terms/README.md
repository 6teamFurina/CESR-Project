# Mixed quadratic response: latest CESR ring

Status: `not_run`.

This study will calculate horizontal, vertical, and mixed second-order
closed-orbit response blocks with SciBmad/GTPSA on the latest ring:
[`Latest_Lattice/latest_cesr_scibmad_repaired.jl`](../../../../Latest_Lattice/latest_cesr_scibmad_repaired.jl).

New outputs should use `gtpsa_results/latest_cesr/`,
`results/latest_cesr/`, and, if retained, `smoke/latest_cesr/`. The control
directions, plane groups, detector labels, and block dimensions must be
discovered from the ring configuration. Report parameter descriptor, RF mode,
direction normalization, seed, amplitudes and units, closure residuals, and
finite-difference validation status.

The first-order response used by the four-sign experiment defaults to
`--response-method=gtpsa`; its paired caches are method-scoped below
`reference/latest_cesr/gtpsa/`. Select
`--response-method=central-difference` explicitly for the independent bounded
BatchParam validation, or use `--recompute-response=true` to atomically rebuild
the selected cache pair while retaining the other method's artifacts. GTPSA
parameterizes only the normal H/V steering subset; inactive skew/group
controls remain primitive values.

The previous reports remain at
[`gtpsa_results/GTPSA_RESULTS.md`](gtpsa_results/GTPSA_RESULTS.md),
[`results/MIXED_TERM_RESULTS.md`](results/MIXED_TERM_RESULTS.md), and
[`smoke/MIXED_TERM_RESULTS.md`](smoke/MIXED_TERM_RESULTS.md). They are legacy
outputs and are not latest-ring evidence. The prior narrative is preserved in
[`README_archived.md`](README_archived.md).

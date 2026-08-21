# Orbit reference data: latest CESR ring

Status: `smoke` for the latest-ring SciBmad response caches.

New orbit studies should load
[`Latest_Lattice/latest_cesr_scibmad_repaired.jl`](../../Latest_Lattice/latest_cesr_scibmad_repaired.jl)
and write ring-scoped reference data under:

```text
reference/latest_cesr/
```

The response cache must be generated from the runtime control registry and
the selected detector/observable registry. Its metadata should preserve the
ordered names, units, RF mode, branch, lattice path, and response shape; do
not encode a fixed filename such as `closed_orbit_response_6x119.csv` in the
latest workflow.

The maintained runner default is now `--response-method=gtpsa`. The checked-in
central-difference files are labeled smoke/validation artifacts: explicitly
select central difference to reproduce them, or recompute when response
metadata does not match the requested method. The current caches are documented
in [`latest_cesr/README.md`](latest_cesr/README.md). The maintained GTPSA pair
is under [`latest_cesr/gtpsa/`](latest_cesr/gtpsa/README.md); the root-level
`latest_cesr` pair is the SciBmad BatchParam central-difference validation at
`h=1e-7 rad`. Both contain a `6 x 103` closed-orbit response and a `288 x 103`
detector response. Sidecar TOML files record the ordered labels, ring/lattice
identity, RF state, method, shapes, chunking, closure diagnostic, pair identity,
and SciBmad version.

The files currently in this directory are retained as legacy references:

- `cesr_bmad_compatible.bmad` and its digested file;
- `closed_orbit_response_6x119.csv`.

They belong to the older CESR export and must not be used as the default
response cache for `latest_cesr`. Bmad/Tao files may be retained for an
explicit cross-code comparison, while SciBmad remains the primary source for
latest-ring quantitative results.

The directory-level description of those retained legacy files is preserved in
[`README_archived.md`](README_archived.md).

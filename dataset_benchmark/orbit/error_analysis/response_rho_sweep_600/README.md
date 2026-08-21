# Response-radius sweep: latest CESR ring

Status: `smoke`. A 7-state integration run completed on the latest ring; the
600-trial production grid has not been run.

The latest-ring sweep will use
[`Latest_Lattice/latest_cesr_scibmad_repaired.jl`](../../../../Latest_Lattice/latest_cesr_scibmad_repaired.jl)
and write new artifacts under:

```text
latest_cesr/
```

The result README at
[`latest_cesr/README.md`](latest_cesr/README.md) identifies the control
registry, active control groups, observable/plane layout, direction count and
seed, radius grid, base kick and units, RF state, response-cache provenance,
solver tolerances, convergence and closure statistics, and incomplete-cell
handling. Figures must be generated from the same `latest_cesr` data directory.

The default first-order backend is GTPSA. Its closed-orbit and detector
responses are generated as one atomically published pair under
`reference/latest_cesr/gtpsa/`. Select
`--response-method=central-difference` for the bounded BatchParam validation,
or pass `--recompute-response=true` to rebuild the selected pair while
preserving the other method's cache.

The parallel production launcher passes the response method explicitly and
serializes the first chunk so the shared pair is complete before the remaining
chunks start:

```powershell
.\run_response_rho_sweep_parallel.ps1 -ResponseMethod gtpsa -Trials 600
```

Merge only disjoint chunks from the same latest-ring run. The merger checks
the lattice, dynamic dimensions, solver contract, response method, and shared
`response_pair_id`; it rejects overlapping scenario/rho/trial cells:

```powershell
python merge_response_rho_sweep_chunks.py `
  --root response_rho_sweep_600/latest_cesr/chunks `
  --output-dir response_rho_sweep_600/latest_cesr/merged
```

The renderer accepts the latest-ring TOML sidecar or merged JSON sidecar and
will not discover the unscoped historical result by default:

```powershell
python render_response_rho_sweep_svg.py `
  --summary response_rho_sweep_600/latest_cesr/gtpsa_smoke/rho_sweep_summary.csv
```

The existing `chunks/`, `combined/`, and `figures/` files are retained as
legacy artifacts from the old CESR export. They are not latest-ring results
and should not be overwritten by a latest run. Their directory-level archived
narrative is preserved in [`README_archived.md`](README_archived.md), with the
broader error-analysis history in [`../README_archived.md`](../README_archived.md).

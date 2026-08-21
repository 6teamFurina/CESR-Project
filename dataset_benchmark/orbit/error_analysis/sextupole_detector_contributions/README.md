# Signed normal-sextupole contributions: latest CESR ring

Status: `smoke validated 2026-08-20`; production is intentionally not run in
this checkout. The validated latest-lattice smoke uses two matched directions.

The maintained calculation uses the validated SciBmad lattice
[`Latest_Lattice/latest_cesr_scibmad_repaired.jl`](../../../../Latest_Lattice/latest_cesr_scibmad_repaired.jl)
and discovers the control, detector, element, and active normal-`Kn2` registries
at runtime. The target is the total second-order detector vector
`Q = (Q_x, Q_y)` from the implicit two-parameter RF-on closed-orbit derivative.
The source boundary is every complete element exit,

`g_j = S_exit,j - A_j S_entrance,j`,

so the all-element reconstruction is a numerical closure reference. The
reported normal-sextupole rows retain only active normal sextupole elements and
their signed projected contribution
`eta_j = <C_j,Q>/<Q,Q>`. A signed projection is not a positive error share;
element vectors can interfere.

The latest runner does not include the old internal-exposure runner. Its
linear maps initialize inactive latest-lattice controls as primitive `0.0` and
parameterize only the selected steering controls, avoiding the known
combined-multipole `sqrt(0)` domain at `SEX_14W`. The direction derivative uses
`Descriptor(6,2,2,2)` and keeps the full six-dimensional source. The output
contains no `hh`/`hv`/`vv` block shares or third-order terms.

Outputs are ring-scoped and kept separate from the historical old-ring files:

- `results/latest_cesr/` is the production destination;
- `smoke/latest_cesr/` contains minimal endpoint checks;
- `reconstruction_summary.csv` records all-element and normal-sextupole-only
  vector closure, signed projection, and source-response closure;
- `sextupole_direction_contributions.csv` contains per-direction signed
  x/y/total contribution vectors summarized by norm and projection, plus the
  local source-kick fields consumed by the beta/phase predictor;
- `metadata.toml` records lattice path, SciBmad provenance, RF state, ordered
  labels, active inventory, source boundary, units, seed, and input path.

Validated smoke (run from `CESR Project/`):

```powershell
$env:JULIA_PKG_PRECOMPILE_AUTO='0'
julia --project=. `
  'dataset_benchmark/orbit/error_analysis/sextupole_detector_contributions/run_sextupole_detector_contributions.jl' `
  --ring=latest --trials=2 `
  --output-dir='dataset_benchmark/orbit/error_analysis/sextupole_detector_contributions/smoke/latest_cesr'
```

The 2026-08-20 smoke produced 76 active normal sextupoles, 144 detectors, and
103 selected controls. The reconstruction summary all-element relative vector
closure was `1.378e-14` (the largest per-direction closure was `1.393e-14`);
all signed projections and source-response checks were finite. These are
endpoint/integration checks, not production or paper statistics.

Production command (not run here):

```powershell
$env:JULIA_PKG_PRECOMPILE_AUTO='0'
julia --project=. `
  'dataset_benchmark/orbit/error_analysis/sextupole_detector_contributions/run_sextupole_detector_contributions.jl' `
  --ring=latest --trials=100 --seed=20260804 --base-kick-rad=5e-6 `
  --output-dir='dataset_benchmark/orbit/error_analysis/sextupole_detector_contributions/results/latest_cesr'
```

The existing `results/` files without the `latest_cesr` component remain
historical old-ring artifacts; see [`README_archived.md`](README_archived.md).

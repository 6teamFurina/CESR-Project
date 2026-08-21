# Sextupole source--beta--phase/Green-function predictor: latest CESR ring

Status: `smoke validated 2026-08-20`; production is intentionally not run in
this checkout. The validated latest-lattice integration smoke uses two matched
directions.

The nominal and direction-matched RF-on optics exporters use SciBmad and the
validated latest lattice
[`Latest_Lattice/latest_cesr_scibmad_repaired.jl`](../../../../Latest_Lattice/latest_cesr_scibmad_repaired.jl).
They discover active normal sextupoles and detector markers at runtime and
write ring-scoped outputs under `results/latest_cesr/`. Direction optics use
the same seeded simultaneous horizontal/vertical steering directions as the
signed contribution runner.

The analyzer consumes the exact total-vector contribution table rather than
the old internal-exposure or block-share files. For each plane it compares:

1. the absolute local normal-sextupole source kick;
2. source times the detector `sqrt(beta_s beta_d)/(2|sin(pi Q)|)` envelope;
3. source times the phase-aware Green function
   `sqrt(beta_s beta_d)/(2 sin(pi Q)) cos(2 pi |Δphi| - pi Q)`.

The exact source contribution remains the full coupled six-dimensional
complete-element SciBmad result. The Green-function predictors are same-plane
uncoupled ranking proxies; their residual disagreement can therefore reflect
coupling, finite-length source terms, and non-sextupole sources. The report
labels nominal versus direction-matched optics, uses phase in turns, and keeps
signed source fields separate from positive predictor magnitudes.

Files in a latest run include `nominal_optics_points.csv`,
`direction_optics_points.csv`, `direction_optics_tunes.csv`,
`sextupole_detector_optics.csv`, `direction_element_correlation_data.csv`,
`element_correlation_data.csv`, `correlation_summary.csv`, two SVG panels,
`RESULTS.md`, and `metadata.toml`. Metadata records the lattice path, ring
id, SciBmad engine/version, RF/branch state, detector/sextupole/control and
observable counts, direction/trial count, seed/base kick, source boundary,
phase units, predictor definitions, and the three cross-validated input
metadata paths.

Validated smoke (run from `CESR Project/`):

```powershell
$env:JULIA_PKG_PRECOMPILE_AUTO='0'
julia --project=. `
  'dataset_benchmark/orbit/error_analysis/sextupole_beta_phase_correlation/export_nominal_optics.jl' `
  --ring=latest --output-dir='dataset_benchmark/orbit/error_analysis/sextupole_beta_phase_correlation/smoke/latest_cesr'
julia --project=. `
  'dataset_benchmark/orbit/error_analysis/sextupole_beta_phase_correlation/export_direction_optics.jl' `
  --ring=latest --trials=2 --output-dir='dataset_benchmark/orbit/error_analysis/sextupole_beta_phase_correlation/smoke/latest_cesr'
```

After the detector smoke and the two optics tables exist, run:

```powershell
& 'C:\Users\JoeyN\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' `
  'dataset_benchmark/orbit/error_analysis/sextupole_beta_phase_correlation/analyze_beta_phase_correlation.py' `
  --optics 'dataset_benchmark/orbit/error_analysis/sextupole_beta_phase_correlation/smoke/latest_cesr/nominal_optics_points.csv' `
  --optics-metadata 'dataset_benchmark/orbit/error_analysis/sextupole_beta_phase_correlation/smoke/latest_cesr/nominal_optics_metadata.toml' `
  --direction-optics 'dataset_benchmark/orbit/error_analysis/sextupole_beta_phase_correlation/smoke/latest_cesr/direction_optics_points.csv' `
  --direction-tunes 'dataset_benchmark/orbit/error_analysis/sextupole_beta_phase_correlation/smoke/latest_cesr/direction_optics_tunes.csv' `
  --contributions 'dataset_benchmark/orbit/error_analysis/sextupole_detector_contributions/smoke/latest_cesr/sextupole_direction_contributions.csv' `
  --closure 'dataset_benchmark/orbit/error_analysis/sextupole_detector_contributions/smoke/latest_cesr/direction_closure.csv' `
  --output-dir 'dataset_benchmark/orbit/error_analysis/sextupole_beta_phase_correlation/smoke/latest_cesr'
```

The 2026-08-20 integration smoke produced 76 active normal sextupoles, 144
detectors (288 x/y observables), and 103 selected controls; it has 152
direction-element rows and a maximum all-element vector closure of `1.393e-14`
(`1.378e-14` in the aggregate reconstruction summary). The predictor outputs
are finite and the report includes pooled, per-direction quantile, and
element-level aggregation. Those correlations are integration evidence for
two matched directions only, not paper statistical conclusions.

Production is the same sequence with `--trials=100`,
`results/latest_cesr/`, and the detector production contribution/closure
tables. Do not mix trial tables, lattice versions, or old-ring files.

The existing `results/` files without the `latest_cesr` component remain
historical old-ring artifacts; see [`README_archived.md`](README_archived.md).

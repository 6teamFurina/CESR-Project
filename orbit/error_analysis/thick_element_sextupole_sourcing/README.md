# Complete-element nonlinear-error attribution: latest CESR ring

Status: `latest_cesr` one-direction x/y smoke validated on 2026-08-20.
Production reruns are intentionally not part of this smoke.

The runner uses the maintained SciBmad lattice
[`Latest_Lattice/latest_cesr_scibmad_repaired.jl`](../../../Latest_Lattice/latest_cesr_scibmad_repaired.jl).
At runtime it obtains the selected ring's control, detector, complete-element,
coordinate, and normal-sextupole registries. No element or detector count is
hard-coded.

For each Gaussian horizontal/vertical corrector-direction pair, the leading
second-order nonlinear detector vector is formed from the internally computed
Hessian contractions. The exact complete-element local source is

```text
g_j = S_exit,j - A_j S_entrance,j
```

and the source vectors are propagated around the periodic ring. The primary
outputs are the summed nonlinear target and its complete-element/family
attribution. Internal Hessian blocks are not exported as paper-facing block
shares, and no third-order calculation is performed.

## Outputs

Each run is ring- and plane-scoped, for example
`horizontal_results/latest_cesr/` or `vertical_results/latest_cesr/`.
Smoke outputs are kept separately under `*_results/latest_cesr_smoke/`.

- `element_contribution_summary.csv`: one row for every complete lattice
  element, with runtime element type, signed total projection `eta_total`, and
  vector-magnitude ratio.
- `family_contribution_summary.csv`: mutually exclusive runtime element-family
  totals; signed projections add to the all-element projection.
- `element_direction_contributions.csv` and
  `family_direction_contributions.csv`: direction-level vectors and
  projection numerators.
- `direction_closure.csv`: total target norm, all-element vector closure,
  signed projection, and family partition closure for each direction.
- `reconstruction_summary.csv`: aggregate closure and finite-value checks.
- `metadata.toml`: lattice/ring provenance, dynamic control and detector
  registries, element/family inventory, units, target definition, and closure
  convention.

Run the calculator and report renderer from the `CESR Project` directory:

```powershell
julia --project=. orbit/error_analysis/thick_element_sextupole_sourcing/run_thick_element_sourcing.jl `
  --ring=latest --output-plane=x --trials=100

python orbit/error_analysis/thick_element_sextupole_sourcing/analyze_thick_element_sourcing.py `
  orbit/error_analysis/thick_element_sextupole_sourcing/horizontal_results/latest_cesr
```

The default output directory follows the selected plane and ring. Supply
`--inputs=...` and `--output-dir=...` for a labeled smoke or custom registry.

## Latest-ring smoke

The smoke used `smoke_corrector_samples_2.csv`, one direction, RF-on, seed
`20260820`, and `base_kick_rad=5e-6`. The latest repaired ring reported 103
controls, 144 detectors, 1177 complete elements, 10 runtime source families,
and 76 active normal sextupoles.

| plane | target RMS [m] | all-element vector closure | family partition relative closure | all-element signed projection |
|---|---:|---:|---:|---:|
| x | `9.1417749e-6` | `4.23146e-15` | `7.54499e-15` | `1.0000000000000033` |
| y | `7.7466703e-6` | `1.27610e-14` | `1.56494e-14` | `1.0000000000000018` |

These checks are vector closures of the summed target, not positive error
shares. Signed family projections can cancel, while magnitude ratios are not
additive because family vectors interfere. Existing results in
`horizontal_results/`, `vertical_results/`, and `family_results/` are legacy
artifacts. The previous study description is preserved in
[`README_archived.md`](README_archived.md).

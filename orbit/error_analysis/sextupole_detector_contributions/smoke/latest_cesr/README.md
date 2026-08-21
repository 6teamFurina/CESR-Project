# latest CESR signed-sextupole smoke

Validated 2026-08-20 with SciBmad and the default latest lattice. This is a
two-direction matched smoke (`trials=2`); production was not run.

- Cardinality: 1177 lattice elements, 76 active normal sextupoles, 144
  detectors (288 x/y observables), and 103 selected controls.
- Aggregate all-element relative vector closure:
  `1.378e-14` in `reconstruction_summary.csv`.
- Largest direction-level all-element closure: `1.393e-14` in
  `direction_closure.csv`.
- Signed projections, source-response checks, and all CSV numeric fields are
  finite. The normal-sextupole residual is retained explicitly; signed
  projections are not positive shares.

Key files are `sextupole_direction_contributions.csv`,
`sextupole_contribution_summary.csv`, `direction_closure.csv`,
`reconstruction_summary.csv`, and `metadata.toml`. The combined source has no
`hh`/`hv`/`vv` block-share or third-order output. See the
[`sextupole_detector_contributions` README](../../README.md) for commands and
scope.

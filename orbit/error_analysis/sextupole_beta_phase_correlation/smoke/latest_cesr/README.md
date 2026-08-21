# latest CESR source--beta--phase smoke

Validated 2026-08-20 with SciBmad and the default latest lattice. This is a
two-direction matched integration smoke (`trials=2`); production was not run.

- Cardinality: 76 active normal sextupoles, 144 detectors (288 x/y
  observables), and 103 selected controls.
- The combined analyzer has 152 direction-element rows and 10 finite predictor
  summary rows (five predictors in each plane).
- Maximum direction-level all-element vector closure is `1.393e-14`; the
  aggregate reconstruction summary is `1.378e-14`.
- Pooled, per-direction quantile, and element-level correlations are included
  as integration evidence for these two matched directions only, not as paper
  statistical conclusions.

Key files are `RESULTS.md`, `correlation_summary.csv`,
`direction_element_correlation_data.csv`, `element_correlation_data.csv`, and
`metadata.toml`. The combined metadata includes cross-validated ring/lattice,
SciBmad/RF/branch, cardinality, seed/base-kick, and all three input-metadata
paths. The source-only predictor uses local source kicks that already include
`K2L`; it does not multiply `K2L` a second time. See the
[`sextupole_beta_phase_correlation` README](../../README.md) for commands and
scope.

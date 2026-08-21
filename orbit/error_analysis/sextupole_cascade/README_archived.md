# Vertical cubic sextupole-cascade experiment

This experiment tests whether the confirmed odd cubic
`vertical-only -> detector y` response is produced by cascaded sextupole
dynamics.

All 76 active order-2 multipole fields are multiplied by a common `lambda2`.
The same paired corrector directions and radii are reused at every strength.
For every modified lattice, the nominal closed orbit, closed-orbit response,
and full detector response are recomputed before extracting

```text
odd_nl(rho) = C3 * rho^3 + C5 * rho^5.
```

From `CESR Project`, run:

```console
julia --project=. orbit/error_analysis/sextupole_cascade/run_sextupole_cascade_experiment.jl
python orbit/error_analysis/sextupole_cascade/analyze_sextupole_cascade.py \
  orbit/error_analysis/sextupole_cascade/results/sextupole_cascade_vectors.csv
```

The targeted wiggler four-corner control is:

```console
julia --project=. orbit/error_analysis/sextupole_cascade/run_wiggler_corner_experiment.jl
python orbit/error_analysis/sextupole_cascade/analyze_wiggler_corners.py \
  orbit/error_analysis/sextupole_cascade/wiggler_corner_results/wiggler_corner_vectors.csv
```

The causal signature for a sextupole cascade is vector-level closure under

```text
C3(lambda2) = C3(0) + lambda2^2 * [C3(1) - C3(0)].
```

`C3(0)` is retained explicitly: if it is nonzero, sextupoles are a material
but non-unique source, and the remaining cubic vector must be attributed to
wigglers or other nonlinear lattice models.

## Formal result

The 100-direction global strength scan shows that a full vector model

```text
C3(lambda2) = A0 + A1*lambda2 + A2*lambda2^2
```

closes all five strength points to a maximum residual of `0.2074%` of the
nominal vector norm. The global component norm fractions, with signed
projections onto the nominal vector in parentheses, are `A0 = 1.1949
(1.1454)`, `A1 = 0.5012 (-0.3707)`, and `A2 = 0.3905 (0.2253)`. Thus a
two-sextupole cascade is measurable but is not the dominant signed source;
the constant non-sextupole component is large and the one-sextupole cross
component cancels part of it.

The wiggler four-corner decomposition gives global norm fractions and signed
projections of `0.6117 (0.5789)` for the non-sextupole/non-wiggler residual,
`0.3044 (-0.1056)` for sextupole-only, `0.6046 (0.5666)` for wiggler-only,
and `0.0958 (-0.0399)` for their interaction. The next attribution target is
therefore the residual left with both `K2=0` and the two planar wigglers off.

See [`results/SEXTUPOLE_CASCADE_RESULTS.md`](results/SEXTUPOLE_CASCADE_RESULTS.md)
and
[`wiggler_corner_results/WIGGLER_CORNER_RESULTS.md`](wiggler_corner_results/WIGGLER_CORNER_RESULTS.md).

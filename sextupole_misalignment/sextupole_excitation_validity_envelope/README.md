# Sextupole local-orbit Taylor validity envelope

This experiment estimates, separately for all 76 active CESR sextupoles, how
far the target-sextupole orbit can be represented by a local Taylor map in the
two model-based corrector bump knobs and the target sextupole strength change.
It uses the repaired latest SciBmad lattice and does not use Bmad/Tao data for
any quantitative result.

## Variables and observable

The normalized map inputs are

\[
q_x = b_x/(1\ {\rm mm}),\qquad
q_y = b_y/(1\ {\rm mm}),\qquad
q_k = \Delta K_2/(0.10\ {\rm m}^{-3}).
\]

Here, `b_x` and `b_y` are requested local-orbit bumps made by the complete
model-based corrector vectors in
`../quadrupole_affinity/exact_11_triplet_validation/results/bump_knobs/local_bump_knobs.csv`.
The individual correctors are therefore not treated as independent Taylor
variables. Their fields are linear functions of `b_x` and `b_y`, and the
compact result table reports the largest absolute model corrector-field change
at each accepted bump boundary. These are native SciBmad `Kn0`/`Ks0` field
parameters, not power-supply currents or calibrated machine readbacks.

The Taylor outputs are all six coordinates at both the entrance and exit of
the target sextupole. The validity gate uses the four transverse positions and
four transverse slopes:

- maximum entrance/exit position error no larger than 1 micrometer;
- maximum entrance/exit slope error no larger than 1 microradian;
- exact RF-on closed orbit converged for every sign at that radius.

The 1-micrometer/1-microradian pair is an explicit, configurable working
definition of Taylor validity. The position gate is deliberately small
compared with the tens-of-micrometers magnetic-center error scale in the
maintained inverse studies; it is not a machine-protection threshold.

The reported limit is the last tested radius for which every signed state in a
family passes. `first_fail` is the next tested radius. A missing `first_fail`
means that the last passing value is only a lower bound at the scan cap, not a
located boundary.

## Calculation

`generate_gtpsa_local_orbit_map.jl` constructs a total-order-four
parameter-dependent periodic fixed point and saves the target entrance/exit
Taylor coefficients. Direct high-order GTPSA propagation is stable for only a
target-dependent subset: 60 of 76 targets completed directly, while at the
other 16 targets the GTPSA C layer can terminate on
an `invalid domain sqrt(0)` error. Those targets use a full-rank total-order-four
fit to the inner exact SciBmad scan. The direct subset is retained as a
derivative-level diagnostic rather than silently presented as all-target
coverage.

`generate_exact_local_orbit_validation.jl` independently solves and tracks 651
exact states per target (49,476 total). It includes axes, x-y diagonals,
maintained-ratio x/K2 and y/K2 rays, fixed-K2 bump rays, fixed-bump K2 rays, and
inner three-parameter points. Radii span 0.25 through 20 in the normalized
coordinates. Of the 49,476 states, 49,362 converge. The 114 failures occur only
on large outer rays and are counted as failed validity points; all maintained
protocol states converge.

## All-ring results

The table below gives the minimum last-passing fourth-order limit across the
76 sextupoles. It is a conservative discrete Taylor-validity boundary, not a
machine-operating limit.

| varied family | all-ring last pass | next tested fail | median last pass | limiting target(s) |
|---|---:|---:|---:|---|
| x bump at `abs(delta K2) = 0.10 m^-3` | 3.0 mm | 4.0 mm | 10.0 mm | `SEX_39W` |
| y bump at `abs(delta K2) = 0.10 m^-3` | 4.0 mm | 5.0 mm | 10.0 mm | `SEX_18W` |
| `abs(delta K2)` at `abs(x bump) = 1.5 mm` | 0.60 m^-3 | 0.80 m^-3 | 2.0 m^-3 | `SEX_18W` |
| `abs(delta K2)` at `abs(y bump) = 1.5 mm` | 0.80 m^-3 | 1.00 m^-3 | 2.0 m^-3 | `SEX_18W`, `SEX_24E` |
| bump and K2 scaled together from `(1.5 mm, 0.10 m^-3)` | 2.5 times = `(3.75 mm, 0.25 m^-3)` | 3 times = `(4.5 mm, 0.30 m^-3)` | 4.5 times | `SEX_18W`, `SEX_39W` |

All 608 signed maintained-protocol states pass the fourth-order gate. On the
direct-GTPSA subset, where this is a fully independent exact comparison, the
worst maintained-protocol errors are 0.00753 micrometer in position and
0.00199 microradian in slope. `validate_results.py` additionally fits each
target using only the normalized cube `max(abs(qx), abs(qy), abs(qk)) <= 1`
and holds the maintained `abs(q_bump) = 1.5` protocol outside that fit. The
held-out all-target maxima are 0.00697 micrometer and 0.00330 microradian.

Second order is not a uniform all-ring model at the same gate. It passes only
586 of 608 maintained signed states; the failures occur at `SEX_18W`,
`SEX_25W`, `SEX_27W`, `SEX_39W`, `SEX_28E`, and `SEX_25E`. Its conservative
all-ring joint scale is only 0.5, corresponding to 0.75 mm and 0.05 m^-3.
Fourth order should therefore be used for the maintained excitation or every
second-order target must use its own smaller bound.

For a guarded synthetic test envelope, a common scale of 2.0 (3.0 mm and
0.20 m^-3 along the maintained-ratio rays) remains one sampled step inside the
worst fourth-order boundary. This is a model-validity recommendation only.
The exact model already reaches full-ring maxima of 5.67 mm in x and 8.56 mm
in y at scale 1; at scale 2 those maxima are 11.39 mm and 17.20 mm. This is why
the scale-2 statement must not be converted directly into a CESR machine scan
without a separate orbit/aperture and actuator review.

## Result files

- `results/analysis/per_sextupole_order4_envelope.csv`: compact 76-target
  fourth-order last-pass/first-fail table, including model corrector-field
  demand.
- `results/analysis/per_sextupole_limits.csv`: order-one, order-two, and
  order-four target summaries.
- `results/analysis/family_validity_limits.csv`: every target, Taylor order,
  and scan family.
- `results/analysis/state_taylor_errors.csv`: exact state-by-state residuals.
- `results/analysis/map_diagnostics.csv`: direct/fitted source, fit rank,
  conditioning, and direct-versus-fit coefficient comparison.
- `results/analysis/summary.json`: aggregate statistics and provenance.
- `results/analysis/validation_summary.json`: independent inventory, prefix,
  maintained-protocol, and small-cube cross-validation checks.
- `results/analysis/small_cube_cross_validation.csv`: per-target held-out
  maintained-protocol errors from the independent small-cube fit.

## Reproduction

From the CESR project root:

```powershell
julia --project=. sextupole_misalignment\sextupole_excitation_validity_envelope\generate_exact_local_orbit_validation.jl --targets=all --overwrite=true
julia --project=. sextupole_misalignment\sextupole_excitation_validity_envelope\generate_gtpsa_local_orbit_map.jl --target=SEX_09AW
python sextupole_misalignment\sextupole_excitation_validity_envelope\analyze_validity_envelope.py
python sextupole_misalignment\sextupole_excitation_validity_envelope\validate_results.py
```

Run the direct map generator once per target; a failed direct target is filled
by the exact-scan fit during analysis. `CESR_VALIDITY_LATTICE` can point the map
generator to an equivalent compatible copy of the latest lattice when a
working-tree lattice edit is in progress.

## Interpretation limits

These limits test Taylor truncation plus exact model closed-orbit convergence.
They do not impose CESR corrector or sextupole power-supply ranges, current-to-
field calibration, hysteresis, aperture, lifetime, interlocks, settling, or
operator constraints, and they contain no randomized latent machine errors.
They must be intersected with those machine limits before use in an experiment.

The latest lattice also emits the documented straight-multipole-in-curved-
reference warning. In particular, curved-DQX pitch responses are not exact
Bmad-parity physics in current SciBmad/BeamTracking. No girder pitch is varied
here, but that model qualification remains part of the lattice provenance.

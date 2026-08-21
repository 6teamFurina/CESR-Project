# Targeted bump-by-K2 inversion: initial smoke results

## Scope

These runs validate the end-to-end experiment plumbing for one target magnet,
`SEX_08W`. They are not the planned ten-magnet protocol study and do not
establish experimental precision. No measurement noise or covariance is used.
The inverse uses the nominal completed GTPSA response map with structural
per-observable RMS scaling.

The assigned truth in both cases is

```text
x offset = +350 micrometers
y offset = -250 micrometers
```

Each exact scan contains five target-orbit bump commands (center and `+/-0.5
mm` in each plane) and three target-strength settings (`delta K2 = 0,+/-0.01
m^-3`). All 15 RF-on scalar closed-orbit/Twiss states converged in both cases.

## Bump-knob result

The model-based 119-corrector knobs close the requested 1 mm target
displacement to machine precision in the nominal first-order map. Detector
orbit leakage RMS is `76.2 micrometers` for the horizontal knob and `69.3
micrometers` for the vertical knob. Maximum corrector kicks for the 1 mm design
amplitude are `0.608 mrad` and `0.326 mrad`, respectively. These knobs are
therefore suitable for the synthetic protocol smoke test but are not claimed
to satisfy CESR hardware or operations limits.

## Truth recovery

Selected two-response-mode estimates are:

| background | observable view | estimated x (um) | estimated y (um) | x error (um) | y error (um) | 2D error (um) |
|---|---|---:|---:|---:|---:|---:|
| target only | orbit | 352.111 | -257.335 | +2.111 | -7.335 | 7.633 |
| target only | orbit + phase + coupling + tune | 364.710 | -251.469 | +14.710 | -1.469 | 14.783 |
| other 75 sextupoles at 300 um RMS | orbit | 360.927 | -223.410 | +10.927 | +26.590 | 28.748 |
| other 75 sextupoles at 300 um RMS | orbit + phase + coupling + tune | 264.600 | -123.888 | -85.400 | +126.112 | 152.307 |

For this linear two-column estimator, projecting the slope vector onto the two
alignment response modes and solving in full observation space give identical
point estimates by construction. Keeping both views is an algebraic check;
their robustness differs only after a noise/nuisance model or a different
nonlinear estimator is introduced.

## Interpretation

The target-only case validates data alignment, signs, bump convention, exact
scan generation, response-slope fitting, and direct truth comparison. The
all-sextupole background case is an intentionally small diagnostic with only
one random realization. It already shows that adding more nominal optics
channels without state conditioning or measured covariance can increase bias.
It does not show that optics should be discarded.

The next protocol stage must compare:

1. a fixed nominal response map;
2. a response map relinearized or conditioned on the measured baseline
   orbit/optics;
3. observable covariance whitening; and
4. multiple random background realizations and representative target magnets.

Only those ensemble inverse errors should determine the final observable set
and input representation.

## Family-wise nominal-model mismatch

A follow-up truth-based diagnostic compared each family's measured exact scan
slope with the nominal local response-map prediction at the known offset. The
RMS residual was divided by the RMS predicted truth signal using the same
structural per-observable scales as the smoke inverse.

| background | orbit | phase | beta | alpha | coupling | tune |
|---|---:|---:|---:|---:|---:|---:|
| target only | 0.075 | 1.776 | 1.956 | 1.357 | 0.087 | 0.314 |
| other 75 sextupoles at 300 um RMS | 0.124 | 1.983 | 2.238 | 1.509 | 0.675 | 0.347 |

These ratios explain why the present naive concatenation favors orbit. The
nominal orbit response remains relatively close to the exact bumped scan,
whereas phase/beta/alpha are strongly affected by state-dependent propagation
and by the corrector bump's passage through other nonlinear elements. Coupling
is close in the target-only case but changes materially in the misaligned
background. Tune and coupling alone also give ill-conditioned two-offset
matrices, so their family-only position estimates are not meaningful joint
`x/y` estimators.

The appropriate optics algorithm is therefore not the current row-concatenated
nominal regression. It should reconstruct local normal/skew feeddown using a
source-to-observable matrix conditioned or relinearized at each baseline/bump,
apply measured covariance by observable block, and then locate the common
two-dimensional zero crossing. Direct nonlinear MAP using the complete
bump-by-`K2` forward scan is the comparison estimator.

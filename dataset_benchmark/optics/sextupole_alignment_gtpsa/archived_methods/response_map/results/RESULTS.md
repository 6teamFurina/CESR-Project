# CESR per-sextupole alignment GTPSA result

## Result

The calculation completed for all 76 active normal CESR sextupoles. Each
sextupole was independently parameterized by

```text
delta Kn2, delta x_offset, delta y_offset
```

using `Descriptor(6, 3, 3, 2)`. The output contains responses of 12 quantities
at each of 99 detector markers and the three RF-on ring eigentunes.

- Successful sextupoles: `76/76`.
- Coefficient rows: `90,516`.
- Rows per sextupole: `1,191 = 99 * 12 + 3`.
- Non-finite first- or second-order coefficients: `0`.
- Stable per-sextupole Twiss time after process-local compilation: median
  `6.232 s` (minimum `4.474 s`).
- Sum of the four workers' Twiss times: `2,232.2 s`; the four process-local
  first calls dominate the maximum of `420.8 s`.

The accumulated `phi_1` and `phi_2` responses use `DET_00W` as a fixed phase
reference. This removes the parameter-dependent additive phase gauge in the
raw parameterized Twiss output and matches a difference-phase measurement.

## Independent numerical validation

`SEX_08W` was checked using exact scalar RF-on closed-orbit and Twiss
calculations at the four corners of each mixed derivative. The steps were

```text
delta Kn2 = +/- 0.01 m^-3
delta offset = +/- 0.1 mm
```

The four-corner finite differences were used only for validation, not for the
production coefficients. Relative L2 residuals by individual output channel
range from `6.7e-8` to `1.05e-4`. The largest relative value belongs to the
very weak `orbit_x` response to vertical offset, whose GTPSA norm is only
`5.67e-4` in the derivative units. All phase, beta, alpha, coupling, and tune
channels are within `3.4e-5`; most are within `3e-5`.

The validation establishes the coefficient normalization and signs for the
tested magnet. It is not a substitute for spot checks at additional
sextupoles or for a finite-amplitude validity-domain scan.

## Response geometry

For each sextupole and each same-unit observable family, the two mixed-response
columns

```text
d2 observable / (d Kn2 d x_offset)
d2 observable / (d Kn2 d y_offset)
```

were assembled into a two-column matrix. The following are unweighted
coefficient-space condition numbers:

| observable family | P10 | median | P90 |
|---|---:|---:|---:|
| orbit | 1.289 | 1.875 | 2.595 |
| phase | 7.954 | 19.917 | 74.286 |
| beta | 1.175 | 1.553 | 2.753 |
| alpha | 1.404 | 2.354 | 4.500 |
| coupling | 152.862 | 378.069 | 2136.302 |
| tune | 361.447 | 712.215 | 2706.982 |

The interpretation is physical rather than paradoxical: normal-optics
quantities strongly constrain the normal feed-down associated with horizontal
offset, while coupling strongly responds to the skew feed-down associated with
vertical offset. Coupling or tune alone therefore gives an imbalanced
two-column matrix. The orbit, beta, and alpha families contain two much more
balanced columns in this nominal coupled CESR model.

These condition numbers must not be treated as experimental position
uncertainties. They are computed separately within same-unit families and do
not include BPM/phase/coupling measurement covariance, `Kn2` calibration
uncertainty, nuisance magnet errors, missing channels, or simulator-to-machine
discrepancy. A physically meaningful joint inverse condition number requires
whitening by those measurement and model-error covariances.

## Per-sextupole thin SVD

Each sextupole's complete mixed-response dictionary was also assembled as

```text
R_j = [d2 O / (d Kn2 d x_offset);
       d2 O / (d Kn2 d y_offset)]
```

with shape `2 x 1191`. A thin SVD was calculated for all 76 matrices under two
scalings: the raw mixed-unit coefficients as an algebraic reference, and a
structural `observable_rms` scaling in which each named observable is divided
by one global RMS mixed-response scale calculated across all magnets,
locations, and the two offset rows.

All 76 matrices have numerical rank two under both scalings. For the
`observable_rms` result, the condition-number minimum/median/P90/maximum is
`1.044 / 1.336 / 2.031 / 3.208`; the corresponding median
`sigma_2 / sigma_1` is `0.749`. The largest condition number belongs to
`SEX_39W`, for which `sigma_2 / sigma_1 = 0.312`. Thus no sextupole loses one
of its two local alignment directions in this nominal structural calculation.

The complete thin-SVD reconstruction and orthogonality checks close to machine
precision: the maximum relative reconstruction error is `6.62e-16`, and the
maximum parameter-rotation and mode-weight orthogonality errors are
`1.25e-15` and `1.74e-15`, respectively.

For a consistently scaled measured response vector `y_j`, the two compressed
mode amplitudes are `z_j = Q_j^T y_j`, where
`R_j = P_j Sigma_j Q_j^T`. Recovering experimental position precision still
requires replacement of `observable_rms` by measured covariance whitening and
the addition or marginalization of nuisance-response modes.

A stacked `152 x 1191` SVD is not the global inverse matrix for independent
one-at-a-time sextupole scans. Those scans produce 76 distinct observation
blocks, so their ideal joint design matrix is block structured with 76 times
as many observation rows. A `152 x 1191` stack can later be useful to ask
whether response shapes across magnets share a common output basis, but it is
not required by the present well-conditioned local inversions and should be
revisited after experimental whitening and nuisance columns are defined.

## Interpretation boundary and next test

This result shows that low-dimensional per-magnet GTPSA can generate the local
mixed response needed for a sextupole alignment inverse model. It does not yet
show that all 152 offsets can be recovered from realistic noisy machine data.
The next bounded study should:

1. choose realistic `Kn2` scan amplitudes and offset ranges;
2. evaluate the second-order Taylor prediction at held-out finite offsets;
3. add `Kn2` calibration, BPM, quadrupole, skew/roll, and drift nuisance terms;
4. whiten the response by realistic measurement covariance;
5. compare SVD/Bayesian linear inversion with a learned inverse model; and
6. validate against a subset of archived CESR sextupole scan data if available.

## Artifacts

- `full/alignment_coefficients.csv`: all first- and second-order coefficients.
- `full/mixed_response_identifiability.csv`: per-sextupole, per-family
  two-column response metrics.
- `full/results_summary.json`: machine-readable structural and timing summary.
- `full/local_response_svd_summary.csv`: singular values, rotations, rank, and
  numerical checks for both scalings and all 76 sextupoles.
- `full/local_response_svd_modes.csv`: the two complete observable-space mode
  weight vectors for every sextupole and scaling.
- `full/local_response_svd_scales.csv`: the structural observable RMS scales.
- `full/local_response_svd_summary.json`: aggregate local-SVD statistics and
  interpretation boundaries.
- `validation/mixed_coefficient_validation.csv`: pointwise GTPSA versus
  four-corner comparison.
- `validation/mixed_coefficient_validation_summary.csv`: validation grouped by
  observable and offset plane.

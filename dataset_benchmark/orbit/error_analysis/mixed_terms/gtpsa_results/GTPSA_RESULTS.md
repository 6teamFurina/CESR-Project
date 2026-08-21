# Direction-contracted GTPSA quadratic orbit response

Two dimensionless second-order GTPSA parameters drive each fixed horizontal
and vertical corrector direction pair,

```text
delta_k = (5 microrad) * (a * h_direction + b * v_direction).
```

The descriptor is `Descriptor(6, 2, 2, 2)`: six initial phase-space
variables, total order two, and two corrector parameters retained through
parameter order two.  The RF-on closed-orbit derivatives are obtained by
implicit differentiation of the one-turn fixed-point equation.  This is a
corrector-parameter GTPSA calculation, not the ordinary second-order
phase-space Twiss map.

For detector orbit `u`, the Taylor derivatives give

```text
Q_hh(rho) = rho^2 * (d2u/da2) / 2
Q_hv(rho) = rho^2 *  d2u/dadb
Q_vv(rho) = rho^2 * (d2u/db2) / 2.
```

The experiment uses the same 100 Gaussian unit-RMS H/V direction pairs and
seed as the finite-difference mixed-term experiment.  Each direction is
tracked once; other radii follow exactly from the quadratic `rho^2` scaling.

## Final median and P10--P90 statistics

At `rho = 1.13`, corresponding to `5.65 microrad` RMS in each active
corrector family, the detector RMS values across directions are reported as
`median [P10, P90]`:

| orbit | Q_hh (micrometre) | Q_hv (micrometre) | Q_vv (micrometre) |
|:---:|---:|---:|---:|
| X | 0.19505 [0.09322, 0.38443] | 0.005691 [0.003267, 0.009874] | 0.56188 [0.27609, 1.11208] |
| Y | 0.004118 [0.001922, 0.008302] | 0.68685 [0.35506, 1.08578] | 0.008793 [0.004164, 0.016240] |

For every direction and orbit plane, define the common-denominator
squared-norm shares

```text
f_ab = ||Q_ab||^2 / (||Q_hh||^2 + ||Q_hv||^2 + ||Q_vv||^2).
```

The adopted final share statistic is also `median [P10, P90]`:

| orbit | f_hh | f_hv | f_vv |
|:---:|---:|---:|---:|
| X | 9.14% [1.83%, 38.26%] | 0.0088% [0.0017%, 0.0313%] | 90.83% [61.72%, 98.15%] |
| Y | 0.0040% [0.0007%, 0.0231%] | 99.977% [99.930%, 99.991%] | 0.0163% [0.0044%, 0.0651%] |

Means remain available in `gtpsa_summary.csv`, but are not used as the primary
paper/report statistic because the direction distributions are broad and
skewed. For example, the X `f_vv` mean is `84.35%`, whereas the adopted median
is `90.83%`.

Because all three blocks are exactly quadratic in this calculation, their
shares are invariant with `rho`.  Their RMS values at any listed radius equal
the table above multiplied by `(rho / 1.13)^2`.

## Independent verification against four-sign finite differences

The implicit fixed-point residuals over all 100 directions are

- maximum first-order residual: `4.34e-19`;
- maximum second-order residual: `8.47e-22`.

At `rho = 1.13`, the mean relative differences between the GTPSA and
four-sign finite-difference RMS values are:

| block | X | Y |
|:---:|---:|---:|
| hh | 1.22e-6 | 4.56e-6 |
| hv | 1.91e-5 | 1.26e-5 |
| vv | 1.00e-5 | 1.64e-5 |

At the smallest radius, relative differences of the very small suppressed
blocks are more affected by floating-point cancellation; the largest
direction-level relative difference is `1.48e-3` for `Y vv`, while the
dominant blocks agree much more closely.  The comparison confirms the finite
difference normalization and the pure-block factor of one half.

The 100-direction GTPSA calculation took `120.1 s` for the core tracking and
implicit differentiation on the recorded machine. The GTPSA median
[P10, P90] values above are the adopted final results. The four-sign data are
retained as an independent validation dataset. The underlying direction-level
Q values, statistics, and finite-difference comparison are stored beside this
report.

## Project method convention

For subsequent corrector-response derivatives and analogous response
coefficients, use GTPSA parameterization with implicit differentiation of the
closed-orbit fixed point as the default numerical method. Use signed finite
differences as an independent validation, higher-order contamination check, or
fallback when the requested observable cannot yet propagate GTPSA types.
Direct nonlinear solves remain the reference for determining the amplitude
range over which a truncated response expansion is valid.

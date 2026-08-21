# Horizontal quadratic-block attribution

This directory preserves the corrector-space Hessian study removed from the
paper because its horizontal-imbalance interpretation is no longer followed
by a dedicated supporting experiment there.  It explains why equal-rms horizontal and
vertical corrector directions produce very different horizontal quadratic
orbit responses, and records the limits of the unsigned source-exposure test.

## Corrector-space Hessian context

For each of 100 fixed horizontal/vertical corrector-direction pairs, write

```text
Delta k = k0 (a h + b v),       k0 = 5 microrad.
```

Implicit differentiation of the RF-on closed-orbit fixed point gives

```text
Q_hh = (rho^2 / 2) u_,aa,
Q_hv =  rho^2      u_,ab,
Q_vv = (rho^2 / 2) u_,bb.
```

Midplane reflection makes the horizontal response even and the vertical
response odd under `v -> -v`.  It therefore permits `Q_hh,x` and `Q_vv,x`
while suppressing `Q_hv,x`; for the vertical output it permits `Q_hv,y` while
suppressing `Q_hh,y` and `Q_vv,y`.

The adopted squared-norm shares are

| orbit | `Q_hh` median [P10, P90], % | `Q_hv` median [P10, P90], % | `Q_vv` median [P10, P90], % |
|---|---:|---:|---:|
| `x` | 9.14 [1.83, 38.26] | 0.0088 [0.0017, 0.0313] | 90.83 [61.72, 98.15] |
| `y` | 0.0040 [0.0007, 0.0231] | 99.977 [99.930, 99.991] | 0.0163 [0.0044, 0.0651] |

The independent 6,401-state four-sign calculation agrees with the six GTPSA
block rms values to relative differences between `1.22e-6` and `1.91e-5` and
leaves `0.953%` of the exact vertical residual rms at `rho = 1.13`.  The full
block calculation and validation are
maintained in [`../mixed_terms/`](../mixed_terms/).

## Why the horizontal response is `Q_vv` dominated

Since the horizontal mixed block is negligible, define

```text
R_Q = ||Q_vv,x|| / ||Q_hh,x||.
```

Across the 100 directions,

```text
R_Q = 3.153 [1.271, 7.321].
```

Thus an equal-rms vertical-corrector direction produces a horizontal quadratic
detector vector about 3.15 times larger than an equal-rms horizontal-corrector
direction.  To test whether this originates at the nonlinear elements, the
first-order internal orbits are evaluated at all 76 active normal sextupoles.
The unsigned local-source exposures are

```text
E_h = sum_j |K2L_j| x_h(j)^2,
E_v = sum_j |K2L_j| y_v(j)^2.
```

Their direction-level ratio is

```text
R_E = E_v / E_h = 2.600 [0.716, 7.370].
```

After dividing out this exposure, the remaining response-efficiency ratio is

```text
R_eta = (||Q_vv,x|| / E_v) / (||Q_hh,x|| / E_h)
      = 1.260 [0.654, 2.137].
```

For every direction, `R_Q = R_E R_eta` algebraically, although their medians
need not multiply.  The direction-level association is strong:

```text
corr(log R_Q, log R_E)                  = 0.872,
corr(log R_Q, log orbit-exposure ratio) = 0.898.
```

This shows that unequal internal orbit exposure is the main first-stage
explanation of the 90.83% horizontal `Q_vv` share.

## Equal-source-exposure control

The H and V directions are symmetrically rescaled so that both source
exposures equal `sqrt(E_h E_v)`.  The resulting horizontal shares are

| block | equal-source-exposure median [P10, P90], % |
|---|---:|
| `Q_hh,x` | 38.62 [18.03, 70.04] |
| `Q_hv,x` | 0.0133 [0.0037, 0.0618] |
| `Q_vv,x` | 61.32 [29.95, 81.95] |

Equalizing the scalar exposure reduces the median `Q_vv,x` share from 90.83%
to 61.32%.  Most of the imbalance therefore comes from how equal corrector rms
samples the sextupoles, rather than from an intrinsic vertical coefficient of
a normal sextupole.  The remaining 61%/39% imbalance contains transport,
phase-cancellation, and non-sextupole effects that this scalar control cannot
separate.

## Element-local continuation

The maintained run also writes

```text
e_h,j = |K2L_j| x_h,j^2,
e_v,j = |K2L_j| y_v,j^2
```

for every direction and active normal sextupole.  Summing the 76 element
values closes each direction's `E_h` and `E_v` to a maximum relative residual
of `7.33e-16`.  The ratio of mean vertical to horizontal exposure is
`2.367927`; this differs from the median direction-level ratio `2.600` because
the two statistics aggregate directions differently.

The first four station pairs, `SEX_31`, `SEX_41`, `SEX_33`, and `SEX_25`,
provide 40.01% of the positive mean vertical-minus-horizontal excess.  East
and West totals agree within about one percent of their mean, indicating a
ring-symmetric optics pattern rather than a one-sided outlier.  When the sign
of `K2` is used only to label the alternating sextupole-location classes, the
negative-`K2` class has exposure ratio 5.0142, while the positive-`K2` class
has ratio 0.7464.  The excess therefore arises mainly at the alternating
negative-`K2` optics locations and is partially offset by the other class.

- [`element_results/top15_sextupole_exposure_differences.svg`](element_results/top15_sextupole_exposure_differences.svg)
  compares the 15 largest element-local differences.
- [`element_results/sextupole_exposures_strengths_ring.svg`](element_results/sextupole_exposures_strengths_ring.svg)
  shows exposure excess and signed `K2L` around the ring.
- [`element_results/ELEMENT_EXPOSURE_RESULTS.md`](element_results/ELEMENT_EXPOSURE_RESULTS.md)
  contains the complete element, station-pair, ring-side, and location-class
  tables.

## Evidence boundary

This experiment is an unsigned source-exposure attribution, not an additive
signed decomposition of the detector error.  In particular,

- `|K2L|u^2` removes the sign of the sextupole kick;
- it does not include detector beta functions, tune, phase advance, or
  cancellation between propagated element vectors;
- midpoint sampling approximates the source inside a thick sextupole; and
- the exact detector vector also contains bends, solenoids, fringe effects,
  geometric tracking, and the other nonlinear elements.

The later complete-element and beta/phase-resolved studies supersede this
proxy when ranking propagated detector contributions.  This directory remains
the causal record of why the equal-corrector-rms horizontal Hessian is strongly
`Q_vv` dominated.

## Reproduction

Run from the `CESR Project` environment:

```powershell
julia --project=. orbit/error_analysis/quadratic_x_attribution/run_internal_exposure_attribution.jl
python orbit/error_analysis/quadratic_x_attribution/analyze_element_exposure.py
```

The first-stage numerical report is
[`results/RESULTS.md`](results/RESULTS.md); the element-local report is linked
above.

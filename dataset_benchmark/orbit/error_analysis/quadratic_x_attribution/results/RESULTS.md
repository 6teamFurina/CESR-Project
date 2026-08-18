# X-quadratic internal-exposure result

Result recorded 2026-08-05. All intervals are `median [P10, P90]` over the
same 100 fixed H/V direction pairs used by the adopted second-order GTPSA
calculation.

## What the 90%/10% number means

Since the X `hv` share is negligible, define

```text
R_Q = ||Q_vv,x|| / ||Q_hh,x||,
f_vv,x approximately R_Q^2 / (1 + R_Q^2).
```

The measured `f_vv,x = 90.83 [61.72, 98.15]%` corresponds direction by
direction to

```text
R_Q = 3.153 [1.271, 7.321].
```

Thus the immediate physical question is why an equal-RMS vertical-corrector
direction generates an X quadratic response vector about 3.15 times larger
than an equal-RMS horizontal-corrector direction.

## Internal normal-sextupole exposure

For a normal sextupole, the leading horizontal kick contains opposite-sign
`x^2` and `y^2` terms. The relevant first-order internal orbits were obtained
with one all-corrector first-order GTPSA map and implicit differentiation of
the RF-on closed orbit. At each of 76 active normal-sextupole elements, the
orbit was sampled from the average of its entrance and exit TPSA maps.

The direction-level source-exposure ratio is

```text
R_E = sum_j |K2L_j| y_v(j)^2 / sum_j |K2L_j| x_h(j)^2
    = 2.600 [0.716, 7.370].
```

The remaining response efficiency after dividing out this exposure is

```text
R_eta = (||Q_vv,x||/E_v) / (||Q_hh,x||/E_h)
      = 1.260 [0.654, 2.137].
```

For every direction, `R_Q = R_E R_eta` algebraically. Medians need not
multiply exactly. Across directions, the Pearson correlations are

```text
corr(log R_Q, log R_E)                         = 0.872
corr(log R_Q, log orbit-exposure ratio)        = 0.898
```

The strong direction-level association shows that unequal internal orbit
exposure is the main explanation of the observed imbalance.

The all-corrector response-matrix implementation was checked against the
independent two-parameter, direction-by-direction GTPSA construction for the
first three direction pairs. Their source-exposure ratios agreed with a
maximum relative difference of approximately `1.1e-13`.

## Equal-source-exposure control

After symmetric rescaling to equal `sum |K2L|u^2` exposure, the X shares are

| block | equal-source-exposure share, % |
|---|---:|
| `hh` | 38.62 [18.03, 70.04] |
| `hv` | 0.0133 [0.0037, 0.0618] |
| `vv` | 61.32 [29.95, 81.95] |

Therefore the original 90.83% `vv` share is not primarily an intrinsic
direction coefficient of a normal sextupole. Most of it arises because equal
corrector RMS does not produce equal quadratic orbit exposure at the
nonlinear elements. A material residual remains: even after this scalar
exposure control, the median balance is about 61%/39%, not 50%/50%.

## Evidence boundary and next test

This is a first-stage attribution, not yet an element-level closure:

- `sum |K2L|u^2` intentionally removes the signs of `K2`, detector transport,
  and phase cancellation;
- midpoint sampling approximates the distributed source inside a thick
  element; and
- the exact detector `Q` includes wigglers, solenoid/tilt effects, geometric
  tracking, and all other nonlinear terms in the model.

The subsequent unsigned element-local decomposition is recorded in
[`../element_results/ELEMENT_EXPOSURE_RESULTS.md`](../element_results/ELEMENT_EXPOSURE_RESULTS.md).
It localizes the proxy imbalance but does not change this evidence boundary.
The next decisive calculation is a signed detector-vector reconstruction,
grouped by sextupole family or ring region, followed by recomputed-lattice
`K2` and wiggler scans/ablations. It must reproduce the baseline
`Q_hh,x` and `Q_vv,x` vectors rather than compare scalar norms alone. The
remaining 61%/39% balance should not be assigned to a particular element
family until that reconstruction closes.

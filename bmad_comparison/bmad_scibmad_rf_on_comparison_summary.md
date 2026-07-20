# CESR RF-On Bmad–SciBmad Comparison

## Configuration

- Beam energy: 5.289 GeV
- RF cavities: `RF_W1`, `RF_W2`, `RF_E1`, and `RF_E2`
- Voltage: 1.5 MV per cavity (6 MV total)
- RF frequency: 499.7669603784 MHz
- Bmad reference: `bmad_reference_rf_on_output.tar.gz`
- SciBmad lattice: `cesr.jl`, loaded with `load_cesr()` and enabled with
  `set_cesr_rf!(ring; on=true)`

The element comparison linearizes every SciBmad element with GTPSA about the
corresponding Bmad entrance orbit. This assigns discrepancies to the element
that creates them instead of hiding them in a one-turn product.

Relative matrix differences use the Bmad map as the normalization scale:

```text
relative max difference = max(abs(R_SciBmad - R_Bmad)) / max(abs(R_Bmad))
relative Frobenius difference = norm(R_SciBmad - R_Bmad) / norm(R_Bmad)
```

The reported percentage is `100 * relative difference`. This matrix-level
normalization is preferred to dividing individual entries because a Bmad entry
may be zero or very close to zero.

## Closed Orbit and Tunes

| Coordinate | Bmad | SciBmad | SciBmad - Bmad | Relative difference |
|---|---:|---:|---:|---:|
| x | -1.666982000e-5 | -1.668519860e-5 | -1.537860e-8 | 0.0923% |
| px | 2.389568950e-3 | 2.390112503e-3 | 5.435531e-7 | 0.0227% |
| y | 1.054880000e-6 | 1.054311223e-6 | -5.687774e-10 | 0.0539% |
| py | 1.777040000e-6 | 1.796930342e-6 | 1.989034e-8 | 1.1193% |
| z | -3.989291500e-4 | -3.993611670e-4 | -4.320170e-7 | 0.1083% |
| pz | -7.163370000e-6 | -7.803721195e-6 | -6.403512e-7 | 8.9392% |

The maximum closed-orbit coordinate difference is `6.404e-7` in `pz`. The
componentwise `py` and `pz` percentages are amplified by their small Bmad
reference values. A more representative whole-vector comparison gives a
relative 2-norm difference of `0.0390%` and a relative max-norm difference of
`0.0268%`.

| Tune | Bmad | SciBmad | SciBmad - Bmad | Relative difference |
|---|---:|---:|---:|---:|
| Qx (fractional) | 0.530014110 | 0.529926212 | -8.790e-5 | 0.0166% |
| Qy (fractional) | 0.578460931 | 0.579077486 | 6.166e-4 | 0.1066% |
| Qz (eigenphase magnitude) | 0.051738795 | 0.051738677 | -1.182e-7 | 0.000228% |

The Bmad transverse tunes are taken from the accumulated Twiss phases printed
by Tao. The longitudinal tune is estimated from the product of Tao's printed
per-element matrices. Those matrices contain only seven printed decimal
places, so their product is less accurate than the Bmad internal one-turn map.
The resulting Bmad eigenvalue moduli differ from one by at most `5.65e-6`.

## Element-by-Element Results

- Maximum local matrix discrepancy: `2.476e-4` at `B06E` (`R44`), equal to
  `0.007645%` after normalization by the Bmad local-map scale.
- Maximum cumulative matrix discrepancy: `1.267e-2` at `Q23E`, equal to
  `0.05221%` of the Bmad cumulative-map scale at that location.
- Maximum normalized cumulative matrix discrepancy: `0.16772%` at `Q39W`.
- Maximum isolated-element exit-orbit discrepancy: `2.526e-6` at `WIG_E`,
  equal to `0.02697%` of the largest Bmad orbit coordinate there.
- Maximum element-length discrepancy: `8.338e-7 m`, limited by printed Tao
  longitudinal positions; the maximum normalized length difference is
  `0.0000283%`.
- Bmad printed affine-map consistency: `1.629e-9`.

| Element | Local matrix difference | Local matrix relative | Exit-orbit difference | Exit-orbit relative |
|---|---:|---:|---:|---:|
| Q00W/CLEO_SOL overlap | 9.145e-5 | 0.005842% | 1.935e-7 | 0.003514% |
| WIG_W | 3.951e-5 | 0.001678% | 2.526e-6 | 0.025336% |
| RF_W1 | 4.857e-8 | 0.00000270% | 6.999e-12 | 0.000000064% |
| RF_W2 | 4.857e-8 | 0.00000270% | 1.571e-11 | 0.000000158% |
| RF_E1 | 4.245e-8 | 0.00000236% | 8.176e-12 | 0.000000079% |
| RF_E2 | 4.191e-8 | 0.00000233% | 7.410e-12 | 0.000000065% |
| WIG_E | 3.775e-5 | 0.001603% | 2.526e-6 | 0.026972% |
| Q00E/CLEO_SOL overlap | 9.141e-5 | 0.005840% | 4.189e-7 | 0.017544% |

The RF cavity maps agree to approximately `5e-8` in absolute terms and
`2.3e-6%` to `2.7e-6%` in relative terms, so RF translation is not a
significant source of the present Bmad–SciBmad difference. The dominant local
matrix discrepancies remain in sector bends, followed by the combined
solenoid–quadrupole overlaps and the two wigglers. The largest affine/orbit
difference remains associated with the wiggler model.

## RF-On Versus RF-Off

| Metric | RF off absolute | RF off relative | RF on absolute | RF on relative |
|---|---:|---:|---:|---:|
| Maximum local matrix difference | 2.478e-4 | 0.007651% | 2.476e-4 | 0.007645% |
| Maximum cumulative matrix difference | 1.258e-2 | 0.165858% | 1.267e-2 | 0.167722% |
| Maximum isolated-element orbit difference | 2.526e-6 | 0.026936% | 2.526e-6 | 0.026972% |

Turning on the RF changes the absolute maximum cumulative matrix discrepancy
by less than one percent and changes the normalized maximum from `0.165858%`
to `0.167722%`. It does not change the identity of the dominant local error
sources. It does, however, make the closed-orbit problem six-dimensional and
avoids the four-variable ForwardDiff failure encountered by the default
coasting-beam closed-orbit calculation.

## Reproduction

```bash
julia --project=. test_codes/test_bmad_scibmad.jl \
  --rf-on \
  --reference=bmad_comparison/bmad_reference_rf_on_output.tar.gz \
  --csv=bmad_comparison/bmad_scibmad_rf_on_comparison.csv

julia --project=. bmad_comparison/compare_rf_on_optics.jl
```

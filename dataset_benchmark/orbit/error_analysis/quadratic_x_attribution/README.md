# X-quadratic 10%/90% attribution

This directory tests why the equal-corrector-RMS X-orbit quadratic response is
approximately 9.14% `Q_hh` and 90.83% `Q_vv` in squared detector-vector norm.

The first-stage diagnostic combines:

- the adopted second-order GTPSA detector `Q` values from
  [`../mixed_terms/gtpsa_results/`](../mixed_terms/gtpsa_results/);
- one all-corrector, first-order GTPSA response matrix evaluated at the 76
  active normal-sextupole elements; and
- the same 100 fixed horizontal/vertical corrector directions used by the
  mixed-term study.

For each direction, the normal-sextupole source-exposure proxy is

```text
E_h = sum_j |K2L_j| x_h(j)^2,
E_v = sum_j |K2L_j| y_v(j)^2.
```

The script also symmetrically rescales the H and V directions so that their
source exposures both equal `sqrt(E_h E_v)`. This is a causal diagnostic only;
the paper's physical result remains defined at equal corrector RMS.

Run from the `CESR Project` environment:

```powershell
julia --project=. dataset_benchmark/orbit/error_analysis/quadratic_x_attribution/run_internal_exposure_attribution.jl
python dataset_benchmark/orbit/error_analysis/quadratic_x_attribution/analyze_element_exposure.py
```

See [`results/RESULTS.md`](results/RESULTS.md) for the measured result and its
limitations. The element-local continuation is in
[`element_results/ELEMENT_EXPOSURE_RESULTS.md`](element_results/ELEMENT_EXPOSURE_RESULTS.md).

## Element-local continuation

The maintained run also writes the individual contributions

```text
e_h,j = |K2L_j| x_h,j^2,
e_v,j = |K2L_j| y_v,j^2
```

for every direction and active normal sextupole. The analyzer checks that the
76 element contributions reconstruct each direction's original `E_h` and
`E_v`, then ranks individual elements, East/West station pairs, ring sides,
and the positive-/negative-`K2` location classes. Since the proxy uses
`|K2L|`, the sign classes label alternating optics locations; they are not a
signed detector-response decomposition.

For independent paper layout, the analyzer writes two separate figures:

- [`element_results/top15_sextupole_exposure_differences.svg`](element_results/top15_sextupole_exposure_differences.svg)
  compares horizontal and vertical exposure for the 15 largest differences;
- [`element_results/sextupole_exposures_strengths_ring.svg`](element_results/sextupole_exposures_strengths_ring.svg)
  shows exposure excess and signed `K2L` around the full ring.

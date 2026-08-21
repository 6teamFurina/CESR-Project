# Signed sextupole contributions to the detector quadratic error

This experiment complements the unsigned normal-sextupole exposure analysis.
It asks which active normal sextupoles reinforce or cancel the final
99-detector horizontal second-order orbit-error vector.

For each fixed horizontal/vertical corrector direction pair, the calculation
combines the adopted GTPSA detector vectors `Q_hh,x`, `Q_hv,x`, and `Q_vv,x`
with the first-order internal orbit at all 76 active normal sextupoles and the
periodic detector response to local horizontal and vertical angular kicks.
The thick element is approximated symmetrically by applying half of each local
kick at the entrance and half at the exit, consistent with midpoint orbit
sampling to leading order.

For block `b`, sextupole `j`, direction `t`, and detector vector `C_b,j,t`, the
ensemble signed projection is

```text
eta_b,j = sum_t dot(C_b,j,t, Q_b,t) / sum_t norm(Q_b,t)^2.
```

The total target is `Q_x = Q_hh,x + Q_hv,x + Q_vv,x`. Positive projection
reinforces the target and negative projection cancels it. These percentages
are signed vector projections, not positive shares, and may be negative or
exceed 100% when large vectors cancel.

For `x = a*x_h + b*x_v` and `y = a*y_h + b*y_v`, the leading thin normal-
sextupole source is

```text
Q_hh: dpx = -K2L/2 * (x_h^2 - y_h^2),  dpy = K2L*x_h*y_h
Q_hv: dpx = -K2L   * (x_h*x_v - y_h*y_v),
      dpy =  K2L   * (x_h*y_v + x_v*y_h)
Q_vv: dpx = -K2L/2 * (x_v^2 - y_v^2),  dpy = K2L*x_v*y_v.
```

The script reports blockwise and total vector closure against the full GTPSA
target. Unless that closure is quantitatively small, the result is a leading
thin-kick normal-sextupole reconstruction rather than a complete element
attribution.

Run from the `CESR Project` root:

```powershell
julia --project=. dataset_benchmark/orbit/error_analysis/sextupole_detector_contributions/run_sextupole_detector_contributions.jl
python dataset_benchmark/orbit/error_analysis/sextupole_detector_contributions/analyze_sextupole_detector_contributions.py
```

Use `--trials=3 --output-dir=.../smoke` for a quick integration test.

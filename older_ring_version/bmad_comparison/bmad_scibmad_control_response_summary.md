# CESR Bmad–SciBmad control-response summary

SciBmad uses 119 first-order GTPSA parameters to differentiate the closed orbit at all 99 detector markers with respect to the Bmad-compatible HKICK/VKICK Overlay knobs.

| Mode | Status | Relative Frobenius | Max absolute difference (m/rad) | Correlation | Runtime (s) |
|---|---|---:|---:|---:|---:|
| rf_on | complete | 2.032285710e-03 | 6.747159290e-02 | 0.999997991471 | 44.444 |
| rf_off | complete | 2.012009559e-03 | 6.726326268e-02 | 0.999998031423 | 41.676 |

## Interpretation

All completed matrices agree with Bmad to at most `0.203229%` in relative Frobenius norm; the lowest full-matrix correlation is `0.999997991471`.

- Bmad `HKICK`/`VKICK` are laboratory-frame kicks. For element alignment tilt `t`, SciBmad uses `HKICK -> (Kn0L, Ks0L) = (-cos(t), -sin(t)) HKICK` and `VKICK -> (Kn0L, Ks0L) = (-sin(t), cos(t)) VKICK`.
- RF-on uses the 6D RF-confined closed orbit. RF-off uses a 4D coasting closed orbit with fixed `z=pz=0`.
- The RF-off baseline orbit uses the six-direction ForwardDiff adapter to extract the transverse `4 x 4` Jacobian, followed by SciBmad's native `BatchSolve.newton!`. All 119 control derivatives in both modes are computed simultaneously with first-order GTPSA parameters.

Detailed per-mode summaries and labeled matrices are stored in the corresponding `bmad_control_response_rf_on` and `bmad_control_response_rf_off` directories.

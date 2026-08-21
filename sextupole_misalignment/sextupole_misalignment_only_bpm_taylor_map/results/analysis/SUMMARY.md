# Only-sextupole-misalignment BPM/Taylor-map results

All forward states use the validated latest repaired SciBmad lattice.  The only
machine error is fixed x/y misalignment on all 76 active normal sextupoles.
BPM errors, time drift, corrector/K2 calibration errors, and quadrupole errors
are absent.

- targets / latent realizations: 76 / 1 per target
- exact states: 9500
- bump grid: 25 points, amplitude +/-0.500 mm per plane
- K2 grid: 5 points, delta K2 range -0.02 to 0.02 m^-3
- BPM count: 111

## BPM-predicted local orbit

| quantity | 2D RMSE [um] | median [um] | P90 [um] | maximum [um] |
|---|---:|---:|---:|---:|
| relative_local_orbit_nominal_k2_nonzero_bumps | 0.040180 | 0.007369 | 0.057320 | 0.326648 |
| relative_local_orbit_all_states | 0.039371 | 0.006678 | 0.056163 | 0.328788 |
| absolute_zero_bump_reference_orbit | 0.010091 | 0.002330 | 0.014513 | 0.043684 |

## Center inverse

| method | beam-relative 2D RMSE [um] | median [um] | P90 [um] | absolute-increment 2D RMSE [um] |
|---|---:|---:|---:|---:|
| fd_linear_source_predicted | 6.394706 | 5.000718 | 8.956349 | 6.394942 |
| fd_quartic_source_predicted | 6.394716 | 5.001335 | 8.955497 | 6.394951 |
| quadratic_o_derivative_predicted | 5.128185 | 3.738924 | 7.549159 | 5.127145 |
| chain_rule_o_derivative_predicted | 5.465936 | 3.981950 | 7.771729 | 5.465113 |
| o_taylor_order3_nominal_local | 33.008657 | 27.946649 | 47.192103 | 33.009199 |
| o_taylor_order4_nominal_local | 6.362654 | 4.656776 | 9.719398 | 6.362141 |
| o_taylor_order5_nominal_local | 5.724636 | 4.490103 | 8.573020 | 5.723937 |
| o_taylor_order4_all_state_local | 6.362657 | 4.656773 | 9.719401 | 6.362143 |
| fd_quartic_source_oracle_local | 6.394493 | 4.986177 | 8.944900 | 6.394715 |
| o_taylor_order4_oracle_all_state_local | 6.361597 | 4.659038 | 9.718034 | 6.361100 |

`fd_*_source_predicted` is the maintained physical two-source inverse after a
linear or quartic-in-K2 derivative extraction.  `quadratic_o_derivative` fits
the observable K2 derivative directly as a quadratic function of the
BPM-predicted target orbit.  `chain_rule_o_derivative` instead obtains the
local orbit Jacobian with respect to the commanded bump and transforms the
observable derivatives by the chain rule.  `o_taylor_orderN` fits all retained
raw nonzero-K2 states to a total-order-N polynomial in local x, y, and K2 and
finds the common zero of its analytic K2 derivative.

Oracle-local rows use exact target coordinates only after all machine-facing
fits and are evaluation diagnostics, not deployable methods.

# Only-sextupole-misalignment BPM/Taylor-map results

All forward states use the validated latest repaired SciBmad lattice.  The only
machine error is fixed x/y misalignment on all 76 active normal sextupoles.
BPM errors, time drift, corrector/K2 calibration errors, and quadrupole errors
are absent.

- targets / latent realizations: 1 / 1 per target
- exact states: 125
- bump grid: 25 points, amplitude +/-0.500 mm per plane
- K2 grid: 5 points, delta K2 range -0.02 to 0.02 m^-3
- BPM count: 111

## BPM-predicted local orbit

| quantity | 2D RMSE [um] | median [um] | P90 [um] | maximum [um] |
|---|---:|---:|---:|---:|
| relative_local_orbit_nominal_k2_nonzero_bumps | 0.045097 | 0.029354 | 0.069368 | 0.096139 |
| relative_local_orbit_all_states | 0.044187 | 0.025191 | 0.070196 | 0.097045 |
| absolute_zero_bump_reference_orbit | 0.006151 | 0.006151 | 0.006151 | 0.006151 |

## Center inverse

| method | beam-relative 2D RMSE [um] | median [um] | P90 [um] | absolute-increment 2D RMSE [um] |
|---|---:|---:|---:|---:|
| fd_linear_source_predicted | 3.651035 | 3.651035 | 3.651035 | 3.647528 |
| fd_quartic_source_predicted | 3.651328 | 3.651328 | 3.651328 | 3.647821 |
| quadratic_o_derivative_predicted | 2.451490 | 2.451490 | 2.451490 | 2.456490 |
| chain_rule_o_derivative_predicted | 2.473944 | 2.473944 | 2.473944 | 2.478730 |
| o_taylor_order3_nominal_local | 9.244054 | 9.244054 | 9.244054 | 9.240876 |
| o_taylor_order4_nominal_local | 2.962466 | 2.962466 | 2.962466 | 2.964227 |
| o_taylor_order5_nominal_local | 2.922802 | 2.922802 | 2.922802 | 2.928902 |
| o_taylor_order4_all_state_local | 2.962485 | 2.962485 | 2.962485 | 2.964246 |
| fd_quartic_source_oracle_local | 3.670629 | 3.670629 | 3.670629 | 3.667110 |
| o_taylor_order4_oracle_all_state_local | 2.964625 | 2.964625 | 2.964625 | 2.966395 |

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

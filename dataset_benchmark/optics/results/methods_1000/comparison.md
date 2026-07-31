# CESR chromatic-optics method comparison

SciBmad pointwise `twiss` is the numerical reference. Detector phase origins are removed per sample. Bmad longitudinal `dphi_3/ddelta` and ring `slip_factor` are sign-aligned using the convention documented in the optics README. Constant columns are excluded from correlations.

| method | samples | physics_seconds | samples_per_second | speedup_vs_bmad | speedup_vs_scibmad_pointwise | result_scope | maximum_closure_residual | closure_scope | minimum_column_correlation | median_column_correlation | maximum_absolute_error | maximum_normalized_column_error |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Bmad/Tao | 1000 | 701.528 | 1.42546 | 1 | 0.813611 | exact per sample | 1.98802e-13 | all samples | 0.095562 | 0.999995 | 44.095 | 450.362 |
| SciBmad pointwise `twiss` | 1000 | 570.771 | 1.75202 | 1.22909 | 1 | exact per sample | 2.43508e-12 | all samples | 1 | 1 | 0 | 0 |
| SciBmad prototype `twiss!` | 1000 | 513.273 | 1.94828 | 1.36677 | 1.11202 | exact per sample | 2.43508e-12 | all samples | 1 | 1 | 0 | 0 |
| SciBmad one parameterized `twiss` | 1000 | 73.1303 | 13.6742 | 9.59284 | 7.80485 | first-order control surrogate | 1.28866e-13 | nominal orbit only | 0.997779 | 0.999996 | 45.7543 | 0.104608 |

## Interpretation

- The reusable `twiss!` prototype is numerically identical to the pointwise reference for every compared field; its maximum absolute difference is zero.
- The parameterized method's lowest column correlation is 0.997779 for `xi_2`. Its closure residual is for the nominal orbit only because this method is a local corrector surrogate.
- Bmad's isolated minimum correlation is 0.095562 for `dorbit_z_ddelta`; the median over nonconstant columns is 0.999995. Large normalized errors in columns whose reference maximum is nearly zero should not be interpreted as a global optics error.

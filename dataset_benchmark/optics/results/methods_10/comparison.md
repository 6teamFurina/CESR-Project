# CESR chromatic-optics method comparison

SciBmad pointwise `twiss` is the numerical reference. Detector phase origins are removed per sample. Bmad longitudinal `dphi_3/ddelta` and ring `slip_factor` are sign-aligned using the convention documented in the optics README. Constant columns are excluded from correlations.

| method | samples | physics_seconds | samples_per_second | maximum_closure_residual | minimum_column_correlation | median_column_correlation | maximum_absolute_error | maximum_normalized_column_error |
|---|---|---|---|---|---|---|---|---|
| Bmad/Tao | 1000 | 701.528 | 1.42546 | 1.98802e-13 | 0.0893954 | 0.999996 | 38.4839 | 482.835 |
| SciBmad pointwise `twiss` | 10 | 6.18455 | 1.61693 | 4.55275e-13 | 1 | 1 | 0 | 0 |
| SciBmad prototype `twiss!` | 10 | 5.81505 | 1.71968 | 4.55275e-13 | 1 | 1 | 0 | 0 |
| SciBmad one parameterized `twiss` | 10 | 14.8803 | 0.672031 | 1.28866e-13 | -6.70976e-15 | 0.999965 | 283.901 | 0.204957 |

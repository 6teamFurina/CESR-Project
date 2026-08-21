# CESR chromatic-optics method comparison

SciBmad pointwise `twiss` is the numerical reference. Detector phase origins are removed per sample. Bmad longitudinal `dphi_3/ddelta` and ring `slip_factor` are sign-aligned using the convention documented in the optics README. Constant columns are excluded from correlations.

| method | samples | physics_seconds | samples_per_second | maximum_closure_residual | minimum_column_correlation | median_column_correlation | maximum_absolute_error | maximum_normalized_column_error |
|---|---|---|---|---|---|---|---|---|
| Bmad/Tao | 1000 | 701.528 | 1.42546 | 1.98802e-13 | 0.095562 | 0.999995 | 44.095 | 450.362 |
| SciBmad pointwise `twiss` | 1000 | 570.771 | 1.75202 | 2.43508e-12 | 1 | 1 | 0 | 0 |
| SciBmad prototype `twiss!` | 1000 | 570.771 | 1.75202 | 2.43508e-12 | 1 | 1 | 0 | 0 |
| SciBmad one parameterized `twiss` | 1000 | 73.1303 | 13.6742 | 1.28866e-13 | 0.997779 | 0.999996 | 45.7543 | 0.104608 |

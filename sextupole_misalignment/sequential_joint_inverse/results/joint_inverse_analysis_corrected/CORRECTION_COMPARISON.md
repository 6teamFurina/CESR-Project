# Fixed-baseline orbit-correction protocol comparison

All rows use the same deterministic 16-machine SciBmad ensemble, 10/3/3
machine split, measurement seeds, noise/drift augmentation, and inverse
definitions.  The only physical protocol change between the two quadrupole-
offset columns is whether the BPM-reference baseline correction is applied and
held fixed during every sextupole scan.

| inverse | zero-offset RMSE [um] | uncorrected offset RMSE [um] | corrected offset RMSE [um] | corrected reduction | corrected excess vs zero |
|---|---:|---:|---:|---:|---:|
| physics_gls | 33.923 | 100.970 | 34.181 | 66.15% | +0.76% |
| shared_target_local_ridge | 33.078 | 99.796 | 33.477 | 66.45% | +1.21% |
| shared_joint_ridge | 33.094 | 99.618 | 33.458 | 66.41% | +1.10% |
| shared_joint_random_feature | 33.092 | 99.119 | 33.444 | 66.26% | +1.06% |

The best learned uncorrected result is
`shared_joint_random_feature` at `99.119 um`.  After fixed
baseline correction, the best learned result is
`shared_joint_random_feature` at `33.444 um`, a
`66.26%`
reduction.  It is only
`1.10%` above the
best zero-offset result and removes `99.45%` of the excess RMSE
attributed to the uncorrected quadrupole-drift protocol.

The zero-offset-trained joint ridge gives `765.859 um` when
evaluated on uncorrected drift but `34.180 um` after the baseline
correction.  Thus correction removes most of the operational distribution
shift before the inverse is asked to estimate sextupole centers.

The corrected result still fails the strict tail gate: its best learned P99 is
`84.486 um` and worst-target RMSE is
`59.585 um`, both above 50
micrometers.  The result establishes the workflow benefit, not final CESR
precision or hardware safety.

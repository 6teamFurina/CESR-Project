# GTPSA-ORM plus noisy-reference correction comparison

All four columns use the same deterministic 16-machine latest-lattice SciBmad
ensemble, 10/3/3 machine split, static nuisance draws, sextupole scans, and
inverse definitions. The new protocol uses the zero-quadrupole-offset
SciBmad/GTPSA closed-orbit Jacobian, with fixed BPM and corrector gains, as its
103-control ORM. The stored reference and each current-orbit correction input
receive independent Gaussian BPM-noise means. The solved baseline command is
then frozen during all 76 sextupole scans.

| inverse | zero offset [um] | uncorrected offset [um] | FD ORM, noiseless correction [um] | GTPSA ORM, noisy correction [um] | GTPSA/noisy change vs FD/noiseless |
|---|---:|---:|---:|---:|---:|
| physics_gls | 33.923 | 100.970 | 34.181 | 34.181 | -0.001% |
| shared_target_local_ridge | 33.078 | 99.796 | 33.477 | 33.477 | +0.000% |
| shared_joint_ridge | 33.094 | 99.618 | 33.458 | 33.458 | -0.001% |
| shared_joint_random_feature | 33.092 | 99.119 | 33.444 | 33.416 | -0.082% |

The correction inputs use `5.0 um` RMS
noise per BPM plane and read, averaged over
`3,072` reads. The expected noise of
each mean is only `0.090 um`; the
realized reference-noise RMS is
`0.091 um`.

The GTPSA ORM agrees with the central finite-difference check to maximum
relative L2 difference
`1.879e-08` across the 16
machines; its maximum periodic-response closure norm is
`3.553e-15`. Relative to the finite-
difference/noiseless baseline, the new commands differ by only
`0.0107 urad` RMS,
the corrected BPM orbit by
`0.0767 um` RMS, and the
corrected target orbit by
`0.1248 um` 2D RMS.

The best learned GTPSA/noisy result is
`shared_joint_random_feature` at `33.416 um` 2D RMSE. It is
`66.287%` below the best
uncorrected result, `1.021%`
above the best zero-offset result, and
`-0.082%` relative to
the best finite-difference/noiseless corrected result. The last difference is
numerically negligible and is not evidence that adding noise improves the
inverse.

The strict tail gate still fails: the best GTPSA/noisy P99 is
`84.382 um` and its worst-target RMSE is
`59.714 um`. This experiment
therefore validates compatibility of the requested workflow at the tested
3,072-read averaging level. It does not establish robustness to single-shot
or low-repeat BPM noise, correlated BPM noise, missing channels, response
measurement error, or hardware effects. A repeat-count/covariance sweep is
required before selecting a CESR correction acquisition protocol.

The ORM is scaled with the realized simulated BPM and baseline-corrector gains,
so this is an exact-calibration/model-conditioned response test rather than an
unknown gain-mismatch test. The 103 baseline controls and 62 local-bump
controls also retain separate deterministic gain registries. That convention
is shared by the two corrected protocols and preserves their paired numerical
comparison, but it must be unified at the physical-device level before a
facility-facing calibration conclusion.

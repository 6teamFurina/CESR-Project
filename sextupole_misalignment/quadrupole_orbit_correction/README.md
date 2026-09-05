# Paired quadrupole-offset orbit correction

This study tests whether machine-available correctors can restore the BPM
closed orbit after a fixed quadrupole alignment drift while every other
maintained static random error remains unchanged.

Each latent machine is evaluated as a strict pair on the repaired latest
SciBmad lattice:

1. quadrupole x/y alignment drift disabled, producing the measured reference
   BPM orbit;
2. the same sextupole offsets, corrector gains, K2 gains, BPM gains,
   quadrupole strength errors, and quadrupole rolls, with independent
   50-micrometer RMS x/y quadrupole offsets enabled;
3. corrector-only orbit restoration to the paired reference BPM vector.

The correction solver receives only reference/current BPM readbacks and a
corrector-to-BPM response matrix. It never consumes the latent quadrupole
offsets or target-sextupole orbit. The maintained scan generator now defaults
to one theoretical response calculated from the nominal latest-lattice
SciBmad/GTPSA model and reused for every latent machine. A realized-model
GTPSA response and the historical measured finite-difference responses remain
explicit comparison options.

The correction basis is the dynamically discovered set of 103 horizontal and
vertical steering Overlay controls in the repaired latest lattice. A coupled
horizontal/vertical SVD-ridge inverse proposes corrector commands. Exact
scalar closed-orbit solves and a backtracking line search verify every accepted
update. Central finite differences are no longer a calculation default; they
must be requested explicitly and remain useful as an independent validation.

The default ridge scale is `1e-2` times the largest response singular value.
It deliberately suppresses weak response modes; the setting is a bounded
simulation choice, not a substitute for CESR corrector limits or an
operations-selected regularization rule.

Run from `CESR Project/`:

```powershell
julia --project=. `
  sextupole_misalignment/quadrupole_orbit_correction/run_quadrupole_orbit_correction.jl

python `
  sextupole_misalignment/quadrupole_orbit_correction/validate_orbit_correction.py
```

The generated scientific summary and arrays are written below
`results/orbit_correction_50um/`.

## Production result (2026-08-28)

The deterministic production benchmark contains 16 paired machines, 113
physical quadrupoles, 111 measurable BPMs, 76 target sextupoles, and 103
normal horizontal/vertical steering controls.  With the response matrix saved
in the zero-quadrupole-offset reference state, correction reduces the
aggregate measured BPM-coordinate RMS difference from 860.882 to 126.034
micrometers and the target-sextupole 2D orbit difference from 1,335.740 to
70.709 micrometers.  The command RMS is 20.364 microradians and the largest
absolute command is 137.143 microradians.

Remeasuring the response after applying the quadrupole offsets gives nearly
the same result: 126.003 micrometers at the BPMs and 70.724 micrometers at the
target sextupoles.  Across machines, the current-versus-reference response
relative-L2 difference has 1.996% median and 3.275% maximum; the two final BPM
vectors differ by 1.307 micrometers RMS.  The fraction of beam-relative
sextupole centers outside the 1.5-millimeter excitation radius changes from
zero in the reference, to 28.618% before correction, and back to zero after
either correction.

This is restoration to the stored nonzero reference orbit, not correction to
BPM zero.  The 126-micrometer aggregate residual also means the finite BPM
vector is not reproduced exactly.

## Fixed-baseline scan and inverse result (2026-08-28)

`generate_corrected_joint_machine_scans.jl` now applies each machine's
stored-reference baseline command and holds it fixed while all 76 sextupoles
are scanned one at a time.  The 62 local bump controls are additive on top of
the 103-control baseline.  The resulting tensor contains 36,480 exact RF-on
SciBmad state lanes.  Its zero-bump/zero-`delta K2` states reproduce the saved
corrected reference with maximum BPM and target errors below `5e-14 m`, and
the saved commands and corrected orbits exactly match the standalone orbit-
correction artifact.

With the same 10/3/3 machine split, noise/drift augmentation, and inverse
definitions, the corrected fixed-physics, target-local ridge, joint ridge, and
joint random-feature RMSE values are 34.181, 33.477, 33.458, and 33.444
micrometers.  The best corrected value is 66.259% below the best uncorrected-
offset value of 99.119 micrometers and only 1.104% above the best zero-offset
value of 33.078 micrometers.  Correction removes 99.447% of the excess RMSE
associated with the uncorrected protocol.

The zero-offset-trained joint ridge improves from 765.859 micrometers when
evaluated on uncorrected drift to 34.180 micrometers after correction.  This
supports correction as a physical preprocessing step rather than asking the
inverse to absorb the annual orbit shift.  The strict tail gate still fails:
the best corrected P99 is 84.486 micrometers and worst-target RMSE is 59.585
micrometers.  All-target context changes ridge RMSE by only -0.058% relative
to the target-local model.

Run the corrected workflow from `CESR Project/`:

```powershell
julia --project=. `
  sextupole_misalignment/quadrupole_orbit_correction/generate_corrected_joint_machine_scans.jl `
  --baseline-response-method=finite_difference `
  --gtpsa-response-model=realized `
  --correction-bpm-noise-rms-m=0.0 `
  --correction-measurement-repeats=1 `
  --corrected-case-name=with_quadrupole_misalignment_corrected

python `
  sextupole_misalignment/sequential_joint_inverse/analyze_joint_inverse.py `
  --comparison-case with_quadrupole_misalignment_corrected `
  --output-dir sextupole_misalignment/sequential_joint_inverse/results/joint_inverse_analysis_corrected

python `
  sextupole_misalignment/sequential_joint_inverse/validate_joint_inverse.py `
  --comparison-case with_quadrupole_misalignment_corrected `
  --analysis-dir sextupole_misalignment/sequential_joint_inverse/results/joint_inverse_analysis_corrected `
  --write-report

python `
  sextupole_misalignment/quadrupole_orbit_correction/compare_corrected_inverse.py
```

The main corrected summary is
`../sequential_joint_inverse/results/joint_inverse_analysis_corrected/SUMMARY.md`;
the matched zero/uncorrected/corrected comparison is
`../sequential_joint_inverse/results/joint_inverse_analysis_corrected/CORRECTION_COMPARISON.md`.

This is a synthetic, noise-free correction benchmark. Corrector commands are
the latest-lattice `HKICK`/`VKICK` control variables in radians; no CESR
hardware limit or operator safety claim is made. Finite BPM matching also does
not guarantee exact orbit equality between BPMs, so the orbit residual at all
76 target sextupoles is reported independently.

## GTPSA ORM plus noisy correction result (2026-08-28)

`generate_gtpsa_noisy_corrected_joint_machine_scans.jl` repeats the complete
16-machine, 76-target workflow with the stored-reference ORM calculated from
the first-order SciBmad/GTPSA implicit periodic closed-orbit Jacobian.  The
static BPM and corrector gain errors are applied to the Jacobian.  Independent
Gaussian noise means are added to the stored reference, every current
correction measurement, and a held-out validation readback.  All other error
draws remain identical to the paired production ensemble, and the accepted
baseline command remains fixed throughout the sextupole scans.

The default correction acquisition uses the same 5-micrometer RMS BPM noise
per plane/read and 3,072 repeats as the downstream sextupole study.  Its mean-
noise standard deviation is 0.090 micrometers, and the realized reference-
noise RMS is 0.091 micrometers.  Across all 16 machines, the GTPSA ORM agrees
with the central finite-difference check to maximum relative L2 difference
`1.879e-8`; the maximum periodic-response closure norm is `3.553e-15`.

Compared with the finite-difference/noiseless correction, the resulting
commands differ by only 0.0107 microradians RMS, the corrected BPM orbit by
0.0767 micrometers RMS, and the corrected target orbit by 0.1248 micrometers
2D RMS.  The downstream best held-out inverse changes from 33.444 to 33.416
micrometers, a numerically negligible -0.082%.  It remains 66.287% below the
uncorrected 99.119-micrometer result and 1.021% above the 33.078-micrometer
zero-offset benchmark.  P99 is 84.382 micrometers and worst-target RMSE is
59.714 micrometers, so the strict tail gate still fails.

Run the matched extension from `CESR Project/`:

```powershell
julia --project=. `
  sextupole_misalignment/quadrupole_orbit_correction/generate_gtpsa_noisy_corrected_joint_machine_scans.jl

python `
  sextupole_misalignment/sequential_joint_inverse/analyze_joint_inverse.py `
  --comparison-case with_quadrupole_misalignment_gtpsa_noisy_corrected `
  --output-dir sextupole_misalignment/sequential_joint_inverse/results/joint_inverse_analysis_gtpsa_noisy_corrected

python `
  sextupole_misalignment/sequential_joint_inverse/validate_joint_inverse.py `
  --comparison-case with_quadrupole_misalignment_gtpsa_noisy_corrected `
  --analysis-dir sextupole_misalignment/sequential_joint_inverse/results/joint_inverse_analysis_gtpsa_noisy_corrected `
  --write-report

python `
  sextupole_misalignment/quadrupole_orbit_correction/compare_gtpsa_noisy_inverse.py
```

The four-protocol result is in
`../sequential_joint_inverse/results/joint_inverse_analysis_gtpsa_noisy_corrected/GTPSA_NOISY_COMPARISON.md`.
Because 3,072 reads suppress the assumed white noise by `sqrt(3072)`, this
test does not demonstrate robustness for single-shot or low-repeat orbit
correction.  Correlated BPM errors, empirical ORM uncertainty, outliers,
missing channels, settling, actuator hysteresis, and hardware constraints also
remain outside the present result.

The response rows and columns are scaled with the realized simulated BPM and
baseline-corrector gains.  This is therefore an exact-calibration/model-
conditioned response test, not a test of an unknown gain mismatch between a
nominal GTPSA model and the machine.  Also, the 103 baseline controls and the
62 local-bump controls currently use separate deterministic 1%-RMS gain draws;
the few channels that map to the same physical kicker are not yet assigned one
shared device-level gain.  Both corrected protocols use the same convention,
so their matched comparison remains valid, but a facility-facing study should
unify that registry before drawing calibration conclusions.

## Nominal GTPSA ORM and threaded full-error scan (2026-08-30)

`generate_corrected_joint_machine_scans.jl` now defaults to this stricter
unknown-calibration experiment; the named
`generate_gtpsa_nominal_corrected_joint_machine_scans.jl` wrapper fixes the
complete maintained acquisition protocol.  A single theoretical 222-by-103
ORM is built from the nominal latest-lattice SciBmad/GTPSA periodic
closed-orbit Jacobian and reused for all 16 latent machines.  Production
explicitly disables the central-finite-difference ORM check, and the response
is not scaled by realized BPM gains, corrector gains, quadrupole errors, or
alignment errors.  Those realizations remain embedded only in the simulated
BPM observations.

With 5-micrometer RMS BPM noise per read averaged over 3,072 reads, the nominal
response correction reduces aggregate BPM-coordinate RMS from 860.882 to
126.037 micrometers and target-sextupole 2D orbit RMS from 1,335.740 to 70.945
micrometers.  This is the correction actually held fixed during the subsequent
one-at-a-time 76-sextupole scan.

The scan is distributed across six Julia threads after each fixed machine
realization is prepared.  Every worker owns an independently loaded mutable
ring/model, so no lattice object is shared across targets.  The 16-by-76
production generation took 312.5 seconds.  Rescanning the first production
target serially gave maximum absolute differences of exactly zero for BPM,
drift-BPM, target, and drift-target tensors.

Run the source generation from `CESR Project/`:

```powershell
julia --threads=auto --project=. `
  sextupole_misalignment/quadrupole_orbit_correction/generate_gtpsa_nominal_corrected_joint_machine_scans.jl
```

The source case is
`../sequential_joint_inverse/results/exact_joint_machines/with_all_errors_gtpsa_nominal_corrected`.
Its downstream BPM/GTPSA state-space inverse and validation are documented in
`../finite_bpm_inversion/README.md`.  This remains a synthetic result: the
maintained error inventory does not yet include bad/missing BPM channels,
correlated empirical ORM errors, settling, hysteresis, or hardware limits.

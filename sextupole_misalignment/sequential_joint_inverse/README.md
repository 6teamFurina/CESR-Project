# Sequential excitation and joint sextupole-center inversion

This study tests the first multi-center route that preserves physical
identifiability: keep one latent machine fixed, excite the 76 active normal
sextupoles one at a time, retain the excitation-target axis, and evaluate all
76 two-plane beam-relative magnetic centers together.

New calculations use the repaired latest SciBmad lattice at
`Latest_Lattice/latest_cesr_scibmad_repaired.jl` with RF on.  The study does
not use Bmad/Tao results as substitutes for the primary calculation.

## Paired machine ensemble

`generate_joint_machine_scans.jl` makes the latent machine realization, rather
than one target scan, the atomic dataset unit.  Within each realization the
following quantities remain fixed while all 76 targets are scanned:

- independent 300-micrometer RMS x/y offsets on all sextupoles;
- 1% RMS corrector gain errors;
- 1% RMS target-K2 intervention gain errors;
- independent quadrupole strength errors uniformly bounded by +/-1%;
- 1 mrad RMS quadrupole rolls;
- 1% RMS BPM x/y gains;
- the target-specific direction used for the drift-response secant.

The two generated cases share every one of those draws.  The second case also
adds independent Gaussian x/y displacement to every physical quadrupole with
50 micrometers RMS per plane.  A physical quadrupole's draw is applied
coherently to all of its tracking slices and remains fixed during every target
scan.  This is the maintained interpretation of residual annual alignment
drift relative to the latest measured nominal geometry.  It is not survey
uncertainty resampled during acquisition and it is not a 50-micrometer radial
hard bound.

The facility statement does not by itself determine whether 50 micrometers is
a per-plane RMS, a two-dimensional radial RMS, or a correlated girder motion.
The primary run uses the conventional per-plane independent-error
interpretation above.  A 50-micrometer radial RMS would correspond to about
35.4 micrometers per plane under an isotropic Gaussian model and should be a
separate sensitivity run rather than silently relabeled as this result.

No orbit correction is applied in this first paired experiment.  The saved
reference BPM and target orbits therefore expose whether the 50-micrometer
quadrupole drift first moves the closed orbit outside the useful excitation
domain.  If that effect dominates, the appropriate extension is a recorded
BPM-only orbit-correction step relative to the yearly nominal orbit; the
learned inverse must not be credited with replacing that machine operation.

## Scan tensor

For every latent machine and target, exact SciBmad states cover the axial bump
cross

```text
(-x,0), (0,-y), (0,0), (0,+y), (+x,0)
```

at `x = y = 1.5 mm` and `delta K2 = -0.10, 0, +0.10 m^-3`.  The primary BPM
tensor has axes

```text
machine, excitation_target, bump, delta_K2, BPM, output_plane
```

and shape `M x 76 x 5 x 3 x 111 x 2`.  An equally shaped paired drift tensor
uses the same target bump knobs and a `-5` to `+5 micrometer` state-order
secant.  Raw scan tensors, target-entry orbits, zero-excitation reference
orbits, all latent draws, target/BPM inventories, commands, and provenance are
preserved.

The generator evaluates each target's 15 physical and 15 drift-secant states
in one 30-lane SciBmad `BatchParam` closed-orbit solve.  Before parameterizing
the selected elements, it materializes their current scalar multipole values.
This avoids promoting zero-valued latest-lattice `DefExpr` components to
`BatchParam` while preserving the physical nuisance values already applied.

## Inverses

`analyze_joint_inverse.py` derives the K2-odd, bump-odd full-BPM parity tensor
and compares four estimators:

1. `physics_gls`: fixed nominal full-BPM source templates and a uniform-
   covariance GLS inverse for every target block;
2. `shared_target_local_ridge`: one target-sharing residual model using the
   selected target scan, target identity, and baseline-orbit modes;
3. `shared_joint_ridge`: the same estimator plus compact context derived from
   all 76 target scans in the latent machine;
4. `shared_joint_random_feature`: the joint representation with a fixed tanh
   hidden layer and validation-selected ridge output.

The latter two consume the complete target-indexed machine tensor and return
all 152 coordinates in one evaluation.  They never sum target scans before
inversion.  The response-residual basis, reference-orbit basis, joint-context
basis, regularization, and nonlinear feature scale are fitted using training
machines only.

Machines, not individual targets or noisy augmentations, are split into
training, validation, and test sets.  The default 16-machine run uses 10/3/3
machines.  Held-out measurement augmentation uses 5-micrometer RMS independent
BPM noise per read, 3,072 repeated balanced eight-state cycles, and a
10-micrometer endpoint-RMS random walk propagated through the exact drift
secant.  Exact target-entry orbit and latent centers are evaluation-only.

The preferred reference gate is aggregate 2D RMSE below 30 micrometers.  The
strict gate requires aggregate RMSE, P99, and every target-level RMSE to remain
below 50 micrometers.  Joint context is counted as helpful only when its
held-out RMSE is lower than the matched shared-local model; model complexity is
not itself evidence of improvement.

## Production pilot result (2026-08-23)

The completed paired run contains 16 full latent machines, split 10/3/3 by
machine for training, validation, and held-out testing.  Without quadrupole
alignment drift, the best learned estimator is the shared target-local ridge:
its held-out 2D RMSE is 33.078 micrometers, P99 is 84.159 micrometers, and the
worst target RMSE is 60.873 micrometers.  The all-target joint ridge gives
33.094 micrometers, a 0.046% degradation relative to the matched local model;
the joint random-feature result is 33.139 micrometers.  The aggregate result is
promising but the preferred 30-micrometer and strict 50-micrometer gates both
fail.

With uncorrected 50-micrometer/plane quadrupole alignment drift, the fixed
physics, local ridge, joint ridge, and joint random-feature RMSE values are
100.970, 99.796, 99.618, and 99.119 micrometers, respectively.  The joint
ridge improves on the local model by only 0.178%; the best nonlinear comparator
still has 371.918-micrometer P99 and 354.570-micrometer worst-target RMSE.
Training the joint ridge without quadrupole drift and evaluating it with the
drift enabled gives 765.859 micrometers RMSE, an explicit failure under this
distribution shift.

The aligned case is also a different excitation-domain problem before orbit
correction: its paired reference BPM orbit changes by 860.882 micrometers RMS,
the beam-relative center truth changes by 1,335.740 micrometers 2D RMS, and
28.618% of machine-target truths lie outside the 1.5-millimeter bump radius,
versus none in the no-alignment case.  Therefore this run supports preserving
one-at-a-time excitation and the target axis, but it does not establish a
material advantage from all-target model context.  The next quadrupole-drift
experiment must add and record BPM-only orbit correction before comparing
inverse architectures.  The complete tables and plots are in
`results/joint_inverse_analysis/SUMMARY.md`.

## Paired BPM-reference correction follow-up (2026-08-28)

`../quadrupole_orbit_correction/` now implements the required pre-scan
correction on the same deterministic 16-machine ensemble.  The correction
target is the closed-orbit BPM readback saved while the quadrupoles are at
their calibrated nominal alignment positions; it is neither BPM electrical
zero nor an orbit that has first been corrected to zero.  All other maintained
random errors remain present and unchanged when the 50-micrometer/plane
quadrupole offsets are switched on.

Using the zero-offset state's stored response matrix and the 103 normal
horizontal/vertical steering Overlay controls reduces the aggregate measured
BPM-coordinate RMS difference from 860.882 to 126.034 micrometers and the
target-sextupole 2D orbit difference from 1,335.740 to 70.709 micrometers.  A
response matrix remeasured after the offsets gives 126.003 and 70.724
micrometers, respectively.  The two corrected BPM vectors differ by only
1.307 micrometers RMS.  The current-versus-stored response-matrix relative-L2
difference has 1.996% median and 3.275% maximum across machines, supporting
the stored reference response as the bounded default for this fixed-error
ensemble.

The fraction of beam-relative sextupole centers outside the maintained
1.5-millimeter scan radius falls from 28.618% before correction to zero with
either response matrix.  This restores the intended excitation range but does
not make every BPM reading exactly identical: the remaining aggregate BPM RMS
is about 126 micrometers, and the largest individual residual is about 1.38
millimeters.  The orbit-restoration benchmark is noise-free and asserts no
CESR hardware limit.

The full corrected follow-up is now complete.  For every machine, the
zero-offset reference ORM is measured once, the resulting 103-control baseline
command is applied after the quadrupole offsets, and that baseline remains
fixed while all 76 targets are scanned.  Local bump commands are additive.
The corrected scan tensor contains 36,480 exact SciBmad state lanes and its
zero-bump/zero-`delta K2` states reproduce the independently saved corrected
reference with maximum BPM and target discrepancies below `5e-14 m`.

The corrected fixed-physics, target-local ridge, joint ridge, and joint
random-feature held-out RMSE values are 34.181, 33.477, 33.458, and 33.444
micrometers.  The best result is therefore 66.259% lower than the best
uncorrected-offset result of 99.119 micrometers and only 1.104% above the best
zero-offset result of 33.078 micrometers; 99.447% of the excess RMSE associated
with the uncorrected protocol is removed.  A joint ridge trained only on the
zero-offset case changes from 765.859 micrometers on uncorrected drift to
34.180 micrometers after correction, showing that the operational distribution
shift is largely removed before inversion.

This does not close the precision problem.  The best corrected P99 is 84.486
micrometers and the worst-target RMSE is 59.585 micrometers, so the preferred
30-micrometer aggregate gate and the strict 50-micrometer tail gate still
fail.  Joint context changes ridge RMSE by only -0.058% relative to the local
model, which is not a material all-target advantage.  Complete corrected
results are in `results/joint_inverse_analysis_corrected/SUMMARY.md` and the
matched three-protocol table is in
`results/joint_inverse_analysis_corrected/CORRECTION_COMPARISON.md`.

## GTPSA ORM and noisy-reference follow-up (2026-08-28)

The second corrected production run keeps the complete static-error ensemble
and scan protocol unchanged, but calculates the 103-control stored-reference
ORM from one first-order SciBmad/GTPSA periodic closed-orbit Jacobian per
machine.  Fixed BPM gains scale its rows and fixed corrector gains scale its
columns.  Independent Gaussian noise is added to the stored reference and to
each current-orbit correction readback before the SVD-ridge command is solved;
latent quadrupole offsets and target-local orbit remain unavailable to the
solver.  The resulting command is again fixed during all 76 sextupole scans.

The correction noise matches the maintained scan measurement model:
5 micrometers RMS per BPM plane/read averaged over 3,072 reads.  Thus each
correction mean has only 0.090 micrometers expected standard deviation, and
the realized stored-reference noise is 0.091 micrometers RMS.  The GTPSA ORM
agrees with a central finite-difference validation to at most `1.879e-8`
relative L2 difference across the 16 machines, with maximum periodic-response
closure `3.553e-15`.

Relative to the finite-difference/noiseless correction, the GTPSA/noisy
baseline commands differ by 0.0107 microradians RMS, the corrected BPM orbit
by 0.0767 micrometers RMS, and the corrected target orbit by 0.1248
micrometers 2D RMS.  The fixed-physics, target-local ridge, joint ridge, and
joint random-feature held-out RMSE values are 34.181, 33.477, 33.458, and
33.416 micrometers.  The best value is 66.287% below the uncorrected result,
1.021% above the zero-offset result, and 0.082% below the matched
finite-difference/noiseless value.  That last tiny signed change is numerical
equivalence, not evidence that added noise improves estimation.

The P99 and worst-target RMSE remain 84.382 and 59.714 micrometers, so the
strict tail gate still fails.  This run validates the requested GTPSA-ORM plus
noisy-orbit workflow only at the tested 3,072-read averaging level.  A repeat-
count sweep with measured BPM covariance is still required before choosing a
CESR acquisition protocol.  Full results and the four-protocol comparison are
in `results/joint_inverse_analysis_gtpsa_noisy_corrected/SUMMARY.md` and
`results/joint_inverse_analysis_gtpsa_noisy_corrected/GTPSA_NOISY_COMPARISON.md`.

This is an exact-calibration/model-conditioned GTPSA response: its rows and
columns use the realized simulated BPM and baseline-corrector gains.  It does
not yet test an unknown calibration mismatch.  The 103 baseline-control gains
and 62 local-bump-control gains are also separate deterministic draws from the
same 1%-RMS prior rather than one unified physical-device registry.  This
convention is identical in the two corrected runs and therefore does not spoil
their paired numerical comparison, but it remains a machine-transfer limit.

## Explicit BPM/GTPSA local-orbit follow-up (2026-08-30)

`../finite_bpm_inversion/analyze_sequential_bpm_gtpsa_inverse.py` now applies
the nearest-upstream/downstream order-one GTPSA transport to this exact
corrected scan tensor before fitting the sextupole centers.  The machine-facing
phase receives only BPM observable readbacks, commands, and the nominal model;
it persists all relative and absolute local-orbit and center estimates before
the exact target orbit or latent offsets are loaded for evaluation.

On deterministic static readbacks, relative local-orbit RMSE is 14.306
micrometers.  The resulting beam-relative center RMSE is 21.658 micrometers,
compared with 20.912 micrometers for the evaluation-only exact-local oracle.
Adding the independently reconstructed absolute reference orbit gives 23.346
micrometers absolute sextupole-offset RMSE.  This closes the earlier data-path
gap: the maintained full-error corrected result no longer needs exact target
orbit as an inverse input.

The stochastic sensitivity remains unresolved.  With 32 held-out measurement
realizations, 5-micrometer per-read BPM white noise averaged over 3,072 reads,
and a 10-micrometer endpoint random walk over the unfiltered repeated 15-state
schedule, local-orbit RMSE remains 13.851 micrometers but center RMSE rises to
60.667 micrometers and P99 to 165.486 micrometers.  This isolates the next task
as covariance-aware K2-slope/drift estimation; it is not evidence that the
BPM/GTPSA local transport failed.

## Nominal-ORM state-space follow-up (2026-08-30)

The maintained full-error follow-up is now
`with_all_errors_gtpsa_nominal_corrected`.  Unlike the preceding
exact-calibration comparison, its correction uses one nominal theoretical
SciBmad/GTPSA ORM with no realized gain/error scaling and no production
finite-difference ORM.  The same fixed sextupole, BPM/corrector/K2,
quadrupole-strength/roll/alignment errors remain active across every one-at-a-
time target scan.  Six independent latest-lattice worker models generated all
16 by 76 targets in 312.5 seconds; the recorded thread-versus-serial maximum
difference is zero.

The downstream acquisition replaces the unfiltered 15-state baseline with
3,072 balanced eight-state cycles, periodic same-bump `K2=0` references every
256 cycles, and a two-plane random-walk state-space smoother that marginalizes
finite reference-calibration noise.  The inverse receives only BPM readbacks,
commands, nominal GTPSA transport/response, and stochastic priors.  Exact
target orbits and all error realizations remain outside the machine-facing
estimator.

For the reconstructed-orbit fixed-GTPSA-template estimator, state-space-
filtered beam-relative center RMSE is 27.081 micrometers and absolute
sextupole-offset RMSE is 28.783 micrometers, with absolute P99 75.207
micrometers.  The hidden-state correction reduces BPM time-state error from
2.632 to 0.320 micrometers.  It changes the displayed center RMSE negligibly
because the balanced parity contrast already cancels first-order drift; the
remaining error is dominated by white noise and static source/model mismatch.
The deterministic static absolute-offset RMSE is 23.232 micrometers.  The
filtered aggregate RMSE passes the 30-micrometer gate, but its 75.207-
micrometer absolute P99 still fails the strict 50-micrometer tail gate.

The complete result and independent validation are in
`../finite_bpm_inversion/results/state_space_sequential_bpm_gtpsa_inverse/SUMMARY.md`
and `../finite_bpm_inversion/results/state_space_sequential_bpm_gtpsa_inverse/VALIDATION.json`.

## Reproduction

From `CESR Project/`:

```powershell
$env:JULIA_PKG_PRECOMPILE_AUTO='0'
julia --project=. `
  sextupole_misalignment/sequential_joint_inverse/generate_joint_machine_scans.jl

julia --project=. `
  sextupole_misalignment/quadrupole_orbit_correction/generate_corrected_joint_machine_scans.jl

julia --project=. `
  sextupole_misalignment/quadrupole_orbit_correction/generate_gtpsa_noisy_corrected_joint_machine_scans.jl

wsl.exe -d Ubuntu-Bmad -- `
  /home/joeyfurina/miniforge3/envs/bmad/bin/python `
  '/mnt/d/Ring_Design_Development/CESR Project/sextupole_misalignment/sequential_joint_inverse/analyze_joint_inverse.py'

wsl.exe -d Ubuntu-Bmad -- `
  /home/joeyfurina/miniforge3/envs/bmad/bin/python `
  '/mnt/d/Ring_Design_Development/CESR Project/sextupole_misalignment/sequential_joint_inverse/analyze_joint_inverse.py' `
  --comparison-case with_quadrupole_misalignment_corrected `
  --output-dir '/mnt/d/Ring_Design_Development/CESR Project/sextupole_misalignment/sequential_joint_inverse/results/joint_inverse_analysis_corrected'

wsl.exe -d Ubuntu-Bmad -- `
  /home/joeyfurina/miniforge3/envs/bmad/bin/python `
  '/mnt/d/Ring_Design_Development/CESR Project/sextupole_misalignment/sequential_joint_inverse/validate_joint_inverse.py' `
  --write-report

wsl.exe -d Ubuntu-Bmad -- `
  /home/joeyfurina/miniforge3/envs/bmad/bin/python `
  '/mnt/d/Ring_Design_Development/CESR Project/sextupole_misalignment/sequential_joint_inverse/validate_joint_inverse.py' `
  --comparison-case with_quadrupole_misalignment_corrected `
  --analysis-dir '/mnt/d/Ring_Design_Development/CESR Project/sextupole_misalignment/sequential_joint_inverse/results/joint_inverse_analysis_corrected' `
  --write-report

wsl.exe -d Ubuntu-Bmad -- `
  /home/joeyfurina/miniforge3/envs/bmad/bin/python `
  '/mnt/d/Ring_Design_Development/CESR Project/sextupole_misalignment/sequential_joint_inverse/analyze_joint_inverse.py' `
  --comparison-case with_quadrupole_misalignment_gtpsa_noisy_corrected `
  --output-dir '/mnt/d/Ring_Design_Development/CESR Project/sextupole_misalignment/sequential_joint_inverse/results/joint_inverse_analysis_gtpsa_noisy_corrected'

wsl.exe -d Ubuntu-Bmad -- `
  /home/joeyfurina/miniforge3/envs/bmad/bin/python `
  '/mnt/d/Ring_Design_Development/CESR Project/sextupole_misalignment/sequential_joint_inverse/validate_joint_inverse.py' `
  --comparison-case with_quadrupole_misalignment_gtpsa_noisy_corrected `
  --analysis-dir '/mnt/d/Ring_Design_Development/CESR Project/sextupole_misalignment/sequential_joint_inverse/results/joint_inverse_analysis_gtpsa_noisy_corrected' `
  --write-report
```

The generated scientific summary is
`results/joint_inverse_analysis/SUMMARY.md`.  The validator checks tensor
shapes, finiteness, paired provenance, zero-excitation closure, the disjoint
machine splits, the complete training/evaluation matrix, per-target coverage,
and saved model/prediction artifacts.

## Interpretation boundary

This is a synthetic sensitivity benchmark, not demonstrated CESR alignment
precision.  The quadrupole-drift magnitude comes from the user's facility
input, while the remaining nuisance distributions are maintained assumptions
rather than measured joint CESR priors.  The drift model is one scalar
target-local mode, BPM white noise is independent, and actuator hysteresis,
K2 polarity asymmetry, settling/outlier masks, missing BPMs, changed optics,
and sim-to-real validation remain outside this first test.  The learned
comparators return point estimates; calibrated posterior covariance and OOD
probability remain future requirements.

The latest lattice emits the documented straight-multipole-in-curved-reference
warning.  No girder pitch is varied here, so the known curved-DQX girder-pitch
discrepancy is not an excitation in this study, but it remains part of the
lattice provenance.

The stochastic augmentation uses an independent balanced-cycle random-walk
realization for each target scan.  It does not yet model one continuous drift
trajectory spanning the complete 76-target acquisition sequence; that
cross-target time-series extension is required before using the joint context
as evidence for operational drift tracking.

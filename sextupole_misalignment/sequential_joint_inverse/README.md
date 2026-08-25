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

## Reproduction

From `CESR Project/`:

```powershell
$env:JULIA_PKG_PRECOMPILE_AUTO='0'
julia --project=. `
  sextupole_misalignment/sequential_joint_inverse/generate_joint_machine_scans.jl

wsl.exe -d Ubuntu-Bmad -- `
  /home/joeyfurina/miniforge3/envs/bmad/bin/python `
  '/mnt/d/Ring_Design_Development/CESR Project/sextupole_misalignment/sequential_joint_inverse/analyze_joint_inverse.py'

wsl.exe -d Ubuntu-Bmad -- `
  /home/joeyfurina/miniforge3/envs/bmad/bin/python `
  '/mnt/d/Ring_Design_Development/CESR Project/sextupole_misalignment/sequential_joint_inverse/validate_joint_inverse.py' `
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

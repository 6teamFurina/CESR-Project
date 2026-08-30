# Finite-BPM sextupole-center inversion

This study succeeds the frozen exact-internal-orbit calculation in
`direct_observable_nuisance_ablation/results/all_76_orbit_protocol/`.  Its goal
is to determine what can be inferred when the measured orbit exists only at a
finite set of BPMs.

All new forward calculations use
`Latest_Lattice/latest_cesr_scibmad_repaired.jl` and SciBmad.  The frozen
all-76 tensor is reused only for staged post-processing experiments with the
same latest lattice and latent machines.

## Two different inverse targets

Let `c_s` be the sextupole magnetic-center offset and `z_s0` the local closed
orbit at the zero-bump, nominal-K2 state.

1. **Beam-relative center**, `d_s = c_s - z_s0`.  This is the bump displacement
   that centers the beam in the sextupole.  It is directly useful for
   beam-based alignment and can be inferred from K2-dependent BPM orbit versus
   known bump commands without measuring `z_s0`.
2. **Absolute mechanical offset**, `c_s`.  This additionally requires a
   BPM-to-local-orbit estimate of `z_s0`.  Without an alignment/model anchor,
   BPM offsets and a common orbit/center translation can be non-identifiable.

These targets must never be mixed in one reported RMSE.

## Staged plan

### Stage A: command-space, finite-BPM beam-relative inversion

`analyze_command_space_finite_bpm.py` fits the K2 slope using only selected BPM
closed-orbit channels.  The five commanded bump coordinates replace exact
internal target coordinates in the fit.  `target_orbits.npy` is loaded only
after fitting to construct the beam-relative truth used for scoring.

The initial ablation uses deterministic ring-uniform BPM subsets with
`1, 2, 4, 8, 16, 32, 64, 111` BPMs.  This is a sensor-count baseline, not an
optimized placement claim.  Later selection should use nominal SciBmad
response information or training-only latent machines and must be evaluated on
held-out machines.

### Stage B: reconstruct nominal local orbit

Estimate `z_s0` from BPM orbit and the lattice model, with `target_orbits.npy`
used only as evaluation truth.  Compare at least:

- nearest/two-sided local BPM transport;
- global regularized orbit-state fit;
- a smoother with latent distributed kicks/model discrepancy.

Report local x/y error and covariance before inserting this estimate into the
sextupole-center inverse.

### Stage C: absolute-offset inversion

Combine the Stage-A relative center with the Stage-B nominal local orbit, then
compare a two-stage estimator against a joint MAP estimator.  Preserve the
other 75 sextupole offsets and the independent quadrupole-strength errors up to
1%.  Add BPM noise, offsets, gains, and missing channels only after the
noise-free finite-location problem is understood.

## Leakage rules

- `bpm_orbits.npy` and known commands may be fit inputs.
- `target_orbits.npy`, target offsets, and nuisance arrays are evaluation/audit
  truth only.
- Any BPM ranking learned from simulated responses must use training
  realizations and be scored on disjoint realizations or latent distributions.
- If the same BPM data are used both to reconstruct local orbit and to fit K2
  response, their correlated uncertainty must be propagated or the channels
  must be split for an ablation.

## Run Stage A

From `CESR Project/`:

```powershell
wsl.exe -d Ubuntu-Bmad -- /home/joeyfurina/miniforge3/envs/bmad/bin/python `
  '/mnt/d/Ring_Design_Development/CESR Project/sextupole_misalignment/finite_bpm_inversion/analyze_command_space_finite_bpm.py'
```

Results are written to `results/command_space_uniform_bpm/`.

## Initial Stage-A result

The first all-76 run completed 4,864 fits. With no BPM measurement errors, the
beam-relative two-dimensional RMSE was `13.913 micrometers` using all 111 BPMs.
Uniform subsets of 8, 16, 32, and 64 BPMs gave `14.052`, `14.108`, `14.013`,
and `14.001 micrometers`, respectively. One BPM produced a poorly constrained
tail; two or more reduced RMSE below `15 micrometers`.

Across nonzero bump states, the nominal commanded displacement differed from
the actual SciBmad local displacement by `22.899 micrometers` two-dimensional
RMS (median `13.026`, P90 `37.964 micrometers`). The BPM-count curve therefore
currently saturates on bump-coordinate/model mismatch rather than on a lack of
BPM output channels. The next test is BPM-conditioned estimation of the local
bump response matrix before optimizing target-specific BPM placement.

## Relative local-orbit predictor comparison

The first BPM-conditioned local-orbit study predicts, for each target
sextupole and each nonzero bump, the nominal-K2 relative local orbit

`delta z_s(b) = z_s(b) - z_s(zero bump)`.

It reuses the frozen all-76 tensor, which contains unknown target offsets,
`300 micrometer` RMS offsets on the other 75 sextupoles, and independent
quadrupole strength errors within `+/-1%`. Exact `target_orbits.npy` values are
loaded only after every prediction and the BPM-only MAP ridge selection are
complete.

Three nominal latest-lattice SciBmad baselines were compared:

1. `command_only`: known bump commands propagated to the target;
2. `two_sided_transport`: command prediction corrected by residual x/y orbit
   at the nearest upstream and downstream BPMs;
3. `global_map`: a command-centered regularized effective-corrector fit to all
   111 BPM residuals, with its ridge ratio chosen by held-out-BPM
   cross-validation.

Across 76 targets, eight latent machines per target, and four nonzero bumps,
the relative local-orbit results are:

| method | x RMSE [micrometers] | y RMSE [micrometers] | 2D RMSE [micrometers] | median [micrometers] | P90 [micrometers] | max [micrometers] |
|---|---:|---:|---:|---:|---:|---:|
| command only | 10.302 | 20.451 | 22.899 | 13.026 | 37.964 | 91.689 |
| two-sided transport | 0.237 | 0.185 | 0.301 | 0.049 | 0.250 | 4.385 |
| global effective-corrector MAP | 7.048 | 13.665 | 15.375 | 8.469 | 24.095 | 96.678 |

The two-sided result improves all 76 targets. Its per-target 2D-RMSE median is
`0.112 micrometer`; 72 targets are at or below `0.5 micrometer`, 74 at or
below `1 micrometer`, and all 76 at or below `2 micrometers`. The largest
two-sided transverse momentum-block condition number is `4.290`, so this
result is not caused by a numerically singular neighbor pair.

The present global MAP improves over command-only on 68 targets but remains
inferior to the local transport estimator. This is a limitation of its
62-corrector discrepancy basis, not evidence that using all BPMs is generally
worse: a local bump is designed to have small BPM leakage, and an effective-
corrector fit can reproduce BPM residuals without preserving the correct
unobserved target displacement. A future global comparison should use a
distributed-kick/state-smoother basis before introducing a neural network.

Generate the nominal model cache from `CESR Project/` with:

```powershell
julia --project=. sextupole_misalignment/finite_bpm_inversion/generate_local_orbit_models.jl
```

Run and validate the comparison with:

```powershell
wsl.exe -d Ubuntu-Bmad -- /home/joeyfurina/miniforge3/envs/bmad/bin/python `
  '/mnt/d/Ring_Design_Development/CESR Project/sextupole_misalignment/finite_bpm_inversion/analyze_local_orbit_predictors.py'

wsl.exe -d Ubuntu-Bmad -- /home/joeyfurina/miniforge3/envs/bmad/bin/python `
  '/mnt/d/Ring_Design_Development/CESR Project/sextupole_misalignment/finite_bpm_inversion/validate_local_orbit_predictor_results.py'
```

Maintained outputs are in `results/local_orbit_predictors/`.

## Two-sided-BPM end-to-end center inversion

The two-sided predicted relative coordinates were propagated through the
maintained all-111-BPM K2-slope center inverse for all 76 targets and eight
latent machines per target. All earlier machine-error settings were retained,
and exact target-local orbit remained evaluation-only.

- completed fits: `608`;
- aggregate x/y/2D center RMSE:
  `3.826 / 4.444 / 5.864 micrometers`;
- median/P90/P99/maximum radial error:
  `3.645 / 9.418 / 17.897 / 25.270 micrometers`;
- per-target RMSE median/P90/maximum:
  `4.325 / 7.968 / 17.569 micrometers`.

For the same tensor, command coordinates gave `13.913 micrometers` 2D RMSE
and exact oracle-local coordinates gave `5.870 micrometers`. Relative to the
oracle, the two-sided center-error-vector difference had `0.192 micrometer`
RMS, `0.025 micrometer` median, `0.181 micrometer` P90, and `1.818 micrometer`
maximum. The correlation between per-case local-orbit prediction RMSE and
center error was `-0.0143`; the remaining center error is therefore dominated
by the maintained source-fit/model limitation rather than by the two-sided
local-orbit estimate.

Run and validate from `CESR Project/` with:

```powershell
wsl.exe -d Ubuntu-Bmad -- /home/joeyfurina/miniforge3/envs/bmad/bin/python `
  '/mnt/d/Ring_Design_Development/CESR Project/sextupole_misalignment/finite_bpm_inversion/analyze_two_sided_center_inversion.py'

wsl.exe -d Ubuntu-Bmad -- /home/joeyfurina/miniforge3/envs/bmad/bin/python `
  '/mnt/d/Ring_Design_Development/CESR Project/sextupole_misalignment/finite_bpm_inversion/validate_two_sided_center_inversion.py'
```

Maintained outputs are in `results/two_sided_center_inversion/`.

## Corrected sequential BPM/GTPSA inverse (2026-08-30)

`analyze_sequential_bpm_gtpsa_inverse.py` now inserts the explicit local-orbit
reconstruction into the complete corrected 16-machine, 76-target workflow.
The physical input is the paired case whose baseline orbit was restored with
the noisy-readback SciBmad/GTPSA ORM.  All sextupole offsets, local-corrector
gains, K2 gains, BPM gains, quadrupole strength errors and rolls, and the
50-micrometer-per-plane quadrupole offsets remain active.  The accepted
baseline corrector command remains fixed while the sextupoles are scanned one
at a time.

For every target, the machine-facing phase uses the zero-bump and bumped BPM
observable readbacks, the known commands, the nominal control responses, and
the nearest-upstream/downstream order-one GTPSA transport.  It saves the
relative target-local orbit, the absolute reference orbit, the beam-relative
center, and the absolute sextupole-offset estimate before loading
`target_orbits.npy`, `reference_target_orbits.npy`, or the latent sextupole
offsets.  Those arrays are evaluation-only.  The validator independently
reconstructs both local-orbit products from BPM data and reports exact
agreement with the saved arrays, as well as a passing structural leakage
check.

The deterministic static-readback result is:

| quantity | 2D RMSE [micrometers] | P99 [micrometers] |
|---|---:|---:|
| relative target-local orbit, nonzero bumps | 14.306 | 36.333 |
| absolute target reference orbit | 7.866 | 35.456 |
| beam-relative center, BPM/GTPSA local orbit | 21.658 | 65.609 |
| beam-relative center, exact-local evaluation oracle | 20.912 | 63.617 |
| absolute sextupole offset, BPM/GTPSA local and reference orbit | 23.346 | 72.376 |

Thus replacing the exact local orbit adds only 0.746 micrometer to aggregate
beam-relative center RMSE in the deterministic full-static-error comparison.
It does not pass the tail gate: the remaining center error is dominated by the
K2-response/source fit and target tails rather than by a catastrophic
local-orbit reconstruction error.

A separate held-out stochastic sensitivity uses the last three machines, 32
measurement realizations, 5-micrometer RMS white noise per BPM plane/read
averaged over 3,072 reads, and a 10-micrometer endpoint-RMS random walk over a
repeated 15-state acquisition.  The local-orbit RMSE remains 13.851
micrometers, but the unfiltered K2-slope center fit reaches 60.667 micrometers
beam-relative RMSE and 165.486 micrometers P99.  This is an acquisition-drift
and slope-weighting failure, not a failure of the BPM/GTPSA local transport.
On that same three-machine subset, the static beam-relative center RMSE is
20.695 micrometers, so the 60.667-micrometer stochastic result is a paired
degradation rather than a comparison against the full 16-machine row.
The 15-state schedule is a transparent baseline rather than an optimized CESR
protocol; the next extension must reuse covariance whitening or the maintained
state-space drift treatment before claiming stochastic precision.

Run and validate from `CESR Project/`:

```powershell
julia --project=. `
  sextupole_misalignment/finite_bpm_inversion/generate_local_orbit_models.jl

python `
  sextupole_misalignment/finite_bpm_inversion/analyze_sequential_bpm_gtpsa_inverse.py

python `
  sextupole_misalignment/finite_bpm_inversion/validate_sequential_bpm_gtpsa_inverse.py
```

The maintained report and validation record are
`results/sequential_bpm_gtpsa_inverse/SUMMARY.md` and
`results/sequential_bpm_gtpsa_inverse/VALIDATION.json`.

## Full-error nominal-ORM state-space inverse (2026-08-30)

`analyze_state_space_bpm_gtpsa_inverse.py` closes the stochastic follow-up on
a newly generated 16-machine, 76-target source tensor.  Orbit correction uses
one nominal theoretical SciBmad/GTPSA ORM of shape 222 by 103.  It is neither
remeasured by a central finite difference nor scaled by any realized BPM or
corrector gain.  The accepted correction command is fixed while the 76
sextupoles are excited one at a time.  Six Julia workers own six independently
loaded mutable latest-lattice models; the production thread-versus-serial
check is exactly zero for BPM, drift-BPM, target, and drift-target coordinates.

The Julia forward generator materializes BPM readbacks with all maintained
static errors embedded.  The Python machine-facing inverse process then opens
only those readbacks, commanded bump/K2 states, nominal order-one GTPSA
response/transport, and declared noise priors.  It never receives the realized
sextupole offsets, BPM/corrector/K2 gains, quadrupole strength/roll/alignment
errors, drift direction or trajectory, or exact target orbit.  Exact target
orbits and offsets are opened only after every estimate is persisted.

The acquisition repeats eight balanced signed signal states 3,072 times.  A
same-bump `K2=0` reference block is inserted every 256 cycles and at the
endpoint.  A two-plane local-orbit random walk is conditioned on those
references; the finite 32-read reference-calibration error is marginalized as
a static nuisance.  The batch Gaussian-conditioning implementation is the
RTS-smoother-equivalent form of the state-space model.

| acquisition and inverse | beam-relative RMSE [micrometers] | relative P99 [micrometers] | absolute-offset RMSE [micrometers] | absolute P99 [micrometers] |
|---|---:|---:|---:|---:|
| deterministic static, noise-floor profiled | 19.470 | 58.138 | 21.343 | 64.823 |
| deterministic static, reconstructed-orbit fixed GTPSA template | 21.126 | 63.269 | 23.232 | 66.133 |
| stochastic, unfiltered fixed GTPSA template | 27.081 | 69.184 | 28.783 | 75.207 |
| stochastic, state-space-filtered fixed GTPSA template | 27.081 | 69.184 | 28.783 | 75.207 |

The state-space correction is active: aggregate BPM time-state error falls
from 2.632 to 0.320 micrometers, an 87.8% reduction.  The final fixed-template
center metric changes by less than 0.001 micrometer at the displayed precision
because the balanced `+,-,-,+` parity contrast already rejects first-order
drift.  The filter therefore cleans the observable time state but does not
materially improve this already balanced center statistic; white noise and
static source/model mismatch dominate the remaining 28.783-micrometer absolute
RMSE.  The finite-calibration BPM/GTPSA local-orbit RMSE is 14.407 micrometers,
and absolute reference-orbit RMSE is 7.974 micrometers.  The 28.783-micrometer
absolute aggregate RMSE passes the maintained 30-micrometer aggregate gate,
but the 75.207-micrometer absolute P99 fails the strict 50-micrometer tail
gate.

Run and validate from `CESR Project/`:

```powershell
julia --threads=auto --project=. `
  sextupole_misalignment/quadrupole_orbit_correction/generate_gtpsa_nominal_corrected_joint_machine_scans.jl

python `
  sextupole_misalignment/finite_bpm_inversion/analyze_state_space_bpm_gtpsa_inverse.py

python `
  sextupole_misalignment/finite_bpm_inversion/validate_state_space_bpm_gtpsa_inverse.py
```

The maintained report and validation record are
`results/state_space_sequential_bpm_gtpsa_inverse/SUMMARY.md` and
`results/state_space_sequential_bpm_gtpsa_inverse/VALIDATION.json`.

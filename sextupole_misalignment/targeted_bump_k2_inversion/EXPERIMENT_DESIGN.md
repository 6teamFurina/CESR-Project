# Targeted bump-by-K2 sextupole-center inversion protocol

## Scientific target

Estimate the two-dimensional magnetic center of one target sextupole in a
machine where the other sextupoles may also be offset. The intervention is a
known two-plane local-orbit bump plus a symmetric scan of only the target
sextupole's `K2`. Simulated offset truth is withheld from the estimator and is
used only for evaluation.

The atomic example is one complete target-magnet scan in one latent machine
realization. A single optics state or BPM row is not an independent example.

## Maximal screening protocol

- Ten representative sextupoles spanning local-SVD condition, `K2`, and ring
  position: `SEX_08W`, `SEX_39W`, `SEX_19W`, `SEX_33E`, `SEX_16E`, `SEX_25W`,
  `SEX_23E`, `SEX_31W`, `SEX_44W`, and `SEX_18E`.
- A `3 x 3` requested target-orbit grid in `x/y`.
- Seven symmetric target-strength settings `0, +/-K, +/-2K, +/-3K`.
- The raw detector/ring scan is retained. Smaller five-point/cross bump and
  three-/five-point `K2` protocols are derived as views of the maximal scan.

The present code starts with a bounded exact smoke case: one target, a
five-point cross bump, and `0,+/-K`. This validates plumbing and truth recovery
before expanding the state count.

## Forward stages

1. Build nominal corrector-space bump knobs from the RF-on corrector-to-orbit
   Jacobian. The target displacement is a hard constraint; global detector
   orbit and corrector norm are minimized.
2. Assign known target truth and, in later stages, background offsets to the
   other 75 sextupoles.
3. For every bump and `K2` state, run an exact scalar SciBmad RF-on closed
   orbit and Twiss calculation.
4. Save the full raw scan, achieved target orbit, commands, masks/metadata,
   and hidden truth.

## Initial inverse baseline

At every bump, fit the full-ring response slope `g = dO/dK2`. Use the completed
target response map

```text
A_j = [d2 O/(dK2 d x_offset), d2 O/(dK2 d y_offset)]
```

to reconstruct two local effective feeddown/alignment coordinates. A beam
bump is equivalent to the opposite magnet displacement, so the physical model
is

```text
g_p - g_nominal = A_j (c_j - b_p).
```

The inferred center from one bump is therefore the reconstructed effective
offset plus the achieved beam bump. All bumps are combined by least squares.
This is the first response-map form of the covariance-weighted feeddown
zero-crossing algorithm. A later local-source matrix can split it further into
dipole, normal-quadrupole, and skew-quadrupole kicks.

## Input views compared

The first executable comparison includes:

1. full `K2` slope vectors;
2. the two alignment response-mode amplitudes;
3. orbit only;
4. orbit plus phase and tune;
5. orbit plus phase, coupling, and tune;
6. all saved detector quantities.

The raw scan tensor is always retained. Direct nonlinear MAP on that tensor,
explicit local kick reconstruction, and learned residual models are later
baselines, not prerequisites for the first truth-recovery check.

## Staged backgrounds

1. target offset only, no measurement noise;
2. all sextupoles with independent offsets;
3. family/spatial correlations and sparse large faults;
4. `K2`, corrector, BPM, quadrupole/roll, drift, and missing-channel nuisance;
5. held-out optics/configuration and real CESR scans where available.

One nuisance family is added at a time before testing the final mixture.

## Forward-model gates

When a GTPSA scan surrogate is introduced, compare it with exact SciBmad using
measured-style covariance `C`. Define

```text
r_sur = ||C^(-1/2)(O_GTSPA - O_exact)|| / sqrt(N_observation).
```

The provisional gate is median `< 0.25` and P95 `< 0.5`. This states that
surrogate error is materially below measurement noise; it is not a position-
accuracy claim.

## Inverse-result gates

Let the required horizontal/vertical center accuracies be `tau_x` and `tau_y`,
to be supplied from CESR physics/operations requirements. Report bias, RMSE,
median and P95 absolute error, catastrophic-failure rate, and interval
coverage separately in `x/y`.

Provisional acceptance criteria are:

```text
|bias_p| < 0.25 tau_p
RMSE_p < tau_p
P(|error_p| > 2 tau_p) < 1%
95% interval empirical coverage between 90% and 98%
```

For protocol or representation compression relative to the maximal/raw
reference, require RMSE ratio `< 1.05`, bias increase `< 0.1 tau`, stable
coverage, and no increase in catastrophic failures.

## Decision logic

- If the nominal local estimator passes with all-sextupole offsets, retain the
  local inverse.
- If baseline-conditioned local inversion passes but nominal inversion does
  not, include baseline orbit/optics context in the final data view.
- If bias remains correlated with other-magnet states, calculate important
  off-diagonal response blocks and add a global refinement.
- If two response modes match full slopes, use modes for the first model while
  preserving raw scans in the dataset.
- If compressed views lose accuracy or robustness, retain the bump-dependent
  `K2` slope or full scan axis.

Development, validation, and test splits are made by complete latent machine
realization. Scans of different target magnets in the same realization remain
in the same split.

## Paired physical-algorithm ladder

The first paired implementation compares the following algorithms on exactly
the same saved scan and latent offset realization:

1. `P0`: fixed nominal GTPSA mixed-derivative inverse;
2. `P1`: background-conditioned mixed response about zero target-offset error;
3. `P2a`: reconstruct two local integrated dipole-source amplitudes, then fit
   the target offset in the conditioned local-source response;
4. `P2b`: reconstruct integrated normal/skew dipole and quadrupole source
   amplitudes, then fit the target offset;
5. `P3`: exact SciBmad full-scan nonlinear Gauss--Newton, run only on a small
   diagnostic subset.

In this first diagnostic, P1--P3 are supplied the saved offsets of the other
75 sextupoles. This deliberate oracle conditioning measures how much error is
caused by nominal propagation and local linearization. A later operational
version must infer or marginalize those nuisance states from baseline data.

The current P2 second stage uses the conditioned linear local-source response.
It does not yet constitute the separate thin-/thick-sextupole analytic
quadratic-kick baseline; that nonlinear local-source fit remains a subsequent
ablation.

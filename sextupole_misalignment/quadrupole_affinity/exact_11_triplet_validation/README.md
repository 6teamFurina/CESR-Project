# Exact 11-condition triplet validation

This directory tests whether the repaired-lattice response-level selection of
three quadrupoles per sextupole survives a finite-amplitude SciBmad protocol.
It uses `Latest_Lattice/latest_cesr_scibmad_repaired.jl` and does not use
Bmad/Tao response data.

## Protocol

For each of the 76 active normal sextupoles, the five retained quadrupoles are
used to form 11 one-at-a-time optics conditions: nominal and `+/-0.1% |K1|`
for each candidate. Every condition contains a `3 x 3` target-bump grid with
`+/-0.5 mm` extrema and five target-sextupole settings
`delta K2 = [-2,-1,0,1,2] x 0.01 m^-3`. The target center truth is
`(+350,-250) um`. The scan therefore contains 495 scalar RF-on SciBmad states
per target and 37,620 states in total.

The current exact observable is the 111-BPM two-plane closed orbit. Five-point
`K2` slopes are fitted first. For each K1 condition, each BPM-plane slope is
then fitted as a two-dimensional quadratic surface in the realized target
orbit. The common center equations from a triplet's nominal and six signed K1
conditions are solved by GLS. All `choose(5,3)=10` triplets are retained.

The precision calculation assumes provisional independent `5 um` raw noise
for every BPM plane and K2 state. That covariance is propagated through the
five-point slope and nine-point quadratic fits into a two-parameter center
Fisher matrix. It is a protocol diagnostic, not a machine precision claim.

## Result

All 76 scans, 760 triplet inversions, and saved arrays pass the validator.

- The response-selected triplet is exact-information rank one for 24/76
  targets and exact-precision rank one for 22/76.
- Those rank changes are not practically significant. Relative to the
  exact-optimal triplet, the response selection's maximum deterministic center
  error penalty is only `6.64e-5 um`, its maximum predicted worst-axis sigma
  penalty is `0.0328 um`, and its maximum information log-determinant gap is
  `8.08e-5`.
- The response-selected deterministic center error spans `1.06--20.28 um`
  with median `5.91 um`. The triplet choice changes this error by at most
  `6.64e-5 um`, so the ten candidates are effectively tied in this observable.
- With the provisional `5 um` raw BPM noise, the predicted closed-orbit-only
  worst-axis center sigma spans `3.17--6.31 mm` with median `4.68 mm`. The
  median raw K2-dependent orbit signal is only `5.74 nm` RMS.
- The five-point K2 response is very linear over `+/-0.02 m^-3`: the median
  linear-fit relative L2 residual is `1.20e-5`, and the median inner-versus-
  outer symmetric-slope difference is `1.82e-6`. A linear K2 slope is adequate
  at this amplitude for the simulated closed-orbit channel.

A matched-measurement-budget audit also changes the interpretation of the
earlier affinity gain. Comparing one nominal plus six K1-conditioned blocks
with seven repeated nominal blocks, rather than with a single nominal block,
gives a median best-triplet log-determinant diversity gain of only
`1.34e-4` and a median worst-axis precision factor of `1.0000295`. Thus the
previous approximately twofold precision factor was almost entirely the
expected benefit of seven independent measurement blocks, not useful optics
diversity from the `0.1%` K1 changes.

The present exact layer therefore **does not validate a unique three-of-five
quadrupole choice**. Keep all five candidates provisional and do not start the
three-knob `3^3` production scan from this ranking. First choose a larger safe
K1 excursion or an alternative optics-knob construction using a matched-budget
objective, and include the direct launch-trajectory phase/coupling channels
used by the response dictionary. Only after that nominal gate shows material
diversity should the fixed other-sextupole, calibration, BPM-noise, and drift
nuisance realizations be run.

## Outputs

- `results/aggregate/exact11_validation_by_sextupole.csv`: response versus
  exact information, precision, and recovery comparisons.
- `results/aggregate/exact11_validation_summary.json`: aggregate exact result.
- `results/aggregate/k2_linearity_by_sextupole.csv` and
  `k2_linearity_summary.json`: five-point K2 diagnostics.
- `results/aggregate/matched_nominal_response_triplets.csv` and
  `matched_nominal_response_summary.json`: equal-measurement-budget response
  audit.
- `results/scans/<target>/`: exact arrays, metadata, timing, and ten triplet
  results for each sextupole.

## Reproduction

Build the latest-lattice two-plane bump knobs once, then run the resumable
target scan (optionally sharded), analyze all target directories, aggregate,
and validate:

```powershell
julia --project=. build_latest_bump_knobs.jl
julia --project=. run_exact_11_targets.jl --targets=all
python analyze_quadratic_bump_center.py --scan-dir=results/scans/sex_09aw
python aggregate_exact_11_results.py
python analyze_k2_linearity.py
python analyze_matched_nominal_gain.py
python validate_exact_11_results.py
```

The generated bump knobs are model-based and have not been approved as CESR
machine knobs. The current inverse also uses simulated target-orbit coordinates
for the bump surface; a later machine-facing implementation must replace these
with a calibrated command-to-local-orbit model or an inferred local orbit.

## Follow-up: `+/-1% |K1|`

A latest-lattice SciBmad optics audit subsequently tested all 113 independent
quadrupole knobs at `+/-1% |K1|`. With the maintained limits of `0.01` maximum
tune shift and `20%` maximum detector beta beating, 110 knobs pass. Three are
rejected by tune shift: `Q10W` (`0.010391`), `Q43E` (`0.011031`), and `Q10E`
(`0.011801`). Their maximum beta beating remains below `11.1%`. Among the
original five-candidate pools, 16 targets retain all five, 48 retain four, and
12 retain three safe candidates. The excluded knobs were replaced from the
original 15-candidate response pool before the exact scan, leaving five safe
1% candidates for every target.

The complete 76-target, 37,620-state exact scan was then repeated at 1%. The
symmetric K1 response was decomposed into

```text
odd  = 0.5 (R_plus - R_minus)
even = 0.5 (R_plus + R_minus) - R_nominal.
```

For the five-point-K2 closed-orbit slope, `RMS(even)/RMS(odd)` has minimum,
median, P90, and maximum `0.00262 / 0.00829 / 0.02180 / 0.05046`. No even term
is detectable under the provisional independent `5 um` raw BPM-noise model.
Thus 1% is sufficiently linear for the 110 knobs that pass the optics audit,
although the three-point symmetric test does not constrain cubic K1 terms and
the maximum deterministic curvature remains about 5% for one target-candidate
pair.

The 1% exact inversion has materially larger *relative* K1 diversity than the
0.1% scan, but the absolute closed-orbit-only benefit remains small. Against
seven repeated nominal blocks, the response-selected triplet's median
log-determinant gain is `0.001765` and median worst-axis precision factor is
`1.000448`; the exact-optimal values are `0.002018` and `1.000505`. The maximum
exact-optimal precision factor is `1.002356`. The median predicted worst-axis
sigma changes only from `4.685 mm` at 0.1% to `4.681 mm` at 1%, and the median
deterministic center error remains `5.91 um`.

Decision: use 1% as the next response/inversion amplitude for the 110 audited
safe knobs, while excluding the three tune-limit failures or assigning them a
separately audited smaller amplitude. Preserve plus and minus conditions as
separate exact blocks. Do not interpret this as validation of a unique
three-quadrupole set: closed orbit alone remains too insensitive. The next
discriminating calculation is the 1% direct launch-trajectory phase/coupling
ablation, followed by nuisance realizations only if its matched-budget gain is
material.

The 1% results are under `results/aggregate_k1_1pct`, the full scans under
`results/scans_k1_1pct`, the optics screen under
`results/k1_1pct_optics_screen`, and the replacement five-candidate sets under
`results/k1_1pct_safe_selection`.

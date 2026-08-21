# GTPSA derivative inverse with drift and BPM white noise

This study replaces the noisy, per-BPM self-normalized center fit with a fixed
two-parameter response model.  The model combines the exact local normal-
sextupole `dO/dK2` source polynomial with the latest-lattice order-1
SciBmad/GTPSA cumulative and one-turn maps.  BPM noise and scan-time drift then
enter through their measurement covariance rather than being reabsorbed into
one fitted scale per BPM channel.

The benchmark contains fixed hidden x/y offsets on all 76 active normal
sextupoles.  No quadrupole, corrector-gain, K2-calibration, BPM-calibration, RF,
or other lattice error is active.  The two stochastic additions considered
here are 5 micrometers RMS independent white noise per BPM plane and read, and a
scalar Gaussian random walk along a fixed random two-plane local-orbit mode
with 10 micrometers RMS endpoint change over a complete target scan.  These are
synthetic sensitivity settings, not measured CESR priors.

## Fixed GTPSA derivative model

Let `Pn` and `Ps` be the full-BPM response templates obtained by transporting
the two local normal-sextupole sources

```text
q_n = (x^2 - y^2) / 2,       q_s = x y
```

from a target sextupole through the periodic ring.  The integrated kick signs
and the 0.272 m sextupole length are included in the templates.  For the
beam-relative magnetic center `c = (cx, cy)`, take a symmetric K2 slope and
then the bump-odd part:

```text
S(b) = [O(b,K+) - O(b,K-)] / (K+ - K-)
Gx   = [S(+bx) - S(-bx)] / (2 b)
Gy   = [S(+by) - S(-by)] / (2 b)
```

The center enters linearly:

```text
Gx = -Pn cx - Ps cy
Gy = -Ps cx + Pn cy .
```

After stacking every BPM x/y channel, this is `g = A c`.  The estimator is the
two-parameter generalized least-squares solution

```text
c_hat = (A' C^-1 A)^-1 A' C^-1 g .
```

K2-odd projection removes K2-independent orbit and static BPM offsets;
bump-odd projection removes bump-even terms.  Most importantly, `A` is fixed
by the physics model and is not renormalized from each noisy scan.

This route is distinct from the direct high-order K2--target-offset GTPSA map
that terminates in the GTPSA C layer with `sqrt(0)` at least at `SEX_14W` and
`SEX_44E`.  It uses the stable order-1 transport maps plus the exact
local `dO/dK2` polynomial and therefore covers all 76 targets.  It is not a
claim that the direct high-order-map failure has been repaired.

## Separate treatment of the two stochastic errors

For white BPM noise with per-read RMS `sigma`, `R` reads of each signed state,
K2 span `Delta K`, and bump span `2b`, every stacked gradient channel has

```text
Var(G) = [2 sigma / (sqrt(R) Delta K (2b))]^2 .
```

The code propagates this variance through `(A'A)^-1`.  A measured, nonuniform
or correlated BPM covariance can replace the present scalar covariance without
changing the inverse.

For drift, the eight signed states are acquired in two time-balanced
`+,-,-,+` blocks, one per bump plane.  Constant and linear common additive
drift have zero contrast weight.  The remaining state-dependent random-walk
response is retained.  If `w_t` is the complete fixed center-estimator weight
of read `t`, the covariance is evaluated exactly as

```text
Cov(c_drift) = q sum_j W_j W_j',       W_j = sum_(t>=j) w_t,
```

where `q` is chosen to keep the full-scan endpoint RMS at 10 micrometers.  A
reverse-cumulative closed form makes this calculation independent of `R`; the
validator also constructs a short read sequence explicitly and obtains the
same covariance.

## SciBmad benchmark and result

The exact paired reference uses
`Latest_Lattice/latest_cesr_scibmad_repaired.jl`, 76 targets, 4 hidden
all-sextupole-offset realizations per target, five bump states
`(0,0),(+/-1.5 mm,0),(0,+/-1.5 mm)`, and K2 values
`-0.10, 0, +0.10 m^-3`.  Baseline and drift-secant scans contain 9,120 exact
RF-on scalar SciBmad states in total.

The conservative default uses 4,096 reads per signed state, or 32,768 reads per
target scan.  Across 512 independent stochastic measurement seeds:

| case | 2D RMSE [um] | median [um] | P90 [um] | P99 [um] |
|---|---:|---:|---:|---:|
| clean | 12.761 | 8.323 | 20.627 | 38.946 |
| BPM white noise | 18.728 | 14.427 | 28.859 | 45.681 |
| random-walk drift | 13.978 | 9.836 | 22.332 | 37.966 |
| combined | 19.575 | 15.225 | 30.163 | 47.032 |

The white-noise and drift contributions beyond the clean model error are 13.70
and 5.76 micrometers RMS respectively.  In the combined case, every target's
own RMSE is below 50 micrometers; the worst is 43.55 micrometers.  Combined P99
is 47.03 micrometers and 99.385% of Monte Carlo center errors are below 50
micrometers.  Gaussian noise is unbounded, so a finite sample maximum cannot be
guaranteed below a fixed threshold; the acceptance gate uses aggregate RMSE,
P99, and every target-level RMSE.

The analytic repeat-count tradeoff is:

| reads/state | reads/target | combined RMSE [um] | worst target RMSE [um] |
|---:|---:|---:|---:|
| 64 | 512 | 110.777 | 146.165 |
| 256 | 2,048 | 56.581 | 76.228 |
| 1,024 | 8,192 | 30.765 | 51.819 |
| 1,280 | 10,240 | 28.219 | 49.769 |
| 2,048 | 16,384 | 23.898 | 46.525 |
| 4,096 | 32,768 | 19.584 | 43.639 |

Thus 1,280 reads/state already clears the RMSE-only target gate, while the
4,096-read default supplies the tested P99 margin.  Machine time is not the
same as analysis runtime: the default requires 2,490,368 BPM acquisitions for
all 76 targets, before state-switching overhead.  The measurement duration is
therefore that count divided by the usable BPM acquisition rate.

Exact SciBmad generation took 646.7 s for the paired 9,120 states.  Once those
states and the GTPSA maps exist, the full covariance propagation, 512-seed
Monte Carlo, tables, and acceptance checks take about 1.3 s.  The code never
replays 4,096 raw arrays per state.

## Repeated eight-state time-series inverse

`analyze_time_series_inverse.py` extends the same eight signed signal states to
a true acquisition-order model.  The scalar physical-orbit drift evolves at
every read and is never reset at a state or cycle boundary.  The full-BPM
white-noise fit is written as covariance-matched GLS.  Under the maintained
equal independent 5 um BPM noise this reduces exactly to the original OLS
matched filter, while the implementation also accepts nonuniform diagonal or
full BPM covariance matrices.

The drift inverse retains the eight-state signal cycle as the basic protocol.
Every 256 cycles, and once at the final endpoint, it replaces that cycle with
four same-bump `delta K2 = 0,+,0,-,0` blocks.  Thus each reference cycle still
contains the eight `K+/-` signal states and adds twelve `K2=0` observations.
The four reference baselines are independently calibrated with 32 reads each;
their finite calibration errors are nuisance states in the filter rather than
treated as exact.  Later references update the already accumulated x/y center
error, which is equivalent to smoothing the time history for the final center.

For 3,072 reads per signal state, the core-only scan has 24,576 acquisitions
per target.  The periodic references increase this to 24,732, plus 128
separate calibration reads.  The random-walk step variance is held fixed per
acquisition, so the slightly longer reference scan has 10.032 um rather than
10.000 um endpoint RMS.  Across the same 76 targets, four hidden machines per
target, and 512 center-error draws, the selected result is:

| estimator | 2D RMSE [um] | P99 [um] | worst target RMSE [um] |
|---|---:|---:|---:|
| balanced eight-state, no drift inverse | 21.108 | 50.502 | 44.694 |
| periodic-reference state-space inverse | 20.297 | 48.681 | 44.561 |

The stochastic drift component falls from 5.757 to 0.319 um.  The filtered
result passes both the preferred aggregate 30 um gate and the conservative
gate requiring aggregate RMSE, P99, and every target-level RMSE below 50 um.
The 20.297 um RMSE is 6.766% of the requested 300 um reference scale; this is a
scale-normalized proxy, not a per-example relative error for centers close to
zero.  Including the 128 calibration reads, the selected protocol uses 24,860
acquisitions per target, or 1,889,360 for all 76 targets before state-switching
and settling overhead.
At 1,280 signal reads/state the analytic filtered RMSE is already 27.632 um,
but the 3,072 setting is retained because it also passes the tested P99 gate.
The unfiltered P99 at the same 3,072 setting is just above 50 um, so the strict
gate crossing is specifically due to the time-series drift inverse.

This remains a synthetic sensitivity result.  The reference cadence is
optimal only for the assumed single drift mode, 5 um BPM noise, uniform
cadence, and 10 um core-scan endpoint drift.  A machine cadence must instead be
chosen from fixed-state BPM time series and revalidated with measured
multidirectional drift and full BPM covariance.

## Compound nuisance revalidation excluding quadrupole misalignment

`analyze_compound_nuisances.py` repeats the nuisance test at the current
excitation amplitudes instead of reusing the older 0.5 mm / 0.02 m^-3 scans.
The exact latest-lattice SciBmad tensors cover all 76 targets, four paired
hidden machines per target, signed 1.5 mm bumps, and K2 extrema +/-0.10 m^-3.
The static compound machine uses exactly the same component draws as the five
one-at-a-time cases.  Quadrupole misalignment is identically zero.

| deterministic case | 2D RMSE [um] | paired increment RMS [um] | P99 [um] |
|---|---:|---:|---:|
| clean reference | 12.761 | 0.000 | 38.946 |
| BPM gain, 1% RMS | 12.752 | 0.369 | 38.919 |
| corrector gain, 1% RMS | 14.476 | 6.421 | 38.265 |
| K2 calibration gain, 1% RMS | 13.394 | 3.016 | 39.134 |
| quadrupole strength, independent +/-1% | 23.825 | 19.384 | 69.750 |
| quadrupole roll, 1 mrad RMS | 13.953 | 5.172 | 37.016 |
| paired linear sum | 25.873 | 21.443 | 68.293 |
| exact static compound | 25.897 | 21.458 | 65.552 |

There is no large nonlinear compound surprise: the exact compound-minus-
linear-sum residual is 1.565 um RMS, 4.834 um at P99, and 7.174 um maximum.
It is 7.293% of the actual compound increment.  The dominant static increment
is quadrupole strength, not BPM gain or a cross-nuisance instability.

Applying the repeated eight-state acquisition and state-space drift inverse to
the same compound machines gives:

| time-series case | 2D RMSE [um] | P99 [um] | worst target RMSE [um] |
|---|---:|---:|---:|
| clean reference + white noise + filtered drift | 20.332 | 48.714 | 44.368 |
| all tested nuisances except quadrupole misalignment | 30.334 | 73.738 | 66.560 |

Two targets, `SEX_09AW` and `SEX_38E`, have time-series RMSE at or above 50
um.  The complete case therefore fails both the preferred aggregate 30 um
gate and the aggregate/P99/all-target 50 um gate.  This is a tail and static-
model problem, not failure of the stochastic mitigation: the compound
white-noise component is 15.814 um, filtered drift is 0.324 um, and their
15.817 um total is unchanged to three decimals from the clean-reference
stochastic component.  Consequently, the non-misalignment nuisances cannot
yet be described collectively as solved; quadrupole-strength conditioning or
machine-specific relinearization is the next required model extension.

## What to try next

- For environmental drift, promote the one fixed local-orbit mode to a small
  measured drift basis and solve the centers and time coefficients jointly
  with a random-walk prior (equivalently, colored-noise GLS or a Kalman
  smoother).  Then optimize the signed-state order against that measured
  covariance.  Interspersed K2=0 references are useful only when their drift
  information exceeds the white noise they add.
- For BPM white noise, measure the full BPM covariance, whiten the response,
  discard unstable or contaminated singular modes, and allocate repeats by
  target information rather than uniformly.  The present equal 5-micrometer
  model makes ordinary and generalized least squares identical.
- A reliable higher-order observation map can subtract finite-amplitude bias
  and permit larger bump/K2 products.  Since center variance scales roughly as
  `1 / [R (Delta K b)^2]`, this is the most direct path to fewer acquisitions.
  Until the direct high-order GTPSA failure is fixed, the full-76 alternative is
  the already validated scan-fitted high-order observation Taylor map.

## Run and validate

From this directory:

```powershell
julia --project=. generate_exact_protocol_scans.jl

python analyze_stochastic_inverse.py
python validate_stochastic_inverse.py
python analyze_time_series_inverse.py
python validate_time_series_inverse.py

julia --project=. generate_exact_protocol_scans.jl `
  --cases=corrector_gain,k2_calibration,quadrupole_strength,quadrupole_roll,combined_without_quadrupole_misalignment,combined_without_quadrupole_misalignment_time_drift

python analyze_compound_nuisances.py
python validate_compound_nuisances.py
```

The first command regenerates the exact SciBmad tensors and is normally needed
only when the lattice, scan amplitudes, or latent-machine ensemble changes.
`results/analysis/SUMMARY.md`, `summary.csv`, `protocol_tradeoff.csv`, and
`per_target_summary.csv` contain the maintained outputs.
The time-series extension writes the corresponding artifacts, the explicit
20-slot reference-cycle schedule, and covariance arrays under
`results/time_series_analysis/`.
The compound revalidation writes its paired decomposition, stochastic
covariances, per-target tails, and acceptance flags under
`results/compound_nuisance_analysis/`.

## Limitations

This is a synthetic numerical sensitivity result, not demonstrated real-
machine precision.  The benchmark has four latent machines per target; the
drift is first order in amplitude, restricted to one fixed local-bump mode per
latent machine, and assumes fixed total endpoint RMS as repeat count changes.
A deployment study needs measured BPM cadence/covariance, nonlocal and
multidirectional drift modes, K2 and corrector readbacks, missing/outlier BPMs,
and enough independent machine configurations to quantify sampling
uncertainty.

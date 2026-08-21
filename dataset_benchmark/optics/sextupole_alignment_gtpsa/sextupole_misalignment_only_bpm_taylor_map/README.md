# Sextupole-misalignment-only BPM/Taylor-map benchmark

This benchmark isolates the inverse problem requested on 2026-08-19. Its only
machine error is a fixed x/y misalignment on all 76 active normal sextupoles.
There is no BPM noise, BPM calibration error, time drift, missing channel,
corrector calibration error, target-K2 calibration error, quadrupole error, or
additional RF/lattice error. Every forward state uses SciBmad and the validated
latest lattice `Latest_Lattice/latest_cesr_scibmad_repaired.jl`.

## Scan and estimand

For each of the 76 target sextupoles, one independent latent machine is drawn:
the target x/y increment is uniform in +/-350 um and the other 75 sextupole
offsets are independent 300 um RMS Gaussians. The errors are held fixed through
a 5 x 5 two-plane bump grid (+/-0.5 mm) and five target-K2 levels
(-0.02, -0.01, 0, 0.01, 0.02 m^-3). This gives 125 exact RF-on SciBmad states
per target and 9,500 states in total.

The deployable input is the 111-BPM orbit tensor plus known bump and K2
commands. Exact target-local orbit and target alignment are loaded only after
all machine-facing fits and are used only for evaluation.

The primary truth is the beam-relative magnetic center

```text
c_rel = c_nominal + delta_c_target - z_target(zero bump, zero delta-K2).
```

An absolute incremental offset is reported separately by adding the BPM-based
estimate of the zero-bump target orbit. The two estimands are never mixed in a
single RMSE.

## Inverse methods

Let `O_j(b,k)` be BPM channel `j`, `b=(bx,by)` the commanded bump, `k` the K2
increment, and `z_hat(b,k)` the two-sided BPM prediction of the target-local
orbit.

1. `fd_*_source_predicted` first extracts `dO/dK2` by a linear or quartic K2
   fit, then applies the maintained physical two-source sextupole inverse using
   `z_hat(b,0)`.
2. `quadratic_o_derivative_predicted` directly fits each observable derivative
   as `S_j(z)=a_j+g_j^T z+1/2 z^T H_j z` and solves the scaled stacked system
   `H_j c = -g_j`. `chain_rule_o_derivative_predicted` instead fits in command
   space and transforms the gradient/Hessian with the fitted `dz_hat/db`
   Jacobian.
3. `o_taylor_orderN_*` fits the raw K2-dependent observation surface
   `O_j(x,y,k)` with every identifiable total-degree monomial through order N,
   differentiates the fitted map analytically with respect to K2, and finds the
   common two-dimensional root. These are empirical observation Taylor maps
   fitted to exact SciBmad scan states; they are not direct GTPSA coefficients.
4. `generate_gtpsa_k2_offset_maps.jl` separately builds direct nominal-lattice
   SciBmad/GTPSA derivatives in K2 and target x/y offset. A normal-form-free
   solver obtains the parameter-dependent periodic fixed point order by order
   from the one-turn map before propagation to every BPM.

## Full 76-target result

The two-sided BPM predictor is substantially more accurate than every center
inverse in this idealized benchmark:

| local-orbit quantity | 2D RMSE [um] | median [um] | P90 [um] | maximum [um] |
|---|---:|---:|---:|---:|
| relative, zero K2, nonzero bumps | 0.040180 | 0.007369 | 0.057320 | 0.326648 |
| relative, all 125 states | 0.039371 | 0.006678 | 0.056163 | 0.328788 |
| absolute zero-bump reference | 0.010091 | 0.002330 | 0.014513 | 0.043684 |

| machine-facing inverse | beam-relative 2D RMSE [um] | median [um] | P90 [um] | absolute-increment RMSE [um] |
|---|---:|---:|---:|---:|
| finite-difference, linear K2 source | 6.394706 | 5.000718 | 8.956349 | 6.394942 |
| finite-difference, quartic K2 source | 6.394716 | 5.001335 | 8.955497 | 6.394951 |
| direct quadratic `dO/dK2(z_hat)` | **5.128185** | **3.738924** | **7.549159** | **5.127145** |
| chain-rule `dO/dK2(b)` via `dz_hat/db` | 5.465936 | 3.981950 | 7.771729 | 5.465113 |
| empirical O Taylor, total order 3 | 33.008657 | 27.946649 | 47.192103 | 33.009199 |
| empirical O Taylor, total order 4 | 6.362654 | 4.656776 | 9.719398 | 6.362141 |
| empirical O Taylor, total order 5 | 5.724636 | 4.490103 | 8.573020 | 5.723937 |

The exact-local-orbit oracle gives 6.394493 um for the quartic source inverse
and 6.361597 um for the order-4 Taylor inverse, essentially unchanged from the
BPM-predicted variants. Thus local-orbit inference is not the limiting error in
this no-BPM-error benchmark. The direct observable-derivative fit is the best
of the tested full-inventory methods. Total order 3 is inadequate; order 4
removes most of its truncation bias and order 5 improves further, but neither
beats the quadratic observable-derivative fit on this scan.

These values are a one-realization-per-target synthetic baseline, not an error
bar for the real machine. More latent realizations are required before quoting
sampling uncertainty.

## Direct GTPSA status

For `SEX_09AW`, the new order-by-order fixed-point map and SciBmad's periodic
Twiss/normal-form map agree over all ten saved mixed-derivative blocks with
relative L2 difference `9.568e-14`; the final fixed-point residual is
`1.536e-13`. On that one smoke case, retained offset orders 1/2/3 give center
errors of 229.8795/0.58993/1.12361 um. This demonstrates that the quadratic
offset block is essential, but one target is not an inventory-level accuracy
claim.

The direct high-order GTPSA path still terminates in the GTPSA C layer with
`invalid domain sqrt(0)` for at least `SEX_14W` and `SEX_44E`, including when
the transverse normal-form calculation is bypassed. The exception cannot be
caught in Julia. Therefore direct GTPSA is retained as a validated subset
diagnostic, while the complete 76-target high-order comparison uses the
scan-fitted SciBmad observation Taylor maps. The lower-level target-dependent
tracking/coordinate-transform failure remains unresolved and must be fixed
before claiming a complete direct-GTPSA result.

## Reproduction and validation

From `CESR Project/`:

```powershell
julia --project=. dataset_benchmark/optics/sextupole_alignment_gtpsa/sextupole_misalignment_only_bpm_taylor_map/generate_only_sextupole_scans.jl `
  --targets=all `
  --output-dir=dataset_benchmark/optics/sextupole_alignment_gtpsa/sextupole_misalignment_only_bpm_taylor_map/results/exact_scans
```

Run `analyze_inverses.py` in the validated WSL `bmad` environment, then run:

```powershell
python dataset_benchmark/optics/sextupole_alignment_gtpsa/sextupole_misalignment_only_bpm_taylor_map/validate_results.py
```

The validator independently checks the error declaration, tensor inventory,
truth mapping, finiteness, every saved summary, derivative solve, Taylor design
rank, and the fixed-point-versus-Twiss GTPSA smoke comparison. Detailed results
are in `results/analysis/SUMMARY.md`, `summary.csv`, `per_case_estimates.csv`,
and `fit_diagnostics.csv`.

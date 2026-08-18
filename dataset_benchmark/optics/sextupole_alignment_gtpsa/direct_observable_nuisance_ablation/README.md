# Direct-observable nuisance ablation

This directory compares sextupole-center fits from closed orbit alone with
fits that additionally use directly measurable excitation readbacks.  All
primary calculations use the repaired latest SciBmad CESR lattice.

The atomic example is a complete 9-orbit-bump by 5-K2 scan for one target in
one fixed latent machine.  Each latent machine contains:

- independent unknown x/y offsets on the other 75 active sextupoles (300 µm
  RMS per plane);
- a randomized unknown target offset;
- independent physical strength errors on all 113 active quadrupole knobs,
  uniform within ±1%;
- no target-local-orbit measurement error at this stage.

The inverse does not receive the sextupole or quadrupole nuisance truth.  It
uses the actual target-local orbit reached by each bump and estimates K2 slopes
within the complete tensor.

## Direct readbacks

- phase/beta inputs: BPM trajectory differences from fixed position/angle
  launches, rather than simulator `beta` or `phi` labels;
- coupling: cross-plane components of the same launch-response measurement;
- tune: frequency fitted from a synthesized TBT trajectory;
- dispersion: BPM orbit difference under a symmetric fixed-energy probe;
- chromaticity: TBT tune difference under that energy probe;
- orbit response: exact closed-orbit finite differences of two physical CESR
  correctors.

The lattice RF cavities are harmonic masters. Directly changing their stored
`rf_frequency` does not define an independent equilibrium-energy state in the
present SciBmad model, so the dispersion/chromaticity readback explicitly uses
a fixed beam-energy delta probe. It is not labeled as an RF-frequency scan.

## Reproduce

From `CESR Project/`:

```powershell
$env:JULIA_PKG_PRECOMPILE_AUTO='0'
julia --project=. dataset_benchmark\optics\sextupole_alignment_gtpsa\direct_observable_nuisance_ablation\generate_paired_scan.jl
```

The analysis requires NumPy and Matplotlib. The validated WSL environment can
run it with:

```powershell
wsl.exe -d Ubuntu-Bmad -- /home/joeyfurina/miniforge3/envs/bmad/bin/python `
  '/mnt/d/Ring_Design_Development/CESR Project/dataset_benchmark/optics/sextupole_alignment_gtpsa/direct_observable_nuisance_ablation/analyze_paired_scan.py' `
  --turns=512 --fft-size=32768
```

## Current pilot

The maintained result is
`results/sex_09aw_paired_pilot/`: eight paired latent realizations for
`SEX_09AW`, totaling 360 primary K2/bump states plus the corrector probe solves.
See its `SUMMARY.md`, `fit_summary.csv`, and plots.

The current equal-block fit is deliberately a structural, noise-free pilot.
It is not a machine-precision claim and does not substitute for a measured
covariance. In particular, a centered sextupole still has an intrinsic
chromatic response, so chromaticity cannot be inserted into the same
zero-at-center feed-down equation used for phase, tune, coupling, and ORM.

## Orbit protocol subsampling

`analyze_protocol_subsampling.py` refits strict subsets of the maintained
9-bump by 5-K2 tensors with a common thin-sextupole source model. On the same
eight latent worlds, 3 versus 5 K2 points changes the noise-free RMSE by less
than `0.0001 µm`. A five-point axial cross gives `2.939 µm`, whereas the full
nine-point grid gives `4.249 µm`; this is model-discrepancy behavior, not proof
that fewer measurements are intrinsically better. A five-point corners-plus-
center geometry gives `16.507 µm`, showing that bump geometry matters more
than count under the present fit.

The paired K1 ablation in `results/sex_09aw_k1_orbit_ablation/` uses the same
eight target truths and quadrupole-error realizations, the five-point cross,
three K2 points, and nominal plus ±1% one-at-a-time changes of `QX4D`, `Q18W`,
and `Q24E`. The common-center RMSE changes from `2.940 µm` without K1 scanning
to `2.960 µm` with all seven K1 conditions. Thus K1 scanning adds no useful
orbit-only improvement in this pilot while multiplying the state count by
seven.

## All-76 economical protocol

`generate_all_targets_orbit_protocol.jl` runs all 76 targets in one Julia
session with eight nuisance realizations per target, five axial-cross bumps,
three outer K2 levels `(-2K,0,+2K)`, and nominal K1 commands. The maintained
`results/all_76_orbit_protocol/` calculation contains 9,120 exact SciBmad
states and 608 center fits. Generation took 553.6 s.

The aggregate 2D center RMSE is `5.870 µm`; realization median/P90/P99 are
`3.664/9.395/17.897 µm`, and the maximum is `25.274 µm`. Per-target RMSE has
median `4.300 µm`, P90 `7.933 µm`, and maximum `17.574 µm`. Forty-six targets
are at or below `5 µm`, 72 at or below `10 µm`, and four exceed `10 µm`.
`SEX_13E`, `SEX_24W`, `SEX_28E`, and `SEX_17W` are the four largest-RMSE
targets and should be the first cases for a nine-bump/model-discrepancy follow-
up. See `results/all_76_orbit_protocol/SUMMARY.md` and the per-target CSV.

This result is frozen as the **oracle-local-orbit baseline**.  The inverse
reads the exact SciBmad x/y orbit at the target sextupole; it is therefore a
conditional inverse floor rather than predicted machine accuracy.  The full
frozen protocol, artifact roles, and reuse boundary are documented in
`results/all_76_orbit_protocol/FROZEN_BASELINE.md`.  Finite-BPM work that does
not use internal orbit as a fit input continues in the sibling
`../finite_bpm_inversion/` directory.

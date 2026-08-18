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
  '/mnt/d/Ring_Design_Development/CESR Project/dataset_benchmark/optics/sextupole_alignment_gtpsa/finite_bpm_inversion/analyze_command_space_finite_bpm.py'
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

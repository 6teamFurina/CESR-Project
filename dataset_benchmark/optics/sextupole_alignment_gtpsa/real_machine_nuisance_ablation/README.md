# Real-machine nuisance ablation for finite-BPM sextupole alignment

This paired study starts from a clean reference that retains the central hard
condition: the target alignment is unknown and the other 75 sextupoles carry
independent `300 micrometer` RMS x/y offsets. Each non-reference row adds
exactly one nuisance. Physical lattice/actuator cases were regenerated with
the validated latest repaired SciBmad lattice; BPM gain and BPM noise were
applied only to simulated readbacks. The inverse always receives nominal bump
and K2 commands plus BPM readings, never latent nuisance values or exact
target-local orbit.

The magnitudes below are representative sensitivity-test settings, not measured
CESR calibration distributions. There are 76
targets and 4 paired latent
machines per target, or 304
fits per row.

## Result table

| added machine error | test magnitude | local-orbit 2D RMSE [um] | center 2D RMSE [um] | aggregate change [um] | paired error-vector increment RMS [um] | median [um] | P90 [um] | max [um] |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| Reference: no added nuisance | none | 0.040 | 6.051 | +0.000 | 0.000 | 4.056 | 8.945 | 24.161 |
| BPM gain | 1% RMS per BPM/plane, fixed in scan | 4.767 | 6.774 | +0.723 | 2.897 | 4.469 | 10.241 | 23.703 |
| Corrector gain | 1% RMS per corrector, fixed in scan | 0.057 | 6.080 | +0.029 | 0.476 | 3.974 | 9.037 | 24.082 |
| K2 calibration gain | 1% RMS intervention gain per target scan | 0.040 | 6.051 | -0.000 | 0.001 | 4.056 | 8.946 | 24.161 |
| Quadrupole strength | independent uniform +/-1% | 0.346 | 5.847 | -0.204 | 1.519 | 3.956 | 8.694 | 24.036 |
| Quadrupole roll | 1 mrad RMS | 0.126 | 5.829 | -0.222 | 2.026 | 4.045 | 8.514 | 27.178 |
| Quadrupole misalignment | 100 um RMS in x and y | 0.236 | 1181.981 | +1175.930 | 1181.175 | 224.169 | 1917.196 | 6204.899 |
| Time drift | +/-5 um target-local linear drift across scan | 0.045 | 1985.647 | +1979.596 | 1985.597 | 2025.733 | 2241.208 | 2499.565 |
| BPM noise | 5 um RMS per BPM/plane/state | 9.649 | 948.518 | +942.467 | 948.746 | 445.897 | 1705.802 | 2496.450 |

The paired error-vector increment compares each nuisance realization with its
matched reference and is the most direct one-at-a-time impact measure. It must
not be confused with the difference between aggregate RMSE values: a nuisance
can partially cancel the reference fit error and lower aggregate RMSE while
still changing individual estimates.

## Interpretation checks

- A fixed multiplicative K2 calibration gain has negligible effect here because
  the maintained fit normalizes every BPM K2-slope channel and fits a free
  propagation matrix. This does not cover K2 hysteresis, polarity asymmetry, or
  point-to-point calibration drift.
- Quadrupole strength and roll slightly lower aggregate RMSE in this finite
  paired sample, but their nonzero paired increment quantifies the actual
  estimate change; they are not beneficial corrections.
- The quadrupole-misalignment row is an **uncorrected-orbit stress test**, not a
  transfer-matrix-only result. Its 100-micrometer RMS offsets push 290/304
  beam-relative truths outside the 0.5-mm bump radius and 176/304 outside the
  current +/-1.5-mm-per-plane fit box. A separate orbit-corrected misalignment
  study is needed to isolate residual matrix mismatch.
- Time drift and BPM noise are intentionally passed to the unchanged three-point
  slope estimator without drift regression, repeated-read averaging, or
  covariance weighting. Their large errors diagnose protocol sensitivity, not
  the best achievable calibrated-machine performance.

## Nuisance definitions

- **BPM gain:** independent multiplicative x/y calibration error, fixed for a
  BPM throughout one scan.
- **Corrector gain:** independent multiplicative error on each corrector's bump
  increment, fixed throughout one scan; the base corrector setting is unchanged.
- **K2 calibration:** one multiplicative error on the target's commanded K2
  intervention, shared by every bump and K2 point in the scan.
- **Quadrupole strength:** independent physical Kn1 errors uniformly bounded by
  `+/-1%`.
- **Quadrupole roll:** independent roll added coherently to every tracking slice
  belonging to the same physical quadrupole.
- **Quadrupole misalignment:** independent x/y displacement added coherently to
  every tracking slice belonging to the same physical quadrupole.
- **Time drift:** a random transverse direction with a target-local command that
  changes linearly from `-5` to `+5 micrometers` in acquisition order. It is zero
  at the zero-bump, nominal-K2 reference state and is propagated physically by
  the same local-bump correctors.
- **BPM noise:** independent Gaussian readout noise for every BPM plane and scan
  state.

The reference intentionally omits quadrupole strength error, unlike the earlier
maintained all-76 result. This makes every row a one-at-a-time ablation; its
absolute reference RMSE therefore need not equal the earlier `5.864 micrometer`
mixed-nuisance result.

## Method and provenance

- lattice: `D:\Ring_Design_Development\CESR Project\Latest_Lattice\latest_cesr_scibmad_repaired.jl`
- exact physical SciBmad states per physical case: `4560`
- summed physical-generation wall time: `2053.8 s`
- local orbit: nominal-model command prediction corrected by the nearest
  upstream/downstream BPM pair
- center inverse: all-111-BPM symmetric three-point K2 slope and the maintained
  shared thin-sextupole source fit
- exact target orbit and target alignment: evaluation only

Run from `CESR Project/`:

```powershell
julia --project=. dataset_benchmark/optics/sextupole_alignment_gtpsa/real_machine_nuisance_ablation/generate_physical_nuisance_scans.jl

wsl.exe -d Ubuntu-Bmad -- /home/joeyfurina/miniforge3/envs/bmad/bin/python `
  '/mnt/d/Ring_Design_Development/CESR Project/dataset_benchmark/optics/sextupole_alignment_gtpsa/real_machine_nuisance_ablation/analyze_nuisance_ablation.py'

wsl.exe -d Ubuntu-Bmad -- /home/joeyfurina/miniforge3/envs/bmad/bin/python `
  '/mnt/d/Ring_Design_Development/CESR Project/dataset_benchmark/optics/sextupole_alignment_gtpsa/real_machine_nuisance_ablation/validate_nuisance_ablation.py'
```

## Limitations

This is still synthetic and noise magnitudes are assumed. Each error is tested
alone, so the table is a sensitivity decomposition rather than a prediction of
the fully combined real-machine error. A later combined study should use
measured calibration priors, correlated girder/family errors, interleaved scan
timing, and covariance-aware or joint-nuisance inference.

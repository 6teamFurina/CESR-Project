# Archived direct-observable and K1 pilots

Archived on 2026-09-05 from `direct_observable_nuisance_ablation/`. These two
completed design-stage pilots use the repaired latest SciBmad CESR lattice,
`Latest_Lattice/latest_cesr_scibmad_repaired.jl`. They preserve the evidence
behind the orbit-only observable choice and omission of active K1 scanning.
Their raw arrays, result tables, and generation metadata were moved unchanged.

Each scan holds the target truth, other-sextupole offsets, and quadrupole
strength errors fixed while commanded bump/K2 states are varied. The inverse
fits local sextupole feed-down using exact target-local orbit coordinates.
These are noise-free oracle-local-orbit comparisons, not finite-BPM machine
precision claims. The all-target frozen baseline and shared fitting functions
remain in [`../../direct_observable_nuisance_ablation/`](../../direct_observable_nuisance_ablation/README.md).

## Direct-observable paired pilot

This pilot tests only `SEX_09AW`. Each latent realization uses:

| Input/intervention | Values | Points |
|---|---|---:|
| Target-sextupole `delta K2` | `-0.02, -0.01, 0, +0.01, +0.02 m^-3` | 5 |
| Commanded local orbit bump | A 3×3 grid with x/y values `-0.5, 0, +0.5 mm` | 9 |
| Quadrupole `K1` command | Nominal; no active K1 scan | 1 |

Each realization therefore contains `9 × 5 = 45` primary scan states. Eight
realizations give 360 primary states, in addition to the corrector, launch,
and energy probes required for the direct readbacks.

## Direct readbacks

These quantities are produced through their corresponding measurement
processes. Internal simulator labels such as ideal `beta`, `phase`, or a
coupling coefficient are not inserted directly as machine observations.

| Physical quantity | Directly measured quantity used in this study |
|---|---|
| Phase/beta information | BPM trajectory differences following known position/angle launches |
| Coupling | Cross-plane components of the same launch-response measurement |
| Tune | Spectral peak or fitted frequency from a synthesized turn-by-turn BPM trajectory |
| Dispersion | BPM orbit difference under a symmetric fixed-energy `delta=+/-0.001` probe |
| Chromaticity | Tune shift under the same fixed-energy probe |
| Orbit response matrix | Closed-orbit finite differences of the physical correctors `HKICK_9AW` and `VX6D`, using a `1e-6` probe field |

The RF cavities in the current SciBmad lattice are harmonic masters. Directly
changing their stored `rf_frequency` does not define an independent
equilibrium-energy state. The dispersion and chromaticity measurements
therefore use a fixed beam-energy delta probe and are not labeled as an RF
frequency scan.

## Paired-pilot results

`results/sex_09aw_paired_pilot/` contains eight paired latent realizations for
`SEX_09AW`. The observable blocks currently use equal structural block
weighting because measured covariance has not yet been supplied.

| Fit inputs | 2D RMSE | Median | P90 | Paired wins over orbit only / 8 |
|---|---:|---:|---:|---:|
| Orbit only | 2.829 µm | 2.930 µm | 3.879 µm | — |
| Orbit + feed-down direct readbacks | 6.923 µm | 5.203 µm | 9.802 µm | 2 |
| Orbit + all direct except chromaticity | 10.298 µm | 9.071 µm | 15.877 µm | 2 |
| Orbit + all direct | 185.943 µm | 62.342 µm | 264.801 µm | 0 |

This is not evidence that the direct observables contain no useful
information. It shows that measurements with different physical forward
models cannot be inserted with equal weight into one zero-at-center feed-down
equation without suitable response models and covariance. In particular, a
centered sextupole has an intrinsic chromatic response, so its chromaticity K2
slope does not obey the zero-center relation used by ordinary feed-down
channels.

## Separate K1 ablation

`results/sex_09aw_k1_orbit_ablation/` separately tests nominal strength and
one-at-a-time `+/-1%` changes of `QX4D`, `Q18W`, and `Q24E`, giving seven K1
conditions. Across the same eight `SEX_09AW` realizations, the orbit-only RMSE
changes from `2.940 µm` without a K1 scan to `2.960 µm` with all seven K1
conditions. The economical protocol therefore omits active K1 scanning. This
does not remove the latent quadrupole strength errors.

## Bump and K2 subsampling

`analyze_protocol_subsampling.py` refits strict subsets of the same 9-bump by
5-K2 tensors. Across the present eight noise-free latent worlds, three versus
five K2 points changes RMSE by less than `0.0001 µm`. The five-point axial
cross gives `2.939 µm`, the full nine-point grid gives `4.249 µm`, and five
corners plus center gives `16.507 µm`. The cross-versus-grid reversal is model
discrepancy between the shared thin-source fit and diagonal finite-amplitude
states; it must not be interpreted as evidence that less data are inherently
better.


## Files and reproduction

- `generate_paired_scan.jl` and `analyze_paired_scan.py` operate on
  `results/sex_09aw_paired_pilot/`.
- `generate_k1_orbit_ablation.jl` and `analyze_k1_orbit_ablation.py` operate on
  `results/sex_09aw_k1_orbit_ablation/`.
- The shared subsampling and source-fitting functions stay in
  `../../direct_observable_nuisance_ablation/analyze_protocol_subsampling.py`.
- Both Julia generators include the maintained
  `../../quadrupole_affinity/exact_11_triplet_validation/common.jl` and read its
  nominal `results/bump_knobs/local_bump_knobs.csv`.

From `CESR Project/`, historical forward reproduction uses the pinned Julia
environment:

```bash
julia --startup-file=no --project=. sextupole_misalignment/archived_methods/direct_observable_k1_pilots/generate_paired_scan.jl
julia --startup-file=no --project=. sextupole_misalignment/archived_methods/direct_observable_k1_pilots/generate_k1_orbit_ablation.jl
```

Existing generator outputs are protected unless overwrite is explicitly enabled.
Python analysis requires NumPy, SciPy, and Matplotlib and writes derived results
inside the selected input directory. To reanalyze the preserved tensors:

```bash
python -B sextupole_misalignment/archived_methods/direct_observable_k1_pilots/analyze_paired_scan.py --turns=512 --fft-size=32768
python -B sextupole_misalignment/archived_methods/direct_observable_k1_pilots/analyze_k1_orbit_ablation.py
python -B sextupole_misalignment/direct_observable_nuisance_ablation/analyze_protocol_subsampling.py
```

Original absolute paths in metadata are provenance records. Current code uses
archive-relative result paths and the maintained shared dependencies above.

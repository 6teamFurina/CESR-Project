# Direct-observable nuisance ablation

This directory studies whether the two-dimensional magnetic center of a target
sextupole can be recovered by changing its `K2` strength and the local orbit
bump, then fitting directly measurable beam responses. It also tests how
closed orbit and additional direct observables behave in the presence of
unknown offsets on the other sextupoles and quadrupole strength errors.

All forward lattice calculations use the repaired latest SciBmad CESR lattice,
`Latest_Lattice/latest_cesr_scibmad_repaired.jl`. Python is used only to read
the SciBmad output, perform NumPy/SciPy inversion, and make plots. Bmad, Tao,
and PyTao are not used to generate the physics results in this study.

## Inverse problem

For a target sextupole `s`, the quantity to be recovered is its two-dimensional
magnetic-center offset

```text
c_s = (c_x,s, c_y,s).
```

One inverse example is a complete orbit-bump by K2 scan tensor from one fixed
latent machine realization, rather than one BPM value or one lattice state.
The inverse first estimates the K2 response of every retained measurement
channel. It then fits the normal-sextupole feed-down structure across five or
nine different local orbit positions to recover `c_s`.

The frozen 76-target result uses the exact SciBmad internal closed-orbit x/y
coordinates at the target sextupole as fit inputs. It is therefore an
**oracle-local-orbit baseline**, not an accuracy claim for a machine in which
only BPM orbits are available. The finite-BPM successor study is maintained in
`../finite_bpm_inversion/`.

## Actively varied inputs

Two main protocols must be distinguished.

### A. Direct-observable paired pilot

This pilot tests only `SEX_09AW`. Each latent realization uses:

| Input/intervention | Values | Points |
|---|---|---:|
| Target-sextupole `delta K2` | `-0.02, -0.01, 0, +0.01, +0.02 m^-3` | 5 |
| Commanded local orbit bump | A 3×3 grid with x/y values `-0.5, 0, +0.5 mm` | 9 |
| Quadrupole `K1` command | Nominal; no active K1 scan | 1 |

Each realization therefore contains `9 × 5 = 45` primary scan states. Eight
realizations give 360 primary states, in addition to the corrector, launch,
and energy probes required for the direct readbacks.

### B. All-76 economical orbit-only protocol

This protocol treats each of the 76 active normal sextupoles as the target in
turn:

| Input/intervention | Values | Points |
|---|---|---:|
| Target-sextupole `delta K2` | `-0.02, 0, +0.02 m^-3` | 3 |
| Commanded local orbit bump | `(-x,0),(0,-y),(0,0),(0,+y),(+x,0)`, with 0.5 mm amplitude | 5 |
| Quadrupole `K1` command | Nominal; no active K1 scan | 1 |

Eight independent latent realizations are generated for every target, giving

```text
76 targets × 8 realizations × 5 bumps × 3 K2 = 9,120 SciBmad states.
```

“Nominal K1” means that the quadrupole commands are not actively scanned. The
physical lattice in every realization still contains random quadrupole
strength errors.

### Separate K1 ablation

`results/sex_09aw_k1_orbit_ablation/` separately tests nominal strength and
one-at-a-time `+/-1%` changes of `QX4D`, `Q18W`, and `Q24E`, giving seven K1
conditions. Across the same eight `SEX_09AW` realizations, the orbit-only RMSE
changes from `2.940 µm` without a K1 scan to `2.960 µm` with all seven K1
conditions. The economical protocol therefore omits active K1 scanning. This
does not remove the latent quadrupole strength errors.

## Measured observables

### Closed orbit

Every scan state stores the horizontal and vertical equilibrium closed orbit
at 111 BPMs:

```text
m = (x_1, y_1, ..., x_111, y_111).
```

There are therefore 222 closed-orbit channels per state. The all-76 protocol
uses only these BPM closed-orbit responses as measured fit outputs.

### Additional direct readbacks in the paired pilot

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

## Latent truth and nuisance variables

The following latent quantities remain fixed within one complete scan tensor
and are regenerated for the next realization:

| Category | Count and distribution | Supplied to the inverse? |
|---|---|---|
| Target-sextupole offset | Independent x/y values uniformly distributed over `[-350,+350] µm` | No; this is the prediction target |
| Other sextupole offsets | The other 75 sextupoles; independent Gaussian x/y offsets with `300 µm` RMS per plane | No; hidden nuisance |
| Quadrupole strength errors | All 113 active quadrupole knobs; independent uniform relative errors over `[-1%,+1%]` | No; hidden nuisance |

The other 75 sextupole offsets and all quadrupole strength errors are actually
applied to the SciBmad forward lattice, so they alter the generated orbit and
direct observables. The inverse does not receive their random truth and does
not explicitly recover them as `150 + 113` joint fit parameters. They act as
unknown background nuisance variables against which the target-center fit is
tested.

The following measurement nuisances have not yet been included:

- BPM position noise, offsets, gain/roll errors, and missing BPMs;
- corrector calibration and readback errors;
- target-local-orbit reconstruction error;
- measured covariance and long-term drift.

The present results are therefore noise-free structural benchmarks. In
particular, the frozen all-76 inverse reads the exact target-local x/y orbit in
order to isolate and validate the sextupole-response inverse before adding
local-orbit reconstruction.

## Direct-observable paired-pilot results

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

## Frozen results for all 76 sextupoles

`results/all_76_orbit_protocol/` uses five axial bumps, three outer K2 values,
nominal K1 commands, and orbit-only measured outputs. It contains 9,120 exact
SciBmad states and 608 center inversions. SciBmad forward generation took
553.6 s.

### Aggregate error

| Metric | 2D center error |
|---|---:|
| Aggregate RMSE over all 608 inversions | **5.870 µm** |
| Median over all realizations | **3.664 µm** |
| P90 over all realizations | **9.395 µm** |
| P99 over all realizations | **17.897 µm** |
| Maximum single-realization error | **25.274 µm** |
| Median of the 76 per-target RMSE values | **4.300 µm** |
| P90 of the per-target RMSE values | **7.933 µm** |
| Maximum per-target RMSE | **17.574 µm** |

### Target coverage

| Threshold | Number of targets |
|---|---:|
| Per-target RMSE `<= 5 µm` | **46 / 76** |
| Per-target RMSE `<= 10 µm` | **72 / 76** |
| Per-target RMSE `> 10 µm` | **4 / 76** |

### Four targets with the largest per-target RMSE

| Target | RMSE | P90 | Maximum error |
|---|---:|---:|---:|
| `SEX_13E` | **17.574 µm** | 24.762 µm | 25.274 µm |
| `SEX_24W` | **13.436 µm** | 16.873 µm | 20.086 µm |
| `SEX_28E` | **12.242 µm** | 19.390 µm | 20.087 µm |
| `SEX_17W` | **12.018 µm** | 15.195 µm | 15.877 µm |

The complete per-target table is stored in
`results/all_76_orbit_protocol/per_target_summary.csv`; individual realization
fits are stored in `per_realization_fits.csv`.

## Interpretation and limitations

The `5.870 µm` result is an **oracle baseline** conditional on exact
target-local orbit. It demonstrates that the K2/orbit protocol contains enough
sextupole-response signal for micrometer-scale center recovery when the other
75 sextupole offsets and ring-wide quadrupole strength errors are present. It
does not predict the final accuracy obtainable from a finite set of machine
BPMs.

The frozen protocol, artifact roles, and reuse boundary are documented in
`results/all_76_orbit_protocol/FROZEN_BASELINE.md`. The next stage in
`../finite_bpm_inversion/` forbids internal target orbit as a fit input. It
first recovers the beam-relative center, then reconstructs the nominal local
orbit from finite BPM data, and finally compares two-stage and joint-MAP
absolute-offset inversion.

## Reproduce

Run the SciBmad forward calculation from `CESR Project/`:

```powershell
$env:JULIA_PKG_PRECOMPILE_AUTO='0'
julia --project=. sextupole_misalignment\direct_observable_nuisance_ablation\generate_paired_scan.jl
```

The Python analysis requires NumPy, SciPy, and Matplotlib:

```powershell
python sextupole_misalignment\direct_observable_nuisance_ablation\analyze_paired_scan.py `
  --turns=512 --fft-size=32768
```

The Python step performs inversion and plotting only; it does not run Bmad,
Tao, or PyTao.

## Other protocol comparisons

`analyze_protocol_subsampling.py` refits strict subsets of the same 9-bump by
5-K2 tensors. Across the present eight noise-free latent worlds, three versus
five K2 points changes RMSE by less than `0.0001 µm`. The five-point axial
cross gives `2.939 µm`, the full nine-point grid gives `4.249 µm`, and five
corners plus center gives `16.507 µm`. The cross-versus-grid reversal is model
discrepancy between the shared thin-source fit and diagonal finite-amplitude
states; it must not be interpreted as evidence that less data are inherently
better.


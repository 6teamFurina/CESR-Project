# Frozen all-target direct-orbit baseline

This directory retains the all-76-sextupole orbit-only oracle baseline and the
shared thin-source fitting functions used by related studies. The early
direct-observable selection and K1-scan pilots were archived on 2026-09-05 in
[`../archived_methods/direct_observable_k1_pilots/`](../archived_methods/direct_observable_k1_pilots/README.md).

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

## Scan protocol

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

## Measured observables

### Closed orbit

Every scan state stores the horizontal and vertical equilibrium closed orbit
at 111 BPMs:

```text
m = (x_1, y_1, ..., x_111, y_111).
```

There are therefore 222 closed-orbit channels per state. The all-76 protocol
uses only these BPM closed-orbit responses as measured fit outputs.

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

## Reproduction

Run from `CESR Project/` with the pinned Julia environment and a Python
environment containing NumPy, SciPy, and Matplotlib:

```bash
julia --startup-file=no --project=. sextupole_misalignment/direct_observable_nuisance_ablation/generate_all_targets_orbit_protocol.jl
python -B sextupole_misalignment/direct_observable_nuisance_ablation/analyze_all_targets_orbit_protocol.py
```

The generator protects existing output unless overwrite is explicitly enabled.
The Python analysis reads the saved SciBmad states and performs inversion and
plotting only.

## Archived pilot studies and shared code

The [pilot archive](../archived_methods/direct_observable_k1_pilots/README.md)
contains the direct-observable comparison, K1-scan negative result, their four
experiment scripts, and both saved pilot tensors. It also documents the
bump/K2 subsampling comparison.

`analyze_protocol_subsampling.py` remains here because the all-target analyzer
imports its `fit_center` and `k2_slope` functions. Its command-line default now
reads the relocated pilot tensor; this path change does not alter the shared
fitting functions. The all-target frozen tensor remains at
`results/all_76_orbit_protocol/` for the finite-BPM comparisons.

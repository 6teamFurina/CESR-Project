# CESR matched dataset benchmark

This benchmark tests a specific digital-twin workload: map 119 CESR horizontal
and vertical corrector commands to horizontal and vertical closed orbit at 99
`DET_*` markers. Bmad and SciBmad consume the same deterministic input CSV and
produce the same labeled `1000 x 198` observable table.

It is not designed to prove that one physics engine is universally faster.
The bounded claim is whether the present CESR interfaces can generate this
matched, exact closed-orbit dataset with higher throughput while preserving
numerical agreement.

## Fairness rules

1. Run both engines on the same physical machine whenever possible.
2. Use the same input CSV, RF state, CESR lattice, detector order, and output
   schema.
3. Keep one model process alive for the entire run.
4. Batch all 119 PyTao variable updates so Tao performs one lattice
   recalculation per sample rather than 119 recalculations.
5. Report initialization, compilation/warmup, physics, and file writing
   separately.
6. Require every sample to converge and compare every one of the 198 outputs.
7. Record CPU model, operating system, thread count, Julia, Bmad, Tao, and
   PyTao versions with the final result.

## Generate the shared inputs

From `CESR Project`:

```console
python dataset_benchmark/orbit/generate_control_samples.py
```

Sample 0 is the zero-control baseline. Samples 1 through 999 use independent
Gaussian corrector commands with `sigma = 5e-6 rad`, clipped at three standard
deviations. The generator seed and full control order are written to the
adjacent metadata JSON.

The committed formal input is:

```text
inputs/cesr_corrector_samples_1000.csv
```

## Response-validity sweep in normalized input radius

`run_response_rho_sweep.jl` measures the error of the direct first-order
`198 x 119` detector-orbit response against converged nonlinear SciBmad RF-on
closed orbits. It runs three input scenarios: all 119 correctors, only the 58
horizontal correctors, and only the 61 vertical correctors.

For each scenario, a Gaussian direction is normalized to exact unit RMS over
the active controls. The applied kick is

```text
delta_k = rho * base_kick * normalized_gaussian_direction
```

so every trial at a requested `rho` has exactly that active-control RMS radius.
The inactive controls are zero, and the same random directions are reused at
every rho. By default, `base_kick = 5e-6 rad`, so `rho = 1` means an active-knob
RMS kick of `5 microrad`. The horizontal-only and vertical-only definitions use
RMS over their active plane, not over all 119 columns.

The extended sweep uses 600 trials at each of 20 positive, approximately
log-spaced rho values from `0.1` through `64`. The zero-control nominal state is
computed once and shared by all three scenarios, giving
`1 + 3 x 20 x 600 = 36,001` nonlinear states. The single Julia runner remains
useful for a short range or one chunk:

```powershell
$env:JULIA_PKG_PRECOMPILE_AUTO='0'

julia --threads=4 --project=. `
  dataset_benchmark/orbit/run_response_rho_sweep.jl `
  --rhos=0,0.1,0.14,0.2,0.28 `
  --trials=600
```

For the complete extended range, a single batched Newton solve can be held back by a few
hard high-rho lanes. The maintained Windows workflow therefore divides the 20
positive rho values into five four-rho chunks. Four independent Julia processes
run concurrently with four threads each; every chunk has 7,201 states, just
above the corresponding `1750 x 4 = 7000` BeamTracking threading threshold.
Each chunk is saved immediately and can be rerun independently:

```powershell
powershell -ExecutionPolicy Bypass -File `
  dataset_benchmark/orbit/run_response_rho_sweep_parallel.ps1
```

Merge the five completed chunks before plotting:

```powershell
python dataset_benchmark/orbit/merge_response_rho_sweep_chunks.py `
  --root=dataset_benchmark/orbit/error_analysis/response_rho_sweep_600/chunks `
  --output-dir=dataset_benchmark/orbit/error_analysis/response_rho_sweep_600/combined
```

For a short validation run:

```powershell
julia --threads=auto --project=. `
  dataset_benchmark/orbit/run_response_rho_sweep.jl `
  --rhos=0,0.5,1 `
  --trials=4 `
  --output-dir=dataset_benchmark/orbit/error_analysis/response_rho_sweep_smoke
```

The runner writes per-trial errors, a per-scenario/rho summary, and TOML
metadata. For each plane, a trial error is the RMSE over all 99 detectors. The
summary reports the mean trial RMSE and the maximum trial RMSE; the plot uses
the latter as a one-sided upper error cap. It also retains the maximum absolute
single-detector error separately.

Plot the completed result with:

```powershell
python dataset_benchmark/orbit/render_response_rho_sweep_svg.py
```

By default, the normal and `--normalize-rho-squared` variants are written to
`error_analysis/response_rho_sweep_600/figures` with distinct filenames.

The plotting command uses only the Python standard library and writes SVG, so
it does not require Matplotlib or another plotting package. By default it
looks first for the merged 600-trial result under
`error_analysis/response_rho_sweep_600/combined`, while retaining the earlier
`results/response_rho_sweep_600/combined` and `results/response_rho_sweep`
locations as fallbacks. The figure uses a
colorblind-safe palette together with distinct markers and line styles for
grayscale readability. It includes a quadratic reference guide anchored at
the smallest complete all-corrector point, standard panel labels, sparse log
ticks, and a direct statement of the physical input scale. If the sweep used a
nondefault base kick, pass it explicitly, for example
`--base-kick-urad=2.5`.

Add `--normalize-rho-squared` to plot `RMSE / rho^2`; a flat curve identifies
quadratic truncation error. Points whose exact reference did not converge for
all 600 trials are marked with a high-contrast cross.

### Current 600-trial extended result

The 2026-08-03 run generated 36,001 unique nonlinear states in five saved
chunks. The first four four-thread Julia processes ran concurrently; the fifth
high-rho chunk ran separately with four threads. All three scenarios converged
for every trial through `rho = 12.8`. At `rho = 18.1`, the all-control and
vertical-only scenarios each had one failure out of 600. Across the full range,
35,901/36,001 states converged, 302 lanes requested full-AD fallback, and 100
lanes remained failed after fallback. Horizontal-only converged 600/600 through
`rho = 64`.

For the complete-reference range through `rho = 12.8`, the all-control X and Y
mean-RMSE local log slopes remain close to 2. The horizontal-only X and Y errors
remain nearly quadratic across the entire converged range. Vertical-only Y
changes smoothly from a slope near 2 at small rho toward a slope near 3 by
`rho = 12.8`, showing a genuine higher-order contribution before solver
failures begin. High-rho all-control and vertical-only means are calculated
only from converged trials and must not be interpreted as unbiased full-sample
errors; the plots mark those incomplete points explicitly.

## SciBmad run

```console
julia --threads=auto --project=. \
  dataset_benchmark/orbit/benchmark_scibmad.jl
```

The SciBmad runner represents each corrector as a `BatchParam` containing all
sample values, solves the RF-on closed orbits as one batch, and then tracks the
solved orbits once around the ring to collect detector coordinates.

By default, the runner calculates the nominal closed orbit and the `6 x 119`
closed-orbit control response, uses `z0 + R * delta-k` as each sample's
initial guess, and runs the frozen-nominal-Jacobian solver with automatic
full-AD fallback. The default tolerances are the normal/Bmad-default values
`reltol=1e-8` and `abstol=1e-10`.

To instead solve the zero-control nominal closed orbit and use the same `z0`
for every sample:

```console
julia --threads=1 --project=. \
  dataset_benchmark/orbit/benchmark_scibmad.jl \
  --initial-guess=nominal-z0
```

The runner records the nominal-orbit cost and the Newton iteration statistics
separately. The original full-AD, zero-initial-guess configuration remains
available explicitly:

```console
julia --threads=1 --project=. \
  dataset_benchmark/orbit/benchmark_scibmad.jl \
  --initial-guess=zero \
  --jacobian-mode=full
```

To calculate the first-order `6 x 119` closed-orbit response
`R = dz_closed/dk`, give every sample its own `z0 + R * delta-k` initial
guess, reuse the nominal Jacobian, and automatically send any failed lanes
back through full-AD Newton:

```console
julia --threads=1 --project=. \
  dataset_benchmark/orbit/benchmark_scibmad.jl \
  --initial-guess=response-linear \
  --jacobian-mode=frozen-nominal \
  --reltol=1e-8 \
  --abstol=1e-10
```

The generated response matrix is saved beside the metadata as
`closed_orbit_response_6x119.csv`. Response construction, nominal-orbit
solution, recurring physics, and warmup/compilation are timed separately.
The default run first loads the validated cache at
`reference/closed_orbit_response_6x119.csv`; it checks the `6 x 119` shape,
coordinate labels, control names and exact control-column order. If the cache
is missing, the runner calculates it with GTPSA and writes it to that path.
Force a refresh after changing the lattice, RF configuration, energy, or
control definitions with:

```console
julia --threads=1 --project=. \
  dataset_benchmark/orbit/benchmark_scibmad.jl \
  --recompute-response=true
```

An alternative cache can be selected with
`--response-matrix-cache=path/to/response.csv`.

## Bmad/Tao run

Run this in the Linux environment containing Bmad, Tao, and PyTao:

```console
python dataset_benchmark/orbit/benchmark_bmad.py
```

The runner keeps one Tao instance alive. PyTao's `cmds(...,
suppress_lattice_calc=True)` applies all 119 controls before triggering one
model recalculation for each sample. This avoids an artificial 119-fold
interface penalty.

## Compare

After both output files are on the same system:

```console
python dataset_benchmark/orbit/compare_outputs.py
```

The report includes convergence, RMSE, maximum orbit difference, correlation,
per-sample relative errors, physics throughput, and the SciBmad/Bmad throughput
ratio.

## Readability evaluation

Line count alone is not a reliable readability measure across Julia and
Python/Tao. Use a small blinded maintenance exercise with CESR researchers:

1. add vertical BPM noise with a supplied seed;
2. add one quadrupole calibration-error input;
3. add tune as an output;
4. identify and report a failed closed-orbit sample.

For each anonymized implementation, record time to a correct change, number of
API concepts consulted, test failures before completion, and whether the
result preserves labels, units, and provenance. The side-by-side source should
still be included in the report, but the maintenance exercise provides a more
defensible readability result than subjective preference or raw lines of code.

## Current limitation

The first implementation supports RF-on only. The existing SciBmad RF-off
coasting-orbit adapter is a scalar ForwardDiff workaround and must be extended
or replaced before an equivalent batched RF-off benchmark is claimed.

## Current results

All rows used the same 1000 samples and compared all 198 detector coordinates
per sample. `Maximum output residual` is the largest absolute detector value
difference from the Bmad table over all 198,000 comparisons. `Maximum closure
residual` is the final six-dimensional one-turn closure norm; `--` means that
the older run did not record this diagnostic.

| Result | Machine | Tolerances `(rel, abs)` | Initial guess / Jacobian | Converged | Physics time (s) | Samples/s | Maximum output residual vs Bmad (m) | Correlation vs Bmad | Maximum closure residual |
|---|---|---|---|---:|---:|---:|---:|---:|---:|
| Bmad/Tao/PyTao | `lnx201` Linux | defaults `(1e-8, 1e-10)` | previous orbit; Tao one-turn matrix reuse | 1000/1000 | 67.370 | 14.843 | 0 | 1 | -- |
| Bmad/Tao/PyTao | local WSL2, Ryzen 9 5900HX | defaults `(1e-8, 1e-10)` | previous orbit; Tao one-turn matrix reuse | 1000/1000 | **11.858** | **84.330** | -- [1] | -- [1] | -- |
| SciBmad, high precision | `lnx201` Linux | `(1e-13, 1e-13)` | zero; full batched AD Jacobian | 1000/1000 | 280.486 | 3.565 | `8.138494e-6` | `0.999999966415495` | -- |
| SciBmad, response + frozen + fallback | `lnx201` Linux | `(1e-8, 1e-10)` | per-sample `z0 + R*delta-k`; frozen nominal Jacobian | 1000/1000 | **33.178** | **30.140** | `8.138494e-6` | `0.999999966415499` | `8.104e-11` |
| SciBmad, high precision | local Windows, Ryzen 9 5900HX | `(1e-13, 1e-13)` | zero; full batched AD Jacobian | 1000/1000 | 64.356 | 15.539 | `8.138494e-6` | `0.999999966415495` | -- |
| SciBmad, normal precision | local Windows, Ryzen 9 5900HX | `(1e-8, 1e-10)` | zero; full batched AD Jacobian | 1000/1000 | 26.457 | 37.798 | `8.138494e-6` | `0.999999966415495` | -- |
| SciBmad, frozen + fallback | local Windows, Ryzen 9 5900HX | `(1e-8, 1e-10)` | fixed nominal `z0`; frozen nominal Jacobian | 1000/1000 | 8.163 | 122.506 | `8.138494e-6` | `0.999999966415489` | `9.802e-11` |
| SciBmad, response + frozen + fallback | local Windows, Ryzen 9 5900HX | `(1e-8, 1e-10)` | per-sample `z0 + R*delta-k`; frozen nominal Jacobian | 1000/1000 | **6.855** | **145.885** | `8.138494e-6` | `0.999999966415499` | `8.104e-11` |

[1] The local WSL2 run used the same committed 1,000-sample CSV and lattice.
The archived `lnx201` Bmad output CSV is not present in this checkout, so its
pointwise residual and correlation against the older Bmad run cannot be
recomputed. Timing and convergence come from
`results/formal_1000/bmad/bmad_rf_on_metadata.json` (Bmad/Tao `20260801-1`).

The final response-initialized row is the maintained default. It checks every
lane's closure and reruns only failed lanes with full-AD Newton; the formal
run needed zero fallbacks. Relative to fixed `z0`, the response predictor
reduced median/mean iterations from `3 / 2.994` to `2 / 1.995` and recurring
physics time by `16.0%`.

On `lnx201`, the maintained SciBmad solver used `33.178 s` versus Bmad's
`67.370 s`, a measured `2.031x` same-host physics throughput advantage. Its
cached `6 x 119` response loaded in `0.000730 s`; the original local matrix
construction took `2.389 s` after compilation. Recompute it after changes to
the lattice, RF configuration, energy, or control definitions.

On the local Ryzen 9 5900HX, WSL2 Bmad used `11.858 s` (`84.330` samples/s),
which is `5.681x` faster than the historical `lnx201` Bmad timing. That ratio
is cross-machine and cross-version. Against the existing local Windows
response-initialized SciBmad result (`6.855 s`), SciBmad is `1.730x` faster
than local WSL2 Bmad. These local runs share the physical machine but used
different operating-system runtimes and were not measured back-to-back.

The Bmad and optimized SciBmad measurements share `lnx201`, but were run at
different times on a shared host. The optimized SciBmad process reported only
`49%` average CPU utilization, so a back-to-back repetition is still desirable
for the strongest hardware-normalized claim. Local Windows versus Linux Bmad
timing ratios are cross-machine context, not controlled speedup claims.
Detailed comparisons remain in `results/formal_1000/`.

## Directory layout

```text
dataset_benchmark/orbit/
|-- archive/                 # original transferred Bmad result package
|-- inputs/                  # shared deterministic 1000-sample input
|-- reference/               # Bmad-compatible CESR lattice used remotely
|-- results/
|   |-- preliminary_10/      # Bmad, CPU, threaded CPU, and CUDA checks
|   `-- formal_1000/         # formal Bmad/SciBmad outputs and report
|-- generate_control_samples.py
|-- benchmark_scibmad.jl
|-- benchmark_bmad.py
|-- compare_outputs.py
`-- README.md
```

The exploratory CUDA result is retained only for provenance. GPU support is
not part of the maintained runner or the formal comparison.

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
python dataset_benchmark/generate_control_samples.py
```

Sample 0 is the zero-control baseline. Samples 1 through 999 use independent
Gaussian corrector commands with `sigma = 5e-6 rad`, clipped at three standard
deviations. The generator seed and full control order are written to the
adjacent metadata JSON.

The committed formal input is:

```text
inputs/cesr_corrector_samples_1000.csv
```

## SciBmad run

```console
julia --threads=auto --project=. \
  dataset_benchmark/benchmark_scibmad.jl
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
  dataset_benchmark/benchmark_scibmad.jl \
  --initial-guess=nominal-z0
```

The runner records the nominal-orbit cost and the Newton iteration statistics
separately. The original full-AD, zero-initial-guess configuration remains
available explicitly:

```console
julia --threads=1 --project=. \
  dataset_benchmark/benchmark_scibmad.jl \
  --initial-guess=zero \
  --jacobian-mode=full
```

To calculate the first-order `6 x 119` closed-orbit response
`R = dz_closed/dk`, give every sample its own `z0 + R * delta-k` initial
guess, reuse the nominal Jacobian, and automatically send any failed lanes
back through full-AD Newton:

```console
julia --threads=1 --project=. \
  dataset_benchmark/benchmark_scibmad.jl \
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
  dataset_benchmark/benchmark_scibmad.jl \
  --recompute-response=true
```

An alternative cache can be selected with
`--response-matrix-cache=path/to/response.csv`.

## Bmad/Tao run

Run this in the Linux environment containing Bmad, Tao, and PyTao:

```console
python dataset_benchmark/benchmark_bmad.py
```

The runner keeps one Tao instance alive. PyTao's `cmds(...,
suppress_lattice_calc=True)` applies all 119 controls before triggering one
model recalculation for each sample. This avoids an artificial 119-fold
interface penalty.

## Compare

After both output files are on the same system:

```console
python dataset_benchmark/compare_outputs.py
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
| SciBmad, high precision | `lnx201` Linux | `(1e-13, 1e-13)` | zero; full batched AD Jacobian | 1000/1000 | 280.486 | 3.565 | `8.138494e-6` | `0.999999966415495` | -- |
| SciBmad, response + frozen + fallback | `lnx201` Linux | `(1e-8, 1e-10)` | per-sample `z0 + R*delta-k`; frozen nominal Jacobian | 1000/1000 | **33.178** | **30.140** | `8.138494e-6` | `0.999999966415499` | `8.104e-11` |
| SciBmad, high precision | local Windows, Ryzen 9 5900HX | `(1e-13, 1e-13)` | zero; full batched AD Jacobian | 1000/1000 | 64.356 | 15.539 | `8.138494e-6` | `0.999999966415495` | -- |
| SciBmad, normal precision | local Windows, Ryzen 9 5900HX | `(1e-8, 1e-10)` | zero; full batched AD Jacobian | 1000/1000 | 26.457 | 37.798 | `8.138494e-6` | `0.999999966415495` | -- |
| SciBmad, frozen + fallback | local Windows, Ryzen 9 5900HX | `(1e-8, 1e-10)` | fixed nominal `z0`; frozen nominal Jacobian | 1000/1000 | 8.163 | 122.506 | `8.138494e-6` | `0.999999966415489` | `9.802e-11` |
| SciBmad, response + frozen + fallback | local Windows, Ryzen 9 5900HX | `(1e-8, 1e-10)` | per-sample `z0 + R*delta-k`; frozen nominal Jacobian | 1000/1000 | **6.855** | **145.885** | `8.138494e-6` | `0.999999966415499` | `8.104e-11` |

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

The Bmad and optimized SciBmad measurements share `lnx201`, but were run at
different times on a shared host. The optimized SciBmad process reported only
`49%` average CPU utilization, so a back-to-back repetition is still desirable
for the strongest hardware-normalized claim. Local Windows versus Linux Bmad
timing ratios are cross-machine context, not controlled speedup claims.
Detailed comparisons remain in `results/formal_1000/`.

## Directory layout

```text
dataset_benchmark/
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

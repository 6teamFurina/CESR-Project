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
  --fallback-full-newton=true \
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

All formal runs completed all 1000 samples. On the shared Linux host
`lnx201.classe.cornell.edu`, Bmad/Tao/PyTao used `67.370 s` for timed physics
work (`14.843 samples/s`), while SciBmad with one Julia thread and module-cache
creation disabled used `280.486 s` (`3.565 samples/s`). The observed
same-host result makes Bmad `4.16x` faster in the timed physics region.
However, GNU `time` reports only `49%` average CPU for the SciBmad process, so
the paired runs should be repeated under controlled host load before treating
this as a stable hardware-normalized ratio. See
`results/formal_1000/bmad_scibmad_lnx201_comparison.md`.

For context, the Windows SciBmad run on an AMD Ryzen 9 5900HX used `65.155 s`
(`15.348 samples/s`), but that cross-machine number cannot establish a speedup.
It is retained in `bmad_scibmad_cross_machine_comparison.md`.

All 1000 samples converged in both engines. Across the complete `1000 x 198`
Bmad/SciBmad tables, correlation is `0.999999966415`, global RMSE is
`2.268e-6 m`, and median per-sample relative 2-norm difference is `0.0338%`.
The Linux and Windows SciBmad output tables are exactly equal, with zero
differing numerical entries.

On the local Windows machine, a controlled initial-guess comparison found that
using the nominal closed orbit reduced mean Newton iterations from `4.034` to
`3.101`, but left the maximum at `12`. The measured physics time was
`64.356 s` from a zero guess and `65.709 s` from the nominal orbit, so a shared
`z0` alone did not accelerate the current single batched solve. The two output
tables agree to a maximum absolute difference of `1.574e-13 m`. See
`results/formal_1000/scibmad_initial_guess_comparison.md`.

Aligning SciBmad's tolerance values with Bmad's defaults (`reltol=1e-8`,
`abstol=1e-10`) reduced the local Windows physics time from `64.356 s` to
`26.457 s`, a `2.432x` speedup. All 1000 samples converged in exactly three
Newton iterations. Relative to the original `1e-13` SciBmad result, the maximum
orbit change was only `1.608e-13 m`. Bmad and SciBmad still use different
mathematical stopping rules, and the existing Bmad timing was measured on a
different machine. See
`results/formal_1000/scibmad_tolerance_comparison.md`.

A local experimental solver then reused one nominal `6 x 6` closed-orbit
Jacobian for all 1000 samples and modified-Newton iterations. With the same
nominal-`z0` initial guess and Bmad-default tolerance values, the full-AD solver
used `22.247 s` while frozen Jacobian used `7.535 s`, a controlled `2.952x`
physics speedup. All samples converged; the maximum final one-turn residual
norm was `9.802e-11`, and the maximum detector-orbit difference from full AD
was `5.566e-10 m`. The runner now checks every lane's final closure and sends
only failed lanes through a full-AD Newton sub-batch. The fallback-enabled
1,000-sample repeat used `8.163 s` with zero fallbacks; a forced single-lane
test recovered that lane successfully. See
`results/formal_1000/scibmad_frozen_jacobian_comparison.md`.

Using the GTPSA-derived `6 x 119` closed-orbit control response to form an
independent initial guess for every sample reduced the frozen solver's
median/mean iterations from `3 / 2.994` to `2 / 1.995`. The recurring physics
time fell from `8.163 s` to `6.855 s` (`1.191x`, or `16.0%` less time), with
all 1000 samples converged, zero fallbacks, and maximum closure norm
`8.104e-11`. The response matrix itself took `2.389 s` to build after
compilation, so this is a net benefit when the matrix is reused for later
digital-twin batches; it is not a cold-start speedup claim. The first-process
warmup/compilation cost for the response path was `100.764 s`.

The requested Bmad, Linux SciBmad, and local SciBmad configurations are
collected in one table in
`results/formal_1000/all_benchmark_results_comparison.md`.

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

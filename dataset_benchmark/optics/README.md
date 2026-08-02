# CESR chromatic-optics benchmark

This directory benchmarks Bmad/Tao and three SciBmad evaluation strategies on
the same 1,000 CESR corrector samples and the same 99 zero-length `DET_*`
markers. RF is off, so the result is four-dimensional coasting closed-orbit
optics.

The common samples, Bmad-compatible lattice, and closed-orbit response cache
live under `../orbit/inputs` and `../orbit/reference`. Optics results remain
under this directory's own `results/` tree.

Each engine writes the zero-order periodic optics and their first derivatives
with respect to relative momentum deviation `delta`:

- `phi_1`, `beta_1`, `alpha_1`, `phi_2`, `beta_2`, `alpha_2`;
- longitudinal accumulated phase `phi_3`;
- coupled-optics `gamma_c`, `c11`, `c12`, `c21`, `c22`;
- all six propagated orbit coordinates;
- one `d..._ddelta` column for every quantity above;
- ring tunes, transverse chromaticities, and the coasting slip coefficient.

The orbit derivatives are dispersion-like phase-space derivatives with respect
to `delta`; they are not derivatives with respect to the 119 correctors.

## Current 1,000-sample results

The table reports stable physics time only. Compilation/warmup, model setup,
closed-orbit setup, and CSV writing are recorded separately in metadata. The
reference for correlations and errors is SciBmad pointwise `twiss`.

| Method | Result type | Physics time | Samples/s | Speedup vs Bmad | Speedup vs pointwise SciBmad | Maximum closure residual | Median column correlation | Minimum column correlation |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| Bmad/Tao | Exact per sample | 701.528 s | 1.425 | 1.000x | 0.814x | `1.99e-13` | 0.999995 | 0.095562 [1] |
| SciBmad pointwise `twiss` | Exact per sample | 570.771 s | 1.752 | 1.229x | 1.000x | `2.44e-12` | 1.000000 | 1.000000 |
| SciBmad prototype `twiss!` | Exact per sample | 513.273 s | 1.948 | 1.367x | 1.112x | `2.44e-12` | 1.000000 | 1.000000 |
| SciBmad one parameterized `twiss` | First-order corrector surrogate | 73.130 s | 13.674 | 9.593x | 7.805x | `1.29e-13` nominal only [2] | 0.999996 | 0.997779 |

[1] Bmad's minimum correlation is the longitudinal
`dorbit_z_ddelta` column. Its median over nonconstant columns is 0.999995, so
the isolated minimum should not be read as the accuracy of the transverse
optics as a whole. Bmad and SciBmad also use independently represented CESR
lattices.

[2] A single parameterized model solves only the nominal coasting closed
orbit. It does not perform 1,000 independent closure solves, so its residual is
not an all-sample maximum. Accuracy is instead checked directly against all
pointwise output rows.

The prototype `twiss!` output is bit-for-bit identical to pointwise SciBmad for
all compared numeric fields. For the parameterized method, the worst
correlation is 0.997779 (`xi_2`), the median is 0.999996, and the largest error
normalized by a column's maximum reference magnitude is 10.46% (`orbit_z`).
Selected ring-level errors are:

| Quantity | Correlation | Maximum absolute error |
|---|---:|---:|
| `Qx_fractional` | 0.999498 | `3.98e-5` |
| `Qy_fractional` | 0.999514 | `8.97e-5` |
| `xi_1` | 0.998645 | `1.73e-3` |
| `xi_2` | 0.997779 | `7.97e-3` |
| `slip_factor` | 0.998228 | `2.16e-4` |

The generated complete report and per-column details are in
`results/methods_1000/comparison.md` and
`results/methods_1000/comparison.csv`.

## SciBmad pointwise and reusable `twiss!`

SciBmad constructs the CESR scalar model and `Descriptor(6, 2)` once. The
existing response-linear/frozen-Jacobian RF-on batch orbit supplies initial
guesses for a short RF-off four-dimensional coasting-orbit solve. That orbit is
passed explicitly to Twiss, so the optics routine does not solve it again.

The installed SciBmad release has `twiss` but no public in-place `twiss!` API.
`twiss_reuse.jl` therefore provides a local prototype that reuses:

- detector/step lookup and ordering;
- the identity `DAMap`;
- segment `DAMap` and TPS storage.

It resets the reusable identity map and recomputes the sample-dependent maps,
normal form, and Twiss values. This keeps the exact pointwise calculation while
removing repeated allocation/setup. The 1,000-sample run reduces Twiss physics
time from 570.771 s to 513.273 s (11.2% speedup). The prototype calls internal
SciBmad stages and may need adjustment when SciBmad internals change.

Run either mode from the CESR project root:

```powershell
$env:JULIA_PKG_PRECOMPILE_AUTO='0'

julia --project=. dataset_benchmark/optics/benchmark_scibmad_chromatic_optics.jl `
  --sample-count=1000 `
  --twiss-mode=fresh `
  --output-dir=dataset_benchmark/optics/results/methods_1000/scibmad_pointwise

julia --project=. dataset_benchmark/optics/benchmark_scibmad_chromatic_optics.jl `
  --sample-count=1000 `
  --twiss-mode=reuse `
  --output-dir=dataset_benchmark/optics/results/methods_1000/scibmad_twiss_reuse
```

## One parameterized Twiss calculation

`benchmark_scibmad_parameterized_twiss.jl` creates
`Descriptor(6, 3, 119, 1)`: six phase-space variables through third order and
119 first-order GTPSA corrector parameters. It calls Twiss once at the nominal
machine, extracts the constant and 119 corrector coefficients, and evaluates
all samples with dense matrix multiplication.

Third phase-space order is necessary here. The corrector dependence of
chromaticity is carried by mixed `delta * corrector` terms and would be cut
from tune polynomials by a second-order phase-space descriptor. Corrector
dependence is still first order, making this a fast local surrogate rather
than an exact pointwise replacement.

For 1,000 samples, the single parameterized Twiss call takes 73.096 s and TPS
coefficient evaluation takes 0.035 s. Its first compilation/warmup takes
278.510 s and is excluded consistently with the other methods.

```powershell
$env:JULIA_PKG_PRECOMPILE_AUTO='0'

julia --project=. dataset_benchmark/optics/benchmark_scibmad_parameterized_twiss.jl `
  --sample-count=1000 `
  --phase-space-order=3 `
  --output-dir=dataset_benchmark/optics/results/methods_1000/scibmad_parameterized_order3
```

## Bmad/Tao

The Bmad benchmark uses one persistent Tao instance with RF off. For each
corrector sample it computes periodic optics at `pz=0`, `pz=+h`, and `pz=-h`.
Bmad's native dispersion fields are used when exposed by the installed Tao.
Beta/alpha chromatic derivatives use the same symmetric momentum difference
as the other local phase and coupling derivatives. `ring_general.chrom_a` and
`ring_general.chrom_b` provide the native ring chromaticities.

Run on lnx201:

```bash
ssh jn577@lnx201.classe.cornell.edu
cd ~/cesr_scibmad
source ~/venvs/pytao/bin/activate

python dataset_benchmark/optics/benchmark_bmad_optics.py \
  --sample-count=1000 \
  --delta-step=1e-5 \
  --output-dir=dataset_benchmark/optics/results/bmad_chromatic_fixed_1000
```

The Bmad timed region contains the 119-variable update, three Bmad optics
recalculations, and all output queries per sample. Initialization, warmup, and
file writing are excluded and separately recorded.

## Shared conventions and validation

- RF is off, so `delta` is a fixed momentum offset rather than a synchrotron
  coordinate.
- Closed-orbit closure is evaluated in `(x, px, y, py)`; `z` need not close in
  a coasting ring.
- SciBmad saves values at detector beginnings. Tao returns detector exits;
  because every `DET_*` is a zero-length marker, these are the same location.
- Phases are stored in turns. Phase origins are removed separately for every
  sample before cross-engine comparison.
- RF-off `phi_3` is accumulated longitudinal advance, not a synchrotron phase.
  The raw Bmad/Tao-derived `dphi_3/ddelta` and ring `slip_factor` use the
  opposite sign from the SciBmad output in this benchmark. The comparison
  script flips these two Bmad columns for validation only; stored raw outputs
  are unchanged.
- The comparison script rejects duplicate keys and all NaN or infinite numeric
  values before calculating metrics.

All methods write a detector table with 99,000 rows and 40 columns, a 1,000-row
ring table, and method metadata. Pointwise SciBmad also writes the 1,000 start
closed orbits.

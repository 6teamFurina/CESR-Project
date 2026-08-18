# Corrector-to-optics Jacobian benchmark

## Developer-recommended per-parameter follow-up

The recommended alternative of one `Descriptor(6,2,1,1)` parameterized Twiss
per corrector has now also been measured. Reusing one typed lattice and
descriptor, 119 serial Twiss calls plus coefficient extraction took 137.206 s;
including separately computed periodic start-orbit response columns took
202.539 s. This is respectively 9.623x and 9.352x slower than the corresponding
wide-descriptor SciBmad times, and 13.276x and 19.597x slower than Bmad's
complete 10.335 s finite-difference result.

The narrow and wide SciBmad detector Jacobians agree to `1.24e-15` relative
Frobenius difference, so this is a like-for-like performance result. The P=1
route has a much cheaper one-time warmup but does not recover that advantage
over 119 stable serial Twiss calls. See
`results/scibmad_per_corrector/RESULTS.md` for timing distributions, numerical
checks, and scope.

The follow-up with quadrupole and sextupole strengths is reported in
`results/extended/RESULTS.md`. The 225-parameter
corrector+quadrupole case completed in 19.720 s with Bmad and 43.930 s for the
most favorable SciBmad optics-only timing. In the 301-parameter
corrector+quadrupole+sextupole case, Bmad completed in 26.558 s, while
SciBmad's first parameterized Twiss did not finish within a one-hour bounded
run. Thus the present implementation becomes less competitive as this
Jacobian is widened; the benchmark does not support a near-term dataset
throughput argument based on one very wide parameterized Twiss call.

This experiment compares two ways to compute the RF-off CESR optics Jacobian
with respect to the same 119 orbit correctors:

- Bmad/Tao: 119 symmetric finite differences, requiring 238 exact periodic
  optics recalculations in one persistent Tao process.
- SciBmad: all controls are first-order GTPSA parameters. The periodic
  `6 x 119` closed-orbit response is obtained by implicit differentiation and
  supplied to one parameterized Twiss calculation.

Both methods save the full `1782 x 119` detector Jacobian (99 detectors times
18 ordinary Twiss/coupling/orbit quantities), the `6 x 119` periodic
closed-orbit response, and the ring-tune Jacobian. Compilation/warmup, model
setup, physics, coefficient extraction, and file writing are timed separately.

## Measured result

The benchmark was run sequentially on the same Ryzen 9 5900HX host. Bmad used
Tao `20260801-1` in the `Ubuntu-Bmad` WSL distribution; SciBmad used Julia
1.12.6 on Windows with one Julia thread.

| Method | Stable physics time | Breakdown | Relative time |
|---|---:|---|---:|
| Bmad/Tao symmetric finite difference | **10.335 s** | 238 optics recalculations plus all updates and queries | 1.000x |
| SciBmad/GTPSA | **21.656 s** | implicit closed orbit 7.399 s + one Twiss 14.168 s + extraction 0.089 s | 2.095x |

The tested hypothesis is therefore false for this workload and current
implementations: SciBmad is not much faster. It is `2.095x` slower than the
complete Bmad finite-difference calculation. Even if the separate implicit
closed-orbit response is excluded, the single SciBmad Twiss call alone is
`1.371x` slower than Bmad's complete 238-point calculation.

This does not mean GTPSA has no advantage. Its derivatives have no corrector
finite-difference step, and higher-order corrector coefficients can be
obtained from a higher-order descriptor. It does mean that for one dense
first-order `1782 x 119` Jacobian, the cost of propagating and normal-forming
the large TPS map exceeds the cost of 238 fast scalar Bmad evaluations.

Compilation is not hidden in the stable numbers. SciBmad's separate warmups
were 52.002 s for the implicit response and 218.410 s for Twiss; Bmad's warmup
was 0.139 s. Model/baseline setup and file writing are also excluded from the
stable physics times and are recorded in metadata.

## Numerical cross-check

The independently represented Bmad and SciBmad CESR lattices give the
following matrix agreement:

| Matrix | Relative Frobenius difference | Cosine correlation |
|---|---:|---:|
| Periodic closed orbit, `6 x 119` | 0.1041% | 0.999999968 |
| Transverse detector orbit, `396 x 119` | 0.1980% | 0.999998105 |
| Complete detector optics, `1782 x 119` | 2.0867% | 0.999887104 |
| Ring tunes, `2 x 119` | 1.1424% | 0.999938526 |

The detector `phi_3` row is convention-dependent in an RF-off coasting ring:
Bmad stores zero for the absent periodic longitudinal mode, whereas SciBmad
retains a longitudinal path/slip response. It is preserved in the raw output
and reported separately by the per-quantity comparison.

Machine-readable timings are in
`results/bmad/bmad_corrector_jacobian_metadata.json` and
`results/scibmad/scibmad_corrector_jacobian_metadata.toml`. Full and
per-quantity agreement metrics are in `results/comparison.json`.

SciBmad uses `Descriptor(6, 2, 119, 1)` by default. Second phase-space order is
needed to retain phase-space/corrector mixed terms in ordinary optics. A
chromatic-optics Jacobian would require third phase-space order and is a
different, more expensive benchmark.

Run SciBmad from the project root:

```powershell
$env:JULIA_PKG_PRECOMPILE_AUTO='0'
julia --project=. dataset_benchmark/optics/corrector_jacobian_benchmark/benchmark_scibmad_corrector_jacobian.jl
```

Run Bmad in the WSL/Linux PyTao environment from the project root:

```bash
source ~/venvs/pytao/bin/activate
python dataset_benchmark/optics/corrector_jacobian_benchmark/benchmark_bmad_corrector_jacobian.py
```

Then compare the labeled matrices:

```powershell
python dataset_benchmark/optics/corrector_jacobian_benchmark/compare_corrector_jacobians.py
```

For the completed 225-parameter extended case:

```powershell
python dataset_benchmark/optics/corrector_jacobian_benchmark/compare_extended_jacobians.py
```

For the one-parameter-per-corrector comparison:

```powershell
python dataset_benchmark/optics/corrector_jacobian_benchmark/compare_per_corrector_jacobians.py
```

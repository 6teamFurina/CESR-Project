# SciBmad one-parameter-per-corrector Twiss benchmark

This benchmark implements the SciBmad developer-recommended narrow-parameter
strategy for a first-order ordinary-optics Jacobian. A single
`Descriptor(6, 2, 1, 1)` and typed CESR lattice are reused. Each of the 119
correctors is activated as the sole GTPSA parameter, one parameterized Twiss
calculation produces one Jacobian column, and the corrector is reset before
the next column.

The run is serial with one Julia thread. Model construction and compilation
warmup are excluded from stable physics timing. The same RF-off coasting
configuration, 99 detectors, 18 detector quantities, periodic start orbit,
and ring tunes are used as in the wide-descriptor and Bmad references.

## Timing

| Method | Detector/tune optics time | Including separately exported start-orbit response | Relative to Bmad |
|---|---:|---:|---:|
| Bmad/Tao, 119 central differences | 10.335 s | included | 1.000x |
| SciBmad, one wide `Descriptor(6,2,119,1)` | 14.257 s | 21.656 s | 1.379x / 2.095x |
| SciBmad, 119 narrow `Descriptor(6,2,1,1)` calls | 137.206 s | 202.539 s | 13.276x / 19.597x |

The narrow SciBmad optics time consists of 0.001 s for corrector activation,
137.065 s for 119 Twiss calls, and 0.140 s for coefficient extraction. Mean
Twiss time is 1.152 s per corrector; the median is 1.205 s and the observed
range is 0.747--1.614 s. Separately computing all 119 periodic start-orbit
response columns takes another 65.333 s.

The P=1 warmup is substantially cheaper than the wide-descriptor warmup:
89.352 s total for response, Twiss, and extraction, compared with 270.801 s
for the wide P=119 run. That startup reduction does not offset performing 119
stable Twiss calculations. In stable serial throughput, the narrow optics-only
route is 9.623x slower than the wide descriptor and 13.276x slower than Bmad.

## Numerical equivalence

| Matrix | Narrow vs wide SciBmad relative Frobenius difference | Narrow SciBmad vs Bmad relative Frobenius difference |
|---|---:|---:|
| Detector optics, `1782 x 119` | `1.24e-15` | 2.0867% |
| Periodic start orbit, `6 x 119` | `9.66e-17` | 0.1041% |
| Ring tunes | `4.73e-17` | 1.1424% |

The two SciBmad constructions therefore produce the same first-order
Jacobian to floating-point precision. The timing difference is not caused by
a changed physical observable, corrector mapping, or derivative convention.

These results only establish serial single-thread performance on this host.
The 119 narrow calculations are independent and could be distributed across
processes, but Bmad's finite-difference columns are also independent. A
parallel comparison would be a separate throughput experiment and is not
included here.

Machine-readable timings are in
`scibmad_per_corrector_jacobian_metadata.toml`; full comparisons are in
`comparison.json`.

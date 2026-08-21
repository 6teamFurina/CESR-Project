# Latest CESR closed-orbit results

Status: `smoke` as of 2026-08-20. No production 1,000-sample timing is stored
here.

All results use SciBmad 0.4.1, Julia 1.12.6, RF on, branch 0, and
[`latest_cesr_scibmad_repaired.jl`](../../../../../Latest_Lattice/latest_cesr_scibmad_repaired.jl).
The runtime registry contains 103 selected steering controls, 144 detectors,
and 288 detector observables.

Two complementary two-sample checks are retained:

| Run | Solver/initial guess | Converged | Timed physics | Max closure |
|---|---|---:|---:|---:|
| [`smoke/`](smoke/README.md) | full-AD Newton / zero | 2/2 | 3.477 s | `1.12e-14` |
| [`smoke_fd_cache/`](smoke_fd_cache/README.md) | frozen nominal Jacobian / cached linear response | 2/2 | 3.137 s | `8.14e-11` |

The second row uses the maintained production-style solver at its default
`reltol=1e-8`, `abstol=1e-10`; the first uses the stricter
`1e-12/1e-13` full-AD path. Warmup and model setup are excluded from the timed
physics values. With only two samples, neither value is a throughput claim.

The known warning about straight multipoles in a curved reference system is
emitted by this lattice. It is the documented curved-DQX limitation and must
remain visible in any study involving the corresponding girder pitch physics.

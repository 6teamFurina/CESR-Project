# Full-Newton latest CESR smoke run

Status: `smoke`, generated 2026-08-20.

- Lattice: [`latest_cesr_scibmad_repaired.jl`](../../../../../../Latest_Lattice/latest_cesr_scibmad_repaired.jl)
- Engine: SciBmad 0.4.1; Julia 1.12.6; RF on; branch 0
- Input: two samples, 103 dynamically selected steering controls
- Output: 144 detectors and 288 ordered `x/y` observables
- Solver: full-AD Newton from zero, `reltol=1e-12`, `abstol=1e-13`
- Convergence: 2/2; iterations min/median/mean/max = 1/3.5/3.5/6
- Timings: warmup 140.169 s; model setup 4.855 s; solve 3.457 s;
  tracking 0.020 s; timed physics 3.477 s
- Maximum final closure norm: `1.1204e-14`

Artifacts:

- [`scibmad_rf_on_samples.csv`](scibmad_rf_on_samples.csv)
- [`scibmad_rf_on_metadata.toml`](scibmad_rf_on_metadata.toml)

This run validates exact BatchParam solving and detector tracking. It is not a
performance benchmark because compilation dominates and only two samples were
used. The curved-reference straight-multipole warning is the documented DQX
limitation of the latest lattice.

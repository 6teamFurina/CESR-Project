# Latest CESR GTPSA response cache

Status: `smoke`, generated 2026-08-20 with the maintained default response
backend.

The paired caches use RF-on SciBmad 0.4.1 tracking on branch 0 of
[`latest_cesr_scibmad_repaired.jl`](../../../../Latest_Lattice/latest_cesr_scibmad_repaired.jl).
They contain the runtime steering registry of 103 controls, with 144 detectors
and 288 ordered `x/y` observables.

| Cache | Shape | Method |
|---|---:|---|
| [`closed_orbit_response.csv`](closed_orbit_response.csv) | `6 x 103` | first-order GTPSA |
| [`detector_response.csv`](detector_response.csv) | `288 x 103` | first-order GTPSA |

Both metadata sidecars record response-pair id
`cdaa3f56-a96d-4198-b687-37f9120ffdd1`, `chunk_count=1`,
`controls_per_batch=0`, and `response_step_rad=0`. The nominal closed-orbit
closure norm was `3.5527e-15`.

Only the selected normal horizontal/vertical steering controls are promoted
to GTPSA parameters. Unselected controls remain primitive `Float64` zeros.
This is necessary for the repaired lattice because promoting the zero
quadrupole component of `SEX_14W` would ask the tracking kernel to evaluate a
TPS `sqrt(0)` before taking its zero-field branch. It is an adapter constraint,
not a failure of the nominal closed orbit or of GTPSA for the steering
response.

The sidecars are:

- [`closed_orbit_response.csv.metadata.toml`](closed_orbit_response.csv.metadata.toml)
- [`detector_response.csv.metadata.toml`](detector_response.csv.metadata.toml)

They preserve ordered controls and observables, shapes, lattice/ring identity,
RF state, tolerances, closure diagnostic, response method, pair identity, and
software version. The pair must be consumed together.

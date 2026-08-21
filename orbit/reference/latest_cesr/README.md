# Latest CESR SciBmad response caches

Status: `smoke`, generated 2026-08-20.

The caches use RF-on SciBmad 0.4.1 tracking on branch 0 of
[`latest_cesr_scibmad_repaired.jl`](../../../Latest_Lattice/latest_cesr_scibmad_repaired.jl).
They were generated from the ordered runtime steering registry: 103 controls
(58 horizontal and 45 vertical), 144 detectors, and 288 `x/y` observables.

| Location | Shapes | Method |
|---|---|---|
| [`gtpsa/`](gtpsa/README.md) | `6 x 103`; `288 x 103` | maintained GTPSA default |
| this directory | `6 x 103`; `288 x 103` | FD validation artifact |

The root-level FD step is `1e-7 rad`, split into 13 chunks of at most 8 controls (positive
and negative lanes are solved together). All lanes converged at
`reltol=1e-12`, `abstol=1e-13`; the maximum final closure norm was
`9.9529e-14`.

The sidecars
[`closed_orbit_response.csv.metadata.toml`](closed_orbit_response.csv.metadata.toml)
and
[`detector_response.csv.metadata.toml`](detector_response.csv.metadata.toml)
store ordered labels, shape, lattice/ring identity, RF state, numerical method,
step, chunking, solver tolerances, closure diagnostic, and software version.
The two sidecars share a response-pair id so readers cannot accidentally mix
matrices from different generations. Cache readers reject label, shape, pair,
ring, lattice, engine, branch, RF, numerical-step, tolerance, or version
mismatches when the metadata is available.

The central-difference method records how the root-level validation cache was
generated; the maintained runner default is `--response-method=gtpsa`, whose
new paired smoke cache is under [`gtpsa/`](gtpsa/README.md). Use the root files
as labeled smoke/validation artifacts and explicitly select central difference
when reproducing them. An earlier many-parameter GTPSA
diagnostic stopped at `SEX_14W` because the adapter TPS-ified unselected
controls, causing the lower-level `sqrt(0)` domain error. With
`zero_value=0.0` and only selected GTPSA controls represented as parameters,
the full 1,177-element map succeeds. The reported relative-L2 differences
against these caches are `2.78e-8` for `6 x 103` and `1.68e-8` for `288 x 103`.
The GTPSA and central-difference files are both SciBmad artifacts and do not use
Bmad/Tao.

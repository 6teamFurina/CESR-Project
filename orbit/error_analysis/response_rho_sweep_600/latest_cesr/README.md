# Latest CESR response-radius results

Status: `smoke`. Minimal runs have completed for the maintained
[`GTPSA path`](gtpsa_smoke/README.md) and the
[`central-difference validation path`](smoke/README.md); the production
600-trial grid has not been run.

The result uses the latest repaired SciBmad lattice, the runtime 103-control
steering registry, 144 detectors/288 observables, and the ring-scoped
first-order cache in [`../../../reference/latest_cesr/`](../../../reference/latest_cesr/README.md).
No historical 119-control/99-detector artifact is mixed into this directory.

The checked-in smoke outputs are intentionally separated by response backend:
`gtpsa_smoke/` is the maintained GTPSA result and `smoke/` is the independent
central-difference validation. Each metadata sidecar records the shared
`response_pair_id`; chunk merges and figures must preserve that provenance.
The no-argument renderer discovers only these latest-ring locations and rejects
unscoped legacy summaries.

The post-audit smoke rerun is in [`smoke_recheck/`](smoke_recheck/README.md).
It reused the GTPSA pair, solved all 7 states with no fallback, and reached a
maximum closure norm of `8.8058e-14`; it is still an integration check, not the
600-trial production result.

The same 7-state case was also executed through the paper-level unified entry
point in [`smoke_runner/`](smoke_runner/README.md), including result-local SVG
rendering. This checks orchestration and paths in addition to the underlying
rho solver.

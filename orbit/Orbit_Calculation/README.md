# Closed-orbit calculation: latest CESR ring

Status: `smoke`. The latest-ring adapter, runtime registries, exact batch solve,
detector tracking, bounded response-cache generation, and cache-loaded frozen-
Jacobian path have been exercised. No production 1,000-sample timing has been
run on this ring.

The default SciBmad model is
[`Latest_Lattice/latest_cesr_scibmad_repaired.jl`](../../Latest_Lattice/latest_cesr_scibmad_repaired.jl).
Use the `latest_cesr` artifact paths:

```text
inputs/latest_cesr/
../reference/latest_cesr/
results/latest_cesr/
```

At startup, the runner should discover and persist:

- the complete control registry and any selected steering groups;
- detector/observable names, order, plane metadata, and units;
- the closed-orbit coordinate dimension and response-matrix shape;
- RF mode, branch, lattice source, seed, tolerances, and software versions.

The primary result is a SciBmad closed-orbit dataset. Bmad/Tao may be run only
as an explicitly labeled matched comparison using the same discovered input
and output registries. Report model setup, warmup, physics, I/O, convergence,
and closure diagnostics separately.

The current lattice exposes 124 writable Overlay/Group controls. The default
steering selector uses the dynamically identified 58 `HKICK` and 45 `VKICK`
controls, rather than a fixed column count. It also discovers 144 detectors
and 288 `x/y` outputs. Custom ring configurations may provide their own
control subset, detector selector, coordinate labels, and ring id.

Here, ring-generic means that registry sizes, ordering, selectors, and artifact
paths are discovered rather than fixed. The maintained solver contract is
still a six-coordinate SciBmad `Bunch`, and the default response controls are
a homogeneous angular-steering subset measured in radians. A future mixed-unit
control registry (for example kicks plus RF or magnet-strength knobs) will need
per-control perturbation scales rather than the single `response-step-rad`.

The checked-in first-order validation response was generated with SciBmad
BatchParam central differences in bounded control chunks. The maintained
runner default is now `--response-method=gtpsa`; the FD cache is a labeled
smoke/validation artifact, not the current default. Select central difference
explicitly, or recompute when response metadata does not match the requested
method. The recorded FD settings are `h=1e-7 rad` and eight controls per
chunk, configurable with `--response-step-rad` and
`--response-controls-per-batch`. This is a cache-generation choice, not a
requirement to replace GTPSA. An earlier many-parameter GTPSA diagnostic
stopped at `SEX_14W` because the adapter TPS-ified unselected controls, causing
the lower-level `sqrt(0)` domain error. With `zero_value=0.0` and only selected
GTPSA controls represented as parameters, the full 1,177-element map succeeds;
the reported relative-L2 differences against the central-difference cache are
`2.78e-8` for `6 x 103` and `1.68e-8` for `288 x 103`. The current
central-difference cache remains explicitly a validation artifact and should
not be mistaken for a GTPSA-default production result.

Current smoke artifacts and exact timings are recorded in
[`results/latest_cesr/README.md`](results/latest_cesr/README.md). The canonical
input is [`inputs/latest_cesr/corrector_samples.csv`](inputs/latest_cesr/corrector_samples.csv);
its ordered 103-column registry is validated against the selected model.

The latest result directory contains per-run READMEs with status and
provenance. The recorded two-sample runs are compilation/integration smoke
tests and must not be quoted as sustained-throughput measurements.

The old fixed-dimension description and benchmark commands are preserved in
[`README_archived.md`](README_archived.md). Existing historical result reports
remain in `results/` and are not latest-ring evidence.

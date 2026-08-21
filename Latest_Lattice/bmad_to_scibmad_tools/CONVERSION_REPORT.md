# Latest CESR Bmad-to-SciBmad conversion audit

Audit date: 2026-08-16

## Version and source status

- Source lattice: `CHESS-U_6000MEV_20251020_S1`, 6 GeV positrons.
- Installed validated environment: Bmad `20260801.1`, Tao `20260801-1`,
  PyTao `1.2.1`.
- Latest conda-forge Bmad found and tested in the isolated environment
  `/home/joeyfurina/.cache/codex-bmad-20260814`: Bmad `20260814.0`, Tao
  `20260814-0`, PyTao `1.2.2`.
- Current project packages: SciBmad `0.4.1`, Beamlines `0.9.2`.
- The upstream repository heads inspected during the audit reported SciBmad
  `0.4.2` and Beamlines `0.10.1`; neither source tree defined `Fork`, `Mirror`,
  or `Wiggler` element constructors.

The latest raw output is
`../bmad_reference/raw_exports/latest_cesr_scibmad_bmad_20260814.jl`. It is the
baseline for all further repairs. The older `20260801-1` output is retained
only to show that the latest Bmad writer materially improved control export.

## What the latest writer preserved

The Bmad model has 12 branches and 1,397 total branch elements. Main branch 0
has 1,331 elements, including tracking elements and lord/control elements.

Tao `20260814-0` emitted 248 SciBmad `DefExpr` assignments. These preserve the
nominal and deferred relationships for the 119 Overlay lords and five Group
lords, including horizontal/vertical correctors, skew controls, sextupole tune
knobs, and beta knobs. The Bmad inventory contains 347 lord-to-slave control
relations; they are recorded in
`../bmad_reference/inventory/bmad_control_relations.csv`. Representative
one-slave, multi-slave, skew, split-superlord, and Group relationships were
changed dynamically after construction. The exhaustive tracking-response audit
in `FULL_CONTROL_VALIDATION.md` covers all 124 lords, all 347 relations, and
475 observation points. Its median per-control maximum relative discrepancy is
`3.74e-4`; the worst informative response is `1.84%` with cosine `0.999834`,
consistent with small finite-amplitude thick-element model differences rather
than a missing relation.

This is a significant improvement over Tao `20260801-1`, which emitted 88
`HKICK` and 45 `VKICK` translation warnings and did not preserve those
controls.

## Writer defects and applied repairs

### 1. Wiggler: three tracking slices lose their magnetic field

The Bmad super-lord `ID_S1A` is a nontrivial planar wiggler:

- total length: `2.355 m`;
- `B_MAX = 1.17 T`;
- `L_PERIOD = 0.19625 m`;
- `N_PERIOD = 12`;
- `K1Y = -1.7087583e-3 m^-2` in the Bmad model.

Superposition with markers/forks splits it into three tracking pieces of
length `0.8846 m`, `0.2930 m`, and `1.1774 m`. The latest writer reports all
three as untranslatable and emits plain `LineElement(L=...)` objects. The
length survives, but the oscillatory field, vertical focusing, path-length
effect, and nonlinear map do not.

The repaired file instantiates three phase-continuous segments from the
validated `PlanarWiggler` implementation in `../../wigglers/wiggler.jl`, with
90/30/120 integration steps. A `6.54758669568665e-15 s` reference-time patch
on the existing internal fork marker restores Bmad's design-particle `z=0`
convention. The five-element block then agrees in affine vector and exit orbit
to `7.12e-15`.

The continuous-field block retains an `R12` increment of `5.89e-6 m` relative
to Bmad's `bmad_standard` wiggler matrix. This is a physical-model distinction,
not an integration error: the value is unchanged between 240 and 960 SciBmad
steps, and Bmad itself produces approximately the same increment when its
wiggler is changed to Runge-Kutta field tracking. Bmad's standard approximation
omits the extra horizontal quiver-path contribution.

### 2. DQX integration steps: twelve strong combined-function bends

Tao exports the correct DQX geometry and `Kn1`, but omits Bmad's
`NUM_STEPS=100`. SciBmad's one-step default caused a maximum local matrix error
of `1.43`. The repaired file assigns `Yoshida(order=6, n_steps=100)` to all 12
DQX bends, reducing their maximum local matrix error to `4.29e-13`.

### 3. Photon forks: eleven invalid and forward-referenced constructors

The raw file emits 11 `Fork(to_line=...)` objects. Every target Beamline is
defined later in the same file, so the first include fails immediately with:

```text
UndefVarError: s4b_line not defined in Main
```

More fundamentally, neither the current project packages nor the inspected
upstream SciBmad/Beamlines source defines a `Fork` constructor. Nine of these
objects are main-ring Bmad `Photon_Fork` elements; the other two fork the
S1 photon branch again.

The repair converts all 11 to zero-length markers in the Beamlines objects,
preserving the positron main-ring map. It also exports
`latest_photon_fork_targets`, `latest_photon_branch_for_fork`, and a
branch-local paraxial ray helper. All eleven branch element counts and lengths
match Bmad exactly. Branch spawning remains explicit rather than automatic.

### 4. Photon mirrors: two unsupported constructors

The writer emits `mirror_s7 = Mirror()` and `mirror_s1a = Mirror()`, but no
such current SciBmad/Beamlines constructor exists. The repair represents them
as zero-length branch-local reference markers and records the Bmad metadata:
`REF_TILT=-pi/2`, `GRAZE_ANGLE=0.004 rad`, and 10 keV reference energy. Since
the reflected reference frame is already encoded in each Bmad photon branch,
the reference ray stays on axis and accumulates the exact branch length.
Reflectivity, finite apertures, curvature, and off-axis specular scattering
are not implemented.

### 5. Girders: 12 control/alignment lords are omitted

The Bmad lattice defines `GIRDER_1AB` through `GIRDER_6CD`. Tao reports each
girder as not translatable. Their start/end markers remain, but the coherent
alignment semantics of each girder do not.

`export_bmad_girder_coefficients.py` expands the Bmad geometry response onto
all 150 member tracking elements/slices.
`../essential_supports/latest_girder_support.jl` exposes
`set_latest_girder!`, which applies any of the six Bmad offset/pitch/tilt
parameters coherently and can reset them to zero. The full audit covers 12
girders x 6 parameters and 972 tracking-response observations.

Translations agree closely. The worst pitch response differs by `21.7%`, at
strong DQX combined-function bends. This is a current library limitation:
SciBmad warns that its straight multipoles in a curved reference system do not
satisfy the exact free-space Maxwell geometry. The nominal DQX linear maps are
still matched to `4.29e-13` after restoring 100 integration steps, but exact
misaligned curved-coordinate multipole parity is not available in current
SciBmad/BeamTracking.

## Validated charged-particle status

The raw export remains preserved and intentionally fails to load. The generated
`../latest_cesr_scibmad_repaired.jl` loads as a 1,177-element, 768.437869 m,
6 GeV positron Beamline.

The branch-0 validation against Tao/Bmad `20260814-0` found:

- exact element count, order, and length for all 1,177 tracking elements;
- outside the wiggler block, maximum local matrix mismatch `4.29e-13`;
- complete wiggler-block affine/exit-orbit mismatch `7.12e-15`;
- RF-on starting closed-orbit mismatch `5.25e-14`;
- one-turn eigenphase tune mismatches about `5.9e-14`, `2.4e-9`, and `5.5e-8`;
- all 124 Overlay/Group controls and all 347 control relations exercised by
  one-pass tracking response, with median maximum relative discrepancy
  `3.74e-4` and worst informative discrepancy `1.84%`;
- all 12 girders exposed through a coherent alignment interface, with the
  curved-DQX pitch limitation quantified above;
- all 11 photon branches registered with exact element counts and lengths and
  a validated branch-local reference-ray path.

These results support nominal numerical parity and complete exported-control
coverage for the charged-particle main ring. The conversion is suitable for
ordinary CESR optics, tracking, corrector/knob studies, coherent girder offsets,
and photon branch geometry. It does not claim exact DQX response under girder
pitch or general photon mirror optics; those require capabilities absent from
the current SciBmad/BeamTracking stack.

# DQX girder-pitch response diagnosis

The largest girder discrepancy is the local `y_pitch` response at the exit of
a strong DQX combined-function sector bend. For `GIRDER_4AB -> DQX4B`, the
central-difference derivative with respect to girder pitch is:

| Model | `dy/d(y_pitch)` (m/rad) | `dpy/d(y_pitch)` (rad/rad) |
|---|---:|---:|
| Bmad | 1.664830640663 | -0.488111169262 |
| SciBmad | 1.348687765076 | -0.282345248822 |

The six-vector relative L2 difference is `0.217422241256`. For a `1 urad`
girder pitch, this corresponds locally to differences of about `0.316 um` in
vertical position and `0.206 urad` in vertical canonical momentum/angle.

## Isolation tests

- Setting the DQX quadrupole component to zero in both codes reduces the
  relative difference to `3.0321e-4` and the maximum derivative-entry
  difference to `1.0943e-3`. Thus the Bmad-derived girder geometry mapping and
  the pure sector-bend response agree closely.
- Changing the Bmad element from `Bmad_Standard` to `Runge_Kutta` leaves the
  full-field discrepancy at `0.217429421057`.
- Increasing the SciBmad Yoshida integration from 100 to 200, 400, and 800
  steps leaves `dy`, `dpy`, and the discrepancy unchanged at the displayed
  precision.

The discrepancy is therefore specifically associated with the quadrupole
field of a pitched curved combined-function bend, rather than the finite-
difference step, insufficient integration steps, the Bmad standard map, or an
incorrect coherent-girder geometry coefficient. Current BeamTracking applies
straight multipole kicks in a curved reference system and emits a warning that
this field does not satisfy the exact free-space Maxwell geometry. Under pitch,
the curvature/multipole/alignment cross terms do not reproduce Bmad's field-
coordinate transformation.

The reproducible diagnostic is in `diagnose_girder_dqx_cause.py` and
`diagnose_girder_dqx_cause.jl`; its Bmad reference is
`bmad_girder_reference/dqx4b_pitch_isolation.json`.

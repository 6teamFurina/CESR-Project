# CESR planar-wiggler zero-input experiment

Date: 2026-08-11

## Question

For the continuous-field planar wiggler in `wiggler.jl`, does a particle with
zero six-dimensional canonical input return to zero at the exit? Are the two
wigglers active in the default SciBmad CESR lattice?

## Reproduction

From `CESR Project`:

```powershell
julia --project=. wigglers\experiment_zero_orbit.jl
```

The script constructs isolated experimental lattices at the CESR reference
momentum of `5.2889999753148 GeV/c`. The nominal wiggler has `B_max = 1.17 T`,
period `0.19625 m`, 12 periods, total length `2.355 m`, sixth-order Yoshida
integration, and 16 steps per period.

## Results

The coordinate order is `(x, px, y, py, z, delta)`.

| Case | x (m) | px | y | py | z (m) | delta |
|---|---:|---:|---:|---:|---:|---:|
| Nominal symmetric wiggler | `-5.2965e-18` | `0` | `0` | `0` | `-2.526140873e-6` | `0` |
| Field off (`B_max=0`) | `0` | `0` | `0` | `0` | `0` | `0` |
| Phase `pi/2`, canonical-zero input | `-1.5704e-17` | `0` | `0` | `0` | `-2.526140873e-6` | `0` |
| Truncated 11.5-period field, phase 0 | `1.293967037e-4` | `0` | `0` | `0` | `-2.420885004e-6` | `0` |
| Nominal field with radiation damping | `-1.163246e-9` | `-4.8789e-19` | `0` | `0` | `-2.526168243e-6` | `-1.078757e-5` |

Splitting the nominal continuous field into 48 sequential quarter-period
elements saved 49 boundary states. The maximum internal horizontal
displacement was `1.293967037e-4 m` (`129.397 micrometres`), while the exit
horizontal displacement was `4.06e-19 m`. The canonical `px` remains zero in
this gauge because the horizontal motion is carried by the vector potential;
the corresponding small-amplitude mechanical-angle scale is
`2.0714e-3 rad`.

The usual analytic oscillation amplitude about the oscillation centre is
`64.6983 micrometres`. Starting at `x=0` at a field maximum makes the internal
trajectory range approximately twice that amount away from its entrance
coordinate, consistent with the tracked `129.397 micrometres` maximum.

The paraxial path-length estimate is

```text
L * angle_amplitude^2 / 4 = 2.526134753e-6 m.
```

The tracked longitudinal magnitude is `2.526140873e-6 m`, a ratio of
`1.000002423` to the estimate. The sign is negative because of the SciBmad
longitudinal-coordinate convention.

Integration convergence was:

| Steps per period | Maximum transverse exit residual | z (m) |
|---:|---:|---:|
| 2 | `5.150e-19` | `-2.561734835e-6` |
| 4 | `5.943e-18` | `-2.526140845e-6` |
| 8 | `9.656e-20` | `-2.526140873e-6` |
| 16 | `5.296e-18` | `-2.526140873e-6` |
| 32 | `8.649e-18` | `-2.526140873e-6` |

The transverse closure is set by symmetry to floating-point precision. The
longitudinal term converges by eight steps per period; the production default
uses 16.

With the first four canonical coordinates set to zero, additional longitudinal
input controls gave:

| Input `(z, delta)` | Maximum output magnitude among `(x,px,y,py)` | Output `(z, delta)` |
|---|---:|---:|
| `(1 mm, 0)` | `5.2965e-18` | `(9.974738591e-4 m, 0)` |
| `(0, +1%)` | `5.0551e-18` | `(-2.476149226e-6 m, +1%)` |
| `(0, -1%)` | `5.7937e-18` | `(-2.577655062e-6 m, -1%)` |
| `(1 mm, +1%)` | `5.0551e-18` | `(9.975238508e-4 m, +1%)` |

Thus the ideal magnetostatic, full-period, radiation-off model keeps the first
four coordinates closed even when `z` or `delta` is nonzero. This is a symmetry
property of this ideal element, not a general guarantee for an arbitrary
wiggler, an asymmetric end field, or radiation-on tracking.

## Default CESR lattice state

`load_cesr()` contains exactly two elements of kind `Wiggler`, `wig_w` and
`wig_e`. Both have:

- `B_max = 1.17 T`;
- period `0.19625 m`;
- phase `-37.6991118431 rad = -12 pi`;
- a nonzero field, so both are active by default;
- radiation damping and radiation fluctuations disabled in their tracking
  method by default.

Thus, "wiggler active" and "radiation tracking active" are separate settings.
In the default CESR lattice the wiggler magnetic fields, focusing,
nonlinearities, internal oscillation, and path-length effect are present, but
stochastic radiation and deterministic radiation energy loss are not tracked.

## Conclusion

For the ideal full 12-period CESR wiggler with radiation disabled, a canonical
zero input returns to transverse zero at the exit to numerical precision, but
does not return to six-dimensional zero: `z` changes by about
`-2.52614 micrometres`. The particle is not on the axis internally. A deliberately
incomplete 11.5-period field leaves a `129.397 micrometre` horizontal exit
offset even though its exit canonical momentum is zero, demonstrating why
both angular and positional closure matter.

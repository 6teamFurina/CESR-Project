# Paired quadrupole-offset orbit correction

This latest-lattice SciBmad experiment keeps every maintained static
random error fixed within each paired machine and switches only the
quadrupole x/y alignment offsets from zero to 50 micrometers RMS per plane.
Correction uses BPM observations, corrector commands, and an orbit-response
matrix; latent quadrupole offsets are never supplied to the solver.

## Aggregate result

| response matrix | before BPM RMS [um] | after BPM RMS [um] | reduction | before target 2D RMS [um] | after target 2D RMS [um] | command RMS | max abs command |
|---|---:|---:|---:|---:|---:|---:|---:|
| reference_orm | 860.881808 | 126.034037 | 6.831x | 1335.740202 | 70.709229 | 2.03637e-5 | 0.000137143 |
| current_orm | 860.881808 | 126.003055 | 6.832x | 1335.740202 | 70.724249 | 2.04948e-5 | 0.00013707 |

The reference-ORM row uses the paired zero-quadrupole-offset response
matrix, representing a stored response measurement from the aligned state.
The current-ORM row remeasures the response after the offsets are applied.
Both matrices include the fixed corrector and BPM gains of that machine.

Across paired machines, the current-versus-reference ORM relative-L2 difference has median `0.0199634` and maximum `0.0327474`. The two correction methods' final BPM vectors differ by `1.306819` micrometers RMS.

## Sextupole scan-range restoration

The fraction of beam-relative sextupole centers outside the configured 1.5-millimeter scan radius is `0.0%` in the zero-offset reference and `28.618%` before correction.

- `reference_orm`: `0.0%` outside after correction.
- `current_orm`: `0.0%` outside after correction.

## Method boundary

- The reference is the BPM orbit of the paired machine with quadrupole
  alignment drift disabled; it is not a zero orbit and no oracle target-local
  coordinates enter correction.
- BPM readings are noise-free in this bounded test. Stable 1% RMS BPM gains
  and 1% RMS corrector gains are included in both the observations and the
  measured response matrices.
- Corrector commands are the 103 latest-lattice HKICK/VKICK Overlay
  steering variables in radians. No CESR
  power-supply or operator limit is asserted.
- Matching finite BPM readings does not prove exact trajectory equality
  between BPMs; the target-sextupole orbit residual is reported separately.
- The latest lattice retains its documented straight-multipole-in-curved-
  reference qualification. Girder pitch is not varied in this experiment.

## Reproduction

```powershell
julia --project=. sextupole_misalignment/quadrupole_orbit_correction/run_quadrupole_orbit_correction.jl
python sextupole_misalignment/quadrupole_orbit_correction/validate_orbit_correction.py
```

Machines: `16`; quadrupoles: `113`; BPMs: `111`; correctors: `103`.

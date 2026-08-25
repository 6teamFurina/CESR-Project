# All-sextupole K2/bump orbit cross-response

This study determines whether changing the bump and K2 at one active normal
sextupole can be treated as local when the orbit is observed at all 76 active
normal sextupoles.  New calculations use the repaired latest SciBmad lattice
at `Latest_Lattice/latest_cesr_scibmad_repaired.jl` with RF on.

## Matrices

The raw response bundle contains:

- `periodic_kick_response.npy`: orbit at every sextupole from a unit local
  `px` or `py` kick at every sextupole;
- `bump_response.npy`: orbit at every sextupole from each target's model-based
  x/y corrector bump knob;
- `sextupole_source_response.npy`: the periodic response to the integrated
  normal and skew source polynomials of each target sextupole;
- `alignment_design.npy`: the selected K2--bump--center derivative, with axes
  `target, bump_axis, observation_sextupole, output_plane, center_axis`.

The alignment design maps a target beam-relative center to the K2-odd,
bump-odd orbit gradient at every sextupole.  It is the propagation object
needed to decide whether a first learned inverse may use target-local,
neighbor-sparse, low-rank global, or fully global orbit inputs.

## Result

The propagation is not target-local.  For the model-based corrector bumps,
the median fraction of radial-response energy at the excited sextupole is
`6.024%` for x bumps and `3.383%` for y bumps.  Across the four
K2--bump--center channels, the median target-only fraction is only
`0.097%`--`0.377%`, with a median participation count of
`31.52`--`33.94` observation sextupoles.

The conclusion is also large in absolute rather than only normalized terms.
For a 0.5 mm target command, the median off-target radial-orbit RMS is
0.228 mm for x bumps and 0.309 mm for y bumps.  The largest individual
off-target responses anywhere in the two matrices are 1.871 mm and 2.708 mm,
respectively.  Thus the present corrector-defined bump is a target-constrained
closed-orbit knob, not a spatially local orbit deformation.

The global bump response and the family of alignment templates nevertheless
have compact shared spatial structure.  The `152 x 152` bump matrix has
effective rank `6.251` and needs 6 modes for 90% of its squared singular-value
energy.  Placing all target/center templates in the same 304-channel
coordinate system gives a `304 x 152` shared-template matrix with effective
rank `22.451` and 20 modes for 90% energy.

That shared-template SVD is not the joint inverse.  The physical design for 76
separate one-target-at-a-time scans retains the target axis and is block
diagonal, with 23,104 rows and 152 center columns.  Its singular values can be
assembled from 76 small `304 x 2` blocks without allocating the large matrix:
it has numerical rank 152, effective rank `143.416`, needs 128 modes for 90%
energy, and has unwhitened condition number `2.060`.  Every target block has
rank two and condition number one to numerical precision.  Thus a shared
full-ring response-mode projection is worth testing, but its truncation must
be selected by covariance-whitened center-recovery performance while
preserving target identity; a 20-mode energy cutoff is not by itself adequate
evidence for recovering multiple misalignments.

A paired exact finite-amplitude SciBmad check covers five selected targets,
signed 0.5 mm x/y bumps, `delta K2 = +/-0.02 m^-3`, and an added
`(+350,-250) um` target offset.  The aligned odd/odd gradient is subtracted
from the matched misaligned gradient.  The compact prediction has a `3.267%`
aggregate relative L2 residual, `0.999475` cosine similarity, and a fitted
scale factor of `1.00421`; the largest target/bump-axis block residual is
`5.781%`.  This validates the sign and scale for the selected cases, not for
an all-magnet misaligned background.

The resulting first-input proposal, parity-contrast definition, required
target-axis convention, and three-way representation ablation are specified in
[`FIRST_MODEL_INPUT.md`](FIRST_MODEL_INPUT.md).

## Descriptor policy

A global descriptor with all sextupole K2, bump, and offset parameters would
generate many irrelevant cross-target derivatives.  This study instead uses:

1. `Descriptor(6, 1)` for cumulative and one-turn maps;
2. `Descriptor(6, 2, 62, 1)` once for first derivatives with respect to the
   current 62 independent correctors;
3. the exact local normal-sextupole polynomial to compose only the required
   K2--bump--center derivative.

No global Hessian is formed.  This factorization is the nominal first-order
periodic response.  It uses an integrated thin source at the target entry;
the selected exact finite-amplitude scans above quantify, but do not remove,
the finite-length and relinearization qualification.

The latest lattice emits the documented straight-multipole-in-curved-reference
warning.  No girder pitch is varied here, so the known curved-DQX girder-pitch
discrepancy is not an excitation in this study, but it remains part of the
lattice provenance.

The orbit values at all sextupole entries are SciBmad model-state observations.
They are not directly measured at CESR.  A machine-facing input must use BPM
readbacks or a finite-BPM sextupole-orbit inference with propagated covariance.

## Reproduction

From `CESR Project/`:

```powershell
julia --project=. sextupole_misalignment/sextupole_cross_response/generate_gtpsa_cross_response.jl

wsl.exe -d Ubuntu-Bmad -- /home/joeyfurina/miniforge3/envs/bmad/bin/python `
  '/mnt/d/Ring_Design_Development/CESR Project/sextupole_misalignment/sextupole_cross_response/analyze_cross_response.py'

julia --project=. sextupole_misalignment/sextupole_cross_response/generate_exact_validation.jl

wsl.exe -d Ubuntu-Bmad -- /home/joeyfurina/miniforge3/envs/bmad/bin/python `
  '/mnt/d/Ring_Design_Development/CESR Project/sextupole_misalignment/sextupole_cross_response/analyze_exact_validation.py'

wsl.exe -d Ubuntu-Bmad -- /home/joeyfurina/miniforge3/envs/bmad/bin/python `
  '/mnt/d/Ring_Design_Development/CESR Project/sextupole_misalignment/sextupole_cross_response/validate_cross_response.py'
```

The generated numerical summary is written to `results/analysis/SUMMARY.md`.
The paired finite-amplitude check is written to
`results/exact_validation_analysis/SUMMARY.md`.

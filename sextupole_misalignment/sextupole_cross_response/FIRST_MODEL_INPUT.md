# First model-input recommendation

## Decision

Keep the excitation-target axis and the full-ring observation axis in the
reference representation.  Do not begin with one target-local orbit value or
with a target-summed global vector.  The present corrector bumps and the
alignment signal both propagate around the ring, while the separate-target
joint inverse derives its rank from retaining which sextupole was excited.

## Atomic example

Use one target sextupole in one latent machine realization as the first atomic
example.  For each x/y bump axis, form the K2-odd, bump-odd orbit gradient at
all 76 sextupole locations and both output planes:

```text
G = [O(+b,+k) - O(+b,-k) - O(-b,+k) + O(-b,-k)] / (4 b k)
```

The reference signal tensor therefore has shape
`bump_axis x observation_sextupole x output_plane = 2 x 76 x 2`, or 304
features per target before masks and context are added.  Preserve the raw scan
states as the primary dataset artifact; `G` is a derived view, not a
replacement for the measurements.

Orbit at every sextupole is an oracle model-state representation, not a direct
CESR measurement.  A machine-facing version must either use the available BPM
readbacks directly or propagate the finite-BPM local-orbit posterior and its
covariance to the sextupole locations.  Do not label inferred sextupole orbit
as measured orbit.

Add the following context rather than asking a model to infer it implicitly:

- target identity and longitudinal position;
- commanded and, when modeled, physical/readback bump and K2 amplitudes;
- zero-excitation orbit and channel-validity mask;
- observation uncertainty or whitening operator;
- lattice/optics identifier and any nuisance settings exposed to the model.

The first output is the two-component beam-relative magnetic center of the
excited target.  For multiple-center recovery, retain a set or tensor of these
target-indexed examples.  With separate target scans, the nominal physical
design is block diagonal with 152 center columns; summing target templates
before inference would destroy that structure.

## Representation ablation

Train and evaluate the same estimator class and the same latent machines with
three matched input views:

1. full 304-channel target-indexed parity tensor;
2. covariance-whitened projection onto shared full-ring response modes;
3. target-only and fixed-neighbor subsets as a locality ablation.

The shared unwhitened template family needs 20 modes for 90% response energy
and 68 for 99%, but these are only starting points for an ablation.  Select the
retained dimension by held-out center error and uncertainty calibration after
whitening, not by response energy alone.  The separate-scan joint design needs
128 modes for 90% energy and 150 for 99%, demonstrating why the target axis
cannot be silently folded into the shared-mode truncation.

## First baselines and splits

Use a covariance-whitened linear/MAP inverse as the required physics baseline.
A small residual MLP, set model, or graph model may then learn deviations from
that baseline; a large unrestricted network is not the first diagnostic.

Evaluate on held-out target sextupoles, latent error combinations, error
magnitudes, noise/drift conditions, missing observations, and changed optics.
Include BPM drift/gain, corrector gain, K2 calibration, quadrupole strength,
roll, and misalignment in staged nuisance ablations.  The selected exact
validation in this folder checks nominal factorization sign and scale only; it
does not replace those distribution-shift tests.

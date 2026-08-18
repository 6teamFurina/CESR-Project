# Sextupole--quadrupole affinity results

Full-ring calculation completed 2026-08-14 for 76 target sextupoles and 106
active quadrupoles. The optics pre-screen retained exactly 15 quadrupoles per
target, so the expensive nuisance-marginalized analysis contains 1,140 pairs.
All 106 quadrupoles remain as heatmap columns; 32 appear in at least one target's
retained set and the other columns are blank.

## Requested metrics

Across the retained pairs, nuisance-marginalized log-determinant information
gain ranges from `0.9630` to `4.8463`, with median `1.3681`. The maximum is
`SEX_33E / Q00E`.

Worst-axis center-precision improvement ranges from `1.1167x` to `2.3970x`,
with median `1.6580x`. The maximum is `SEX_34E / Q00W`.

The metrics select materially different quadrupoles:

- Information winner counts: `Q00W` for 37 targets, `Q00E` for 36, `Q01E`
  for 2, and `Q01W` for 1.
- Precision winner counts: `Q01E` for 42 targets, `Q01W` for 32, `Q00E` for
  1, and `Q00W` for 1.
- Only 4 of 76 targets have the same winner under both definitions.

This split is physically useful rather than a plotting artifact. Log
determinant rewards improvement in the area of the two-dimensional center
error ellipse, while the requested precision ratio protects the worse of the
coordinate-axis uncertainties. If the operational objective is a robust
two-plane center estimate, the precision winner is the more direct primary
choice; the information winner should remain in the next exact-validation set
because it may add complementary two-dimensional information.

The simple beta/phase pre-screen is useful but not sufficient by itself. Its
rank-1 quadrupole is the final information winner for only 7 targets and the
final precision winner for 34 targets; median winning ranks are 3 and 2,
respectively. The nuisance-marginalized scoring stage should therefore not be
replaced by raw optics proximity.

## Independent numerical check

The nominal Bmad finite-difference target responses were compared with the
previously saved SciBmad/GTPSA mixed `K2-offset` derivatives for all 76 targets
and all 10 retained observable families. Across 760 family comparisons, the
minimum cosine similarity is `0.997884`. The largest family-relative L2
difference is `0.06570` for the weak `SEX_41E` vertical-orbit response; its
cosine similarity remains `0.997884`. The full check is saved in
`nominal_bmad_vs_scibmad.csv`.

## Interpretation limits

The values use provisional diagonal measurement noise and independent `300 um`
RMS priors for the other 75 sextupole offsets. Their 150 target-`K2` cross-
response columns are explicitly Schur-marginalized, but are evaluated at
nominal quadrupole optics and reused for the `K1 +/-` blocks. The information
gain also compares three measurement blocks (`nominal`, `q+`, and `q-`) with
one nominal block, exactly as requested; its absolute value therefore includes
the extra measurement count, although target-by-target quadrupole ranking is
still meaningful.

No orbit-bump replication is included in this screen. The natural next stage
is to carry, for every sextupole, the best-information and best-precision
quadrupoles (plus any desired runner-up) into the exact `3 x 3` bump,
five-point `K2`, and three-point `K1` protocol. These heatmaps are a candidate
design result, not a CESR position-precision claim.

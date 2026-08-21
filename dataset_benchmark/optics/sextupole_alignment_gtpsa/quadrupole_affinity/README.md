# Sextupole--quadrupole affinity screen

This study ranks which quadrupole changes are most useful for estimating the
two-dimensional magnetic center of each active normal sextupole. The maintained
calculation now uses the repaired CHESS-U 6 GeV SciBmad lattice in
`Latest_Lattice/latest_cesr_scibmad_repaired.jl`.

The new lattice contains 76 active nonzero-`Kn2` sextupoles and 113 independent
active quadrupole knobs. There are 114 nonzero-`Kn1` tracking elements because
`Q49` is split into two tracking slices; both slices are changed together and
appear as one physical quadrupole in the screen. Measurements use the 111
`DET_*` markers whose label identifies a physical BPM. Auxiliary `DET_*`
markers are excluded.

## Calculation

The inexpensive first pass changes each independent active quadrupole by

```text
Kn1 = Kn1_0 +/- 0.1% |Kn1_0|
```

and evaluates at every target sextupole

```text
sqrt((delta log beta_x)^2 + (delta log beta_y)^2
     + (delta phi_x)^2 + (delta phi_y)^2).
```

The phase terms are fixed-BPM difference phases and are expressed in radians
inside the screening norm. Tune changes above `0.01` or BPM beta beating above
`20%` are rejected. The 15 largest allowed optics-leverage values are retained
for every target. Unselected pairs remain blank in the final heatmaps.

For every target, batched SciBmad/GTPSA maps calculate

```text
A_s = d^2 O / (d Kn2_s d [x_s, y_s])
```

at nominal optics and at the retained quadrupole's positive and negative
conditions. The measurement vector has 1110 directly measurable entries. At
each of 111 BPMs it contains the horizontal and vertical trajectory for a fixed
reference launch, followed by the two trajectory components produced by four
small launch probes:

```text
x = 1 mm, px = 0.1 mrad, y = 1 mm, py = 0.1 mrad.
```

These eight differential-trajectory readings are the measured transfer-map
columns, expressed as BPM displacements rather than abstract matrix entries.
Same-plane columns carry beta/phase information and cross-plane columns carry
coupling information. This representation is more direct experimentally and
avoids differentiating the gauge-dependent Sagan--Rubin coupled-Twiss
decomposition at exactly zero coupling.

The nominal GTPSA calculation obtains every sextupole's own mixed `Kn2-offset`
response. For a target sextupole, the two columns from each of the other 75
sextupoles form 150 explicit nuisance directions. The saved coefficients are
the requested mixed Hessian entries, not finite differences.

The analysis assumes provisional independent `5 um` BPM noise for each fixed or
probed trajectory reading. It is converted to slope noise for the five-point scan
`delta Kn2 = [-2, -1, 0, 1, 2] x 0.01 m^-3`. Other-sextupole offsets receive
independent `300 um` RMS Gaussian priors.

Some sextupoles carry zero-valued quadrupole-control fields in the lattice.
When their `Kn2` becomes a TPS parameter, SciBmadStandard's scalar branch would
form `sqrt(0 TPS)` before falling back to DriftKick. The generator selects that
same DriftKick branch explicitly for these elements; the scalar nominal map is
unchanged while the parameterized map remains analytic.

For target-center columns `A`, nuisance columns `N`, and whitened weights
`C^-1/2`, retained center information is the Schur complement

```text
F_center = F_cc - F_cn (F_nn + Sigma_nuisance^-1)^-1 F_nc.
```

The two reported affinities are

```text
information gain = log det(F_[nominal,q+,q-]) - log det(F_nominal)

precision improvement = max(sigma_x,0, sigma_y,0)
                        / max(sigma_x,q, sigma_y,q).
```

## Reproduction

Use the pinned Julia environment under `CESR Project`. Run the three response
stages as separate Julia processes so each stage owns one GTPSA descriptor:

```powershell
$env:JULIA_PKG_PRECOMPILE_AUTO='0'

julia --startup-file=no --project='CESR Project' `
  'CESR Project/dataset_benchmark/optics/sextupole_alignment_gtpsa/quadrupole_affinity/generate_scibmad_affinity_responses.jl' `
  --stage=screen

julia --startup-file=no --project='CESR Project' `
  'CESR Project/dataset_benchmark/optics/sextupole_alignment_gtpsa/quadrupole_affinity/generate_scibmad_affinity_responses.jl' `
  --stage=nominal

julia --startup-file=no --project='CESR Project' `
  'CESR Project/dataset_benchmark/optics/sextupole_alignment_gtpsa/quadrupole_affinity/generate_scibmad_affinity_responses.jl' `
  --stage=candidates
```

Then calculate the two affinities and render the heatmaps:

```powershell
python 'CESR Project/dataset_benchmark/optics/sextupole_alignment_gtpsa/quadrupole_affinity/analyze_affinity.py'
python 'CESR Project/dataset_benchmark/optics/sextupole_alignment_gtpsa/quadrupole_affinity/render_affinity_heatmaps.py'
python 'CESR Project/dataset_benchmark/optics/sextupole_alignment_gtpsa/quadrupole_affinity/select_greedy_quadrupole_sets.py'
python 'CESR Project/dataset_benchmark/optics/sextupole_alignment_gtpsa/quadrupole_affinity/validate_greedy_quadrupole_selection.py'
python 'CESR Project/dataset_benchmark/optics/sextupole_alignment_gtpsa/quadrupole_affinity/select_best_quadrupole_triplets.py'
python 'CESR Project/dataset_benchmark/optics/sextupole_alignment_gtpsa/quadrupole_affinity/validate_quadrupole_triplet_selection.py'
```

Response bundles are resumable by target and candidate condition. Results are
kept under `results/scibmad_latest`, separate from the archived Bmad/Tao
prototype. Both static and interactive heatmaps use target sextupole on the
x-axis and quadrupole on the y-axis.

## Repaired-lattice result

The completed run screened all `76 x 113 = 8588` pairs and evaluated the 15
retained candidates per target (`1140` pairs). The union contains 47
quadrupoles, so only those 47 appear as heatmap rows; non-retained pair cells
are blank.

Across the retained pairs, nuisance-marginalized information gain spans
`1.39795--2.14188` (median `1.85276`) and worst-axis precision improvement
spans `1.45573--1.71227` (median `1.60334`). `Q10W` is the best-information
choice for 38 targets and the best-precision choice for 35. Both metrics choose
the same quadrupole for 62 of 76 targets. The scalar optics screen's rank-one
candidate wins the final information metric for only 9 targets and the
precision metric for only 10, so the response-level calculation remains
necessary.

All 76 nominal bundles, all 2280 candidate-sign responses, and all 1140 score
rows pass the shape and finiteness validator. Explicit DriftKick selection for
the zero-quadrupole-control sextupoles was also checked against the scalar
SciBmadStandard map: both the saved trajectory and Jacobian differences were
exactly zero.

## Greedy five-candidate selection

The next selection stage uses the 15 response-level candidates already saved
for each target. Starting from the nominal block, it adds one quadrupole at a
time. Every trial stacks nominal plus the symmetric `K1 +/-` blocks of all
quadrupoles selected so far, recomputes the other-sextupole Schur
marginalization, and maximizes cumulative log-determinant information. The
worst-axis precision improvement is retained as the secondary diagnostic and
deterministic tie-breaker.

Five complementary quadrupoles are retained per target. The first three are
marked as the provisional operational set; candidates four and five remain for
the exact finite-amplitude and observable-ablation study. The repaired-lattice
run produced 380 selection rows. The three-candidate sets use 25 distinct
quadrupoles across all targets, while the full five-candidate sets use 34.

At three candidates, the median cumulative information gain is `3.12216` and
the median worst-axis precision improvement is `2.20396`. At five candidates,
the corresponding medians are `3.72965` and `2.57317`. The first three retain a
median `83.54%` of the five-candidate log-determinant gain; this quantifies the
cost of reducing the eventual protocol from five knobs to three, but it is not
yet an exact bump-grid performance result.

The detailed ranked rows, one-line-per-target sets, metadata, and summary are
under `results/scibmad_latest/selection`. Step one agrees with the original
single-quadrupole information winner for all 76 targets, providing a closure
check on the greedy implementation.

## Exhaustive five-to-three selection

For every target, the retained five candidates define 11 available optics
conditions: nominal and one-at-a-time `K1 +/-` for each candidate. All ten
three-of-five subsets are evaluated. A triplet uses nominal plus its six
relevant candidate-sign blocks, with the nuisance Schur complement recomputed
for that seven-condition subset. The maximum-information triplet is the primary
selection; the maximum-precision triplet and both rankings are saved as
diagnostics.

For all 76 targets, exhaustive subset selection reproduces the first three
greedy choices exactly. The selected triplets use 25 distinct quadrupoles. The
median information gain is `3.12216`, and the median worst-axis precision
improvement is `2.20396`. Information- and precision-optimal triplets differ for
32 targets, so the precision alternative remains important in the subsequent
finite-amplitude inverse test.

The 760 evaluated triplets are in `triplet_combinations.csv`; the selected
information triplet and best-precision alternative for every target are in
`best_triplets_by_sextupole.csv`. This closes the response-level 11-condition
selection but does not replace the planned exact bump-by-`K2` validation or the
combined three-knob `3^3` grid.

## Interpretation boundary

The affinity heatmap is a nominal-orbit, single-quadrupole design screen. The
greedy follow-on selects a complementary multi-quadrupole set from those saved
responses, but neither stage yet replicates the planned `3 x 3` bump grid or
uses an archived CESR measurement covariance. The other-sextupole nuisance
dictionary is evaluated at nominal quadrupole optics and reused in the
candidate `+/-` blocks. These are candidates for the next exact protocol
experiment, not real-machine position-precision claims.

The greedy set calculation advances the multi-quadrupole design but retains the
same nominal-launch dictionary and reused nominal nuisance responses. Its first
three candidates therefore remain provisional until the exact `3 x 3` bump,
five-point `K2`, finite-`K1`, direct-observable protocol is evaluated.

## Exact 11-condition audit

The subsequent repaired-lattice audit generated all `76 x 495 = 37,620`
closed-orbit states and evaluated all 760 three-of-five triplets. It did not
validate a unique triplet at `+/-0.1% |K1|`: the response selection matched the
exact closed-orbit information winner for 24 targets and the precision winner
for 22, but changing to the exact winner improved deterministic center error by
at most `6.64e-5 um`. The ten triplets are effectively tied.

More importantly, a fair seven-block comparison against seven repeated nominal
blocks gives a median best-triplet precision factor of only `1.0000295`. The
earlier `~2.2` factor compared seven measurements with one and therefore mostly
measured repetition, not K1-generated optics diversity. The five candidates
remain useful as a starting pool, but the present first-three ranking must not
be frozen or used to launch the combined `3^3` scan. The next design iteration
must use a larger safe K1 excursion or stronger optics knobs, a matched-budget
objective, and direct launch-trajectory phase/coupling observables before adding
the nuisance ensemble. Full protocol, results, and reproduction commands are in
`exact_11_triplet_validation/README.md`.

The follow-up `+/-1% |K1|` SciBmad audit finds 110/113 knobs within the
maintained tune and beta-beating limits; `Q10W`, `Q43E`, and `Q10E` exceed the
`0.01` tune-shift limit. After replacing those knobs from the retained
15-candidate pools, all 76 targets were rerun. The symmetric K1 even/odd
response ratio has median `0.00829`, P90 `0.02180`, and maximum `0.05046`, so
1% is adopted for the audited-safe knobs in the next observable study.
Nevertheless, the exact closed-orbit-only matched-budget precision factor is
only `1.000505` at the median even for the exact-optimal triplet. This does not
restore a unique three-quadrupole selection; direct launch-trajectory
phase/coupling remains the next required ablation.

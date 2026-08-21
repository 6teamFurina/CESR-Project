# Greedy quadrupole-set selection

This repaired-lattice SciBmad result selects five complementary quadrupoles for
each of the 76 active sextupoles. The first three are the provisional
operational set; the remaining two are retained for the finite-amplitude and
observable-ablation stage.

The selector greedily maximizes cumulative nuisance-marginalized
log-determinant information. At each step it stacks the nominal response and
the `K1 +/-` response blocks for every selected quadrupole, then recomputes the
Schur complement for the 150 other-sextupole offset nuisance directions.
Worst-axis precision remains an explicit diagnostic and tie-breaker.

## Result

- Targets: 76.
- Retained selections: 380, exactly five unique quadrupoles per target.
- Provisional operational selections: 228, exactly three per target.
- Distinct quadrupoles used by the three-candidate sets: 25.
- Distinct quadrupoles used by the five-candidate sets: 34.
- Median three-candidate information gain: 3.12216.
- Median five-candidate information gain: 3.72965.
- Median fraction of five-candidate information retained by the first three:
  83.54% (range 81.67%--85.48%).
- Median worst-axis precision improvement: 2.20396 at three candidates and
  2.57317 at five.
- The first greedy choice reproduces the original best single-quadrupole
  information choice for all 76 targets.

## Exhaustive five-to-three check

The five retained candidates provide 11 response-level optics conditions:
nominal and `K1 +/-` for each candidate. All ten possible three-candidate
subsets were evaluated for every sextupole. Each subset uses nominal plus its
six relevant candidate-sign blocks and recomputes the nuisance Schur
complement.

- Evaluated triplets: 760.
- Targets whose exhaustive information-optimal triplet differs from the first
  three greedy candidates: 0 of 76.
- Information-optimal versus precision-optimal triplet disagreements: 32 of
  76.
- Selected information gain: 2.23389--3.75650, median 3.12216.
- Selected worst-axis precision improvement: 1.80735--2.57149, median 2.20396.
- Distinct quadrupoles in the selected triplets: 25.

The maximum-information triplet is the provisional primary selection. The
maximum-precision alternative is retained because the exact finite-amplitude
offset-recovery test may prefer it for some of the 32 disagreement targets.

The complete one-line-per-target list is in
`quadrupole_sets_by_sextupole.csv`. The step-by-step selections and cumulative
metrics are in `greedy_quadrupole_selection.csv`; machine-readable aggregate
statistics and assumptions are in `selection_summary.json` and
`selection_metadata.json`.

The exhaustive triplet table is in `triplet_combinations.csv`; the selected
triplet and precision alternative per target are in
`best_triplets_by_sextupole.csv`. Their aggregate assumptions and statistics
are in `triplet_selection_metadata.json` and
`triplet_selection_summary.json`.

## Boundary

This selection uses the existing nominal-launch SciBmad/GTPSA response bundles.
The other-sextupole nuisance dictionary is still evaluated at nominal optics
and reused in the candidate blocks. The first three choices are therefore
provisional until the `3 x 3` orbit-bump, five-point `K2`, finite-`K1`, and
direct-observable protocol is run.

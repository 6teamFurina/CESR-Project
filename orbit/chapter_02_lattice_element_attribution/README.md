# Chapter 2: complete lattice-element attribution

Status: `production_complete` as of 2026-08-22.

This chapter applies the paper's exact complete-element Hessian-source method
to the validated latest CESR SciBmad lattice. For each fixed paired horizontal
and vertical steering direction, it forms the leading second-order nonlinear
detector target and reconstructs it from every complete element boundary,

```text
g_j = S_exit,j - A_j S_entrance,j.
```

The production ensemble uses 100 directions for each detector output plane.
Signed family projections are additive; family magnitude ratios are not
additive because propagated element vectors interfere.

Run after Chapter 1:

```powershell
.\run_production.ps1
```

Outputs are written directly below this chapter:

```text
results/latest_cesr/horizontal/  detector-x raw and analyzed tables
results/latest_cesr/vertical/    detector-y raw and analyzed tables
tables/                          paired family-attribution CSV/LaTeX tables
figures/                         complete-element and normal-sextupole SVG/PDF panels
RESULTS.md                       generated result summary
```

This is a SciBmad result. Bmad/Tao data are not substituted for any production
quantity. The latest-lattice curved-DQX girder-pitch limitation must remain
visible if these results are later interpreted as a girder-pitch study.

The completed horizontal and vertical 100-direction ensembles reconstruct the
complete-element targets with relative closures `1.31e-14` in both planes.
The 76 active normal sextupoles supply signed projections of `96.508%` in x
and `98.647%` in y. See [`RESULTS.md`](RESULTS.md) for the complete provenance,
interpretation boundary, and paper artifacts.

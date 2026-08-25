# Chapter 1: nonlinear closed-orbit response-radius sweep

Status: `production_complete` as of 2026-08-22.

This chapter reproduces the paper's nonlinear RF-on closed-orbit response-
error method on the validated latest CESR SciBmad lattice. Gaussian steering
directions are normalized to unit RMS over the active controls, reused at
every normalized radius `rho`, and evaluated for the all-, horizontal-, and
vertical-control subspaces. The exact nonlinear detector orbit is compared
with the nominal first-order GTPSA response.

The default production run uses 600 directions at the 20 positive radii used
by the paper workflow, plus one shared zero-input reference. Its expected
state count is therefore `1 + 3 x 20 x 600 = 36,001`.

Run from this directory or from any PowerShell working directory:

```powershell
.\run_production.ps1
```

The runner writes directly below this chapter:

```text
results/latest_cesr/chunks/   independently calculated rho groups
results/latest_cesr/merged/   unified production tables and metadata
tables/                       paper-width CSV and LaTeX tables
figures/                      stacked SVG and PDF paper figures
RESULTS.md                    generated result summary
```

All quantitative outputs must identify SciBmad, RF-on branch 0, the latest
repaired lattice, the dynamic control/detector registries, the response method,
the direction seed, solver tolerances, convergence, fallback, and closure.
The known curved-DQX girder-pitch limitation remains applicable if girder
pitch physics is interpreted; this response study does not vary girder pitch.

The completed run contains 36,001 unique exact states. It converged 36,000,
with one failed all-control direction at `rho=51.2`; 62 of 63 fallback lanes
were recovered. The maximum closure among converged states was approximately
`1.0e-13`. See [`RESULTS.md`](RESULTS.md) for the generated summary and links
to the paper artifacts.

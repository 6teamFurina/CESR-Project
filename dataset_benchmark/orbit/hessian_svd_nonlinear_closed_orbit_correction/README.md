# Hessian-SVD Nonlinear Closed-Orbit Correction in a CESR Model Based on SciBmad

## Status

Working JACoW-style scaffold for the follow-on orbit-correction paper. The
abstract, introduction, problem definition, Hessian-SVD formulation, correction
algorithm/experiments, and conclusion are explicit placeholders. They must not
be read as completed scientific results.

The only substantive result currently transferred into this paper is the
vertical-only signed-parity study. It is positioned as a model-order and
trust-region boundary for future Hessian-based correction experiments.

## Main files

- `hessian_svd_nonlinear_closed_orbit_correction.tex` - LaTeX source.
- `hessian_svd_nonlinear_closed_orbit_correction.pdf` - compiled working paper.
- `jacow.cls` - local JACoW class copied from the companion calculation paper.

## Figure assets

The source SVGs remain independent files and the PDFs are their LaTeX-ready
counterparts:

- `figures/vertical_parity_crossover_loglog.svg`
- `figures/vertical_parity_crossover_loglog.pdf`
- `figures/vertical_parity_odd_fraction_percentiles.svg`
- `figures/vertical_parity_odd_fraction_percentiles.pdf`
- `figures/vertical_parity_growth_linear.svg`
- `figures/vertical_parity_growth_linear.pdf`

The linear-growth pair is retained as an unused diagnostic and is not included
in the current LaTeX scaffold.

## Established evidence included

The transferred experiment used 100 fixed vertical-only directions, 27 radii,
and 5401 states including the nominal state. All states converged without
fallback. The fitted vertical even and odd residual trends crossed at
`rho = 1.31`; at `rho = 6.4`, the mean vertical odd/even ratio was `4.91`,
while the horizontal control ratio was `3.05e-4`.

These observations establish a channel-dependent quadratic-model boundary.
They do **not** yet demonstrate Hessian-SVD correction performance, select a
correction algorithm, or identify the physical source of the cubic term.

## Build

From this directory, run:

```powershell
& 'C:\texlive\2024\bin\windows\latexmk.exe' -lualatex -interaction=nonstopmode -halt-on-error hessian_svd_nonlinear_closed_orbit_correction.tex
```

The checked-in figure PDFs make SVG conversion unnecessary for a normal build.

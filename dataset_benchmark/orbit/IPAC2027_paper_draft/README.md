# IPAC'27 orbit-response paper draft

This folder contains a working JACoW-style LaTeX draft based only on the
currently completed CESR orbit and optics benchmarks, the orbit response-radius
sweep, and the signed vertical-parity experiment.

## Files

- `ipac2027_orbit_response_draft.tex`: main manuscript source.
- `jacow.cls`: official JACoW class v3.01, downloaded from the maintained
  JACoW template repository.
- `figures/*.svg`: exact copies of the three maintained project figures.
- `figures/*.pdf`: LaTeX-ready conversions of those SVG sources.
- `svg_to_pdf.py`: local deterministic converter used for the three simple
  project SVGs.

## Intentional placeholders

The following sections deliberately contain no claimed result yet:

1. Abstract;
2. Introduction;
3. the four-sign mixed-Hessian attribution of the all-corrector second-order
   vertical residual;
4. the sextupole/octupole attribution of the vertical-only cubic response.

The placeholder boxes should be removed only after the corresponding text or
experiment is complete.

## Build

From this directory:

```console
python svg_to_pdf.py
latexmk -lualatex ipac2027_orbit_response_draft.tex
```

The current source compiles successfully to a three-page PDF with JACoW class
v3.01. The IPAC'27 contributed-paper limit is three pages plus an optional
references-only page. This is a scientific working draft; final page fitting,
author metadata, paper code, abstract, and JACoW Cat Scan validation remain to
be done after the missing experiments are inserted.

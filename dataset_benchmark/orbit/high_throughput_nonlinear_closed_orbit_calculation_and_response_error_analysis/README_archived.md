# High-Throughput Nonlinear Closed-Orbit Calculation and Response-Error Analysis

This folder contains a calculation-centered JACoW-style LaTeX draft based on
the matched nonlinear-rho SciBmad--Bmad benchmark, the 600-direction orbit
response-radius sweep, the exact complete-element Hessian-source
decomposition, the normal-sextupole physical source--beta--phase predictor,
and the archived internal-exposure attribution.
The paper positions these calculations as numerical and sampling foundations
for a future CESR orbit dataset; projected-Hessian orbit correction is reserved
for separate work.

## Files

- `high_throughput_nonlinear_closed_orbit_calculation_and_response_error_analysis.tex`:
  main manuscript source, named after the paper title.
- `jacow.cls`: official JACoW class v3.01, downloaded from the maintained
  JACoW template repository.
- `../error_analysis/response_rho_sweep_600/figures/`: the maintained raw
  horizontal/vertical nonlinear-response plot used at the top of page 2.  The
  $\rho^2$-normalized companion is retained as a source diagnostic but is not
  included in the manuscript.
- `figures/`: PDF build assets for the general response sweep and the paired
  element-attribution figure.  The maintained source SVGs remain under
  `../error_analysis/thick_element_sextupole_sourcing/horizontal_results/` and
  `../error_analysis/thick_element_sextupole_sourcing/vertical_results/`, with
  the input-exposure source under
  `../error_analysis/quadratic_x_attribution/element_results/`.
- `svg_to_pdf.py`: local deterministic converter for the maintained SVG
  sources used by the manuscript.

## Remaining draft metadata

The scientific narrative, Abstract, and Introduction are populated.  The
detailed signed-parity/cubic analysis has been transferred conceptually to the
companion Hessian-SVD orbit-correction paper; this manuscript retains only one
brief observation in the general amplitude-dependence subsection.  The new
second-order narrative proceeds from general amplitude dependence to exact
horizontal and vertical complete-element source attribution, then to the
normal-sextupole source--beta--phase predictor.  The removed corrector-space
Hessian and horizontal-imbalance discussion is preserved in the
`quadratic_x_attribution` README.  Temporary
author name, corresponding-author email, acknowledgements, paper code, and
conference metadata must still be replaced before circulation or submission.
The complete-element source sums close the adopted horizontal and vertical
quadratic detector targets to numerical precision.  The complete single-column
family table and paired side-by-side horizontal/vertical normal-sextupole
contribution panels are output-side attributions; the unsigned exposure study
is retained as a distinct input-side diagnostic outside the manuscript.
Recomputed
lattice ablations and third-order source attribution remain outside the
completion requirement for this calculation paper.

## Build

From this directory:

```console
python svg_to_pdf.py
latexmk -lualatex high_throughput_nonlinear_closed_orbit_calculation_and_response_error_analysis.tex
```

The restored source presently exceeds the final page limit while a
content-preserving compression strategy is developed.  The IPAC'27
contributed-paper limit is
three pages plus an optional references-only page; compression and figure
selection will be performed after the full draft is assembled.  Final author
metadata, paper code, acknowledgements, conference checks, and JACoW Cat Scan
validation also remain.

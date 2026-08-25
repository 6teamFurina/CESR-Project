# Latest-CESR calculation-paper manuscript

This folder is a format-preserving latest-lattice revision of the JACoW-style
working paper in
`../high_throughput_nonlinear_closed_orbit_calculation_and_response_error_analysis/`.
The historical manuscript remains unchanged.

## Primary artifact

- `high_throughput_nonlinear_closed_orbit_calculation_and_response_error_analysis.tex`
  retains the original JACoW class, two-column structure, section hierarchy,
  figure placements, table styles, bibliography, and temporary author block.
- The compiled PDF has the same basename.

Build status: LuaLaTeX compilation succeeds in four pages (three body pages
plus the dedicated references page), with no overfull boxes or unresolved
citations/cross-references. Every rendered page was visually inspected.

## Latest-CESR evidence used

All numerical conclusions and plotted data in this revision are SciBmad
production results generated from
`../../Latest_Lattice/latest_cesr_scibmad_repaired.jl` with RF on and branch 0.

- Chapter 1: `../chapter_01_nonlinear_response_rho_sweep/`
  - 36,001 unique nonlinear states over 20 radii and three steering subspaces.
  - 36,000 converged states, 63 recorded fallback attempts, and 62 recoveries.
  - 103 controls (58 horizontal, 45 vertical), 144 detectors, and 288 x/y
    observables.
- Chapter 2: `../chapter_02_lattice_element_attribution/`
  - 100 directions per output plane, 1,177 complete elements, 10 runtime
    families, and 76 active normal sextupoles.
  - All-element relative closure of `1.308639e-14` in x and `1.311835e-14`
    in y.
- Matched latest-ring timing benchmark:
  - 9,001/9,001 shared states returned in both implementations.
  - Reusable-model 16-thread SciBmad required 17.989 s for physics and
    22.842 s including setup.
  - The repeated 16-thread Bmad/Tao reference required 136.797 s for the
    9,000 nonzero states used in the paper table.

The `data/` directory contains the compact CSV tables used to audit the
manuscript numbers. The `figures/` directory contains PDF and SVG copies of
the production figures under the original manuscript figure basenames.

## Deliberately deferred claims

The historical Bmad/Tao timing and cross-code response numbers are not reused;
the manuscript uses the identical-input latest-ring timing benchmark instead.
The historical source--beta--phase correlation coefficients are also omitted;
the latest-ring direction-resolved family percentiles replace that table.

This charged-particle study does not vary girder pitch or photon branches, so
the documented curved-DQX pitch and general mirror-optics limitations are not
exercised by the reported calculations. The results are synthetic numerical
evidence and are not a real-machine precision or operating-limit claim.

# Extended optics-Jacobian benchmark

The extended experiment retains the same RF-off CESR data product as the
119-corrector benchmark: 99 detectors times 18 ordinary optics quantities,
the periodic start-orbit response, and ring tunes. Physical controls are
grouped by CESR base name, including superposition slices:

- 119 corrector kicks;
- 106 active quadrupole strengths;
- 76 active sextupole strengths.

Bmad/Tao uses symmetric scalar finite differences. SciBmad assigns all
strengths first-order GTPSA parameter variables in one
`Descriptor(6, 2, P, 1)` model. Corrector, quadrupole, and sextupole steps in
the Bmad reference are respectively `1e-6 rad`, `1e-6 m^-2`, and `1e-4 m^-3`.

## Timing result

Stable physics timings were measured sequentially on the same Ryzen 9 5900HX
host. Bmad used Tao `20260801-1` in WSL; SciBmad used Julia 1.12.6 on Windows
with one Julia thread.

| Parameters | Bmad symmetric finite difference | SciBmad optics only | SciBmad including separate start-orbit response | Outcome |
|---|---:|---:|---:|---|
| 119 correctors | 10.335 s (238 solves) | 14.257 s (1.379x Bmad) | 21.656 s (2.095x Bmad) | completed |
| + 106 quadrupoles = 225 | 19.720 s (450 solves) | 43.930 s (2.228x Bmad) | 61.399 s (3.114x Bmad) | completed |
| + 76 sextupoles = 301 | 26.558 s (602 solves) | not available | not available | first parameterized Twiss exceeded 1 h |

“SciBmad optics only” is the parameterized Twiss call plus coefficient
extraction. It excludes the separately computed periodic start-orbit response
to give SciBmad the most favorable fair comparison for the detector/tune
Jacobian. Bmad's timed region includes all lattice updates, periodic optics
recalculations, and output queries.

Bmad remains almost linear at about `0.088 s` per parameter (one `+/-` pair).
From 119 to 225 parameters, its time grows by `1.908x`, close to the parameter
ratio. SciBmad's optics-only time grows by `3.081x`; consequently its
disadvantage widens from `1.379x` to `2.228x` even before charging the separate
closed-orbit response.

For 225 parameters, the stable SciBmad Twiss call allocated 106.25 GB
cumulatively and took 43.833 s; the separate implicit response allocated
52.80 GB cumulatively and took 17.469 s. These allocation totals are not peak
resident memory. First-call warmups were 779.494 s for Twiss and 47.996 s for
the response.

For 301 parameters, the final bounded run completed model setup in 1.339 s and
the response warmup in 62.339 s, then did not finish its first parameterized
Twiss before the 3600 s process limit. The process was terminated at the
limit, and no SciBmad matrix or stable timing is claimed. The corresponding
Bmad result completed in 26.558 s.

## Numerical cross-check for 225 parameters

The 225-parameter matrices use identical labeled columns in both engines.
Agreement is reported separately by family so that the larger quadrupole
responses do not hide corrector behavior.

| Matrix | Corrector relative Frobenius difference | Quadrupole relative Frobenius difference | Quadrupole cosine correlation |
|---|---:|---:|---:|
| Complete detector optics | 2.0867% | 1.0623% | 0.999967594 |
| Periodic start orbit | 0.1041% | 0.0444% | 0.999999901 |
| Ring tunes | 1.1424% | 0.5173% | 0.999990321 |

These are code-to-code comparisons between independently represented CESR
lattices, not real-machine validation. The close labeled agreement shows that
the timing result is not explained by a gross quadrupole mapping error.

Machine-readable details are in each case's `bmad` and `scibmad` directories.
The 225-parameter family- and quantity-resolved comparison is in
`correctors_quads/comparison.json`; the bounded 301-parameter failure record is
in `correctors_quads_sextupoles/scibmad/scibmad_extended_jacobian_timeout.toml`.

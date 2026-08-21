# Direct-observable nuisance ablation

This paired SciBmad pilot uses 1 complete scan tensors for `SEX_09AW`.
Every tensor contains independent unknown offsets on the other 75 sextupoles
(Gaussian RMS 300 µm per plane) and
independent, fixed quadrupole strength errors uniformly bounded by
±1.0%. The inverse never receives either
nuisance truth. Target local-orbit coordinates are treated as exact.

The direct readbacks are generated from measurement processes: fixed
launch/angle BPM trajectory differences (the raw inputs to phase, beta and
cross-plane coupling reconstruction), TBT tune spectra,
a fixed beam-energy ±delta probe (dispersion orbit difference and tune shift), and
finite differences of two actual correctors (ORM). The energy probe uses
`delta=0.001` rather than retuning the harmon-master RF frequency.

| model | 2D RMSE [µm] | median [µm] | P90 [µm] | paired wins / 1 |
|---|---:|---:|---:|---:|
| orbit only | 3.830 | 3.830 | 3.830 | — |
| orbit + feed-down direct readbacks | 6.538 | 6.538 | 6.538 | 0 |
| orbit + all direct except chromaticity | 15.828 | 15.828 | 15.828 | 0 |
| orbit + all direct | 60.597 | 60.597 | 60.597 | 0 |

These are noise-free, structurally block-normalized physics fits, not predicted
machine precision. Each observable block receives equal total weight because a
measured covariance has not yet been supplied. Chromaticity is shown separately
because a centered sextupole has an intrinsic chromatic response, so its K2 slope
does not obey the same zero-at-center relation as ordinary linear feed-down.

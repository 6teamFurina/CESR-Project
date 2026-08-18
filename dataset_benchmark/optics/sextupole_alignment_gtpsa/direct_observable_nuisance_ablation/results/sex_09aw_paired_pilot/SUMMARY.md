# Direct-observable nuisance ablation

This paired SciBmad pilot uses 8 complete scan tensors for `SEX_09AW`.
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

| model | 2D RMSE [µm] | median [µm] | P90 [µm] | paired wins / 8 |
|---|---:|---:|---:|---:|
| orbit only | 2.829 | 2.930 | 3.879 | — |
| orbit + feed-down direct readbacks | 6.923 | 5.203 | 9.802 | 2 |
| orbit + all direct except chromaticity | 10.298 | 9.071 | 15.877 | 2 |
| orbit + all direct | 185.943 | 62.342 | 264.801 | 0 |

These are noise-free, structurally block-normalized physics fits, not predicted
machine precision. Each observable block receives equal total weight because a
measured covariance has not yet been supplied. Chromaticity is shown separately
because a centered sextupole has an intrinsic chromatic response, so its K2 slope
does not obey the same zero-at-center relation as ordinary linear feed-down.

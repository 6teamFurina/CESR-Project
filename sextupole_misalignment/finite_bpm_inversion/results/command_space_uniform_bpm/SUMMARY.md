# Command-space finite-BPM baseline

This experiment fits the beam-relative sextupole center using only selected
BPM closed-orbit channels and the five known bump commands. Exact internal
target orbit is excluded from every fit and from BPM selection. It is loaded
only afterward to evaluate the relative truth `c_s - z_s0`.

- targets / latent realizations: 76 / 8 per target
- total fits: 4864
- K2 protocol: three outer points
- BPM selection: deterministic ring-uniform subsets
- BPM noise/offset/gain errors: none
- internal target orbit used by inverse: **no**
- absolute mechanical offset estimated: **no**
- 111-BPM beam-relative 2D RMSE:
  **13.913 micrometers**
- nominal-command versus actual local-bump displacement 2D
  RMS / median / P90: **22.899 / 13.026 / 37.964 micrometers**

| retained BPMs | RMSE [micrometers] | median [micrometers] | P90 [micrometers] |
|---:|---:|---:|---:|
| 1 | 25.273 | 9.019 | 24.747 |
| 2 | 14.787 | 8.345 | 23.379 |
| 4 | 14.601 | 8.699 | 23.871 |
| 8 | 14.052 | 8.300 | 22.772 |
| 16 | 14.108 | 8.496 | 22.435 |
| 32 | 14.013 | 8.380 | 22.162 |
| 64 | 14.001 | 8.335 | 22.579 |
| 111 | 13.913 | 8.276 | 22.566 |

The result measures command-space recovery of the beam-centering displacement.
It does not yet include reconstruction of the nominal local orbit required to
recover an absolute mechanical offset. The nearly flat result from 8 through
111 BPMs, together with the local-bump mapping diagnostic above, makes bump
calibration/model mismatch the next variable to test before optimizing BPM
placement. See `summary.csv` for the complete BPM-count ablation.

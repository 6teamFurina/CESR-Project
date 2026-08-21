# Interleaved and repeated BPM acquisition protocol

This paired study tests whether repeated per-point averaging and interleaved
`0,+,0,-,0` K2 acquisition rescue the finite-BPM sextupole-center inverse under
uncorrelated BPM readout noise and correlated random scan drift. Every fit uses
one of 76 target sextupoles in 4 latent machines;
all 76 sextupole x/y offsets are fixed but hidden during a scan. Quadrupole
strength, roll, and alignment errors are absent in this bounded test.

The clean exact-state reference has beam-relative 2D RMSE
`6.054 micrometers` (median
`4.056`, P90 `8.942`).

## Acquisition and error model

- **Blocked means:** acquire `0` N times, `+` N times, then `-` N times at each
  bump and apply the unchanged symmetric three-point slope.
- **Interleaved direct:** repeat `0,+,0,-` N times, append a final `0`, average
  by K2 state, and apply the unchanged slope.
- **Interleaved 0-reference:** use the same reads, linearly interpolate the
  two adjacent zero-state readings around every `+` and `-` read, subtract
  those references, then average the paired symmetric slopes.
- BPM noise is independent Gaussian noise with `5.0 micrometer` RMS per BPM
  plane and acquisition.
- Drift is a scalar Gaussian random walk along each latent machine's fixed
  random two-plane local-bump direction. Its expected end-to-end RMS change is
  held at `10.0 micrometers` for every complete
  scan, so repeat-count comparisons hold total drift severity fixed rather
  than assuming each extra read lengthens the scan.
- The nominal states and state-specific drift secants come from paired exact
  RF-on latest-lattice SciBmad scans. Replaying arbitrary acquisition histories
  is a local linear interpolation of those physical states, not a new exact
  closed-orbit solve for every repeated read.

## Results

| nuisance | order | slope estimator | repeats N | acquisitions/scan | 2D RMSE [um] | median [um] | P90 [um] | paired increment vs clean [um] | fit at boundary |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|
| bpm_noise | blocked | direct | 1 | 15 | 1022.196 | 619.266 | 1720.788 | 1021.837 | 22.70% |
| bpm_noise | interleaved | direct | 1 | 25 | 985.424 | 516.505 | 1716.940 | 984.743 | 19.41% |
| bpm_noise | interleaved | reference_interpolated | 1 | 25 | 995.585 | 552.998 | 1749.917 | 995.154 | 20.72% |
| bpm_noise | blocked | direct | 4 | 60 | 966.507 | 504.230 | 1720.418 | 966.901 | 20.72% |
| bpm_noise | interleaved | direct | 4 | 85 | 1009.911 | 600.520 | 1743.078 | 1010.042 | 19.74% |
| bpm_noise | interleaved | reference_interpolated | 4 | 85 | 991.925 | 554.027 | 1734.859 | 992.004 | 19.41% |
| bpm_noise | blocked | direct | 16 | 240 | 955.648 | 469.709 | 1769.087 | 956.203 | 18.75% |
| bpm_noise | interleaved | direct | 16 | 325 | 975.714 | 490.499 | 1723.755 | 975.677 | 18.75% |
| bpm_noise | interleaved | reference_interpolated | 16 | 325 | 964.403 | 457.947 | 1742.626 | 964.243 | 17.11% |
| bpm_noise | blocked | direct | 64 | 960 | 987.027 | 599.312 | 1701.116 | 986.991 | 19.08% |
| bpm_noise | interleaved | direct | 64 | 1285 | 961.501 | 534.789 | 1714.548 | 961.607 | 18.09% |
| bpm_noise | interleaved | reference_interpolated | 64 | 1285 | 960.128 | 542.671 | 1722.473 | 960.296 | 17.11% |
| random_walk_drift | blocked | direct | 1 | 15 | 816.306 | 427.866 | 1532.433 | 816.221 | 12.83% |
| random_walk_drift | interleaved | direct | 1 | 25 | 779.642 | 374.154 | 1450.972 | 780.028 | 10.20% |
| random_walk_drift | interleaved | reference_interpolated | 1 | 25 | 749.912 | 400.606 | 1365.614 | 749.800 | 8.88% |
| random_walk_drift | blocked | direct | 4 | 60 | 776.767 | 359.080 | 1516.752 | 776.675 | 12.83% |
| random_walk_drift | interleaved | direct | 4 | 85 | 800.554 | 299.453 | 1594.378 | 800.906 | 11.51% |
| random_walk_drift | interleaved | reference_interpolated | 4 | 85 | 677.367 | 288.882 | 1281.451 | 677.639 | 7.89% |
| random_walk_drift | blocked | direct | 16 | 240 | 847.053 | 396.201 | 1608.444 | 846.981 | 13.82% |
| random_walk_drift | interleaved | direct | 16 | 325 | 562.682 | 215.379 | 973.323 | 562.267 | 5.59% |
| random_walk_drift | interleaved | reference_interpolated | 16 | 325 | 579.779 | 192.843 | 935.562 | 579.830 | 6.25% |
| random_walk_drift | blocked | direct | 64 | 960 | 702.550 | 368.806 | 1167.560 | 702.746 | 6.91% |
| random_walk_drift | interleaved | direct | 64 | 1285 | 441.881 | 90.085 | 418.710 | 441.850 | 3.95% |
| random_walk_drift | interleaved | reference_interpolated | 64 | 1285 | 337.492 | 66.816 | 238.163 | 337.292 | 1.64% |
| combined | blocked | direct | 1 | 15 | 1052.444 | 646.749 | 1792.584 | 1052.292 | 22.70% |
| combined | interleaved | direct | 1 | 25 | 998.758 | 508.567 | 1741.234 | 998.501 | 19.41% |
| combined | interleaved | reference_interpolated | 1 | 25 | 1010.217 | 563.455 | 1747.605 | 1009.742 | 22.04% |
| combined | blocked | direct | 4 | 60 | 942.187 | 481.671 | 1713.199 | 941.540 | 18.42% |
| combined | interleaved | direct | 4 | 85 | 1048.314 | 635.721 | 1795.364 | 1048.312 | 23.68% |
| combined | interleaved | reference_interpolated | 4 | 85 | 976.273 | 537.107 | 1702.657 | 976.338 | 18.42% |
| combined | blocked | direct | 16 | 240 | 952.009 | 480.791 | 1735.703 | 952.554 | 18.75% |
| combined | interleaved | direct | 16 | 325 | 974.223 | 496.661 | 1726.998 | 974.239 | 18.42% |
| combined | interleaved | reference_interpolated | 16 | 325 | 971.460 | 474.465 | 1703.462 | 971.189 | 16.45% |
| combined | blocked | direct | 64 | 960 | 947.084 | 555.416 | 1688.729 | 947.469 | 17.76% |
| combined | interleaved | direct | 64 | 1285 | 985.675 | 566.662 | 1739.956 | 985.667 | 18.75% |
| combined | interleaved | reference_interpolated | 64 | 1285 | 990.697 | 582.576 | 1734.519 | 990.731 | 19.08% |

## Interpretation

- **Interleaving helps correlated drift, but does not restore the clean
  inverse.** At `N=64`, the
  drift-only RMSE falls from `702.550`
  micrometers for blocked acquisition to
  `337.492` micrometers for the
  interleaved reference estimator. Its median falls from
  `368.806` to
  `66.816` micrometers. The remaining
  long tail is still far above the clean result.
- **Simple averaging is insufficient at the assumed BPM noise.** The typical
  clean outer-K2 half-difference is only `4.246 nm` RMS
  over the retained BPM channels, while one read has `5000 nm` RMS noise.
  Approximately `1,386,416` independent reads would
  be needed merely to reduce the raw mean noise to that typical signal scale;
  this is not a claim that the final center precision would then be 6
  micrometers. At `N=64`, the
  noise-only interleaved-reference RMSE remains
  `960.128` micrometers and
  `17.11%` of fits
  reach the +/-1.5-mm search boundary.
- **The combined case is noise dominated.** At the largest tested repeat
  count, blocked means give `947.084`
  micrometers and the interleaved reference estimator gives
  `990.697` micrometers. Adjacent zero
  references add their own white-noise variance, so drift cancellation cannot
  compensate while BPM noise dominates.
- The maintained inverse normalizes every BPM slope channel by its own
  across-bump RMS. When raw readout noise dominates a channel, this
  self-normalization largely removes the amplitude benefit of averaging until
  the averaged noise approaches the nanometer-scale K2 signal. The next
  bounded estimator comparison should therefore use known measurement
  covariance or fixed model/clean structural scaling, together with explicit
  drift regression; acquisition order alone is not enough.

The full curves are in `results/protocol_rmse_vs_repeats.png`.

## Integrity and provenance checks

- Lattice: `D:\Ring_Design_Development\CESR Project\Latest_Lattice\latest_cesr_scibmad_repaired.jl`.
- Baseline quadrupole-offset maximum absolute value:
  `0.000000 micrometers`.
- Latent sextupole-offset RMS over the retained physical source tensor:
  `298.143 micrometers`.
- At a 5-micrometer drift displacement, the state-specific exact-scan secants
  differ from the nominal latest-lattice bump map by
  `44.753 nm` RMS over BPM channels and
  states (relative L2 `0.01997`). This
  difference is retained by the replay rather than replaced with the nominal
  map.
- The batched inverse globally profiles the same two-source, per-channel-
  normalized objective on a 2D multiresolution grid. Against the maintained
  six-start local solver on clean data, its center-vector difference has RMS
  `0.585 micrometers`, median
  `0.004`, P90
  `0.006`, and maximum
  `10.201`. The aggregate clean RMSE
  changes only from `6.051` to `6.054 micrometers`; the larger maximum flags
  an alternate profiled minimum rather
  than target-truth leakage.
- Exact target orbit and sextupole alignment enter only after fitting, for
  evaluation.

## Run and validate

From `CESR Project/` with the validated `Ubuntu-Bmad` Python environment:

```powershell
wsl.exe -d Ubuntu-Bmad -- /home/joeyfurina/miniforge3/envs/bmad/bin/python `
  '/mnt/d/Ring_Design_Development/CESR Project/dataset_benchmark/optics/sextupole_alignment_gtpsa/interleaved_measurement_protocol/analyze_interleaved_protocol.py'

wsl.exe -d Ubuntu-Bmad -- /home/joeyfurina/miniforge3/envs/bmad/bin/python `
  '/mnt/d/Ring_Design_Development/CESR Project/dataset_benchmark/optics/sextupole_alignment_gtpsa/interleaved_measurement_protocol/validate_interleaved_protocol.py'
```

## Limitations

The 5-micrometer BPM noise and 10-micrometer random-walk span are sensitivity
settings, not measured CESR priors. The drift is restricted to the local-bump
orbit mode and one fixed direction per latent machine. Its interpolation is
first order in drift amplitude. The fixed-total-span timing assumption favors
a protocol comparison; a machine deployment study must use measured sampling
cadence, BPM covariance, nonlocal drift modes, outliers, and missing channels.

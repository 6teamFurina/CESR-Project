# GTPSA-derivative stochastic inverse result

The full latest-lattice benchmark uses all 76 active normal sextupoles and
4 hidden all-sextupole-misalignment realizations per target.  The protocol
uses delta-K2 extrema -0.100/0.100 m^-3, signed
local bumps +/-1.500 mm, and
4096 reads per signed state.  Every read has
5.0 um independent BPM white noise; the physical
orbit drift is a random walk with 10.0
um endpoint-change RMS.

| case | 2D RMSE [um] | median [um] | P90 [um] | P99 [um] | maximum [um] |
|---|---:|---:|---:|---:|---:|
| clean | 12.761 | 8.323 | 20.627 | 38.946 | 40.779 |
| bpm_white_noise | 18.728 | 14.427 | 28.859 | 45.681 | 82.837 |
| random_walk_drift | 13.978 | 9.836 | 22.332 | 37.966 | 57.198 |
| combined | 19.575 | 15.225 | 30.163 | 47.032 | 85.947 |

- combined worst target-level RMSE: 43.547 um
- combined draws below 50.0 um: 99.385%
- required threshold: 50.0 um
- acceptance gate: PASS
- exact SciBmad generation: 646.7 s for
  9120 paired states
- stochastic inverse and 512 measurement seeds:
  1.302 s

The estimator does not fit a separate propagation vector for every noisy BPM
channel.  It fixes the two local sextupole source templates using the validated
SciBmad/GTPSA transport maps, takes K2-odd/bump-odd contrasts in a time-balanced
`+,-,-,+` order, and solves only for x/y center with the known white-noise
covariance.  Even and bump-independent terms cancel from this contrast.  The
random-walk covariance is propagated with an exact reverse-cumulative closed
form, so runtime is independent of the repeat count.

This is a synthetic sensitivity benchmark.  The 5 um BPM noise and 10 um drift
span are assumed settings rather than measured CESR priors.  Gaussian white
noise is unbounded, so the threshold is imposed on aggregate RMSE, P99, and
every target-level RMSE rather than on the single largest Monte Carlo draw.

# Unknown-background linear versus nonlinear P2b

## Question

Does adding a nonlinear target-offset fit in four-source space improve the
estimate when the other 75 sextupoles have fixed but unknown offset errors?

All inverses use the same maintained `SEX_08W` exact scan with target truth
`(+350,-250) micrometers`, a fixed `300 micrometer` RMS background realization,
five cross bumps, and `delta K2 = 0,+/-0.01 m^-3`. The inverse dictionaries are
built from the nominal lattice: none of the saved background offsets are read
by P1 or P2b.

## Compared algorithms

- Nominal-bump-conditioned P1: an exact finite-difference mixed response is
  calculated at every bump with all non-target offset errors set to zero.
- Linear P2b: the same response is factorized through four local integrated
  sources and mapped linearly to target offset.
- Nonlinear P2b: the four reconstructed source slopes are fitted to a per-bump
  quadratic polynomial in target `(x,y)` offset. The polynomial is calibrated
  from a nominal `3 x 3` target-offset grid spanning `+/-0.5 mm`.

This is a first nonlinear second-stage test using the existing three-point K2
slopes. It is not yet the full per-finite-K2 thick-sextupole source fit.

## Results

| background seen by inverse | algorithm | estimated x (um) | estimated y (um) | 2D error (um) |
|---|---|---:|---:|---:|
| unknown, nominal model | P1 combined | 251.375 | -123.775 | 160.186 |
| unknown, nominal model | linear P2b | 251.830 | -123.703 | 159.964 |
| unknown, nominal model | nonlinear P2b | 251.972 | -123.438 | 160.086 |
| no background errors (closure) | linear P2b | 351.829 | -251.093 | 2.130 |
| no background errors (closure) | nonlinear P2b | 350.310 | -249.657 | 0.462 |

The nominal nonlinear source calibration itself closes accurately: the
per-bump source-space polynomial residuals are between `1.30e-10` and
`1.75e-10` in the saved local-source units.

## Interpretation

The nonlinear second stage works when the forward background matches the
nominal inverse: it reduces the target-only closure error from `2.130` to
`0.462 micrometers`. It does not repair the unknown-background case; the error
changes only from `159.964` to `160.086 micrometers`.

For this realization, the dominant error is therefore upstream of the
nonlinear center fit. Unknown offsets alter the source-to-observable propagation
and bias the four reconstructed sources. Once those biased sources have been
formed, a more accurate local sextupole polynomial fits the wrong sources and
cannot recover the true center.

This is one noiseless target/background realization, not an ensemble result.
The next experiment should compare reconstructed sources directly with exact
local-source truth under multiple backgrounds, then test baseline-measurement
conditioning or nuisance-augmented source reconstruction before adding more
complexity to the nonlinear second stage.

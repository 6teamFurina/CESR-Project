# Paired P0--P3 physical-inversion benchmark

## Scope

All methods consume the same saved `SEX_08W` exact scan:

- target truth: `(x, y) = (+350, -250) micrometers`;
- other 75 sextupoles: one fixed independent `300 micrometer` RMS realization;
- five achieved two-plane bump states;
- `delta K2 = 0, +/-0.01 m^-3`;
- no added measurement noise or missing channels.

P0 uses the nominal response map. P1--P3 know the saved other-sextupole
realization while keeping the target offset hidden from initialization. This
is an oracle-background diagnostic, not a claim of real-machine precision.

## Results

| method | observable view | estimated x (um) | estimated y (um) | 2D error (um) |
|---|---|---:|---:|---:|
| P0 nominal mixed GTPSA | orbit | 360.927 | -223.410 | 28.748 |
| P1 background-conditioned mixed response | orbit | 350.714 | -257.327 | 7.362 |
| P2a two local dipole sources | orbit | 350.714 | -257.327 | 7.362 |
| P0 nominal mixed GTPSA | orbit + phase + coupling + tune | 264.600 | -123.888 | 152.307 |
| P1 background-conditioned mixed response | orbit + phase + coupling + tune | 351.867 | -251.067 | 2.150 |
| P2b four local sources | orbit + phase + coupling + tune | 351.879 | -251.064 | 2.160 |
| P3 exact full-scan inverse, initial P2b point | same combined view | 351.879 | -251.064 | 2.160 |
| P3 after one exact Gauss--Newton update | same combined view | 350.001 | -250.000 | 0.000568 |
| P3 after two updates | same combined view | 350.000 | -250.000 | 0.000004 |

## Runtime on the maintained case

The common 15-state measured/synthetic scan is excluded from all inverse
timings. Julia/SciBmad compilation warm-up is reported separately.

| method stage | warm exact states | warm time (s) | note |
|---|---:|---:|---|
| P1 conditioned mixed-response construction | 50 | 53.384 | five bumps, nested central differences in K2 and target x/y offset |
| P2 four-source construction, incremental beyond P1 | 40 | 45.404 | central differences of Kn0L, Ks0L, Kn1L, and Ks1L at five bumps |
| P2 total model preparation including P1 | 90 | 98.788 | required by the present P2 implementation |
| P0--P2 Python solve and output, combined wall time | 0 | 0.953 | includes interpreter startup, CSV I/O, all P0/P1/P2 cases, and output |
| P3 exact nonlinear inverse, three updates | 150 | 139.259 | ten complete 15-state forward-scan evaluations |
| P3 Julia/SciBmad warm-up | 1 | 36.578 | excluded from the warm inverse time |
| P3 complete external wall time | -- | 190.5 | load, warm-up, inverse, and output |

These are serial timings from one run on the current workstation, not a
throughput benchmark. P1/P2 central differences can be replaced by GTPSA
parameter derivatives, and P3 candidate/Jacobian evaluations can be cached,
parallelized, or evaluated with a validated surrogate.

## Interpretation

The nominal optics failure is primarily a response-conditioning failure in
this realization. Adding phase/coupling/tune to the fixed nominal response
increases the error from `28.748` to `152.307 micrometers`; using the same
observables with a response conditioned on the actual background reduces the
error to `2.150 micrometers`. Therefore this case does not support discarding
optics measurements.

P2a and P2b match their corresponding P1 estimates closely because the first
P2 implementation is a covariance-aware factorization through the local
two-/four-source subspace, followed by the same conditioned linear offset
model. This validates the local-source reconstruction and shows that it does
not lose useful information here; it does not yet demonstrate an advantage
from a separate nonlinear sextupole-kick equation.

P3 removes the remaining local-linearization error and recovers the generating
truth to numerical precision. This near-zero error is expected because P3 uses
the same exact SciBmad model, the true other-sextupole background, and noiseless
synthetic observations. It is an inverse-consistency check and an upper-bound
oracle, not a realistic accuracy forecast.

## Next decision experiment

The next useful expansion is an ensemble of target truths and background
realizations with three levels of information:

1. oracle background conditioning, as here;
2. conditioning inferred only from baseline orbit/optics;
3. nominal model with explicit nuisance columns and priors.

Add measured-style covariance only after the noiseless algorithmic bias is
characterized. The nonlinear analytic local-kick fit should then be compared
with the current linear P2 and exact P3.

# SciBmad closed-orbit initial-guess comparison

Date: 2026-07-30  
Machine: local Windows PC, AMD Ryzen 9 5900HX  
Julia: 1.12.6, one Julia thread  
Dataset: the same 1,000 CESR corrector samples, 119 controls, RF on  
Solver tolerances: `reltol = abstol = 1e-13`

## Result

| Initial guess | Converged | Closed-orbit solve | Solve + detector tracking | Throughput | Newton iterations (min / median / mean / max) |
|---|---:|---:|---:|---:|---:|
| Six-dimensional zero | 1000/1000 | 63.786 s | 64.356 s | 15.539 samples/s | 3 / 3 / 4.034 / 12 |
| Nominal closed orbit `z0` | 1000/1000 | 65.124 s | 65.709 s | 15.219 samples/s | 0 / 3 / 3.101 / 12 |

Computing `z0` separately required 0.151 s of model setup and 0.053 s of
closed-orbit solving (0.204 s total). This cost is reported separately and is
not included in the 65.709 s physics time.

Using `z0` reduced the mean iteration count by 23.1%, but the measured physics
time was 2.1% slower in this single paired run. This small timing difference
should be treated as ordinary run-to-run variation, not evidence that `z0`
intrinsically makes the calculation slower.

The important observation is that both batches still required a maximum of 12
Newton iterations. SciBmad's current batched Newton loop continues evaluating
the full batch until the slowest lane has converged. It freezes updates for
already-converged lanes, but their lower average iteration count does not
shorten the number of full-batch residual/Jacobian evaluations. Consequently,
using one shared `z0` alone does not improve this 1,000-sample batch runtime.

## Numerical equivalence

The two runs used different initial guesses but converged to the same physical
solutions within solver precision:

- compared observable entries: 198,000
- maximum absolute difference: `1.574e-13 m`
- RMSE: `1.773e-14 m`
- convergence flags: identical

## Nominal closed orbit

In SciBmad coordinate order `[x, px, y, py, z, pz]`,

```text
[-1.6685198602033575e-05,
  2.3901125031495305e-03,
  1.0543112225980547e-06,
  1.7969303419546103e-06,
 -3.9936116702395509e-04,
 -7.8037211957000425e-06]
```

The next useful experiment is a sample-specific first-order predictor
`z_initial = z0 + (dz/dk) * delta_k`. It must reduce the maximum iteration
count, not only the mean, to accelerate the current single 1,000-sample batch.

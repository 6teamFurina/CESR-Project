# Nonlinear-rho orbit benchmark results

The comparison contains 9,000 nonzero shared inputs (600 for every scenario/rho cell), plus one shared baseline.

| Metric | SciBmad | Bmad/Tao |
|---|---:|---:|
| Converged nonzero inputs | 9000/9000 | 9000/9000 |
| Paired converged inputs | 9000/9000 | 9000/9000 |
| Summed physics time | 28.204 s | 102.725 s |
| Throughput | 319.104 samples/s | 87.613 samples/s |
| SciBmad speedup, physics only | 3.642x | 1x |
| SciBmad setup + physics time | 30.690 s | n/a |
| SciBmad speedup including per-group model setup | 3.347x | 1x |
| SciBmad initial-guess + model setup + physics time | 30.907 s | n/a |
| SciBmad speedup including all runtime setup | 3.324x | 1x |

## Zero-input baseline agreement

| Plane | RMSE [m] | Maximum absolute difference [m] | Correlation |
|---|---:|---:|---:|
| x | 3.1053668e-06 | 7.4018708e-06 | 0.99999997 |
| y | 5.3155981e-08 | 2.4343402e-07 | 1 |

## Per-cell results

The response RMSE compares each engine's orbit after subtracting that engine's own zero-input baseline.

| Scenario | rho | Sci conv. | Bmad conv. | response x RMSE [m] | x relative RMSE | response y RMSE [m] | y relative RMSE | Sci/Bmad speedup (physics) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| all | 1.13 | 600/600 | 600/600 | 1.2565763e-07 | 0.0493% | 8.4084862e-07 | 0.2458% | 4.1712151x |
| all | 3.2 | 600/600 | 600/600 | 3.5991844e-07 | 0.0499% | 2.4118219e-06 | 0.2490% | 3.9047185x |
| all | 4.53 | 600/600 | 600/600 | 5.1736089e-07 | 0.0507% | 3.4674589e-06 | 0.2529% | 3.1970929x |
| all | 6.4 | 600/600 | 600/600 | 7.5406613e-07 | 0.0523% | 5.0562011e-06 | 0.2610% | 3.2862863x |
| all | 9.05 | 600/600 | 600/600 | 1.1333396e-06 | 0.0555% | 7.6325421e-06 | 0.2785% | 2.7729236x |
| horizontal | 1.13 | 600/600 | 600/600 | 1.2190868e-07 | 0.0478% | 7.8127903e-09 | 0.1869% | 5.0589126x |
| horizontal | 3.2 | 600/600 | 600/600 | 3.4898964e-07 | 0.0483% | 2.2033952e-08 | 0.1861% | 5.2549482x |
| horizontal | 4.53 | 600/600 | 600/600 | 4.9862009e-07 | 0.0488% | 3.1110611e-08 | 0.1856% | 4.243857x |
| horizontal | 6.4 | 600/600 | 600/600 | 7.1568502e-07 | 0.0496% | 4.3794821e-08 | 0.1849% | 4.3330451x |
| horizontal | 9.05 | 600/600 | 600/600 | 1.0404607e-06 | 0.0510% | 6.1629334e-08 | 0.1840% | 4.3043048x |
| vertical | 1.13 | 600/600 | 600/600 | 1.0334603e-08 | 0.3583% | 8.3307343e-07 | 0.2388% | 4.8097325x |
| vertical | 3.2 | 600/600 | 600/600 | 5.59537e-08 | 0.5549% | 2.3872227e-06 | 0.2416% | 3.4830219x |
| vertical | 4.53 | 600/600 | 600/600 | 1.0734724e-07 | 0.6363% | 3.4256999e-06 | 0.2449% | 2.9410601x |
| vertical | 6.4 | 600/600 | 600/600 | 2.1011947e-07 | 0.7044% | 4.9730169e-06 | 0.2516% | 2.8939268x |
| vertical | 9.05 | 600/600 | 600/600 | 4.1957905e-07 | 0.7572% | 7.4275511e-06 | 0.2657% | 2.6334148x |

## Timing interpretation

SciBmad was run in one Julia process, with 600 simultaneous TPSA lanes per cell. Its physics-only number includes frozen-Jacobian iterations, explicit closure checks, tracking, and any full-AD fallback. The end-to-end variant additionally includes construction of each 600-lane model, but excludes compilation warmup, CSV I/O, and the one-time first-order initial-guess preparation reported in metadata.

Bmad was run sequentially in one persistent Tao/PyTao process in Ubuntu-Bmad. Its timed region includes corrector updates, one Tao model recalculation, and observable reads. Its convergence flag is based on Tao good_model flags and finite data; this path does not expose the explicit one-turn closure norm used by SciBmad.

The two engines ran on the same physical machine but in different host runtimes (SciBmad on Windows Julia and Bmad inside WSL Ubuntu), so timing is an application-level comparison rather than a microarchitectural kernel benchmark.

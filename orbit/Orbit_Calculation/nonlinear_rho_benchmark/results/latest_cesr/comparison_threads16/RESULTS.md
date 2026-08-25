# Nonlinear-rho orbit benchmark results

The comparison contains 9,000 nonzero shared inputs (up to 600 for any scenario/rho cell), plus one shared baseline.

| Metric | SciBmad | Bmad/Tao |
|---|---:|---:|
| Converged nonzero inputs | 9000/9000 | 9000/9000 |
| Paired converged inputs | 9000/9000 | 9000/9000 |
| Julia threads | 16 | n/a |
| CPU multithreaded tracking | True | native Tao path |
| Summed physics time | 47.746 s | 166.608 s |
| Throughput | 188.497 samples/s | 54.019 samples/s |
| SciBmad speedup, physics only | 3.489x | 1x |
| SciBmad setup + physics time | 168.423 s | n/a |
| SciBmad speedup including per-group model setup | 0.989x | 1x |
| SciBmad initial-guess + model setup + physics time | 174.171 s | n/a |
| SciBmad speedup including all runtime setup | 0.957x | 1x |

## Zero-input baseline agreement

| Plane | RMSE [m] | Maximum absolute difference [m] | Correlation |
|---|---:|---:|---:|
| x | 2.2307821e-15 | 7.0241015e-15 | 0.5972088 |
| y | 3.1305895e-16 | 7.2688577e-16 | -0.28381093 |

## Per-cell results

The response RMSE compares each engine's orbit after subtracting that engine's own zero-input baseline.

| Scenario | rho | Sci conv. | Bmad conv. | response x RMSE [m] | x relative RMSE | response y RMSE [m] | y relative RMSE | Sci/Bmad speedup (physics) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| all | 1.13 | 600/600 | 600/600 | 8.0387074e-09 | 0.0036% | 6.8648042e-09 | 0.0019% | 3.5544296x |
| all | 3.2 | 600/600 | 600/600 | 6.4560839e-08 | 0.0101% | 6.0554324e-08 | 0.0060% | 3.5785172x |
| all | 4.53 | 600/600 | 600/600 | 1.2964387e-07 | 0.0144% | 1.3265229e-07 | 0.0094% | 3.5054269x |
| all | 6.4 | 600/600 | 600/600 | 2.59859e-07 | 0.0204% | 3.0455733e-07 | 0.0152% | 3.2905199x |
| all | 9.05 | 600/600 | 600/600 | 5.2422215e-07 | 0.0290% | 7.43733e-07 | 0.0262% | 3.3134574x |
| horizontal | 1.13 | 600/600 | 600/600 | 3.148776e-09 | 0.0014% | 1.1491192e-14 | 100.0000% | 3.7194866x |
| horizontal | 3.2 | 600/600 | 600/600 | 2.5346661e-08 | 0.0041% | 3.1304577e-16 | 99.9957% | 3.610068x |
| horizontal | 4.53 | 600/600 | 600/600 | 5.0976542e-08 | 0.0058% | 3.1303085e-16 | 99.9894% | 3.5189792x |
| horizontal | 6.4 | 600/600 | 600/600 | 1.0238414e-07 | 0.0082% | 3.1304654e-16 | 99.9917% | 3.5226703x |
| horizontal | 9.05 | 600/600 | 600/600 | 2.0700327e-07 | 0.0118% | 3.1296671e-16 | 99.9693% | 3.4639955x |
| vertical | 1.13 | 600/600 | 600/600 | 6.9836043e-09 | 0.5530% | 1.080615e-09 | 0.0003% | 3.5580328x |
| vertical | 3.2 | 600/600 | 600/600 | 5.6048589e-08 | 0.5533% | 2.4510638e-08 | 0.0025% | 3.3236009x |
| vertical | 4.53 | 600/600 | 600/600 | 1.1245906e-07 | 0.5539% | 6.9580434e-08 | 0.0049% | 3.3873962x |
| vertical | 6.4 | 600/600 | 600/600 | 2.2506105e-07 | 0.5552% | 1.964963e-07 | 0.0099% | 3.4138908x |
| vertical | 9.05 | 600/600 | 600/600 | 4.5271759e-07 | 0.5581% | 5.5721015e-07 | 0.0198% | 3.6551604x |

## Timing interpretation

SciBmad was run in one Julia process, with the configured simultaneous TPSA lanes per cell. Its physics-only number includes frozen-Jacobian iterations, explicit closure checks, tracking, and any full-AD fallback. The end-to-end variant additionally includes construction of each batch model, but excludes compilation warmup, CSV I/O, and the one-time first-order initial-guess preparation reported in metadata.

Bmad was run sequentially in one persistent Tao/PyTao process in Ubuntu-Bmad. Its timed region includes corrector updates, one Tao model recalculation, and observable reads. Its convergence flag is based on Tao good_model flags and finite data; this path does not expose the explicit one-turn closure norm used by SciBmad.

The two engines ran on the same physical machine but in different host runtimes (SciBmad on Windows Julia and Bmad inside WSL Ubuntu), so timing is an application-level comparison rather than a microarchitectural kernel benchmark. The recorded Julia thread count and CPU-multithreading flag define the SciBmad hardware usage; the Bmad reference uses Tao's native execution path rather than a manually threaded sample loop.

# Nonlinear-rho orbit benchmark results

The comparison contains 9,000 nonzero shared inputs (up to 600 for any scenario/rho cell), plus one shared baseline.

| Metric | SciBmad | Bmad/Tao |
|---|---:|---:|
| Converged nonzero inputs | 9000/9000 | 9000/9000 |
| Paired converged inputs | 9000/9000 | 9000/9000 |
| Julia threads | 1 | n/a |
| CPU multithreaded tracking | False | native Tao path |
| Summed physics time | 105.703 s | 166.608 s |
| Throughput | 85.145 samples/s | 54.019 samples/s |
| SciBmad speedup, physics only | 1.576x | 1x |
| SciBmad setup + physics time | 273.372 s | n/a |
| SciBmad speedup including per-group model setup | 0.609x | 1x |
| SciBmad initial-guess + model setup + physics time | 281.879 s | n/a |
| SciBmad speedup including all runtime setup | 0.591x | 1x |

## Zero-input baseline agreement

| Plane | RMSE [m] | Maximum absolute difference [m] | Correlation |
|---|---:|---:|---:|
| x | 2.2307821e-15 | 7.0241015e-15 | 0.5972088 |
| y | 3.1305895e-16 | 7.2688577e-16 | -0.28381093 |

## Per-cell results

The response RMSE compares each engine's orbit after subtracting that engine's own zero-input baseline.

| Scenario | rho | Sci conv. | Bmad conv. | response x RMSE [m] | x relative RMSE | response y RMSE [m] | y relative RMSE | Sci/Bmad speedup (physics) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| all | 1.13 | 600/600 | 600/600 | 8.0387074e-09 | 0.0036% | 6.8648042e-09 | 0.0019% | 1.8314206x |
| all | 3.2 | 600/600 | 600/600 | 6.4560839e-08 | 0.0101% | 6.0554324e-08 | 0.0060% | 1.6537727x |
| all | 4.53 | 600/600 | 600/600 | 1.2964387e-07 | 0.0144% | 1.3265229e-07 | 0.0094% | 1.4081648x |
| all | 6.4 | 600/600 | 600/600 | 2.59859e-07 | 0.0204% | 3.0455733e-07 | 0.0152% | 1.4480552x |
| all | 9.05 | 600/600 | 600/600 | 5.2422216e-07 | 0.0290% | 7.4373304e-07 | 0.0262% | 1.3722474x |
| horizontal | 1.13 | 600/600 | 600/600 | 3.1487759e-09 | 0.0014% | 1.1491194e-14 | 100.0001% | 1.8398261x |
| horizontal | 3.2 | 600/600 | 600/600 | 2.5346661e-08 | 0.0041% | 3.1304921e-16 | 99.9968% | 1.7679107x |
| horizontal | 4.53 | 600/600 | 600/600 | 5.0976542e-08 | 0.0058% | 3.1303283e-16 | 99.9901% | 1.6226038x |
| horizontal | 6.4 | 600/600 | 600/600 | 1.0238414e-07 | 0.0082% | 3.1304625e-16 | 99.9916% | 1.6601433x |
| horizontal | 9.05 | 600/600 | 600/600 | 2.0700327e-07 | 0.0118% | 3.1298673e-16 | 99.9757% | 1.6140968x |
| vertical | 1.13 | 600/600 | 600/600 | 6.9836043e-09 | 0.5530% | 1.080615e-09 | 0.0003% | 1.6170631x |
| vertical | 3.2 | 600/600 | 600/600 | 5.6048589e-08 | 0.5533% | 2.4510638e-08 | 0.0025% | 1.510946x |
| vertical | 4.53 | 600/600 | 600/600 | 1.1245906e-07 | 0.5539% | 6.9580434e-08 | 0.0049% | 1.5290186x |
| vertical | 6.4 | 600/600 | 600/600 | 2.2506105e-07 | 0.5552% | 1.964963e-07 | 0.0099% | 1.4906771x |
| vertical | 9.05 | 600/600 | 600/600 | 4.5271759e-07 | 0.5581% | 5.5721015e-07 | 0.0198% | 1.5055954x |

## Timing interpretation

SciBmad was run in one Julia process, with the configured simultaneous TPSA lanes per cell. Its physics-only number includes frozen-Jacobian iterations, explicit closure checks, tracking, and any full-AD fallback. The end-to-end variant additionally includes construction of each batch model, but excludes compilation warmup, CSV I/O, and the one-time first-order initial-guess preparation reported in metadata.

Bmad was run sequentially in one persistent Tao/PyTao process in Ubuntu-Bmad. Its timed region includes corrector updates, one Tao model recalculation, and observable reads. Its convergence flag is based on Tao good_model flags and finite data; this path does not expose the explicit one-turn closure norm used by SciBmad.

The two engines ran on the same physical machine but in different host runtimes (SciBmad on Windows Julia and Bmad inside WSL Ubuntu), so timing is an application-level comparison rather than a microarchitectural kernel benchmark. The recorded Julia thread count and CPU-multithreading flag define the SciBmad hardware usage; the Bmad reference uses Tao's native execution path rather than a manually threaded sample loop.

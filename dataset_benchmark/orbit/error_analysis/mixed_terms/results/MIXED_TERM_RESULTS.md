# All-corrector horizontal--vertical mixed-term experiment

The four-sign finite difference directly separates the pure horizontal
`Q_hh`, pure vertical `Q_vv`, and mixed `Q_hv` response vectors. All
reconstructions and remainders are computed as vectors before detector RMS.

## Mean decomposition

| rho | H/V kick RMS (urad) | plane | Q_hh (um) | Q_vv (um) | Q_hv (um) | Q_hv slope | mixed energy share | full remainder / exact |
|---:|---:|:---:|---:|---:|---:|---:|---:|---:|
| 0.1 | 0.5 | X | 0.00165635 | 0.00534217 | 5.05141e-05 | n/a | 0.0146% | 0.0283% |
| 0.1 | 0.5 | Y | 3.69664e-05 | 7.69898e-05 | 0.00559752 | n/a | 99.9627% | 0.0844% |
| 0.14 | 0.7 | X | 0.00324645 | 0.0104707 | 9.90076e-05 | 2.0000 | 0.0146% | 0.0396% |
| 0.14 | 0.7 | Y | 7.24541e-05 | 0.000150899 | 0.0109711 | 2.0000 | 99.9627% | 0.1181% |
| 0.2 | 1 | X | 0.0066254 | 0.0213687 | 0.000202056 | 2.0000 | 0.0146% | 0.0565% |
| 0.2 | 1 | Y | 0.000147866 | 0.000307957 | 0.0223901 | 2.0000 | 99.9627% | 0.1688% |
| 0.28 | 1.4 | X | 0.0129858 | 0.0418826 | 0.000396031 | 2.0000 | 0.0146% | 0.0791% |
| 0.28 | 1.4 | Y | 0.000289817 | 0.000603596 | 0.0438846 | 2.0000 | 99.9627% | 0.2363% |
| 0.4 | 2 | X | 0.0265016 | 0.0854748 | 0.000808227 | 2.0000 | 0.0146% | 0.1130% |
| 0.4 | 2 | Y | 0.000591463 | 0.00123183 | 0.0895604 | 2.0000 | 99.9627% | 0.3375% |
| 0.57 | 2.85 | X | 0.0538149 | 0.173567 | 0.00164121 | 2.0000 | 0.0146% | 0.1611% |
| 0.57 | 2.85 | Y | 0.00120104 | 0.00250138 | 0.181864 | 2.0000 | 99.9627% | 0.4809% |
| 0.8 | 4 | X | 0.106007 | 0.3419 | 0.00323293 | 2.0000 | 0.0146% | 0.2261% |
| 0.8 | 4 | Y | 0.00236585 | 0.00492731 | 0.358242 | 2.0000 | 99.9627% | 0.6750% |
| 1.13 | 5.65 | X | 0.2115 | 0.682144 | 0.00645024 | 2.0000 | 0.0146% | 0.3193% |
| 1.13 | 5.65 | Y | 0.00472024 | 0.00983074 | 0.71475 | 2.0000 | 99.9627% | 0.9533% |

## Direction-resolved result at the largest fitted radius

At `rho = 1.13` (5.65 urad RMS in each family):

| plane | mixed/pure P10 | median | P90 | mixed-energy share P10 | median | P90 | reconstruction improvement P10 | median | P90 |
|:---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| X | 0.004129 | 0.009872 | 0.01781 | 0.002% | 0.010% | 0.032% | 58.561% | 92.888% | 98.074% |
| Y | 38.52 | 64.73 | 129.2 | 99.933% | 99.976% | 99.994% | 99.957% | 99.996% | 99.999% |
| X+Y | 0.659 | 1.075 | 1.581 | 30.280% | 53.620% | 71.429% | 99.955% | 99.995% | 99.998% |

## Compact checks

- Mean Y `Q_hv/rho^2` changes by `0.0004233%` over the fitted interval.
- At the smallest radius, mean Y mixed energy share is `99.9627%`.
- At the largest radius, mean Y mixed energy share is `99.9627%`.
- At the largest radius, the full signed reconstruction leaves `0.953287%` of the exact Y residual RMS.
- At the largest radius, adding `Q_hv` reduces the mean squared Y reconstruction error by `99.982%` relative to the pure-only reconstruction.
- For the combined 198-dimensional X+Y vector at the largest radius, the pooled mixed squared-norm share is `47.0412%`; its direction-resolved P10/median/P90 shares are `30.28% / 53.62% / 71.43%`.

The cross-term hypothesis is supported only if `Q_hv` remains quadratic,
dominates the Y direction distribution, and the signed reconstruction closes
the four exact joint responses with a much smaller higher-order remainder.

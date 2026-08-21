# All-corrector horizontal--vertical mixed-term experiment

The four-sign finite difference directly separates the pure horizontal
`Q_hh`, pure vertical `Q_vv`, and mixed `Q_hv` response vectors. All
reconstructions and remainders are computed as vectors before detector RMS.

## Mean decomposition

| rho | H/V kick RMS (urad) | plane | Q_hh (um) | Q_vv (um) | Q_hv (um) | Q_hv slope | mixed energy share | full remainder / exact |
|---:|---:|:---:|---:|---:|---:|---:|---:|---:|
| 0.1 | 0.5 | X | 0.00185604 | 0.00559683 | 3.58685e-05 | n/a | 0.0084% | 0.0172% |
| 0.1 | 0.5 | Y | 3.04903e-05 | 6.40431e-05 | 0.00502417 | n/a | 99.9867% | 0.0703% |
| 0.2 | 1 | X | 0.00742412 | 0.0223873 | 0.000143483 | 2.0001 | 0.0084% | 0.0342% |
| 0.2 | 1 | Y | 0.000121961 | 0.000256107 | 0.0200967 | 2.0000 | 99.9867% | 0.1398% |

## Direction-resolved result at the largest fitted radius

At `rho = 0.2` (1 urad RMS in each family):

| plane | mixed/pure P10 | median | P90 | mixed-energy share P10 | median | P90 | reconstruction improvement P10 | median | P90 |
|:---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| X | 0.005365 | 0.00838 | 0.0114 | 0.003% | 0.008% | 0.013% | 99.330% | 99.603% | 99.876% |
| Y | 85.91 | 86.66 | 87.41 | 99.986% | 99.987% | 99.987% | 100.000% | 100.000% | 100.000% |

## Compact checks

- Mean Y `Q_hv/rho^2` changes by `2.1441e-05%` over the fitted interval.
- At the smallest radius, mean Y mixed energy share is `99.9867%`.
- At the largest radius, mean Y mixed energy share is `99.9867%`.
- At the largest radius, the full signed reconstruction leaves `0.139828%` of the exact Y residual RMS.
- At the largest radius, adding `Q_hv` reduces the mean squared Y reconstruction error by `99.9997%` relative to the pure-only reconstruction.

The cross-term hypothesis is supported only if `Q_hv` remains quadratic,
dominates the Y direction distribution, and the signed reconstruction closes
the four exact joint responses with a much smaller higher-order remainder.

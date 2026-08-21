# K2/wiggler cubic-vector corner attribution

The four-corner inclusion-exclusion decomposition is

`C11 = C00 + (C10-C00) + (C01-C00) + (C11-C10-C01+C00)`.

| Component | Global norm / nominal | Global signed projection | Direction median norm [P10, P90] |
|---|---:|---:|---:|
| residual | 0.611670 | 0.578853 | 0.373856 [0.088196, 0.704985] |
| sextupole | 0.304397 | -0.105594 | 0.676804 [0.128797, 1.051059] |
| wiggler | 0.604632 | 0.566636 | 0.461420 [0.032578, 0.839816] |
| interaction | 0.095775 | -0.039894 | 0.077226 [0.047331, 0.181803] |

Norm fractions do not add because the signed detector vectors interfere. Signed projections add to one (up to numerical extraction error) and quantify alignment with the nominal response.

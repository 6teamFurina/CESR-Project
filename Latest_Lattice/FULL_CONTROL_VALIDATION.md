# Full latest-CESR control tracking validation

- Bmad/Tao reference: `20260814-0`, branch 0 changed to open geometry for a one-pass map.
- Controls tested: `124` (all 119 Overlay and 5 Group lords).
- Bmad lord-to-slave relationships covered: `347`.
- Response observations compared: `475`.
- Fixed nonzero particle start: `[0.001, 0.0002, -0.0008, 0.00015, 0.0005, 0.0002]`.
- Maximum absolute six-vector entry difference: `6.640276e-03`.
- Maximum relative L2 difference for informative responses: `1.836114e-02`.

## Largest relative discrepancies

| Control | Type | Observation | Bmad norm | Relative L2 | Max abs | Cosine |
|---|---|---:|---:|---:|---:|---:|
| `SK_Q18W` | `Overlay` | 295 | 6.336929e-05 | 1.836114e-02 | 1.078849e-06 | 0.999834149 |
| `RAW_XQUNEING_1` | `Group` | 753 | 2.874278e-06 | 1.072408e-02 | 3.071388e-08 | 0.999982766 |
| `H_CANT_TRIM_S4` | `Overlay` | 1169 | 7.776289e-01 | 8.675617e-03 | 6.640276e-03 | 0.999999258 |
| `H_CANT_TRIM_S4` | `Overlay` | 1177 | 8.438008e-01 | 7.886746e-03 | 6.566005e-03 | 0.999999558 |
| `RAW_XQUNEING_1` | `Group` | 699 | 4.707242e-06 | 5.270536e-03 | 2.466310e-08 | 0.999986135 |
| `RAW_XQUNEING_2` | `Group` | 254 | 2.903222e-07 | 4.945409e-03 | 1.310795e-09 | 0.999994842 |
| `RAW_XQUNEING_1` | `Group` | 667 | 4.837005e-06 | 4.517558e-03 | 2.167392e-08 | 0.999993961 |
| `SK_Q24E` | `Overlay` | 831 | 3.064674e-04 | 3.819730e-03 | 1.138386e-06 | 0.999993039 |
| `RAW_XQUNEING_1` | `Group` | 651 | 5.579498e-06 | 3.781578e-03 | 2.097643e-08 | 0.999994962 |
| `RAW_XQUNEING_1` | `Group` | 479 | 3.104089e-06 | 3.691304e-03 | 1.128086e-08 | 0.999996669 |
| `RAW_XQUNEING_1` | `Group` | 787 | 5.145996e-06 | 3.569272e-03 | 1.825176e-08 | 0.999997863 |
| `RAW_XQUNEING_2` | `Group` | 865 | 4.533460e-06 | 3.503396e-03 | 1.511721e-08 | 0.999993947 |
| `RAW_XQUNEING_1` | `Group` | 737 | 6.338509e-06 | 3.397623e-03 | 2.139921e-08 | 0.999994922 |
| `RAW_XQUNEING_2` | `Group` | 815 | 4.307633e-06 | 3.348985e-03 | 1.414898e-08 | 0.999997281 |
| `RAW_XQUNEING_1` | `Group` | 823 | 9.601671e-06 | 3.299347e-03 | 3.151111e-08 | 0.999999156 |
| `SK_Q24E` | `Overlay` | 1177 | 8.114139e-04 | 3.284926e-03 | 2.629991e-06 | 0.999999555 |
| `RAW_XQUNEING_1` | `Group` | 407 | 2.296663e-06 | 3.194732e-03 | 7.156523e-09 | 0.999998991 |
| `RAW_BETASING_8` | `Group` | 1177 | 7.798254e-05 | 3.188419e-03 | 1.874639e-07 | 0.999995729 |
| `RAW_XQUNEING_1` | `Group` | 743 | 6.030168e-06 | 2.933502e-03 | 1.757743e-08 | 0.999995854 |
| `H26W` | `Overlay` | 1177 | 1.330786e+00 | 2.753620e-03 | 3.458538e-03 | 0.999996304 |

The CSV contains both six-component response vectors for every observation.

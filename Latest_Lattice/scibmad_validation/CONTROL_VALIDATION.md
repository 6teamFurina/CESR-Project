# Latest CESR control validation

Representative Overlay and Group variables were changed after lattice construction to verify that `DefExpr` relationships remain live.

| Control variable | Element | Attribute | Expected slope | Measured slope | Abs. error |
|---|---|---|---:|---:|---:|
| `H11W_hkick` | `b11w` | `Kn0` | -1.521083279005e-01 | -1.521083279006e-01 | 7.658e-14 |
| `H12W_hkick` | `b12w` | `Kn0` | -7.605416395027e-02 | -7.605416395030e-02 | 3.829e-14 |
| `H12W_hkick` | `b13w` | `Kn0` | -7.604231968009e-02 | -7.604231967995e-02 | 1.393e-13 |
| `V09AW_vkick` | `sex_09aw` | `Ks0` | 3.676470588235e+00 | 3.676470588235e+00 | 4.441e-16 |
| `SK_Q14W_k1` | `sex_14w` | `Ks1` | -1.000000000000e+00 | -1.000000000000e+00 | 0.000e+00 |
| `H_CANT_S3_hkick` | `hs3a` | `Kn0` | -5.914545454545e+00 | -5.914545454545e+00 | 8.882e-16 |
| `H_CANT_S3_hkick` | `hs3c` | `Kn0` | -5.914545454545e+00 | -5.914545454545e+00 | 8.882e-16 |
| `H_CANT_S3_hkick` | `hs3b!s1` | `Kn0` | 2.000000000000e+01 | 2.000000000000e+01 | 0.000e+00 |
| `H_CANT_S3_hkick` | `hs3b!s2` | `Kn0` | 2.000000000000e+01 | 2.000000000000e+01 | 0.000e+00 |
| `H48W_hkick` | `b48w!s1` | `Kn0` | -3.395223735058e-01 | -3.395223735058e-01 | 4.247e-14 |
| `H48W_hkick` | `b48w!s2` | `Kn0` | -3.395223735058e-01 | -3.395223735058e-01 | 4.247e-14 |
| `RAW_XQUNEING_1_command` | `sex_12w` | `Kn2` | -8.718000000000e-03 | -8.717999999908e-03 | 9.237e-14 |
| `RAW_XQUNEING_1_command` | `sex_27w` | `Kn2` | -5.174220000000e-01 | -5.174220000002e-01 | 1.561e-13 |
| `RAW_XQUNEING_2_command` | `sex_47e` | `Kn2` | 2.996750000000e-01 | 2.996749999999e-01 | 5.629e-14 |

Maximum tested slope error: `1.561e-13`.

The tests cover a one-slave bend overlay, a two-slave overlay, vertical and skew correctors, the repaired split-superlord cases, and both sextupole Group commands.

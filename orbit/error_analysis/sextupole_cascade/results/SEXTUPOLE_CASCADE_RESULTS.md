# Sextupole-cascade strength-scan result

The scan uses 100 fixed vertical-corrector directions. Every lattice variant has its own nominal closed orbit and linear detector response.

## Global result

- With all order-2 multipoles removed, the C3 norm fraction has median 0.8414 and P10--P90 [0.1537, 1.3822] across directions.
- Its signed projection onto the nominal C3 vector has median 0.6685 and P10--P90 [0.0681, 1.2353].
- The largest global vector residual of the anchored model C3(lambda)=C3(0)+lambda^2[C3(1)-C3(0)] over interior scan points is 12.5362% of the nominal norm.
- A full vector fit C3(lambda)=A0+A1*lambda+A2*lambda^2 has maximum scan-point residual 0.2074% of the nominal norm.
- The fitted component norm fractions (signed projections onto nominal in parentheses) are: A0 1.1949 (1.1454), A1 0.5012 (-0.3707), and A2 0.3905 (0.2253).

A material A2 component supports a two-sextupole cascade. A material A1 component instead indicates one sextupole interaction combined with a fixed nonlinear source or strength-dependent feed-down. A nonzero A0 proves that sextupoles are not the only source. The next discriminating control is the nonlinear wiggler model.

## Lambda scan

| lambda2 | measured norm / nominal | pure lambda2^2 residual | constant + lambda2^2 residual | full quadratic residual | full-fit direction median [P10, P90] |
|---:|---:|---:|---:|---:|---:|
| 0.00 | 1.194982 | 1.194982 | 0.000000 | 0.001034 | 0.002090 [0.000332, 0.005458] |
| 0.25 | 1.105105 | 1.044904 | 0.094505 | 0.002074 | 0.004181 [0.000662, 0.010904] |
| 0.50 | 1.039896 | 0.797290 | 0.125362 | 0.000057 | 0.000119 [0.000017, 0.000351] |
| 0.75 | 1.002085 | 0.449482 | 0.093659 | 0.002053 | 0.004176 [0.000663, 0.010943] |
| 1.00 | 1.000000 | 0.000000 | 0.000000 | 0.001029 | 0.002090 [0.000333, 0.005472] |

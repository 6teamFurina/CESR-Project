# Exact finite-amplitude cross-response validation

The paired exact SciBmad check uses 5 selected targets, aligned and
misaligned machine states, signed 0.500 mm bumps, and
`delta K2 = +/-0.020 m^-3`.  The aligned K2-odd/bump-odd gradient is
subtracted from the matched misaligned gradient before comparison with the
nominal GTPSA alignment design.

- aggregate relative L2 residual: `3.267478e-02`;
- cosine similarity: `0.999474810`;
- fitted scale multiplying the GTPSA prediction:
  `1.004207573`;
- residual in the raw `+delta K2` minus `-delta K2` orbit difference at one
  signed bump, RMS:
  `0.598099 nm`;
- the same two-state residual, P90:
  `0.979322 nm`;
- the same two-state residual, maximum:
  `2.099392 nm`;
- full four-corner odd/odd contrast residual RMS:
  `1.196199 nm`;
- maximum target/axis block relative L2 residual:
  `5.780579e-02`.

This validates the compact source factorization only for the selected
finite-amplitude, single-target-offset cases.  It does not validate a
misaligned all-magnet background or measured-machine covariance.

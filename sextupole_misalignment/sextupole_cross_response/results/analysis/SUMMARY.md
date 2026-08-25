# All-sextupole GTPSA cross-response result

The latest repaired SciBmad lattice was linearized with GTPSA at the nominal
RF-on closed orbit.  The calculation covers all 76 active normal sextupoles
as both excitation and observation locations.

## Compact derivative construction

No descriptor containing all sextupole K2, bump, and offset parameters was
formed.  The response uses an order-1 GTPSA periodic local-kick map, one
first-derivative corrector calculation, and the exact local normal-sextupole
polynomial.  The resulting selected K2--bump--center derivative is saved as
`alignment_design.npy` with axes
`target, bump_axis, observation_sextupole, output_plane, center_axis`.

## Locality

- x-bump median target-only energy fraction:
  `6.024%`;
  median participation count:
  `15.04` sextupoles.
- y-bump median target-only energy fraction:
  `3.383%`;
  median participation count:
  `23.57` sextupoles.
- Per unit target bump, the median off-target radial-orbit RMS is
  `0.456` for x commands and
  `0.617` for y commands.  For a
  0.5 mm command these are `0.228`
  and `0.309 mm`; the largest
  individual off-target responses in the complete matrix are
  `1.871` and
  `2.708 mm`.
- Across the four K2--bump--center channels, the median target-only energy
  fraction ranges from `0.097%` to
  `0.377%`; the median participation count
  ranges from `31.52` to
  `33.94` sextupoles.

These matrices are not block-local under a target-only definition.  Small
absolute orbit at a distant location must not be interpreted as negligible
until the response is whitened by the intended measurement covariance.

## Effective rank

| matrix | numerical rank | effective rank | modes for 90% energy | modes for 99% energy |
|---|---:|---:|---:|---:|
| `bump_152x152` | 14 | 6.251 | 6 | 7 |
| `periodic_kick_152x152` | 152 | 12.071 | 11 | 43 |
| `sextupole_source_152x152` | 152 | 12.071 | 11 | 43 |
| `shared_alignment_template_304x152` | 152 | 22.451 | 20 | 68 |
| `separate_scan_block_design_23104x152` | 152 | 143.416 | 128 | 150 |


The shared `304 x 152` template matrix puts every target/center template into
one common channel coordinate system.  Its compact spectrum is evidence for a
candidate shared observation basis, not for a 20-dimensional joint inverse.
The physical design for separate one-target-at-a-time scans retains the target
axis and is block diagonal: it has full column rank 152, effective
rank `143.416`, and unwhitened condition number
`2.060`.  Every individual `304 x 2`
target block has rank two; its condition-number median / maximum is
`1.000000000000 / 1.000000000000`.
These Euclidean nominal values are structural checks, not experimental
position precision.  A truncated shared basis must be tested against
covariance-whitened center-recovery error rather than selected by response
energy alone.

## Interpretation boundary

The source factorization is a nominal, first-order periodic-response result.
It uses an integrated thin normal-sextupole source placed at the element entry.
Finite-length effects, finite bump/K2 amplitudes, misaligned-background
relinearization, BPM noise, and machine-operating limits are outside this
matrix and require separate exact SciBmad validation.

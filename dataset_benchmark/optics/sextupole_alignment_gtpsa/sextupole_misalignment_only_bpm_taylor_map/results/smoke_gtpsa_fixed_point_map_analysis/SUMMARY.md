# Direct SciBmad/GTPSA K2--offset map inverse

The map is generated directly by SciBmad/GTPSA on the nominal validated latest
lattice.  Measured quartic-in-K2 BPM slopes and BPM-predicted target-local
coordinates come from the common sextupole-misalignment-only exact scan.  The
zero-bump slope is subtracted from both measurement and map, so the inverse is
driven by bump-dependent K2 response rather than an absolute nominal baseline.

- targets / realizations: 1 / 1 per target
- descriptor: `Descriptor(6, 4, 3, 4)`
- maximum saved offset order: 3
- BPM channels: 222

| method | beam-relative 2D RMSE [um] | median [um] | P90 [um] | maximum [um] |
|---|---:|---:|---:|---:|
| gtpsa_k2_offset_order1_predicted_local | 229.879504 | 229.879504 | 229.879504 | 229.879504 |
| gtpsa_k2_offset_order2_predicted_local | 0.589930 | 0.589930 | 0.589930 | 0.589930 |
| gtpsa_k2_offset_order3_predicted_local | 1.123612 | 1.123612 | 1.123612 | 1.123612 |

These direct-map methods are more model-dependent than the scan-profiled
physical-source and empirical Taylor-surface methods.  Other 75 sextupole
offsets are hidden and are not used to condition the nominal GTPSA map.

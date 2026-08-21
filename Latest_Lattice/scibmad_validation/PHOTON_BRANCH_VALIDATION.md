# Latest CESR photon branch registry validation

All eleven Bmad photon fork targets are registered as independently queryable Beamlines branches and runnable branch-local paraxial ray lines.

| Bmad branch | Name | Elements | Length (m) | Length error (m) |
|---:|---|---:|---:|---:|
| 1 | `S4B_LINE` | 4 | 8.369 | 0.000e+00 |
| 2 | `S7A_LINE` | 9 | 32.795 | 0.000e+00 |
| 3 | `S7B1_S7B2_LINE` | 9 | 34.495 | 0.000e+00 |
| 4 | `S1A1_S1A2_S1A3_LINE` | 6 | 22.123 | 0.000e+00 |
| 5 | `S1A1_LINE` | 4 | 1 | 0.000e+00 |
| 6 | `S1A2_S1A3_LINE` | 3 | 1 | 0.000e+00 |
| 7 | `S2A_LINE` | 4 | 9.743 | 0.000e+00 |
| 8 | `S2B_LINE` | 4 | 8.369 | 0.000e+00 |
| 9 | `S3A_LINE` | 4 | 9.743 | 0.000e+00 |
| 10 | `S3B_LINE` | 4 | 8.369 | 0.000e+00 |
| 11 | `S4A_LINE` | 4 | 9.743 | 0.000e+00 |

`latest_photon_branch(name)` and `latest_photon_branch_for_fork(name)` provide lookup, while `track_latest_photon_branch(name; ray=...)` propagates a paraxial ray through the branch drifts. The reference ray remains on axis and accumulates the exact Bmad branch length.

The archived flat mirrors have `REF_TILT=-pi/2`, `GRAZE_ANGLE=0.004 rad`, and 10 keV reference energy. They remain identity markers in branch-local coordinates because the Bmad branch already defines the reflected reference frame. This interface does not model reflectivity, apertures, curvature, or off-axis mirror scattering; SciBmad/Beamlines currently has no photon `Fork`/`Mirror` tracker.

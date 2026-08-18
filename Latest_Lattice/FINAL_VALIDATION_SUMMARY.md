# Latest CESR SciBmad acceptance summary

Audit date: 2026-08-16

Source: archived `CHESS-U_6000MEV_20251020_S1` 6 GeV positron lattice.
Authoritative writer: Bmad/Tao `20260814-0`. Generated model:
`latest_cesr_scibmad_repaired.jl`.

| Capability | Result | Acceptance |
|---|---|---|
| Main-ring structure | 1,177 tracking elements; 768.4378690000005 m; exact order and lengths | Pass |
| Non-wiggler local linear maps | Maximum matrix difference `4.29e-13` | Pass |
| `ID_S1A` wiggler | Full-block affine/exit difference `7.12e-15`; residual `R12=5.89e-6 m` matches Bmad Runge-Kutta field tracking | Pass with documented model choice |
| RF-on closed orbit | Maximum starting-orbit difference `5.25e-14` | Pass |
| One-turn eigenphases | Magnitude differences about `5.9e-14`, `2.4e-9`, `5.5e-8` | Pass |
| Overlay/Group controls | All 124 lords, 347 relations, and 475 observations covered; median maximum relative difference `3.74e-4`, worst informative `1.84%` | Pass |
| Girders | All 12 girders, six parameters, 150 members/slices, and 972 observations covered | Pass for coherent geometry and offsets; qualified for pitch |
| Photon branches | All 11 fork targets registered; exact element counts and lengths; branch-local reference rays validated | Pass for branch geometry/reference rays |

## Known limits

1. Girder pitch through the strong curved DQX combined-function bends differs
   from Bmad by up to `21.7%`. Current SciBmad/BeamTracking uses straight
   multipoles in a curved reference system and does not provide Bmad's exact
   curved-coordinate misalignment field. Nominal DQX maps and coherent offsets
   are not affected by this qualification.
2. Photon fork selection is explicit. The branch-local ray helper does not
   model mirror reflectivity, apertures, curvature, or off-axis specular
   scattering because current SciBmad/Beamlines has no photon `Fork`/`Mirror`
   tracker. The reflected reference geometry and path lengths are preserved.

Within these two explicit library limits, the generated lattice is accepted as
functionally and numerically equivalent for normal CESR charged-particle
optics, tracking, correction/knob studies, and coherent girder-offset work.
By project policy recorded 2026-08-16, it is the default lattice for all new
CESR and ring-analysis work; the maintained `cesr.jl` is retained only for
historical reproduction unless the user explicitly requests it.

# Latest CESR nominal RF-on ring comparison

- Bmad/Tao reference: `20260814-0`, branch 0.
- SciBmad lattice: `latest_cesr_scibmad_repaired.jl`.
- Bmad starting closed orbit: `[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]`.
- SciBmad starting closed orbit: `[1.1091870886232994e-15, -4.611685348661217e-17, 3.539647287700555e-30, 2.5153541705864106e-29, -3.4658307095783547e-15, 5.250405540231082e-14]`.
- Maximum starting closed-orbit difference: `5.250405540231e-14`.
- Bmad transverse tunes from Twiss phase: `[0.5560211358176979, 0.6370675534129386]`.
- Bmad one-turn eigenphase magnitudes: `[0.03395450500571177, 0.3629324465870806, 0.4439788641849802]`.
- SciBmad signed tunes: `[-0.44397880948323115, -0.36293244415462744, -0.033954505005652794]`.
- SciBmad eigenphase magnitudes: `[0.033954505005652794, 0.36293244415462744, 0.44397880948323115]`.

The remaining nominal discrepancy is dominated by the continuous-field wiggler map versus Bmad's faster standard-matrix approximation; see `compare_wiggler_block.jl`.

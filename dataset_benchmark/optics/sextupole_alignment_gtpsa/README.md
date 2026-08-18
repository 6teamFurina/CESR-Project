# CESR sextupole magnetic-center studies

This directory separates two related but scientifically distinct experiments.

## `response_map/`

The completed per-sextupole GTPSA response-map study. It contains the three-
parameter `Kn2`, `x_offset`, and `y_offset` calculation for all 76 active
normal sextupoles, the independent finite-difference check, and the local
`2 x 1191` thin-SVD analysis.

## `targeted_bump_k2_inversion/`

The follow-up protocol experiment. For a target sextupole with known simulated
truth, it combines known two-plane local-orbit bumps with a symmetric target-
`K2` scan, reconstructs the target magnetic center, and compares the inferred
`x/y` offsets directly with truth. Its purpose is to decide the final scan
protocol, observable set, input representation, and whether a locally
conditioned inverse remains unbiased when the other sextupoles are also
misaligned.

The completed response-map coefficients are inputs to the new inverse
experiment; they are not themselves treated as independent training samples.

## `quadrupole_affinity/`

The nuisance-marginalized quadrupole-selection study. Its maintained version
uses the repaired CHESS-U 6 GeV SciBmad lattice, keeps 15 of 113 independent
active quadrupole knobs per target sextupole using exact scalar local
beta/phase leverage, and calculates the target mixed `Kn2-offset` response
under nominal and `Kn1 +/-` optics with batched GTPSA. The response dictionary
uses fixed and four-probe launch trajectories at 111 BPMs (1110 direct BPM
readings), rather than differentiated coupled-Twiss gauge quantities. The
information and precision heatmaps include 150 response columns from the
offsets of the other 75 sextupoles as explicit nuisance directions. Heatmap
columns are target sextupoles and rows are quadrupoles. This is a nominal-launch
pre-screen for the later bump-grid protocol, not a final machine-precision
result.

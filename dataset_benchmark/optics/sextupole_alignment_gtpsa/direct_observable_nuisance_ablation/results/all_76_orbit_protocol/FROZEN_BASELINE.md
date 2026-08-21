# Frozen oracle-local-orbit baseline

Status: **frozen on 2026-08-17**.

This directory is the maintained reference result for the idealized case in
which the transverse closed-orbit coordinates at the target sextupole are
known exactly.  It must be cited as the **oracle-local-orbit baseline**, not as
predicted CESR machine accuracy.

## Scientific question

Conditional on exact target-local orbit, how accurately can the common
thin-sextupole source fit recover the target magnetic center from an economical
closed-orbit scan when the other sextupoles and quadrupoles are unknown?

## Fixed protocol

- lattice: `Latest_Lattice/latest_cesr_scibmad_repaired.jl`;
- engine: exact scalar RF-on SciBmad closed orbit and tracking;
- targets: all 76 active normal sextupoles;
- latent realizations: 8 independently generated machines per target;
- observations: 111 BPM x/y closed-orbit coordinates;
- orbit intervention: five-point axial cross at nominal amplitude 0.5 mm;
- sextupole intervention: three outer K2 levels, `(-0.02, 0, +0.02) m^-3`;
- K1 commands: nominal only;
- target offset: independent uniform x/y values within +/-350 micrometers;
- other 75 sextupoles: independent 300 micrometer RMS x/y offsets;
- all 113 quadrupoles: independent physical strength errors within +/-1%;
- BPM noise, BPM offsets, and missing BPMs: absent;
- local-orbit uncertainty: absent.

There are 9,120 SciBmad states and 608 center fits.  SciBmad generation took
553.6 s in the recorded run.

## Frozen result

- aggregate two-dimensional RMSE: **5.870 micrometers**;
- realization median / P90 / P99: **3.664 / 9.395 / 17.897 micrometers**;
- maximum realization error: **25.274 micrometers**;
- per-target RMSE median / P90 / maximum:
  **4.300 / 7.933 / 17.574 micrometers**.

## Oracle boundary

The fit in `analyze_all_targets_orbit_protocol.py` reads
`target_orbits.npy` and uses the exact target x/y coordinates at nominal K2 as
the five bump coordinates.  A real machine does not directly measure these
coordinates.  It measures the BPM tensor in `bpm_orbits.npy` plus the applied
magnet commands.  Therefore the 5.870 micrometer result isolates the
sextupole-response inverse after removing local-orbit reconstruction error.

The saved `target_orbits.npy` may be used as ground truth to evaluate a future
BPM-to-local-orbit estimator.  It must not be passed to a finite-BPM fit as an
input feature.

## Artifact roles

- `bpm_orbits.npy`: observable closed-orbit tensor;
- `bump_points.csv`: known commanded scan coordinates;
- `target_orbits.npy`: oracle/evaluation-only internal coordinates;
- `target_truth.npy`: target offset label, evaluation only;
- `latent_sextupole_offsets.npy`: hidden nuisance truth, audit only;
- `latent_quadrupole_relative_errors.npy`: hidden nuisance truth, audit only;
- `per_realization_fits.csv` and `per_target_summary.csv`: frozen oracle fits;
- `scan_metadata.toml`: complete generation settings and provenance.

## Successor study

New finite-BPM work belongs in the sibling directory
`finite_bpm_inversion/`.  It must report beam-relative center recovery
separately from absolute mechanical-offset recovery and must always identify
whether internal target orbit was used for fitting or only for evaluation.


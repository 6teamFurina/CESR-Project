# Nonlinear-rho RF-off optics benchmark

This experiment evaluates the same 9,001 corrector states used by the
nonlinear-rho closed-orbit benchmark with the same optics data product as the
earlier 1,000-sample chromatic-optics study.

The shared inputs are not regenerated here. They remain under
`../../orbit/Orbit_Calculation/nonlinear_rho_benchmark/shared_input/`:

- one zero-control baseline;
- scenarios `all`, `horizontal`, and `vertical`;
- radii `1.13`, `3.2`, `4.53`, `6.4`, and `9.05`;
- 600 fixed unit-RMS Gaussian directions per scenario and radius;
- base kick `5e-6 rad`, giving active-control RMS kicks from `5.65e-6` to
  `4.525e-5 rad`.

## Observable definition

The calculation is the same RF-off, four-dimensional coasting periodic optics
used in `../README.md`. At all 99 zero-length `DET_*` markers it saves normal
mode phase, beta and alpha; longitudinal phase; coupled-optics fields; all six
orbit coordinates; and first derivatives with respect to relative momentum
deviation. Ring outputs contain the two tunes, chromaticities, and coasting
slip coefficient.

The SciBmad runner separately records RF-on initial-orbit convergence, RF-off
coasting-orbit convergence, and Twiss convergence. The Bmad runner records a
paired RF-off optics status and transverse closure norm. Both runners
checkpoint every original `(scenario, rho)` cell, so an unstable sample cannot
discard earlier cells.

## Run

From the CESR project root on Windows:

```powershell
$env:JULIA_PKG_PRECOMPILE_AUTO='0'
julia --project=. dataset_benchmark/optics/nonlinear_rho_benchmark/run_scibmad_nonlinear_rho_optics.jl `
  --max-groups=16 `
  --output-dir=dataset_benchmark/optics/nonlinear_rho_benchmark/results/full_9001/scibmad
```

Run Bmad in the local `Ubuntu-Bmad` WSL environment:

```powershell
wsl.exe -d Ubuntu-Bmad -- /home/joeyfurina/miniforge3/envs/bmad/bin/python `
  "/mnt/d/Ring_Design_Development/CESR Project/dataset_benchmark/optics/nonlinear_rho_benchmark/run_bmad_nonlinear_rho_optics.py" `
  --max-groups=16 `
  --output-dir "/mnt/d/Ring_Design_Development/CESR Project/dataset_benchmark/optics/nonlinear_rho_benchmark/results/full_9001/bmad"
```

After both runs finish, generate the convergence and per-category accuracy
report:

```powershell
wsl.exe -d Ubuntu-Bmad -- /home/joeyfurina/miniforge3/envs/bmad/bin/python `
  "/mnt/d/Ring_Design_Development/CESR Project/dataset_benchmark/optics/nonlinear_rho_benchmark/compare_nonlinear_rho_optics.py" `
  --scibmad-dir "/mnt/d/Ring_Design_Development/CESR Project/dataset_benchmark/optics/nonlinear_rho_benchmark/results/full_9001/scibmad" `
  --bmad-dir "/mnt/d/Ring_Design_Development/CESR Project/dataset_benchmark/optics/nonlinear_rho_benchmark/results/full_9001/bmad" `
  --output-dir "/mnt/d/Ring_Design_Development/CESR Project/dataset_benchmark/optics/nonlinear_rho_benchmark/results/full_9001/comparison"
```

The comparison removes detector phase origins per sample, sign-aligns the two
documented longitudinal Bmad conventions, and reports both absolute
cross-engine differences and input responses after subtracting each engine's
own zero-control baseline. Accuracy is kept separate for ordinary Twiss,
coupled optics, chromatic derivatives, orbit-derived fields, ring tunes, and
ring chromatic quantities.

## Recorded result (2026-08-10)

Both engines produced valid RF-off optics for all 9,001 states. SciBmad's
RF-on initial orbit, RF-off coasting orbit, and exact pointwise Twiss stages
were each `9001/9001`; Bmad/Tao was `9001/9001`. The largest SciBmad RF-on and
coasting closure norms were `9.999e-11` and `9.906e-11`, respectively. Bmad's
largest transverse closure norm was `4.141e-12`.

Every successful sample has exactly 99 unique detector rows and one ring row,
and all compared numeric fields are finite. The complete result tree is about
1.205 GB. Stable calculation times are recorded in metadata, but the two full
runs overlapped for part of their execution and must not be used as a formal
speed comparison.

After subtracting each engine's own zero-input baseline, the median response
NRMSE across all 15 nonzero cells is:

| Category | Median response NRMSE | Important qualification |
|---|---:|---|
| ordinary Twiss | 2.115% | worst 11.086% for weak vertical-only `phi_2` response at rho 9.05 |
| coupled optics | 1.803% | horizontal-only coupling responses are weak and give 17--19% relative errors |
| chromatic derivatives | 2.064% | worst 7.632% for vertical-only `dphi_2/ddelta` at rho 9.05 |
| orbit-derived | 0.203% | dominated in the tail by the known anomalous `dorbit_z_ddelta` comparison |
| ring tunes | 2.061% | vertical-only `Qy` reaches 9.840% at rho 9.05 |
| ring chromatic quantities | 1.929% | vertical-only `xi_2` reaches 5.895% at rho 9.05 |

For `all` inputs, the ordinary-Twiss median response NRMSE grows from 1.884%
at rho 1.13 to 2.936% at rho 9.05; the chromatic-derivative median grows from
1.795% to 2.493%. Correlations for the regular non-orbit quantities remain
high even at the largest radius. The longitudinal `dorbit_z_ddelta` field is a
separate exception: it has approximately 101--213% response NRMSE depending
on cell and can have near-zero correlation. This is consistent with the
isolated longitudinal derivative discrepancy already seen in the 1,000-sample
benchmark and must not be used to characterize transverse optics accuracy.

These results establish solver survival over the tested nonlinear-rho domain,
not agreement with the real machine. The two codes use independently
represented CESR lattices, and relative errors become large for weak response
channels. A dataset should therefore retain per-field scales, engine
provenance, closure diagnostics, and field-specific validation thresholds.
No source-attribution study is part of this experiment.

The full generated report is
`results/full_9001/comparison/RESULTS.md`; machine-readable convergence,
per-cell/per-quantity, and category summaries are stored beside it.

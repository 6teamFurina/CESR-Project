# Targeted bump-by-K2 sextupole offset inversion

This experiment performs a closed synthetic inverse test with known truth:

```text
known sextupole offset
  -> corrector-generated local orbit bumps
  -> exact target-K2 scans
  -> full-ring orbit/optics slopes
  -> physical response-map/feeddown inversion
  -> direct truth comparison
```

The full protocol and pass/fail criteria are in `EXPERIMENT_DESIGN.md`.

## Exact smoke workflow

Run from the `CESR Project` directory.

Build the nominal local-bump knobs:

```powershell
julia --project=. dataset_benchmark/optics/sextupole_alignment_gtpsa/targeted_bump_k2_inversion/build_local_bump_knobs.jl `
  --targets=SEX_08W
```

Generate a five-bump, three-`K2` exact scan with known truth:

```powershell
julia --project=. dataset_benchmark/optics/sextupole_alignment_gtpsa/targeted_bump_k2_inversion/run_exact_scan.jl `
  --target=SEX_08W `
  --true-x-offset-m=3.5e-4 `
  --true-y-offset-m=-2.5e-4
```

Run the response-map/feeddown inverse and compare with truth:

```powershell
python dataset_benchmark/optics/sextupole_alignment_gtpsa/targeted_bump_k2_inversion/invert_scan.py
```

The default smoke case is deliberately noise-free and has no background
sextupole errors. It tests the end-to-end convention, signs, observation
alignment, bump realization, and inverse plumbing. Passing it does not
establish experimental precision. The next runs add all-sextupole offset
backgrounds and measured-style covariance before expanding to the ten-magnet
protocol screen.

Initial truth-recovery smoke results are recorded in `results/RESULTS.md`.

Validate both maintained smoke result sets with:

```powershell
python dataset_benchmark/optics/sextupole_alignment_gtpsa/targeted_bump_k2_inversion/validate_smoke_results.py
```

## Paired P0--P3 benchmark

The maintained paired comparison uses the same saved `SEX_08W` scan and the
same realization of the other 75 sextupole offsets for every inverse. Generate
the background-conditioned P1 and local-source P2 dictionaries, run P0--P2,
then run the small exact P3 inverse:

```powershell
julia --project=. dataset_benchmark/optics/sextupole_alignment_gtpsa/targeted_bump_k2_inversion/generate_conditioned_inverse_models.jl

& 'C:\Users\JoeyN\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' `
  dataset_benchmark/optics/sextupole_alignment_gtpsa/targeted_bump_k2_inversion/benchmark_physical_inverses.py

julia --project=. dataset_benchmark/optics/sextupole_alignment_gtpsa/targeted_bump_k2_inversion/run_exact_p3_inverse.jl

& 'C:\Users\JoeyN\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' `
  dataset_benchmark/optics/sextupole_alignment_gtpsa/targeted_bump_k2_inversion/validate_paired_benchmark.py
```

P1--P3 currently know the saved offsets of the other 75 sextupoles. They are
oracle-background diagnostics that isolate response conditioning and nonlinear
truncation; they are not yet operational estimators for unknown nuisance
states. Results and interpretation are in `results/PAIRED_P0_P3_RESULTS.md`.

The follow-up unknown-background test builds nominal P1/P2b dictionaries while
hiding the saved offsets of the other 75 sextupoles, and adds a quadratic
target-offset fit in the reconstructed four-source space:

```powershell
julia --project=. dataset_benchmark/optics/sextupole_alignment_gtpsa/targeted_bump_k2_inversion/generate_conditioned_inverse_models.jl `
  --background-mode=nominal `
  --nonlinear-calibration=true `
  --output-dir=dataset_benchmark/optics/sextupole_alignment_gtpsa/targeted_bump_k2_inversion/results/unknown_background_benchmark/nominal_models

& 'C:\Users\JoeyN\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' `
  dataset_benchmark/optics/sextupole_alignment_gtpsa/targeted_bump_k2_inversion/benchmark_nonlinear_p2b.py
```

The nonlinear second stage improves the matched-background closure case but
does not repair the unknown-background source bias. See
`results/UNKNOWN_BACKGROUND_NONLINEAR_P2B_RESULTS.md`.

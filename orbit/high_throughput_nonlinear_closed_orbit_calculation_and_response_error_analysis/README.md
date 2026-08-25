# High-throughput nonlinear closed-orbit paper: latest CESR ring

Status: `two_of_three_production_threads_complete`; the latest-ring rho sweep
and complete-element-family attribution were completed on 2026-08-22. The
normal-sextupole physical beta/phase production thread and manuscript rebuild
have not been run.

## Implemented latest-ring experiment threads

All three retained SciBmad threads now run against the latest repaired CESR
lattice without fixed control, detector, element, family, or sextupole counts.
The checked integration artifacts establish code-path and reconstruction
closure only:

- the rho workflow converged 7/7 states with no fallback for
  `rho={0,0.1}` and two directions per nonzero all/H/V scenario;
- complete-element sourcing covered all 1,177 elements and closed the summed
  x/y nonlinear targets to `4.23e-15` and `1.28e-14` in one-direction smokes;
- the two-direction normal-sextupole workflow covered 76 active normal
  sextupoles and closed the all-element x/y target to `1.38e-14` before
  selecting the sextupole rows;
- matching nominal and direction-dependent SciBmad optics were exported for
  the sextupole source--beta--phase/Green-function ranking analysis.

The completed paper-scoped production results are under
[`chapter_01_nonlinear_response_rho_sweep`](../chapter_01_nonlinear_response_rho_sweep/README.md)
and
[`chapter_02_lattice_element_attribution`](../chapter_02_lattice_element_attribution/README.md).
The rho run contains 36,001 unique states with 36,000 converged; the complete-
element x/y ensembles contain 100 directions each and close to `1.31e-14`.
The maintained experiment descriptions and checked smoke artifacts remain under
[`response_rho_sweep_600`](../error_analysis/response_rho_sweep_600/README.md),
[`thick_element_sextupole_sourcing`](../error_analysis/thick_element_sextupole_sourcing/README.md),
[`sextupole_detector_contributions`](../error_analysis/sextupole_detector_contributions/README.md),
and
[`sextupole_beta_phase_correlation`](../error_analysis/sextupole_beta_phase_correlation/README.md).
Run any individual thread, or all three, through the ring-scoped entry point:

```powershell
.\run_latest_cesr_experiments.ps1 -Mode smoke -Experiment all -Ring latest
.\run_latest_cesr_experiments.ps1 -Mode production -Experiment all -Ring latest
```

The launcher defaults to GTPSA for the first-order response and never invokes
Bmad/Tao. `central-difference` remains an explicit bounded SciBmad validation
choice. Ring-specific selectors and artifact paths live in the shared adapter
and configuration; experiment dimensions are discovered from the selected
model at runtime.

The read-only cross-thread validator checks finite tables, latest-lattice and
SciBmad provenance, dynamic registry agreement, element/family closure, and
the absence of paper-facing block-share/third-order fields:

```powershell
& 'C:\Users\JoeyN\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' `
  .\validate_latest_cesr_results.py
```

Its current expected result remains `OVERALL: MISSING_PRODUCTION` because the
normal-sextupole contribution and beta/phase production destinations are still
absent. The rho and complete-element chapter destinations must validate as
production, and all smoke checks must pass. Add `--require-production` only
after the remaining two destinations are generated.

## Frozen latest-ring paper scope

The required SciBmad production results are now limited to:

1. nonlinear orbit-response error versus normalized corrector amplitude
   `rho`;
2. complete-element-family attribution of the nonlinear error;
3. normal-sextupole element and physical attribution of that error.

Detailed `hh`/`hv`/`vv` response-block shares and their input-space
interpretation are not required. The vertical signed-parity/cubic analysis is
also outside this paper. Those retained studies are follow-up material rather
than completion gates. The attribution may use a direction-contracted Hessian
internally in the verified quadratic interval, but the manuscript should
report the summed nonlinear target and its family/sextupole attribution rather
than the omitted block decomposition.

The default model for a new calculation is
[`Latest_Lattice/latest_cesr_scibmad_repaired.jl`](../../Latest_Lattice/latest_cesr_scibmad_repaired.jl).
Latest data, generated figures, and any rebuilt paper assets should be kept in
ring-scoped locations such as `figures/latest_cesr/` and a linked
`results/latest_cesr/` directory. Every quantitative table and figure must
identify the lattice, branch, control/observable registries, RF state, solver
configuration, and SciBmad provenance.

The checked-in `.tex`, PDF, and `figures/` files currently describe the old
CESR export. They are not silently reclassified as latest results. Before a
new paper build, regenerate all numbers and figures with SciBmad, then record
the status and paths here. Bmad/Tao may appear only as an explicitly labeled
cross-code validation and may be added after the SciBmad result threads are
complete. Until then, cross-code agreement and speedup claims remain pending.

The previous manuscript workflow and result description are preserved in
[`README_archived.md`](README_archived.md). The latest-lattice curved-DQX
girder-pitch limitation must be stated if the paper exercises girder pitch.

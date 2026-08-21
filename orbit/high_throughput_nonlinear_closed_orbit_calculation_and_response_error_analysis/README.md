# High-throughput nonlinear closed-orbit paper: latest CESR ring

Status: `implementation_smoke_complete`; the latest-ring production statistics
and manuscript rebuild have not been run.

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

These small runs are not the 600-direction rho result or the 100-direction
attribution ensembles, so their family percentages and correlations are not
paper conclusions. The maintained experiment descriptions and checked smoke
artifacts are under
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

Its current expected result is `OVERALL: MISSING_PRODUCTION` with every smoke
check passing. Add `--require-production` in a paper-build gate; until all
four ring-scoped production destinations exist, that mode returns exit code 2.

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

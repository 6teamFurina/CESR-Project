# Sextupole workspace cleanup audit

Date: 2026-09-05. Status: the first cleanup and the first three archive
recommendations were authorized by the user and completed. The isolated
Taylor-map benchmark, K1 triplet/exact-11 study, and earlier sequential/ORM
comparisons remain in place.

## Finding

Most historical method directories also contain useful scientific controls or
shared production code. The completed first cleanup removed duplicate historical
methods, superseded development
outputs, and regenerable caches: **303 files, 98.73 MiB of logical file content**
across 24 explicitly listed paths. Allocated disk space may differ.

The original reviewed paths, file counts, byte totals, and reasons are in
[`cleanup_candidates_2026-09-05.tsv`](cleanup_candidates_2026-09-05.tsv).
The per-path completion status and UTC timestamps are in
[`cleanup_completed_2026-09-05.tsv`](cleanup_completed_2026-09-05.tsv).

Immediately before removal, all 24 paths still matched the audited counts and
sizes, and the two historical copies were compared with their archives again by
direct byte equality. Deletion was restricted to those exact paths. Immediately
afterward, every other existing file in this study still existed with unchanged
size and modification time. Both canonical archives and both protected GTPSA
smoke fixtures were retained.

The subsequent three-study archive relocation is recorded separately in
[`archived_methods/archive_moves_2026-09-05.tsv`](archived_methods/archive_moves_2026-09-05.tsv).
It moved 139 files in 11 source/destination paths; it did not delete their
historical data. See the [archive index](archived_methods/README.md).

## Evidence and scope

- Read the project-memory research context, particularly the sextupole archive,
  K1-scan retirement, finite-BPM, nominal-GTPSA, and burst-acquisition decisions.
  Checked the project/study READMEs against source code and saved metadata.
- Inspected the sextupole workspace and searched workspace source/documentation
  for imports, Julia includes, default input paths, and candidate references.
  Inspected metadata references between development result directories as well.
- Compared the two historical directory pairs by relative filename and direct
  byte equality. No cryptographic checksums were calculated.
- The initial audit used static dependencies and provenance. After deletion,
  four existing validators passed: full-error state-space inversion, burst-size
  sweep, the formal excitation-validity envelope, and the isolated sextupole
  Taylor-map benchmark including the fixed-point-versus-Twiss smoke check. No
  new lattice scan was generated. Storage totals are filesystem inventory
  measurements, not accelerator results.
- The bundled Python lacked SciPy/Matplotlib for the two finite-BPM validators;
  those passed using the existing `cesr-analysis` environment at
  `/Users/furinaovo/miniforge3/envs/cesr-analysis/bin/python`. The envelope and
  isolated Taylor-map validators passed with the bundled Python. All runs used
  `-B` to avoid recreating the deleted bytecode caches. The two finite-BPM
  `VALIDATION.json` reports were refreshed; the envelope validation report was
  written outside the study at `/tmp/cesr-cleanup-envelope-validation.json`.
- No Git repository was detected at the workspace root or its parents. Do not
  assume that removing these files can be undone with Git.
- A missing static reference is not proof that an external notebook or manual
  command never uses a path. Recommendations below are scoped to the inspected
  workspace and the currently recorded production workflow.

## Completed first cleanup

The following table records the pre-removal evidence and the actions that have
now been completed. The original candidate paths no longer exist.

| Paths relative to this directory | Size (MiB) | Evidence and disposition |
| --- | ---: | --- |
| `response_map/` | 73.82 | All 29 filenames also exist under `archived_methods/response_map/`; 27 files are byte-identical. The differences are the archive notice/commands in README and the extra parent-directory traversal in the Julia entry point. Remove the top-level copy; retain the archive. |
| `targeted_bump_k2_inversion/` | 9.19 | All 63 filenames also exist in the archive; 61 files are byte-identical. Only README and `common.jl` differ, reflecting the archive location. Remove the top-level copy; retain the archive. |
| `sequential_joint_inverse/results/{smoke,full_target_smoke,analysis_smoke,analysis_smoke_output,analysis_smoke_output2}/` | 10.24 | Small development datasets and two analysis iterations. Their metadata describe 1-machine/1-target, 1-machine/76-target, or 4-machine/4-target runs. Maintained analysis uses `exact_joint_machines/`. The two analysis-output metadata files refer to `analysis_smoke`; remove this group together. |
| `sextupole_excitation_validity_envelope/results/pilot_*` and `results/smoke_*`, limited to the 11 manifest entries | 4.16 | Single-target development grids and intermediate map/analysis outputs. The maintained generator/analyzer/validator use `exact_validation/`, `gtpsa_maps/`, and `analysis/`. No resolved production dependency on these pilot/smoke paths was found. |
| `direct_observable_nuisance_ablation/results/smoke/` | 0.90 | One-realization preliminary direct-observable scan. The documented paired pilot has eight realizations; the all-target frozen baseline is separate. No maintained consumer of this smoke directory was found. |
| `.DS_Store` and the four listed `__pycache__/` directories | 0.43 | Finder metadata and regenerable Python bytecode. Preserve all corresponding source files. |

The top-level historical copies had **88 byte-identical files out of 92**.
Every numerical result file in those copies matched its archived counterpart.
The remaining four differences are relocation/documentation changes, not unique
scientific results. The archived copy is already the designated canonical
historical location in [project memory](../PROJECT_MEMORY.md).

The one external historical response consumer found,
[`validate_nominal_responses.py`](archived_methods/bmad_quadrupole_affinity/validate_nominal_responses.py),
already defaults to `archived_methods/response_map/results/full/alignment_coefficients.csv`.
It does not require the redundant top-level response directory.

## Retired production methods: archive selectively, retain the evidence

| Method or experiment | Current role | Recommendation |
| --- | --- | --- |
| Historical P0–P3/P1–P2 response-dictionary inverses | Older-lattice provenance and conditioning diagnostics | Keep one copy under `archived_methods/`; delete only its redundant top-level copies in the first cleanup. |
| Bmad/Tao quadrupole-affinity screen | Explicit historical cross-code reference; superseded as the primary affinity result by latest-lattice SciBmad | The old `quadrupole_affinity/results/responses/` and `results/affinity/` total 94.75 MiB. Consider a separate reference archive with `generate_bmad_affinity_responses.py` and `validate_nominal_responses.py`. Do not remove these results while expecting that comparison validator to work. Retain `results/scibmad_latest/`. |
| Greedy K1 selection, exact 11-condition triplets, and the seven-condition K1 scan | Design-stage experiments; K1 scanning was excluded from the current orbit-only production protocol | Keep the equal-measurement-budget comparison and negative-result summaries. The old scan outputs can be moved to a research archive after their consumers are redirected. Do not delete the entire affinity tree: it contains shared Julia code and current bump knobs. |
| Exact target-local-orbit estimator in `direct_observable_nuisance_ablation/` | Frozen oracle baseline; superseded as a machine-facing inverse by finite-BPM reconstruction | Retain the baseline and raw scan tensor. Earlier finite-BPM comparisons still reuse it, and it isolates transport error from center-model error. |
| Unfiltered 15-state sequential BPM inverse and error-conditioned/noisy-ORM correction comparisons | Superseded production settings, but useful failure/control experiments | Retain comparison summaries and the code currently imported by the state-space implementation. Do not remove the sequential inverse or generic correction generator wholesale. |
| Blocked/interleaved acquisition and isolated/compound nuisance studies | Earlier mechanism and ablation controls | Keep the decisive results explaining parity/drift handling and nuisance sensitivity. The modern full-error result does not replace these controlled comparisons. |
| Single-turn visits (`B=1`) | Historical control and regression baseline for the burst study | Keep the baseline arrays and validation. Project memory explicitly retires it as the operational candidate while retaining its control role. |
| Direct high-order GTPSA offset maps | Limited validated subset plus an unresolved failure boundary | Keep compact success/failure evidence and regression fixtures. A failed method can still establish a capability limitation that the current low-order approach does not resolve. |

The initial table above records the pre-migration assessment. The user later
approved the first three archive groups below; project memory now records that
organization decision. The physics conclusions and remaining study scopes are
unchanged.

## Archive decisions and remaining candidates

Items 1–3 are complete. Items 4–6 were not relocated. Archive migrations retain
code, metadata, decisive results, and validators together, with relative paths
and README links repaired.

1. **Completed: old Bmad/Tao affinity reference.** The two historical generator
   and comparison scripts and `quadrupole_affinity/results/{responses,affinity}/`
   now live in `archived_methods/bmad_quadrupole_affinity/`. The latest SciBmad
   response generator, analysis/selection utilities, and `results/scibmad_latest/`
   remain in their existing locations.
2. **Completed: early direct-observable and K1 pilots.** Moved
   `generate_paired_scan.jl`, `analyze_paired_scan.py`,
   `generate_k1_orbit_ablation.jl`, `analyze_k1_orbit_ablation.py`, and
   `results/{sex_09aw_paired_pilot,sex_09aw_k1_orbit_ablation}/` from
   `direct_observable_nuisance_ablation/` into
   `archived_methods/direct_observable_k1_pilots/`. The observable-selection and
   K1 negative-result narrative moved into its README.
   `analyze_protocol_subsampling.py` remains in the maintained directory because
   the all-target analyzer imports its fitting functions; its CLI default now
   reads the archived pilot tensor. The all-target frozen tensor and its
   generator/analyzer remain in place for finite-BPM baselines.
3. **Completed: blocked/interleaved protocol comparison.** The whole study now
   lives in `archived_methods/interleaved_measurement_protocol/`. Its imports and
   input paths resolve to the maintained `finite_bpm_inversion/` and
   `real_machine_nuisance_ablation/` studies. The result tables and drift/averaging
   evidence were preserved, and its existing validator passes from the archive.
4. **Isolated sextupole Taylor-map comparison.** The whole
   `sextupole_misalignment_only_bpm_taylor_map/` study can be organized as a
   retained reference benchmark rather than an active production branch. No
   incoming production-code reference was found. It has outgoing dependencies
   on finite-BPM helpers, the nominal local-orbit model, and affinity bump knobs;
   update those paths on migration. Keep all its compact GTPSA regression
   fixtures and its high-order-map limitation evidence together.
5. **K1 triplet selection and exact-grid studies, after splitting shared code.**
   Archive the selection/11-condition experiment runners, analyses, and
   `exact_11_triplet_validation/results/{scans,scans_k1_1pct,aggregate,aggregate_k1_1pct,k1_1pct_optics_screen,k1_1pct_safe_selection,logs}/`
   as one reproducible design-stage package. Preserve `common.jl`, the shared
   latest-lattice generator, and `results/bump_knobs/` as active dependencies, or
   relocate them into a common module first. This is not a whole-directory move.
6. **Earlier finite-BPM and ORM comparison outputs, after redirecting consumers.**
   The unfiltered sequential result and the finite-difference/error-conditioned
   correction comparisons can be separated from the current nominal-GTPSA
   production result. Their historical scan cases and summaries remain inputs
   to dedicated validators/comparison scripts. Shared inverse and correction
   modules must stay active or be extracted before their experiment wrappers
   are archived. This is a lower-priority refactor, not part of this cleanup.

## Important deletion traps

1. **Some smoke data are regression fixtures.**
   [`sextupole_misalignment_only_bpm_taylor_map/validate_results.py`](sextupole_misalignment_only_bpm_taylor_map/validate_results.py)
   checks `smoke_gtpsa_maps/` against `smoke_gtpsa_fixed_point_maps/` at lines
   189–207. The check is conditional on both directories existing: deleting
   either silently removes that coverage rather than necessarily failing the
   validator. Keep both and their compact supporting diagnostic outputs. The
   entire Taylor-map smoke family is excluded from the first-cleanup manifest.
2. **Affinity is also a shared library location.**
   `finite_bpm_inversion/generate_local_orbit_models.jl`,
   `sequential_joint_inverse/generate_joint_machine_scans.jl`, and
   `real_machine_nuisance_ablation/generate_physical_nuisance_scans.jl` include
   `quadrupole_affinity/exact_11_triplet_validation/common.jl`. That file includes
   `quadrupole_affinity/generate_scibmad_affinity_responses.jl`. The state-space
   inverse also reads `exact_11_triplet_validation/results/bump_knobs/`.
3. **Earlier inverse scripts supply current functions.**
   `finite_bpm_inversion/analyze_state_space_bpm_gtpsa_inverse.py` imports
   `analyze_sequential_bpm_gtpsa_inverse.py` and
   `gtpsa_derivative_stochastic_inverse/analyze_stochastic_inverse.py`.
   The sequential module imports `build_two_sided_maps` from
   `analyze_local_orbit_predictors.py`. Retiring their original experiments
   does not make these modules unused.
4. **The correction entry points share implementation.**
   `generate_gtpsa_nominal_corrected_joint_machine_scans.jl` includes
   `generate_corrected_joint_machine_scans.jl`; the latter includes
   `run_quadrupole_orbit_correction.jl` and `gtpsa_noisy_response.jl`.
   The generic generator now defaults to nominal GTPSA. Its historical name
   and optional finite-difference branch are not grounds for deleting it.

## Keep in the active workspace

- Nominal-GTPSA correction code and the full-error source case
  `sequential_joint_inverse/results/exact_joint_machines/with_all_errors_gtpsa_nominal_corrected/`,
  including its shared machine metadata and paired latent evaluation truth.
- `finite_bpm_inversion/results/local_orbit_model/`, the state-space inverse,
  `results/state_space_sequential_bpm_gtpsa_inverse/`, and
  `results/burst_size_sweep/`, with their validators.
- `finite_bpm_inversion/results/runtime_breakdown_20260905/`: a current runtime
  investigation, not an obsolete folder merely because it contains probes/logs.
- `sextupole_cross_response/`: relevant to the recorded next multiplexing study.
- The maintained `sextupole_excitation_validity_envelope/` results and validator:
  `exact_validation/`, `gtpsa_maps/`, and `analysis/`.
- Shared GTPSA templates, two-sided transport code, latest-lattice bump knobs,
  and retained scientific comparison/negative-result evidence described above.

The larger organizational issue is that shared helpers live inside experiment
directories. Before archiving entire retired methods, extract those helpers
into a common module, redirect imports/data paths, and run the affected existing
validators. The first-cleanup manifest avoids that refactor and leaves all
existing validators and production source code in place except for the two
redundant copies whose canonical archived versions are preserved.

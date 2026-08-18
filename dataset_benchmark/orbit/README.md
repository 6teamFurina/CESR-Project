# CESR orbit dataset studies

The orbit work is organized into calculation, analysis, and manuscript work
packages with shared references:

- [`Orbit_Calculation/`](Orbit_Calculation/): corrector input generation,
  SciBmad and Bmad/Tao orbit-dataset runners, comparison code, archived inputs,
  and benchmark results.
- [`error_analysis/`](error_analysis/): response-radius sweeps, signed-parity
  experiments, scaling analysis, figures, and their results.
- [`high_throughput_nonlinear_closed_orbit_calculation_and_response_error_analysis/`](high_throughput_nonlinear_closed_orbit_calculation_and_response_error_analysis/):
  the calculation and second-order physical-attribution paper, named after its
  title rather than a conference year.
- [`hessian_svd_nonlinear_closed_orbit_correction/`](hessian_svd_nonlinear_closed_orbit_correction/):
  the companion correction-paper scaffold.  It currently contains the
  transferred higher-order parity/cubic validity analysis and placeholders for
  the Hessian-SVD method and correction experiments.
- [`reference/`](reference/): the shared Bmad-compatible lattice and validated
  closed-orbit response cache used by both work packages.

Each work package has its own README with runnable commands and current
technical results.

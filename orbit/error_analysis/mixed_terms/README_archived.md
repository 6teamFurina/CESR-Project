# GTPSA all-corrector horizontal--vertical mixed-term response

The final quadratic-block values are calculated with two second-order GTPSA
corrector parameters per fixed H/V direction pair and implicit
differentiation of the RF-on closed-orbit fixed point. Results are reported as
`median [P10, P90]` across directions.

The earlier four-sign experiment is retained as an independent validation of
the GTPSA normalization, block signs, and vector reconstruction.

For each fixed pair of horizontal and vertical Gaussian directions, each
normalized to exact unit RMS within its own corrector family, the runner uses
the same directions at every radius and evaluates pure `+/-h`, pure `+/-v`, and
all four joint sign states. It extracts `Q_hh`, `Q_vv`, and `Q_hv` with signed
finite differences and reconstructs the four nonlinear joint residual vectors
before taking detector RMS norms.

To reproduce the adopted result from `CESR Project`, run:

```console
julia --project=. orbit/error_analysis/mixed_terms/run_mixed_term_gtpsa.jl
```

To reproduce the four-sign validation, run:

```powershell
julia --project=. orbit/error_analysis/mixed_terms/run_mixed_term_experiment.jl
python orbit/error_analysis/mixed_terms/analyze_mixed_terms.py
```

The default production experiment uses 100 direction pairs and 8 radii from
`ρ=0.1` to `ρ=1.13`, where `ρ=1` is 5 microrad RMS in each corrector
family. The simultaneous all-corrector vector also has exactly this global RMS.

The main quantitative tests are:

- quadratic scaling of `Q_hv`;
- direction-resolved mixed-to-pure ratios and mixed squared-norm share;
- vector-level reconstruction of all four joint signs;
- error reduction when `Q_hv` is added to a pure-block reconstruction;
- closed-orbit convergence and closure norms relative to the measured signal.

## Adopted GTPSA direction contraction

`run_mixed_term_gtpsa.jl` replaces the corrector finite differences with two
second-order GTPSA parameters for every fixed H/V direction pair.  It obtains
the RF-on closed-orbit derivatives by implicit differentiation of the
one-turn fixed-point equation and directly returns `Q_hh`, `Q_hv`, and
`Q_vv`.

The committed 100-direction result and the finite-difference comparison are
in [`gtpsa_results/GTPSA_RESULTS.md`](gtpsa_results/GTPSA_RESULTS.md).
The adopted reporting convention is `median [P10, P90]`; means remain in the
CSV for secondary diagnostics. For X, the final `f_vv` result is
`90.83% [61.72%, 98.15%]`.

For each orbit component, the reported `hh`, `hv`, and `vv` squared-norm
shares use a common three-block denominator and are summarized as the median
with a P10--P90 interval across the 100 fixed directions. The adopted GTPSA
statistics are in
[`gtpsa_results/gtpsa_summary.csv`](gtpsa_results/gtpsa_summary.csv); the
four-sign reconstruction diagnostics remain in
[`results/MIXED_TERM_RESULTS.md`](results/MIXED_TERM_RESULTS.md).

The results characterize this RF-on SciBmad lattice and are not a
machine-validated CESR error budget.

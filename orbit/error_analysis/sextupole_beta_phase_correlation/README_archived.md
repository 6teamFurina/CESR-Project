# Sextupole contribution correlation with beta and phase

This experiment tests whether nominal beta functions and phase advance explain
the relative magnitude of the maintained per-sextupole quadratic detector
contributions. It is deliberately limited to the 76 active normal sextupoles.

The experiment reuses:

- the 100 direction pairs and local `x_h`, `y_v` orbit responses from
  `quadratic_x_attribution/element_results/element_exposure_directions.csv`;
- the saved horizontal and vertical thick-element contribution norms from
  `thick_element_sextupole_sourcing/{horizontal,vertical}_results/`; and
- the saved total detector-vector norm for each direction.

RF-on Twiss functions and phase advances are newly computed for both the
nominal lattice and each of the same 100 simultaneous `h+v` direction states.
The thick-element source is inserted at the complete-element exit, so the
sextupole optics point is evaluated at that same boundary. SciBmad phases are
stored in turns.

Run from the `CESR Project` root:

```powershell
julia --project=. orbit/error_analysis/sextupole_beta_phase_correlation/export_nominal_optics.jl
julia --project=. orbit/error_analysis/sextupole_beta_phase_correlation/export_direction_optics.jl
python orbit/error_analysis/sextupole_beta_phase_correlation/analyze_beta_phase_correlation.py
```

The outputs in `results/` are:

- `nominal_optics_points.csv`: nominal Twiss data at all sextupole source
  boundaries and detectors;
- `direction_optics_points.csv`: the same 175 points for all 100 direction
  states, with `direction_optics_tunes.csv` retaining the matched tunes;
- `sextupole_detector_optics.csv`: all 76 x 99 beta products, phase advances,
  closed-orbit response envelopes, and phase-aware response factors;
- `direction_element_correlation_data.csv`: the joined 100 x 76 x 2-plane
  analysis table;
- `direction_sextupole_transport_factors.csv`: the direction-matched
  detector-vector beta-envelope and beta/phase transport factors;
- `element_correlation_data.csv`: direction-aggregated element values;
- `correlation_summary.csv`: pooled, direction-level, and element-level
  Spearman and log-Pearson correlations; and
- `direction_optics_variation_summary.csv`: beta, tune, and aggregated
  transport-factor changes from the nominal optics; and
- `RESULTS.md` plus two SVG scatter plots.

The envelope predictors omit the cosine phase factor; the phase-aware
predictors include it before taking the detector-vector L2 norm. Nominal and
direction-matched variants are reported separately. All remain approximate
uncoupled-style predictors; the actual contribution is the full coupled
six-dimensional thick-element result.

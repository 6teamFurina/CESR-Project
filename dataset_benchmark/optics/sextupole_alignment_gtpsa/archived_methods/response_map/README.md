# CESR sextupole-alignment GTPSA study

> Historical archive: this study uses the older `cesr_model.jl` lattice. The
> commands below are retained for reproduction only; new work must use the
> latest repaired SciBmad lattice.

This study parameterizes each of the 76 active normal CESR sextupoles
independently with

```text
delta Kn2, delta x_offset, delta y_offset
```

and extracts first- and second-order responses of RF-on periodic detector
orbit and optics, plus the three ring eigentunes.

The required descriptor is `Descriptor(6, 3, 3, 2)`: parameter order two
retains the `Kn2-offset` mixed derivatives, while total/phase-space order three
is necessary to retain those two parameters multiplied by a phase-space
variable in the linear-optics feed-down terms. `Descriptor(6, 2, 3, 2)` is
therefore insufficient for phase, coupling, beta, and tune mixed responses.

The saved `d2_*` values are true derivatives returned by `GTPSA.hessian`.
For a diagonal term, the corresponding Taylor polynomial coefficient is
`d2/2`.

The accumulated `phi_1` and `phi_2` series are gauge-normalized by subtracting
their value at the first detector, `DET_00W`, at every TPS order. This matches
the fixed-reference difference-phase observable and removes the
parameter-dependent additive phase origin used internally by parameterized
Twiss.

Run one contiguous inventory block with:

```powershell
julia --project=. dataset_benchmark/optics/sextupole_alignment_gtpsa/archived_methods/response_map/run_sextupole_alignment_gtpsa.jl `
  --start=1 --stop=19 --part-label=1 `
  --output-dir=dataset_benchmark/optics/sextupole_alignment_gtpsa/archived_methods/response_map/results/full
```

Validate one magnet against exact four-corner numerical differences with:

```powershell
julia --project=. dataset_benchmark/optics/sextupole_alignment_gtpsa/archived_methods/response_map/validate_sextupole_alignment_gtpsa.jl
```

The numerical differences are validation only; the production coefficient
calculation does not use finite differences.

After all four parts finish, merge and validate them with:

```powershell
python dataset_benchmark/optics/sextupole_alignment_gtpsa/archived_methods/response_map/aggregate_results.py
```

Compute the thin SVD of each sextupole's local `2 x 1191` mixed-response
matrix with:

```powershell
python dataset_benchmark/optics/sextupole_alignment_gtpsa/archived_methods/response_map/analyze_local_response_svd.py
```

The script retains an unscaled algebraic reference and an
`observable_rms`-standardized structural analysis. The latter equalizes named
observable types using scales calculated across the full 76-magnet response
dictionary. It is not an experimental uncertainty estimate: that requires
whitening with measured orbit, phase, coupling, and tune covariance.

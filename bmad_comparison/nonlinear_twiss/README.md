# SciBmad nonlinear Twiss outputs

This directory computes first momentum derivatives of the CESR normal-mode
Twiss parameters with a second-order GTPSA descriptor and the patched
four-dimensional RF-off coasting closed orbit.

Run the calculation from the CESR project root:

```console
julia --project=. bmad_comparison/nonlinear_twiss/compute_nonlinear_twiss.jl
```

The `RF_off/` directory contains:

- `nonlinear_twiss.csv`: closed orbit, tunes, global chromaticities,
  accumulated chromatic phase derivatives, beta/alpha functions, and their
  first derivatives with respect to relative momentum deviation `delta` at
  every saved longitudinal position;
- `chromaticity_along_ring.svg`: accumulated `d(phi_i)/d(delta)`, whose ring-end
  value is the corresponding global chromaticity;
- `beta_derivative_along_ring.svg`: `d(beta_i)/d(delta)` for modes 1 and 2;
- `alpha_derivative_along_ring.svg`: `d(alpha_i)/d(delta)` for modes 1 and 2.

The phase `phi_i` is expressed in turns, so its momentum derivative is directly
comparable with `d(Q_i)/d(delta)` at the end of the ring.

## Why the RF cavities are off

With the RF cavities off, `pz` (equivalently the relative momentum deviation
`delta`) is a coasting-beam parameter. It can therefore be used as a controllable
independent variable when evaluating `dQ/d(delta)`, `d(beta)/d(delta)`, and
`d(alpha)/d(delta)`.

With the RF cavities on, `pz` (`delta`) becomes a dynamical variable of the
six-dimensional synchrotron motion rather than an external parameter. It cannot
be used as a controllable independent variable for these derivatives. Therefore,
this nonlinear Twiss study uses the RF-cavity-off lattice, consistent with the
conventional chromatic-Twiss definition.

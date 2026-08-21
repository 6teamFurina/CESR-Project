# Essential SciBmad lattice support

This directory contains files loaded at runtime by
`../latest_cesr_scibmad_repaired.jl`. They are part of the executable SciBmad
lattice, not part of the Bmad-to-SciBmad conversion workspace.

- `latest_lattice_support.jl` defines the phase-continuous latest-CESR wiggler
  segments and loads the girder interface.
- `latest_girder_support.jl` defines `set_latest_girder!`.
- `alignment_coefficients.csv` contains the Bmad-linearized geometry
  coefficients read by the girder interface when the lattice is loaded.

The wiggler support also loads the shared implementation at
`../../wigglers/wiggler.jl`. Therefore a distributable CESR project bundle
must retain that file and the relative directory layout. Conversion,
diagnostic, and validation programs are kept separately in
`../bmad_to_scibmad_tools/` and are not needed for an ordinary lattice load.

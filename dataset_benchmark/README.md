# CESR dataset benchmarks

The benchmark suite is split by generated dataset:

- [`orbit/`](orbit/README.md): RF-on closed-orbit generation, shared corrector
  samples, response-matrix cache, Bmad-compatible lattice, and orbit results.
- [`optics/`](optics/README.md): RF-off coasting Twiss/chromatic-optics
  generation and Bmad/SciBmad method comparisons.

`optics` intentionally reads the common corrector samples and reference files
from `orbit/inputs` and `orbit/reference`; these inputs are not duplicated.
Generated results remain inside the corresponding dataset directory.

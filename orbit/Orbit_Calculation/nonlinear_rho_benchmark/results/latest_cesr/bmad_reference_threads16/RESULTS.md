# Explicit 16-Thread Bmad/Tao Benchmark

The latest-CESR nonlinear-rho reference was rerun in one persistent Tao/PyTao
process with `OMP_NUM_THREADS=16` and Tao `global n_threads=16`. The 9,001 input
states were submitted sequentially; Tao's native OpenMP regions were allowed to
use up to 16 threads.

All `9,001/9,001` states returned usable finite outputs. The complete timed
region was `137.075 s`; the 9,000 nonzero states used in the matched paper table
took `136.797 s`, or `65.791 states/s`. The timed region includes corrector
updates, one Tao model recalculation per state, and observable reads. Tao
initialization (`0.638 s`), warmup (`0.048 s`), and file output are excluded.

The explicit-thread output file is text-identical to the retained first
explicit-thread run over all 9,002 CSV lines; its paper interval is 1.436%
longer. Relative to the reusable-model 16-thread SciBmad calculation, the
application-level wall-clock ratios are `7.605x` for physics and `5.989x`
including reusable-model SciBmad setup.

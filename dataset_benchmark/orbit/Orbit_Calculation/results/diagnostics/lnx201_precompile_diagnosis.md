# lnx201 SciBmad startup diagnosis

Exported: 2026-07-30 23:49 EDT

## Finding

The benchmark was not blocked by the response-matrix cache, memory, or disk
space. The attempted fallback run used `--compiled-modules=no` and eventually
failed inside BeamTracking:

```text
UndefVarError: reference_momentum_shift! not defined in BeamTracking
Suggestion: define the const at top-level before running function that uses it
(stricter Julia v1.12+ rule).
```

The error originated in BeamTracking's generated kernel while SciBmad was
forming the nominal closed-orbit Jacobian. This mode must not be used for the
benchmark on Julia 1.12.

## Exported state

- active benchmark/precompile Julia process: none
- stale lock:
  `/home/jn577/.julia/compiled/v1.12/SciBmad/3SAY4_WKTnQ.ji.pidfile`
- memory: 125 GiB total, 107 GiB available
- home storage: 10 GiB total, 8.3 GiB available
- failed command wall time: 10:05.91
- CPU utilization: 49%
- maximum resident memory: 1,215,108 KiB
- exit status: 1

## Recovery

With no Julia process active, remove only the stale SciBmad pidfile and load
the environment normally with one precompile task:

```bash
rm -f ~/.julia/compiled/v1.12/SciBmad/3SAY4_WKTnQ.ji.pidfile

cd ~/cesr_scibmad
export JULIA_PKG_PRECOMPILE_AUTO=0
export JULIA_NUM_PRECOMPILE_TASKS=1
export JULIA_NUM_THREADS=1

/usr/bin/time -v ~/.juliaup/bin/julia --project=. \
  -e 'using Beamlines, GTPSA, SciBmad; println("environment-load-ok")' \
  2>&1 | tee scibmad_environment_load.log
```

After `environment-load-ok`, run the packaged script without
`--compiled-modules=no`:

```bash
bash run_lnx201_latest.sh
```

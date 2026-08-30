#!/usr/bin/env julia

"""Generate the unknown-error correction/scanning case from a nominal GTPSA ORM.

The orbit-response matrix is calculated once from the nominal latest SciBmad
lattice.  No central-difference ORM is applied and no latent BPM, corrector,
sextupole, quadrupole, roll, or alignment realization is supplied to that
response model.  Observable BPM readbacks still contain their fixed gains and
measurement noise, and physical corrector commands still contain their hidden
gains.  The accepted baseline command is held fixed during every target scan.
"""

include(joinpath(@__DIR__, "generate_corrected_joint_machine_scans.jl"))

function main_gtpsa_nominal_corrected_scans(args=ARGS)
    protocol = [
        "--baseline-response-method=gtpsa",
        "--gtpsa-response-model=nominal",
        "--validate-gtpsa-with-finite-difference=false",
        "--correction-bpm-noise-rms-m=5.0e-6",
        "--correction-measurement-repeats=3072",
        "--correction-noise-seed=20261123",
        "--corrected-case-name=$GTPSA_NOMINAL_CORRECTED_CASE",
        "--target-parallelism=threads",
        "--scan-thread-count=0",
    ]
    return main_corrected_scans(vcat(protocol, args))
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main_gtpsa_nominal_corrected_scans())
end

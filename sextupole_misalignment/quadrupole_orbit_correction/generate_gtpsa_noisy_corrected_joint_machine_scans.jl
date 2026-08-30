#!/usr/bin/env julia

"""Generate fixed-baseline scans using a GTPSA ORM and noisy BPM correction.

The static latent machine is identical to the maintained paired ensemble.  A
first-order SciBmad/GTPSA ORM is calculated in the zero-quadrupole-offset
reference state, including the machine's fixed BPM and corrector gains.  The
stored reference orbit and every current correction measurement use distinct
Gaussian BPM-mean noise draws from the same 5-micrometer-per-read model.  The
default 3,072 repeats match the maintained sextupole inverse acquisition model.
The resulting baseline command is held fixed during every target scan.
"""

include(joinpath(@__DIR__, "generate_corrected_joint_machine_scans.jl"))

function main_gtpsa_noisy_corrected_scans(args=ARGS)
    protocol = [
        "--baseline-response-method=gtpsa",
        "--validate-gtpsa-with-finite-difference=true",
        "--correction-bpm-noise-rms-m=5.0e-6",
        "--correction-measurement-repeats=3072",
        "--correction-noise-seed=20261123",
        "--corrected-case-name=$GTPSA_NOISY_CORRECTED_CASE",
    ]
    return main_corrected_scans(vcat(protocol, args))
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main_gtpsa_noisy_corrected_scans())
end

#!/usr/bin/env julia

"""Generate the paired exact SciBmad states for the enhanced GTPSA protocol.

Only the clean and time-drift secant tensors are needed.  BPM white noise and
random acquisition histories are applied analytically in the fast inverse, so
they never require additional closed-orbit solves.
"""

include(joinpath(
    @__DIR__, "..", "real_machine_nuisance_ablation",
    "generate_physical_nuisance_scans.jl",
))

function protocol_main(args=ARGS)
    defaults = [
        "--cases=baseline,time_drift",
        "--realizations=4",
        "--target-limit=0",
        "--bump-amplitude-m=0.0015",
        "--k2-step-m3=0.05",
        "--output-root=$(joinpath(@__DIR__, "results", "exact_k5_b3"))",
    ]
    return main(vcat(defaults, args))
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(protocol_main())
end

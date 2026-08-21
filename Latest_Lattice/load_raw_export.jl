"""
Smoke-test the unmodified SciBmad file emitted by Bmad/Tao 20260814-0.

This script is expected to fail until the unsupported/forward-referenced
photon-branch elements in the raw export are repaired.
"""

include(joinpath(@__DIR__, "latest_cesr_scibmad_bmad_20260814.jl"))
println("raw_export_loaded=true")
println("cesr_type=", typeof(cesr))
println("cesr_elements=", length(cesr))

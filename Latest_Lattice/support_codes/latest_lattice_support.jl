using SciBmad

include(joinpath(@__DIR__, "..", "..", "wigglers", "wiggler.jl"))
using .WigglerModels: planar_wiggler_four_potential

include(joinpath(@__DIR__, "latest_girder_support.jl"))

const LATEST_WIGGLER_B_MAX = 1.17
const LATEST_WIGGLER_PERIOD = 0.19625
const LATEST_WIGGLER_TOTAL_LENGTH = 2.355
const LATEST_WIGGLER_KW = 2pi / LATEST_WIGGLER_PERIOD
const LATEST_WIGGLER_PHASE0 = -LATEST_WIGGLER_KW * LATEST_WIGGLER_TOTAL_LENGTH / 2
# Bmad advances the reference time through its planar wiggler so the design
# particle exits at z=0. The continuous SciBmad field otherwise reports the
# 1.962917 um quiver-path delay as an affine z offset.
const LATEST_WIGGLER_REFERENCE_DT = 6.54758669568665e-15

"""
    LatestWigglerSegment(; L, s_offset, n_steps, kwargs...)

Construct one phase-continuous piece of the Bmad `ID_S1A` planar wiggler.
`s_offset` is the piece entrance measured from the upstream end of the full
2.355 m element. The explicit step counts reproduce Bmad's 90/30/120 split
after the superimposed sector marker and photon fork divide the 240-step
wiggler into three tracking slaves.
"""
function LatestWigglerSegment(;
    L,
    s_offset,
    n_steps::Integer,
    order::Integer=6,
    radiation_damping_on::Bool=false,
    radiation_fluctuations_on::Bool=false,
    kwargs...,
)
    L > 0 || throw(ArgumentError("L must be positive"))
    0 <= s_offset < LATEST_WIGGLER_TOTAL_LENGTH ||
        throw(ArgumentError("s_offset is outside the full ID_S1A wiggler"))
    n_steps > 0 || throw(ArgumentError("n_steps must be positive"))
    order in (2, 4, 6, 8) ||
        throw(ArgumentError("Yoshida order must be one of 2, 4, 6, or 8"))

    phase = LATEST_WIGGLER_PHASE0 + LATEST_WIGGLER_KW * s_offset
    return LineElement(;
        kind="Wiggler",
        L,
        four_potential=planar_wiggler_four_potential,
        four_potential_params=(LATEST_WIGGLER_B_MAX, LATEST_WIGGLER_KW, phase),
        four_potential_normalized=false,
        tracking_method=Yoshida(;
            order,
            n_steps,
            radiation_damping_on,
            radiation_fluctuations_on,
        ),
        kwargs...,
    )
end

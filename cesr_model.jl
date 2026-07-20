using Beamlines

if !isdefined(@__MODULE__, :load_cesr)
    include(joinpath(@__DIR__, "cesr.jl"))
end

if !isdefined(@__MODULE__, :attach_cesr_controls!)
    include(joinpath(@__DIR__, "cesr_controls.jl"))
end

"""
    load_cesr_model(; zero_value=0.0, rf_on=nothing, rf_voltage=CESR_RF_VOLTAGE)

Load an independent copy of the static CESR ring and attach the complete
DefExpr-based control layer.

`rf_on=nothing` preserves the static RF state. Set `rf_on=true` to enable all
four cavities, or `rf_on=false` to explicitly set their voltage to zero.

Returns a named tuple with fields `ring` and `controls`.
"""
function load_cesr_model(;
    zero_value=0.0,
    rf_on::Union{Nothing,Bool}=nothing,
    rf_voltage::Real=CESR_RF_VOLTAGE,
)
    ring = load_cesr()
    if !isnothing(rf_on)
        set_cesr_rf!(ring; on=rf_on, voltage=rf_voltage)
    end
    controls = attach_cesr_controls!(ring; zero_value)
    return (; ring, controls)
end

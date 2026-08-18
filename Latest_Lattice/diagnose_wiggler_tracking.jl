#!/usr/bin/env julia

"""Check convergence of the continuous-field ID_S1A map against Bmad."""

using LinearAlgebra
using Printf
using GTPSA
import Beamlines

include(joinpath(@__DIR__, "compare_local_maps.jl"))

function bmad_block_map(reference, first_index, last_index)
    R = Matrix{Float64}(I, 6, 6)
    vec0 = zeros(6)
    for index in first_index:last_index
        element = reference[index + 1]
        vec0 = element.R * vec0 + element.vec0
        R = element.R * R
    end
    return R, vec0
end

function scibmad_block_map(first_index, last_index, entrance, order, step_scale)
    elements = [
        Beamlines.deepcopy_no_beamline(element)
        for element in cesr.line[first_index:last_index]
    ]
    for element in elements
        element_kind(element) == "Wiggler" || continue
        element.tracking_method = Yoshida(
            order=order,
            n_steps=step_scale * element.tracking_method.n_steps,
        )
    end
    line = Beamlines.Beamline(
        elements;
        p_over_q_ref=cesr.p_over_q_ref,
        species_ref=cesr.species_ref,
    )
    map = WigglerModels.gtpsa_transport_map(line; v0=entrance, order=1)
    R = Matrix(GTPSA.jacobian(map.v))
    output = GTPSA.scalar.(map.v)
    return R, output - R * entrance
end

function main()
    reference = read_reference(REFERENCE)
    first_index, last_index = 997, 1001
    entrance = reference[first_index].orbit_out
    R_bmad, vec0_bmad = bmad_block_map(reference, first_index, last_index)
    println("order,step_scale,total_steps,max_abs_dR,relative_fro_dR,max_abs_dvec0")
    for order in (4, 6, 8), step_scale in (1, 2, 4)
        R, vec0 = scibmad_block_map(
            first_index, last_index, entrance, order, step_scale,
        )
        @printf(
            "%d,%d,%d,%.12e,%.12e,%.12e\n",
            order, step_scale, 240 * step_scale,
            maximum(abs, R - R_bmad), norm(R - R_bmad) / norm(R_bmad),
            maximum(abs, vec0 - vec0_bmad),
        )
        if order == 8 && step_scale == 4
            println("delta_R=")
            display(R - R_bmad)
            println("delta_vec0=$(vec0 - vec0_bmad)")
        end
    end
end

main()

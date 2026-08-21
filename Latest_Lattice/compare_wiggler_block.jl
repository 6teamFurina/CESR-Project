#!/usr/bin/env julia

"""Compare the complete split ID_S1A block rather than gauge-dependent slices."""

using LinearAlgebra
using Printf
using GTPSA
import Beamlines

include(joinpath(@__DIR__, "compare_local_maps.jl"))

function main()
    reference = read_reference(REFERENCE)
    first_index = 997
    last_index = 1001
    elements = [
        Beamlines.deepcopy_no_beamline(element)
        for element in cesr.line[first_index:last_index]
    ]
    reference_time_shift = isempty(ARGS) ? Float64(elements[4].dt) : parse(Float64, only(ARGS))
    # PatchParams cannot coexist with the four-potential on a Wiggler, so use
    # the existing zero-length photon-fork marker inside the split block.
    elements[4].dt = reference_time_shift
    line = Beamlines.Beamline(
        elements;
        p_over_q_ref=cesr.p_over_q_ref,
        species_ref=cesr.species_ref,
    )
    entrance = reference[first_index].orbit_out
    map = WigglerModels.gtpsa_transport_map(line; v0=entrance, order=1)
    R_scibmad = Matrix(GTPSA.jacobian(map.v))
    output_scibmad = GTPSA.scalar.(map.v)
    vec0_scibmad = output_scibmad - R_scibmad * entrance

    R_bmad = Matrix{Float64}(I, 6, 6)
    vec0_bmad = zeros(6)
    for index in first_index:last_index
        element = reference[index + 1]
        vec0_bmad = element.R * vec0_bmad + element.vec0
        R_bmad = element.R * R_bmad
    end
    output_bmad = R_bmad * entrance + vec0_bmad

    names = join((reference[index + 1].name for index in first_index:last_index), ", ")
    println("indices=$first_index:$last_index")
    println("elements=$names")
    @printf("downstream_reference_time_shift=%.16e s\n", reference_time_shift)
    @printf("length=%.16g m\n", sum(element.L for element in elements))
    @printf("max_abs_dR=%.12e\n", maximum(abs, R_scibmad - R_bmad))
    @printf("relative_fro_dR=%.12e\n", norm(R_scibmad - R_bmad) / norm(R_bmad))
    @printf("max_abs_dvec0=%.12e\n", maximum(abs, vec0_scibmad - vec0_bmad))
    @printf("max_abs_exit_orbit=%.12e\n", maximum(abs, output_scibmad - output_bmad))
    println("SciBmad vec0=$(vec0_scibmad)")
    println("Bmad vec0=$(vec0_bmad)")
end

main()

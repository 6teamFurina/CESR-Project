#!/usr/bin/env julia

"""Sweep SciBmad integration settings for the first DQX combined-function bend."""

using LinearAlgebra
using Printf

include(joinpath(@__DIR__, "compare_local_maps.jl"))

function main()
    reference = read_reference(REFERENCE)
    index = 22
    element = cesr.line[index]
    bmad = reference[index + 1]
    entrance = reference[index]
    original_method = element.tracking_method

    println("element=$(element.name) Bmad=$(bmad.name) L=$(element.L)")
    println("order,n_steps,max_abs_dR,relative_fro,R11,R12,R21,R22")
    try
        for order in (2, 4, 6, 8), n_steps in (1, 2, 5, 10, 25, 50, 100, 200)
            element.tracking_method = Yoshida(; order, n_steps)
            R, _, _ = local_map(element, entrance.orbit_out)
            delta = R - bmad.R
            @printf(
                "%d,%d,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e\n",
                order, n_steps, maximum(abs, delta), norm(delta) / norm(bmad.R),
                R[1, 1], R[1, 2], R[2, 1], R[2, 2],
            )
        end
    finally
        element.tracking_method = original_method
    end
end

main()

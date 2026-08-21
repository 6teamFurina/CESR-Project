#!/usr/bin/env julia

"""Compare nominal RF-on closed orbit and tunes for the latest CESR ring."""

using LinearAlgebra
using Printf
using GTPSA
using SciBmad

include(joinpath(@__DIR__, "compare_local_maps.jl"))

const OPTICS_SUMMARY = joinpath(@__DIR__, "RING_OPTICS_COMPARISON.md")
const ELEMENT_INDEX = joinpath(
    @__DIR__, "bmad_reference_branch0", "element_index.csv",
)

function final_bmad_twiss()
    lines = readlines(ELEMENT_INDEX)
    header = split(first(lines), ',')
    values = split(last(lines), ',')
    columns = Dict(name => values[index] for (index, name) in enumerate(header))
    return (
        beta_a=parse(Float64, columns["beta_a"]),
        alpha_a=parse(Float64, columns["alpha_a"]),
        phi_a=parse(Float64, columns["phi_a"]),
        beta_b=parse(Float64, columns["beta_b"]),
        alpha_b=parse(Float64, columns["alpha_b"]),
        phi_b=parse(Float64, columns["phi_b"]),
    )
end

function main()
    reference = read_reference(REFERENCE)
    R_bmad = Matrix{Float64}(I, 6, 6)
    for element in @view reference[2:end]
        R_bmad = element.R * R_bmad
    end
    bmad_mode_tunes = sort(
        abs.(angle.(filter(value -> imag(value) > 0, eigvals(R_bmad)))) ./ (2pi),
    )
    final_twiss = final_bmad_twiss()
    bmad_transverse_tunes = mod.([final_twiss.phi_a, final_twiss.phi_b] ./ (2pi), 1)

    closed = find_closed_orbit(cesr)
    optics = twiss(
        cesr;
        GTPSA_descriptor=GTPSA.Descriptor(6, 1),
        v0=closed.v0,
        v0_and_coast=(closed.v0, closed.coasting_beam),
    )
    scibmad_tunes = GTPSA.scalar.(optics.tunes)
    scibmad_mode_tunes = sort(abs.(scibmad_tunes))
    bmad_orbit = reference[1].orbit_out
    scibmad_orbit = vec(closed.v0)
    orbit_delta = scibmad_orbit - bmad_orbit

    open(OPTICS_SUMMARY, "w") do io
        println(io, "# Latest CESR nominal RF-on ring comparison")
        println(io)
        println(io, "- Bmad/Tao reference: `20260814-0`, branch 0.")
        println(io, "- SciBmad lattice: `latest_cesr_scibmad_repaired.jl`.")
        println(io, "- Bmad starting closed orbit: `$(bmad_orbit)`." )
        println(io, "- SciBmad starting closed orbit: `$(scibmad_orbit)`." )
        @printf(io, "- Maximum starting closed-orbit difference: `%.12e`.\n", maximum(abs, orbit_delta))
        println(io, "- Bmad transverse tunes from Twiss phase: `$(bmad_transverse_tunes)`." )
        println(io, "- Bmad one-turn eigenphase magnitudes: `$(bmad_mode_tunes)`." )
        println(io, "- SciBmad signed tunes: `$(scibmad_tunes)`." )
        println(io, "- SciBmad eigenphase magnitudes: `$(scibmad_mode_tunes)`." )
        println(io)
        println(io, "The remaining nominal discrepancy is dominated by the continuous-field " *
                    "wiggler map versus Bmad's faster standard-matrix approximation; " *
                    "see `compare_wiggler_block.jl`.")
    end

    println("Bmad starting orbit:    $bmad_orbit")
    println("SciBmad starting orbit: $scibmad_orbit")
    @printf("Maximum orbit difference: %.12e\n", maximum(abs, orbit_delta))
    println("Bmad transverse tunes:  $bmad_transverse_tunes")
    println("Bmad eigenphase tunes:  $bmad_mode_tunes")
    println("SciBmad signed tunes:   $scibmad_tunes")
    println("SciBmad magnitude tunes: $scibmad_mode_tunes")
    println("Wrote $OPTICS_SUMMARY")
end

main()

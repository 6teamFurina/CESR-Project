#!/usr/bin/env julia

"""Match the Bmad DQX/girder isolation test in SciBmad."""

using LinearAlgebra
using Printf
using SciBmad

const LATTICE_DIR = normpath(joinpath(@__DIR__, ".."))

include(joinpath(LATTICE_DIR, "latest_cesr_scibmad_repaired.jl"))

const REFERENCE = joinpath(LATTICE_DIR, "bmad_reference", "girder", "dqx4b_pitch_isolation.json")
const STEP = 1.0e-6

function json_vector(text, key)
    match_result = match(Regex("\\\"$key\\\"\\s*:\\s*\\[([^]]+)\\]"), text)
    isnothing(match_result) && error("Missing $key in $REFERENCE")
    return parse.(Float64, strip.(split(match_result.captures[1], ',')))
end

function response()
    function track_to_dqx()
        bunch = Bunch(v=zeros(1, 6))
        SciBmad.BTBL.check_bl_bunch!(bunch, cesr, false)
        for index in 1:22
            track!(bunch, cesr.line[index])
        end
        return vec(Float64.(bunch.coords.v))
    end
    set_latest_girder!(cesr, "GIRDER_4AB"; y_pitch=STEP)
    plus = track_to_dqx()
    set_latest_girder!(cesr, "GIRDER_4AB"; y_pitch=-STEP)
    minus = track_to_dqx()
    set_latest_girder!(cesr, "GIRDER_4AB")
    return (plus - minus) / (2STEP)
end

function report(label, bmad, scibmad)
    difference = scibmad - bmad
    @printf("%s\n", label)
    @printf("  Bmad:    %s\n", bmad)
    @printf("  SciBmad: %s\n", scibmad)
    @printf("  relative L2: %.12e\n", norm(difference) / norm(bmad))
    @printf("  max abs:     %.12e\n", maximum(abs, difference))
end

function main()
    text = read(REFERENCE, String)
    bmad_nominal = json_vector(text, "nominal")
    bmad_runge_kutta = json_vector(text, "runge_kutta")
    bmad_no_k1 = json_vector(text, "no_k1")
    original_kn1 = dqx4b.Kn1
    original_method = dqx4b.tracking_method
    scibmad_nominal = response()
    dqx4b.Kn1 = 0.0
    scibmad_no_k1 = response()
    dqx4b.Kn1 = original_kn1
    report("nominal DQX4B", bmad_nominal, scibmad_nominal)
    report("Bmad DQX4B Runge-Kutta versus SciBmad Yoshida", bmad_runge_kutta, scibmad_nominal)
    report("DQX4B with K1=0", bmad_no_k1, scibmad_no_k1)
    println("step convergence with nominal K1")
    for n_steps in (100, 200, 400, 800)
        dqx4b.tracking_method = Yoshida(order=6, n_steps=n_steps)
        candidate = response()
        relative = norm(candidate - bmad_nominal) / norm(bmad_nominal)
        @printf("  %4d steps: dy=% .12e dpy=% .12e relative_L2=%.12e\n",
                n_steps, candidate[3], candidate[4], relative)
    end
    dqx4b.tracking_method = original_method
end

main()

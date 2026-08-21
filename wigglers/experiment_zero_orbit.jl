#!/usr/bin/env julia

"""
Numerical checks for the CESR planar-wiggler zero-input orbit.

Run from the CESR project environment with:

    julia --project=. wigglers/experiment_zero_orbit.jl

The experiment distinguishes transverse closure from the longitudinal path-
length term and inspects whether the two wigglers in `cesr.jl` are active in
the default lattice.
"""

using Beamlines
using Printf
using SciBmad

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(PROJECT_ROOT, "cesr.jl"))

const B_MAX = 1.17
const L_PERIOD = 0.19625
const N_PERIOD = 12
const WIGGLER_LENGTH = N_PERIOD * L_PERIOD
const P0C = 5.2889999753148e9
const ZERO6 = zeros(6)

function one_element_line(element)
    return Beamlines.Beamline(
        [element];
        species_ref=Species("electron"),
        E_ref=P0C,
    )
end

function track_exit(element; v0=ZERO6)
    tracked = track(
        one_element_line(element);
        v0=reshape(collect(v0), 1, 6),
        use_KA=false,
        use_explicit_SIMD=false,
    )
    return vec(tracked.v[1, :, end])
end

function print_vector(label, vector)
    @printf("%-30s", label)
    for value in vector
        @printf(" % .12e", value)
    end
    println()
end

function nominal_wiggler(; slices_per_period=16, phase=nothing,
                         radiation_damping_on=false)
    return PlanarWiggler(
        alias="WIG_EXPERIMENT",
        B_max=B_MAX,
        L_period=L_PERIOD,
        N_period=N_PERIOD,
        L=WIGGLER_LENGTH,
        phase=phase,
        slices_per_period=slices_per_period,
        order=6,
        radiation_damping_on=radiation_damping_on,
        radiation_fluctuations_on=false,
    )
end

"""Construct a deliberately incomplete periodic field for a closure control."""
function truncated_wiggler(; periods=11.5, phase=0.0, slices_per_period=16)
    length = periods * L_PERIOD
    k_w = 2pi / L_PERIOD
    return LineElement(
        kind="TruncatedWiggler",
        L=length,
        four_potential=WigglerModels.planar_wiggler_four_potential,
        four_potential_params=(B_MAX, k_w, phase),
        four_potential_normalized=false,
        tracking_method=Yoshida(
            order=6,
            n_steps=round(Int, periods * slices_per_period),
            radiation_damping_on=false,
            radiation_fluctuations_on=false,
        ),
    )
end

"""Split the same continuous field into short elements so internal positions are saved."""
function sampled_wiggler_line(; samples_per_period=4, steps_per_sample=4)
    n_samples = N_PERIOD * samples_per_period
    sample_length = WIGGLER_LENGTH / n_samples
    k_w = 2pi / L_PERIOD
    phase0 = -k_w * WIGGLER_LENGTH / 2
    elements = map(1:n_samples) do index
        s0 = (index - 1) * sample_length
        LineElement(
            kind="WigglerSample",
            L=sample_length,
            four_potential=WigglerModels.planar_wiggler_four_potential,
            four_potential_params=(B_MAX, k_w, phase0 + k_w * s0),
            four_potential_normalized=false,
            tracking_method=Yoshida(
                order=6,
                n_steps=steps_per_sample,
                radiation_damping_on=false,
                radiation_fluctuations_on=false,
            ),
        )
    end
    return Beamlines.Beamline(
        elements;
        species_ref=Species("electron"),
        E_ref=P0C,
    )
end

function run_internal_orbit_check()
    line = sampled_wiggler_line()
    states = Vector{Vector{Float64}}()
    state = zeros(6)
    push!(states, copy(state))
    for element in line.line
        isolated_element = Beamlines.deepcopy_no_beamline(element)
        state = track_exit(isolated_element; v0=state)
        push!(states, copy(state))
    end
    trajectory = reduce(hcat, states)
    x_values = vec(trajectory[1, :])
    px_values = vec(trajectory[2, :])
    @printf("saved longitudinal stations   %d\n", length(x_values))
    @printf("max internal |x| [m]          %.12e\n", maximum(abs, x_values))
    @printf("max internal canonical |px|   %.12e\n", maximum(abs, px_values))
    print_vector("sampled-lattice exit", vec(trajectory[:, end]))
end

function run_convergence_check()
    println("slices/period     max transverse residual          longitudinal coordinate")
    for slices_per_period in (2, 4, 8, 16, 32)
        output = track_exit(nominal_wiggler(; slices_per_period))
        transverse_residual = maximum(abs, output[1:4])
        @printf("%13d     %.12e             % .12e\n",
                slices_per_period, transverse_residual, output[5])
    end
end

function run_longitudinal_input_check()
    element = nominal_wiggler()
    cases = (
        ("z = +1 mm", [0.0, 0.0, 0.0, 0.0, 1e-3, 0.0]),
        ("delta = +1%", [0.0, 0.0, 0.0, 0.0, 0.0, 1e-2]),
        ("delta = -1%", [0.0, 0.0, 0.0, 0.0, 0.0, -1e-2]),
        ("z=1 mm, delta=+1%", [0.0, 0.0, 0.0, 0.0, 1e-3, 1e-2]),
    )
    for (label, input) in cases
        print_vector(label, track_exit(element; v0=input))
    end
end

function run_controls()
    nominal = track_exit(nominal_wiggler())
    field_off = track_exit(PlanarWiggler(
        B_max=0.0,
        L_period=L_PERIOD,
        N_period=N_PERIOD,
        L=WIGGLER_LENGTH,
        slices_per_period=16,
        order=6,
    ))
    phase_shifted = track_exit(nominal_wiggler(; phase=pi / 2))
    truncated = track_exit(truncated_wiggler())
    radiation_on = track_exit(nominal_wiggler(; radiation_damping_on=true))

    print_vector("nominal symmetric", nominal)
    print_vector("field off", field_off)
    print_vector("phase = pi/2 (canonical 0)", phase_shifted)
    print_vector("11.5 periods, phase = 0", truncated)
    print_vector("radiation damping on", radiation_on)

    scales = WigglerModels.planar_wiggler_scales(
        B_max=B_MAX,
        L_period=L_PERIOD,
        p0c=P0C,
    )
    paraxial_path_length = WIGGLER_LENGTH * scales.angle_amplitude^2 / 4
    @printf("paraxial extra path [m]       %.12e\n", paraxial_path_length)
    @printf("tracked |z| / estimate        %.12f\n",
            abs(nominal[5]) / paraxial_path_length)
end

function inspect_cesr_wigglers()
    ring = load_cesr()
    wigglers = filter(element -> uppercase(String(element.kind)) == "WIGGLER", ring.line)
    @printf("CESR wiggler count             %d\n", length(wigglers))
    for element in wigglers
        b_max, k_w, phase = element.four_potential_params
        @printf("%-12s B_max=%g T  period=%.12g m  phase=%.12g rad  active=%s\n",
                String(element.name), b_max, 2pi / k_w, phase, string(!iszero(b_max)))
        method = element.tracking_method
        @printf("             radiation damping=%s  fluctuations=%s\n",
                string(method.radiation_damping_on),
                string(method.radiation_fluctuations_on))
        print_vector("  isolated zero-input exit", track_exit(element))
    end
end

println("coordinate order: x, px, y, py, z, delta")
println("\n--- zero-input controls ---")
run_controls()
println("\n--- internal sampled orbit ---")
run_internal_orbit_check()
println("\n--- integration convergence ---")
run_convergence_check()
println("\n--- zero transverse input with longitudinal offsets ---")
run_longitudinal_input_check()
println("\n--- default CESR lattice inspection ---")
inspect_cesr_wigglers()

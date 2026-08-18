#!/usr/bin/env julia

"""Reconstruct signed normal-sextupole contributions to detector Q_x."""

const SEXT_DETECTOR_HERE = @__DIR__
const QUADRATIC_ATTRIBUTION_HERE = normpath(joinpath(
    SEXT_DETECTOR_HERE, "..", "quadratic_x_attribution",
))
include(joinpath(QUADRATIC_ATTRIBUTION_HERE, "run_internal_exposure_attribution.jl"))

using Dates
using LinearAlgebra
using Printf
using TOML

function parse_sext_detector_args(args)
    options = Dict{String,String}(
        "trials" => "100",
        "seed" => "20260804",
        "base-kick-rad" => "5e-6",
        "output-dir" => joinpath(SEXT_DETECTOR_HERE, "results"),
    )
    for argument in args
        startswith(argument, "--") || error("Arguments must have --name=value form: $argument")
        fields = split(argument[3:end], '='; limit=2)
        length(fields) == 2 || error("Missing value in argument: $argument")
        haskey(options, fields[1]) || error("Unknown option: --$(fields[1])")
        options[fields[1]] = fields[2]
    end
    return options
end

function detector_positions(ring, detectors)
    lookup = Dict(name => index for (index, name) in enumerate(detectors))
    positions = fill(NaN, length(detectors))
    s_m = 0.0
    for element in ring.line
        s_m += Float64(GTPSA.scalar(Beamlines.deval(element.L)))
        name = uppercase(String(element.name))
        haskey(lookup, name) || continue
        positions[lookup[name]] = s_m
    end
    all(isfinite, positions) || error("Not all detector positions were found")
    return positions
end

"""Periodic detector response to split entrance/exit kicks at sextupoles."""
function build_local_kick_detector_response(closed_orbit, inventory, detectors)
    element_indices = sort!(collect(keys(inventory)))
    source_row = Dict(index => row for (row, index) in enumerate(element_indices))
    n_sources = length(element_indices)
    descriptor = Descriptor(6, 1, 2 * n_sources, 1)
    variables = vars(descriptor)
    parameters = params(descriptor)
    model = load_cesr_model(zero_value=zero(parameters[1]), rf_on=true)
    input_map = [closed_orbit[index] + copy(variables[index]) for index in 1:6]
    bunch = Bunch(v=reshape(input_map, 1, 6))
    SciBmad.BTBL.check_bl_bunch!(bunch, model.ring, false)

    detector_index = Dict(name => index for (index, name) in enumerate(detectors))
    detector_maps = Vector{typeof(copy.(input_map))}(undef, length(detectors))
    found = falses(length(detectors))
    for (element_index, element) in enumerate(model.ring.line)
        row = get(source_row, element_index, 0)
        if !iszero(row)
            # A half kick on each side is the second-order symmetric thin-kick
            # approximation to a source sampled at the element midpoint.
            bunch.coords.v[1, 2] += parameters[row] / 2
            bunch.coords.v[1, 4] += parameters[n_sources + row] / 2
        end
        track!(bunch, element)
        if !iszero(row)
            bunch.coords.v[1, 2] += parameters[row] / 2
            bunch.coords.v[1, 4] += parameters[n_sources + row] / 2
        end
        name = uppercase(String(element.name))
        haskey(detector_index, name) || continue
        detector = detector_index[name]
        found[detector] && error("Detector $name occurs more than once")
        detector_maps[detector] = copy.(vec(bunch.coords.v))
        found[detector] = true
    end
    all(found) || error("Missing detector maps: $(detectors[.!found])")
    one_turn_map = copy.(vec(bunch.coords.v))

    jacobian = Matrix(GTPSA.jacobian(one_turn_map; include_params=true))
    expected_columns = 6 + 2 * n_sources
    size(jacobian) == (6, expected_columns) ||
        error("Unexpected local-kick one-turn Jacobian size $(size(jacobian))")
    A = jacobian[:, 1:6]
    B = jacobian[:, 7:end]
    fixed_point_matrix = I - A
    closed_response = fixed_point_matrix \ B
    lifted = vcat(closed_response, Matrix{Float64}(I, 2 * n_sources, 2 * n_sources))

    response_x_px = zeros(length(detectors), n_sources)
    response_x_py = similar(response_x_px)
    for detector in eachindex(detector_maps)
        detector_jacobian = Matrix(GTPSA.jacobian(
            detector_maps[detector]; include_params=true,
        ))
        response = detector_jacobian * lifted
        response_x_px[detector, :] .= view(response, 1, 1:n_sources)
        response_x_py[detector, :] .= view(response, 1, n_sources + 1:2 * n_sources)
    end

    one_turn_response = jacobian * lifted
    fixed_point_closure = maximum(abs, one_turn_response - closed_response)
    solve_closure = maximum(abs, fixed_point_matrix * closed_response - B)
    names = [String(inventory[index].name) for index in element_indices]
    k2l = [Float64(inventory[index].k2l) for index in element_indices]
    s_m = [Float64(inventory[index].s_m) for index in element_indices]
    return (;
        response_x_px, response_x_py, element_indices, names, k2l, s_m,
        solve_closure, fixed_point_closure,
    )
end

function local_block_sources(k2l, x_h, x_v, y_h, y_v)
    hh_px = @. -0.5 * k2l * (x_h^2 - y_h^2)
    hh_py = @. k2l * x_h * y_h
    hv_px = @. -k2l * (x_h * x_v - y_h * y_v)
    hv_py = @. k2l * (x_h * y_v + x_v * y_h)
    vv_px = @. -0.5 * k2l * (x_v^2 - y_v^2)
    vv_py = @. k2l * x_v * y_v
    return Dict(:hh => (hh_px, hh_py), :hv => (hv_px, hv_py), :vv => (vv_px, vv_py))
end

function propagated_contributions(response, px, py)
    return response.response_x_px .* reshape(px, 1, :) .+
           response.response_x_py .* reshape(py, 1, :)
end

function relative_vector_residual(target, reconstructed)
    denominator = norm(target)
    return denominator == 0 ? NaN : norm(target - reconstructed) / denominator
end

function main_sextupole_detector_contributions(args=ARGS)
    options = parse_sext_detector_args(args)
    trials = parse(Int, options["trials"])
    trials >= 1 || error("--trials must be positive")
    seed = parse(Int, options["seed"])
    base_kick = parse(Float64, options["base-kick-rad"])
    base_kick > 0 || error("--base-kick-rad must be positive")
    output_dir = abspath(options["output-dir"])

    input_reference = read_samples(joinpath(
        MIXED_CALCULATION_DIR, "inputs", "cesr_corrector_samples_1000.csv",
    ))
    names = input_reference.names
    float_model = load_cesr_model(zero_value=0.0, rf_on=true)
    detectors = detector_names(float_model.ring)
    detector_s_m = detector_positions(float_model.ring, detectors)
    closed_orbit = nominal_closed_orbit(float_model.ring)
    inventory = normal_sextupole_inventory(float_model.ring)
    samples = generate_mixed_samples(names, [1.0], trials, seed, base_kick)

    @printf("Building 119-corrector internal response at %d sextupoles...\n", length(inventory))
    internal = build_internal_response(names, closed_orbit, base_kick, inventory)
    @printf("Building %d-parameter local-kick detector response...\n", 2 * length(inventory))
    kick_response = build_local_kick_detector_response(closed_orbit, inventory, detectors)
    internal.element_names == kick_response.names || error("Sextupole order mismatch")
    maximum(abs, internal.k2l - kick_response.k2l) < 1e-14 || error("K2L mismatch")

    n_sources = length(kick_response.names)
    blocks = (:hh, :hv, :vv)
    numerator = Dict(block => zeros(n_sources) for block in blocks)
    magnitude_sq = Dict(block => zeros(n_sources) for block in blocks)
    denominator = Dict(block => 0.0 for block in blocks)
    residual_sq = Dict(block => 0.0 for block in blocks)
    reconstruction_dot = Dict(block => 0.0 for block in blocks)
    total_numerator = zeros(n_sources)
    total_magnitude_sq = zeros(n_sources)
    total_denominator = 0.0
    total_residual_sq = 0.0
    total_reconstruction_dot = 0.0
    exposure_h_sum = zeros(n_sources)
    exposure_v_sum = zeros(n_sources)
    direction_rows = NamedTuple[]
    element_direction_rows = NamedTuple[]

    solve_seconds = @elapsed begin
        for trial in 1:trials
            @printf("Direction %d/%d\n", trial, trials)
            h_direction = view(samples.horizontal_directions, trial, :)
            v_direction = view(samples.vertical_directions, trial, :)
            exact = direction_gtpsa_q(names, detectors, closed_orbit,
                h_direction, v_direction, base_kick)
            x_h = internal.response_x * h_direction
            x_v = internal.response_x * v_direction
            y_h = internal.response_y * h_direction
            y_v = internal.response_y * v_direction
            exposure_h = abs.(kick_response.k2l) .* x_h .^ 2
            exposure_v = abs.(kick_response.k2l) .* y_v .^ 2
            exposure_h_sum .+= exposure_h
            exposure_v_sum .+= exposure_v
            sources = local_block_sources(kick_response.k2l, x_h, x_v, y_h, y_v)

            targets = Dict(block => exact.q[Symbol("x_$(block)")] for block in blocks)
            contributions = Dict{Symbol,Matrix{Float64}}()
            reconstructions = Dict{Symbol,Vector{Float64}}()
            for block in blocks
                px, py = sources[block]
                contribution = propagated_contributions(kick_response, px, py)
                reconstruction = vec(sum(contribution; dims=2))
                contributions[block] = contribution
                reconstructions[block] = reconstruction
                target = targets[block]
                denominator[block] += sum(abs2, target)
                residual_sq[block] += sum(abs2, target - reconstruction)
                reconstruction_dot[block] += dot(reconstruction, target)
                for source in 1:n_sources
                    vector = view(contribution, :, source)
                    numerator[block][source] += dot(vector, target)
                    magnitude_sq[block][source] += sum(abs2, vector)
                end
            end

            total_target = targets[:hh] + targets[:hv] + targets[:vv]
            total_contribution = contributions[:hh] + contributions[:hv] + contributions[:vv]
            total_reconstruction = vec(sum(total_contribution; dims=2))
            total_denominator += sum(abs2, total_target)
            total_residual_sq += sum(abs2, total_target - total_reconstruction)
            total_reconstruction_dot += dot(total_reconstruction, total_target)
            for source in 1:n_sources
                vector = view(total_contribution, :, source)
                total_numerator[source] += dot(vector, total_target)
                total_magnitude_sq[source] += sum(abs2, vector)
                push!(element_direction_rows, (; trial, element_order=source,
                    element_name=kick_response.names[source], s_m=kick_response.s_m[source],
                    k2l_m2=kick_response.k2l[source], source_exposure_h=exposure_h[source],
                    source_exposure_v=exposure_v[source],
                    total_projection_numerator=dot(vector, total_target),
                    total_contribution_norm_m=norm(vector)))
            end

            push!(direction_rows, (; trial,
                q_hh_norm_m=norm(targets[:hh]), q_hv_norm_m=norm(targets[:hv]),
                q_vv_norm_m=norm(targets[:vv]), q_total_norm_m=norm(total_target),
                hh_relative_closure=relative_vector_residual(targets[:hh], reconstructions[:hh]),
                hv_relative_closure=relative_vector_residual(targets[:hv], reconstructions[:hv]),
                vv_relative_closure=relative_vector_residual(targets[:vv], reconstructions[:vv]),
                total_relative_closure=relative_vector_residual(total_target, total_reconstruction),
                total_signed_projection=dot(total_reconstruction, total_target) / sum(abs2, total_target)))
        end
    end

    element_rows = NamedTuple[]
    for source in 1:n_sources
        push!(element_rows, (; element_order=source,
            element_name=kick_response.names[source], s_m=kick_response.s_m[source],
            k2l_m2=kick_response.k2l[source],
            mean_exposure_h=exposure_h_sum[source] / trials,
            mean_exposure_v=exposure_v_sum[source] / trials,
            mean_exposure_excess_v_minus_h=(exposure_v_sum[source] - exposure_h_sum[source]) / trials,
            eta_hh=numerator[:hh][source] / denominator[:hh],
            eta_hv=numerator[:hv][source] / denominator[:hv],
            eta_vv=numerator[:vv][source] / denominator[:vv],
            eta_total=total_numerator[source] / total_denominator,
            magnitude_hh=sqrt(magnitude_sq[:hh][source] / denominator[:hh]),
            magnitude_hv=sqrt(magnitude_sq[:hv][source] / denominator[:hv]),
            magnitude_vv=sqrt(magnitude_sq[:vv][source] / denominator[:vv]),
            magnitude_total=sqrt(total_magnitude_sq[source] / total_denominator)))
    end

    summary_row = (; trials, active_normal_sextupoles=n_sources, detectors=length(detectors),
        internal_response_closure=internal.first_residual,
        local_kick_solve_closure=kick_response.solve_closure,
        local_kick_fixed_point_closure=kick_response.fixed_point_closure,
        hh_concatenated_relative_closure=sqrt(residual_sq[:hh] / denominator[:hh]),
        hv_concatenated_relative_closure=sqrt(residual_sq[:hv] / denominator[:hv]),
        vv_concatenated_relative_closure=sqrt(residual_sq[:vv] / denominator[:vv]),
        total_concatenated_relative_closure=sqrt(total_residual_sq / total_denominator),
        hh_reconstruction_signed_projection=reconstruction_dot[:hh] / denominator[:hh],
        hv_reconstruction_signed_projection=reconstruction_dot[:hv] / denominator[:hv],
        vv_reconstruction_signed_projection=reconstruction_dot[:vv] / denominator[:vv],
        total_reconstruction_signed_projection=total_reconstruction_dot / total_denominator,
        sum_eta_hh=sum(numerator[:hh]) / denominator[:hh],
        sum_eta_hv=sum(numerator[:hv]) / denominator[:hv],
        sum_eta_vv=sum(numerator[:vv]) / denominator[:vv],
        sum_eta_total=sum(total_numerator) / total_denominator, solve_seconds)

    mkpath(output_dir)
    element_path = write_namedtuple_csv(joinpath(output_dir, "sextupole_contribution_summary.csv"), element_rows)
    direction_path = write_namedtuple_csv(joinpath(output_dir, "direction_closure.csv"), direction_rows)
    element_direction_path = write_namedtuple_csv(joinpath(output_dir, "sextupole_direction_contributions.csv"), element_direction_rows)
    summary_path = write_namedtuple_csv(joinpath(output_dir, "reconstruction_summary.csv"), [summary_row])
    detector_path = write_namedtuple_csv(joinpath(output_dir, "detectors.csv"), [
        (; detector_order=index, detector_name=detectors[index], s_m=detector_s_m[index])
        for index in eachindex(detectors)
    ])
    metadata_path = joinpath(output_dir, "metadata.toml")
    open(metadata_path, "w") do io
        TOML.print(io, Dict(
            "format" => "cesr-signed-sextupole-detector-contributions-v1",
            "date" => string(Dates.today()), "trials" => trials, "seed" => seed,
            "base_kick_rad" => base_kick,
            "target" => "Q_x = Q_hh,x + Q_hv,x + Q_vv,x at 99 detectors",
            "local_source" => "signed thin normal-sextupole quadratic kick split equally between element entrance and exit",
            "projection" => "ensemble dot(C_j,Q)/ensemble norm(Q)^2",
        ); sorted=true)
    end

    @printf("Solve time: %.3f s\n", solve_seconds)
    println("Elements:           $element_path")
    println("Directions:         $direction_path")
    println("Element directions: $element_direction_path")
    println("Summary:            $summary_path")
    println("Detectors:          $detector_path")
    println("Metadata:           $metadata_path")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main_sextupole_detector_contributions())
end

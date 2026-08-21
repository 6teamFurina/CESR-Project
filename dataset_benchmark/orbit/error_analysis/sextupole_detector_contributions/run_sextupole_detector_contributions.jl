#!/usr/bin/env julia

"""Latest-ring signed normal-sextupole sourcing of the nonlinear orbit error.

The maintained target is the total second-order detector vector obtained from
the two-parameter implicit closed-orbit derivative.  Local sources are the
exact thick-element Hessian sources

    g_j = S_exit,j - A_j S_entrance,j,

at every complete SciBmad element.  The reported element contribution keeps
only the active normal-sextupole elements; the all-element reconstruction is
retained as a closure check.

This file deliberately does not include or call the historical internal
exposure runner.  The latest lattice is initialized with primitive inactive
controls and only selected steering controls receive GTPSA parameters.  This
avoids promoting inactive combined-multipole controls at their zero operating
point.
"""

const SEXT_DETECTOR_HERE = @__DIR__
const MIXED_TERM_HERE = normpath(joinpath(SEXT_DETECTOR_HERE, "..", "mixed_terms"))
include(joinpath(MIXED_TERM_HERE, "run_mixed_term_gtpsa.jl"))

using Dates
using LinearAlgebra
using Printf
using TOML

function parse_sext_detector_args(args)
    options = Dict{String,String}(
        "ring" => "latest",
        "trials" => "100",
        "seed" => "20260804",
        "base-kick-rad" => "5e-6",
        "output-dir" => "",
        "inputs" => "",
    )
    for argument in args
        startswith(argument, "--") || error("Arguments must have --name=value form: $argument")
        fields = split(argument[3:end], '='; limit=2)
        length(fields) == 2 || error("Missing value in argument: $argument")
        haskey(options, fields[1]) || error("Unknown option: --$(fields[1])")
        options[fields[1]] = fields[2]
    end
    ring = Symbol(lowercase(options["ring"]))
    ring in (:latest, :latest_cesr, :repaired_latest, :legacy, :legacy_cesr, :historical) ||
        error("--ring must be latest or legacy")
    isempty(options["inputs"]) && (options["inputs"] = default_ring_paths(; ring).inputs)
    isempty(options["output-dir"]) &&
        (options["output-dir"] = joinpath(SEXT_DETECTOR_HERE, "results", ring_artifact_id(ring)))
    return options
end

constant_term(value) = Float64(GTPSA.scalar(Beamlines.deval(value)))

"""Discover active normal Kn2/Kn2L elements without a fixed CESR count."""
function latest_normal_sextupole_inventory(ring)
    inventory = Dict{Int,NamedTuple}()
    s_start = 0.0
    for (index, element) in enumerate(ring.line)
        length_m = constant_term(element.L)
        kn2 = constant_term(element.Kn2)
        kn2l = constant_term(element.Kn2L)
        k2l = iszero(kn2) ? kn2l : kn2 * length_m
        if !iszero(k2l)
            inventory[index] = (;
                name=String(element.name),
                k2l,
                s_m=s_start + length_m,
            )
        end
        s_start += length_m
    end
    isempty(inventory) && error("No active normal sextupoles were found")
    return inventory
end

"""Build first-order internal orbit responses with only selected controls typed."""
function build_selected_internal_response(
    names,
    closed_orbit,
    base_kick,
    inventory;
    model_factory=load_ring_model,
    config=nothing,
)
    state_dimension = length(closed_orbit)
    descriptor = Descriptor(state_dimension, 1, length(names), 1)
    variables = vars(descriptor)
    parameters = params(descriptor)
    # Keep all inactive Overlay/Group controls primitive.  Only selected
    # steering controls below carry GTPSA parameters.
    model = configured_model(model_factory, config; zero_value=0.0, rf_on=true)
    for (index, name) in enumerate(names)
        model.controls[name] = base_kick * parameters[index]
    end
    input_map = [closed_orbit[index] + copy(variables[index]) for index in 1:state_dimension]
    bunch = Bunch(v=reshape(input_map, 1, state_dimension))
    SciBmad.BTBL.check_bl_bunch!(bunch, model.ring, false)
    element_indices = sort!(collect(keys(inventory)))
    source_row = Dict(index => row for (row, index) in enumerate(element_indices))
    maps = NamedTuple[]
    for (element_index, element) in enumerate(model.ring.line)
        active = get(source_row, element_index, 0)
        entrance = active == 0 ? nothing : copy.(vec(bunch.coords.v))
        track!(bunch, element)
        if active != 0
            exit_map = copy.(vec(bunch.coords.v))
            midpoint = [(entrance[index] + exit_map[index]) / 2 for index in 1:state_dimension]
            push!(maps, (; element_index, map=midpoint))
        end
    end
    one_turn_map = copy.(vec(bunch.coords.v))
    jacobian = Matrix(GTPSA.jacobian(one_turn_map; include_params=true))
    expected_columns = state_dimension + length(names)
    size(jacobian) == (state_dimension, expected_columns) ||
        error("Unexpected selected-control internal Jacobian size $(size(jacobian))")
    A = jacobian[:, 1:state_dimension]
    B = jacobian[:, state_dimension + 1:end]
    first = (I - A) \ B
    lifted = vcat(first, Matrix{Float64}(I, length(names), length(names)))
    first_residual = norm((I - A) * first - B, Inf)
    transverse = transverse_coordinate_indices(config; state_dimension)
    response_x = zeros(length(maps), length(names))
    response_y = similar(response_x)
    for (row, item) in enumerate(maps)
        item_jacobian = Matrix(GTPSA.jacobian(item.map; include_params=true))
        size(item_jacobian) == (state_dimension, expected_columns) ||
            error("Unexpected selected-control internal map Jacobian size")
        response_x[row, :] .= vec(view(item_jacobian, transverse.x, :)' * lifted)
        response_y[row, :] .= vec(view(item_jacobian, transverse.y, :)' * lifted)
    end
    return (; response_x, response_y, first_residual)
end

"""Build cumulative and local linear maps with primitive lattice controls."""
function build_linear_lattice_data(
    closed_orbit,
    detectors;
    model_factory=load_ring_model,
    config=nothing,
)
    state_dimension = length(closed_orbit)
    descriptor = Descriptor(state_dimension, 1)
    variables = vars(descriptor)
    model = configured_model(model_factory, config; zero_value=0.0, rf_on=true)
    input_map = [closed_orbit[index] + copy(variables[index]) for index in 1:state_dimension]
    bunch = Bunch(v=reshape(input_map, 1, state_dimension))
    SciBmad.BTBL.check_bl_bunch!(bunch, model.ring, false)

    n_elements = length(model.ring.line)
    local_maps = Vector{Matrix{Float64}}(undef, n_elements)
    cumulative_maps = Vector{Matrix{Float64}}(undef, n_elements)
    element_names = Vector{String}(undef, n_elements)
    element_types = Vector{String}(undef, n_elements)
    element_s_m = zeros(n_elements)
    detector_lookup = Dict(name => row for (row, name) in enumerate(detectors))
    detector_element_indices = zeros(Int, length(detectors))
    detector_maps = Vector{Matrix{Float64}}(undef, length(detectors))
    found = falses(length(detectors))
    previous = Matrix{Float64}(I, state_dimension, state_dimension)
    s_m = 0.0
    for (element_index, element) in enumerate(model.ring.line)
        track!(bunch, element)
        current = Matrix(GTPSA.jacobian(copy.(vec(bunch.coords.v)); include_params=true))
        size(current) == (state_dimension, state_dimension) ||
            error("Unexpected cumulative map size $(size(current))")
        local_maps[element_index] = current / previous
        cumulative_maps[element_index] = current
        previous = current
        length_m = constant_term(element.L)
        s_m += length_m
        element_names[element_index] = String(element.name)
        element_types[element_index] = string(element.kind)
        element_s_m[element_index] = s_m
        name = uppercase(String(element.name))
        haskey(detector_lookup, name) || continue
        detector = detector_lookup[name]
        found[detector] && error("Detector $name occurs more than once")
        detector_element_indices[detector] = element_index
        detector_maps[detector] = current
        found[detector] = true
    end
    all(found) || error("Missing detector maps: $(detectors[.!found])")
    return (; local_maps, cumulative_maps, one_turn_map=previous,
        element_names, element_types, element_s_m,
        detector_element_indices, detector_maps)
end

"""Periodic response from one element-exit source to detector coordinates."""
function build_element_source_response(linear, element_indices, output_coordinate)
    fixed_point = I - linear.one_turn_map
    responses = Matrix{Float64}[]
    for element_index in element_indices
        exit_map = linear.cumulative_maps[element_index]
        downstream_to_end = linear.one_turn_map / exit_map
        closed_start = fixed_point \ downstream_to_end
        response = zeros(length(linear.detector_maps), size(linear.one_turn_map, 1))
        for detector in eachindex(linear.detector_maps)
            detector_map = linear.detector_maps[detector]
            total = detector_map * closed_start
            if linear.detector_element_indices[detector] >= element_index
                total += detector_map / exit_map
            end
            response[detector, :] .= view(total, output_coordinate, :)
        end
        push!(responses, response)
    end
    return responses
end

"""Reconstruct a detector vector from local source states."""
function reconstruct_sources(linear, sources, output_coordinate, mask=nothing)
    state_dimension = size(linear.one_turn_map, 1)
    size(sources) == (state_dimension, length(linear.local_maps)) ||
        error("Source matrix size mismatch")
    isnothing(mask) || length(mask) == length(linear.local_maps) ||
        error("Source mask size mismatch")
    end_source = zeros(state_dimension)
    for element_index in eachindex(linear.local_maps)
        end_source = linear.local_maps[element_index] * end_source
        if isnothing(mask) || mask[element_index]
            end_source += view(sources, :, element_index)
        end
    end
    initial = (I - linear.one_turn_map) \ end_source
    detector_values = zeros(length(linear.detector_maps))
    detector_by_element = Dict(
        index => detector for (detector, index) in enumerate(linear.detector_element_indices)
    )
    state = initial
    for element_index in eachindex(linear.local_maps)
        state = linear.local_maps[element_index] * state
        if isnothing(mask) || mask[element_index]
            state += view(sources, :, element_index)
        end
        detector = get(detector_by_element, element_index, 0)
        iszero(detector) || (detector_values[detector] = state[output_coordinate])
    end
    return detector_values, norm(state - initial, Inf)
end

"""First- and second-order derivatives after closed-orbit substitution."""
function closed_map_derivatives(map, closed_derivatives)
    state_dimension = size(closed_derivatives.lifted, 1) - 2
    jacobian, hessians = map_derivatives(map, state_dimension)
    first = jacobian * closed_derivatives.lifted
    second = zeros(state_dimension, 2, 2)
    for coordinate in 1:state_dimension, left in 1:2, right in left:2
        value = dot(
            view(jacobian, coordinate, 1:state_dimension),
            view(closed_derivatives.second, :, left, right),
        ) + dot(
            view(closed_derivatives.lifted, :, left),
            hessians[coordinate] * view(closed_derivatives.lifted, :, right),
        )
        second[coordinate, left, right] = value
        second[coordinate, right, left] = value
    end
    return (; first, second)
end

function one_turn_direction_map(ring, input_map)
    bunch = Bunch(v=reshape(copy.(input_map), 1, length(input_map)))
    SciBmad.BTBL.check_bl_bunch!(bunch, ring, false)
    for element in ring.line
        track!(bunch, element)
    end
    return copy.(vec(bunch.coords.v))
end

"""Exact all-element source states and x/y targets for one direction pair."""
function direction_thick_sources_both(
    names,
    detectors,
    closed_orbit,
    h_direction,
    v_direction,
    base_kick,
    linear;
    model_factory=load_ring_model,
    config=nothing,
)
    state_dimension = length(closed_orbit)
    setup = direction_parameterized_model(
        names, h_direction, v_direction, base_kick, state_dimension;
        model_factory, config,
    )
    input_map = [closed_orbit[index] + copy(setup.variables[index]) for index in 1:state_dimension]
    one_turn_map = one_turn_direction_map(setup.model.ring, input_map)
    closed_derivatives = implicit_closed_orbit_derivatives(one_turn_map, state_dimension)
    bunch = Bunch(v=reshape(copy.(input_map), 1, state_dimension))
    SciBmad.BTBL.check_bl_bunch!(bunch, setup.model.ring, false)
    entrance = closed_map_derivatives(input_map, closed_derivatives)
    n_elements = length(setup.model.ring.line)
    sources = zeros(state_dimension, n_elements)
    target_x = zeros(length(detectors))
    target_y = zeros(length(detectors))
    detector_lookup = Dict(name => row for (row, name) in enumerate(detectors))
    transverse = transverse_coordinate_indices(config; state_dimension)
    found = falses(length(detectors))

    for (element_index, element) in enumerate(setup.model.ring.line)
        track!(bunch, element)
        exit_map = copy.(vec(bunch.coords.v))
        exit_derivatives = closed_map_derivatives(exit_map, closed_derivatives)
        local_matrix = linear.local_maps[element_index]
        hh = 0.5 .* (
            view(exit_derivatives.second, :, 1, 1) -
            local_matrix * view(entrance.second, :, 1, 1)
        )
        hv = (
            view(exit_derivatives.second, :, 1, 2) -
            local_matrix * view(entrance.second, :, 1, 2)
        )
        vv = 0.5 .* (
            view(exit_derivatives.second, :, 2, 2) -
            local_matrix * view(entrance.second, :, 2, 2)
        )
        sources[:, element_index] .= hh + hv + vv

        name = uppercase(String(element.name))
        if haskey(detector_lookup, name)
            detector = detector_lookup[name]
            found[detector] && error("Detector $name occurs more than once")
            target_x[detector] =
                0.5 * exit_derivatives.second[transverse.x, 1, 1] +
                exit_derivatives.second[transverse.x, 1, 2] +
                0.5 * exit_derivatives.second[transverse.x, 2, 2]
            target_y[detector] =
                0.5 * exit_derivatives.second[transverse.y, 1, 1] +
                exit_derivatives.second[transverse.y, 1, 2] +
                0.5 * exit_derivatives.second[transverse.y, 2, 2]
            found[detector] = true
        end
        entrance = exit_derivatives
    end
    all(found) || error("Missing detector targets: $(detectors[.!found])")
    return (; sources, target_x, target_y,
        first_closure_residual=closed_derivatives.first_residual,
        second_closure_residual=closed_derivatives.second_residual)
end

function relative_residual(target, reconstruction)
    denominator = norm(target)
    return denominator == 0 ? NaN : norm(target - reconstruction) / denominator
end

function local_thin_source_kicks(k2l, x_h, x_v, y_h, y_v)
    px = @. -0.5 * k2l * (x_h^2 - y_h^2) - k2l * (x_h * x_v - y_h * y_v) -
        0.5 * k2l * (x_v^2 - y_v^2)
    py = @. k2l * x_h * y_h + k2l * (x_h * y_v + x_v * y_h) + k2l * x_v * y_v
    return px, py
end

function main_sextupole_detector_contributions(args=ARGS; model_factory=nothing, config=nothing)
    options = parse_sext_detector_args(args)
    ring = Symbol(lowercase(options["ring"]))
    model_factory, config = resolve_ring_model_factory(model_factory, config; ring)
    trials = parse(Int, options["trials"])
    trials >= 1 || error("--trials must be positive")
    seed = parse(Int, options["seed"])
    base_kick = parse(Float64, options["base-kick-rad"])
    base_kick > 0 || error("--base-kick-rad must be positive")
    output_dir = abspath(options["output-dir"])
    input_path = abspath(options["inputs"])
    input_reference = read_samples(input_path)
    names = input_reference.names
    validate_control_names(names, config)

    float_model = configured_model(model_factory, config; zero_value=0.0, rf_on=true)
    detectors = configured_detector_names(float_model, config)
    closed_orbit = nominal_closed_orbit(float_model.ring)
    state_dimension = length(closed_orbit)
    inventory = latest_normal_sextupole_inventory(float_model.ring)
    sextupole_indices = sort!(collect(keys(inventory)))
    samples = generate_mixed_samples(
        names, [1.0], trials, seed, base_kick;
        config, model=float_model,
    )
    transverse = transverse_coordinate_indices(config; state_dimension)
    momenta = transverse_momentum_indices(config; state_dimension)

    @printf("Building selected-control internal map at %d active normal sextupoles...\n", length(inventory))
    internal = build_selected_internal_response(
        names, closed_orbit, base_kick, inventory;
        model_factory, config,
    )
    @printf("Building primitive-control linear maps for %d elements and %d detectors...\n",
        length(float_model.ring.line), length(detectors))
    linear = build_linear_lattice_data(closed_orbit, detectors; model_factory, config)
    response_x = build_element_source_response(linear, sextupole_indices, transverse.x)
    response_y = build_element_source_response(linear, sextupole_indices, transverse.y)

    n_sources = length(sextupole_indices)
    n_detectors = length(detectors)
    total_denominator = 0.0
    x_denominator = 0.0
    y_denominator = 0.0
    total_residual_sq = 0.0
    x_residual_sq = 0.0
    y_residual_sq = 0.0
    total_sext_residual_sq = 0.0
    x_sext_residual_sq = 0.0
    y_sext_residual_sq = 0.0
    total_reconstruction_dot = 0.0
    x_reconstruction_dot = 0.0
    y_reconstruction_dot = 0.0
    total_numerator = zeros(n_sources)
    x_numerator = zeros(n_sources)
    y_numerator = zeros(n_sources)
    total_magnitude_sq = zeros(n_sources)
    x_magnitude_sq = zeros(n_sources)
    y_magnitude_sq = zeros(n_sources)
    direction_rows = NamedTuple[]
    element_direction_rows = NamedTuple[]
    max_periodic_source_closure = 0.0
    max_first_closure = 0.0
    max_second_closure = 0.0
    max_sext_response_check = 0.0

    solve_seconds = @elapsed begin
        for trial in 1:trials
            @printf("Latest signed-source direction %d/%d\n", trial, trials)
            h_direction = view(samples.horizontal_directions, trial, :)
            v_direction = view(samples.vertical_directions, trial, :)
            result = direction_thick_sources_both(
                names, detectors, closed_orbit,
                h_direction, v_direction, base_kick, linear;
                model_factory, config,
            )
            max_first_closure = max(max_first_closure, result.first_closure_residual)
            max_second_closure = max(max_second_closure, result.second_closure_residual)

            all_x, closure_x = reconstruct_sources(linear, result.sources, transverse.x)
            all_y, closure_y = reconstruct_sources(linear, result.sources, transverse.y)
            max_periodic_source_closure = max(max_periodic_source_closure, closure_x, closure_y)
            target = vcat(result.target_x, result.target_y)
            all_reconstruction = vcat(all_x, all_y)

            source_mask = falses(length(linear.local_maps))
            source_mask[sextupole_indices] .= true
            sext_x, _ = reconstruct_sources(linear, result.sources, transverse.x, source_mask)
            sext_y, _ = reconstruct_sources(linear, result.sources, transverse.y, source_mask)
            sext_reconstruction = vcat(sext_x, sext_y)
            expected_sext_x = zeros(n_detectors)
            expected_sext_y = zeros(n_detectors)
            for source in 1:n_sources
                element_index = sextupole_indices[source]
                expected_sext_x .+= response_x[source] * view(result.sources, :, element_index)
                expected_sext_y .+= response_y[source] * view(result.sources, :, element_index)
            end
            max_sext_response_check = max(
                max_sext_response_check,
                maximum(abs, sext_x - expected_sext_x),
                maximum(abs, sext_y - expected_sext_y),
            )

            x_h = internal.response_x * h_direction
            x_v = internal.response_x * v_direction
            y_h = internal.response_y * h_direction
            y_v = internal.response_y * v_direction
            px_source, py_source = local_thin_source_kicks(
                [inventory[index].k2l for index in sextupole_indices],
                x_h, x_v, y_h, y_v,
            )

            total_denominator += sum(abs2, target)
            x_denominator += sum(abs2, result.target_x)
            y_denominator += sum(abs2, result.target_y)
            total_residual_sq += sum(abs2, target - all_reconstruction)
            x_residual_sq += sum(abs2, result.target_x - all_x)
            y_residual_sq += sum(abs2, result.target_y - all_y)
            total_sext_residual_sq += sum(abs2, target - sext_reconstruction)
            x_sext_residual_sq += sum(abs2, result.target_x - sext_x)
            y_sext_residual_sq += sum(abs2, result.target_y - sext_y)
            total_reconstruction_dot += dot(all_reconstruction, target)
            x_reconstruction_dot += dot(all_x, result.target_x)
            y_reconstruction_dot += dot(all_y, result.target_y)

            for source in 1:n_sources
                element_index = sextupole_indices[source]
                local_source = view(result.sources, :, element_index)
                contribution_x = response_x[source] * local_source
                contribution_y = response_y[source] * local_source
                contribution = vcat(contribution_x, contribution_y)
                projection = dot(contribution, target)
                total_numerator[source] += projection
                x_numerator[source] += dot(contribution_x, result.target_x)
                y_numerator[source] += dot(contribution_y, result.target_y)
                total_magnitude_sq[source] += sum(abs2, contribution)
                x_magnitude_sq[source] += sum(abs2, contribution_x)
                y_magnitude_sq[source] += sum(abs2, contribution_y)
                push!(element_direction_rows, (; trial,
                    element_order=source,
                    element_index,
                    element_name=String(inventory[element_index].name),
                    s_m=Float64(inventory[element_index].s_m),
                    k2l_m2=Float64(inventory[element_index].k2l),
                    source_kick_x_rad=Float64(px_source[source]),
                    source_kick_y_rad=Float64(py_source[source]),
                    source_kick_norm_rad=norm([local_source[momenta.px], local_source[momenta.py]]),
                    x_projection_numerator=dot(contribution_x, result.target_x),
                    y_projection_numerator=dot(contribution_y, result.target_y),
                    total_projection_numerator=projection,
                    x_contribution_norm_m=norm(contribution_x),
                    y_contribution_norm_m=norm(contribution_y),
                    total_contribution_norm_m=norm(contribution)))
            end

            push!(direction_rows, (; trial,
                q_x_norm_m=norm(result.target_x),
                q_y_norm_m=norm(result.target_y),
                q_total_norm_m=norm(target),
                all_element_relative_closure=relative_residual(target, all_reconstruction),
                normal_sextupole_relative_closure=relative_residual(target, sext_reconstruction),
                all_element_signed_projection=dot(all_reconstruction, target) / sum(abs2, target),
                normal_sextupole_signed_projection=dot(sext_reconstruction, target) / sum(abs2, target),
                first_fixed_point_closure=result.first_closure_residual,
                second_fixed_point_closure=result.second_closure_residual,
                periodic_source_closure_x=closure_x,
                periodic_source_closure_y=closure_y))
        end
    end

    element_rows = NamedTuple[]
    for source in 1:n_sources
        element_index = sextupole_indices[source]
        push!(element_rows, (; element_order=source,
            element_index,
            element_name=String(inventory[element_index].name),
            s_m=Float64(inventory[element_index].s_m),
            k2l_m2=Float64(inventory[element_index].k2l),
            eta_x=x_numerator[source] / x_denominator,
            eta_y=y_numerator[source] / y_denominator,
            eta_total=total_numerator[source] / total_denominator,
            magnitude_x=sqrt(x_magnitude_sq[source] / x_denominator),
            magnitude_y=sqrt(y_magnitude_sq[source] / y_denominator),
            magnitude_total=sqrt(total_magnitude_sq[source] / total_denominator)))
    end

    summary = (; trials,
        elements=length(linear.local_maps),
        active_normal_sextupoles=n_sources,
        detectors=n_detectors,
        state_dimension,
        control_count=length(names),
        internal_response_closure=internal.first_residual,
        first_fixed_point_closure_max=max_first_closure,
        second_fixed_point_closure_max=max_second_closure,
        periodic_source_closure_max=max_periodic_source_closure,
        sextupole_response_check_max=max_sext_response_check,
        all_element_relative_closure=sqrt(total_residual_sq / total_denominator),
        normal_sextupole_relative_closure=sqrt(total_sext_residual_sq / total_denominator),
        x_all_element_relative_closure=sqrt(x_residual_sq / x_denominator),
        y_all_element_relative_closure=sqrt(y_residual_sq / y_denominator),
        x_normal_sextupole_relative_closure=sqrt(x_sext_residual_sq / x_denominator),
        y_normal_sextupole_relative_closure=sqrt(y_sext_residual_sq / y_denominator),
        all_element_signed_projection=total_reconstruction_dot / total_denominator,
        x_all_element_signed_projection=x_reconstruction_dot / x_denominator,
        y_all_element_signed_projection=y_reconstruction_dot / y_denominator,
        normal_sextupole_signed_projection=sum(total_numerator) / total_denominator,
        x_normal_sextupole_signed_projection=sum(x_numerator) / x_denominator,
        y_normal_sextupole_signed_projection=sum(y_numerator) / y_denominator,
        solve_seconds)

    mkpath(output_dir)
    element_path = write_namedtuple_csv(joinpath(output_dir, "sextupole_contribution_summary.csv"), element_rows)
    direction_path = write_namedtuple_csv(joinpath(output_dir, "direction_closure.csv"), direction_rows)
    element_direction_path = write_namedtuple_csv(joinpath(output_dir, "sextupole_direction_contributions.csv"), element_direction_rows)
    summary_path = write_namedtuple_csv(joinpath(output_dir, "reconstruction_summary.csv"), [summary])
    detector_s = zeros(n_detectors)
    s_accum = 0.0
    detector_lookup = Dict(name => index for (index, name) in enumerate(detectors))
    for element in float_model.ring.line
        s_accum += constant_term(element.L)
        name = uppercase(String(element.name))
        haskey(detector_lookup, name) && (detector_s[detector_lookup[name]] = s_accum)
    end
    detector_path = write_namedtuple_csv(joinpath(output_dir, "detectors.csv"), [
        (; detector_order=index, detector_name=detectors[index], s_m=detector_s[index])
        for index in eachindex(detectors)
    ])
    metadata_path = joinpath(output_dir, "metadata.toml")
    open(metadata_path, "w") do io
        metadata = merge(ring_metadata(config; ring), Dict(
            "format" => "cesr-latest-signed-sextupole-detector-contributions-v2",
            "date" => string(Dates.today()),
            "trials" => trials,
            "seed" => seed,
            "base_kick_rad" => base_kick,
            "input_csv" => input_path,
            "control_count" => length(names),
            "control_names" => names,
            "state_dimension" => state_dimension,
            "element_count" => length(linear.local_maps),
            "detector_count" => n_detectors,
            "detector_names" => detectors,
            "observable_labels" => vcat(["$(name):x" for name in detectors], ["$(name):y" for name in detectors]),
            "active_normal_sextupoles" => n_sources,
            "active_normal_sextupole_names" => [String(inventory[index].name) for index in sextupole_indices],
            "target" => "total second-order nonlinear detector vector Q=(Q_x,Q_y)",
            "source_boundary" => "exact complete-element exit boundary g_j = S_exit,j - A_j*S_entrance,j",
            "source_convention" => "normal-sextupole elements retained from exact all-element Hessian source; no hh/hv/vv outputs",
            "descriptor" => "Descriptor($(state_dimension), 2, 2, 2) for direction derivatives; primitive controls for linear maps",
            "phase_space_coordinates" => ["x", "px", "y", "py", "z", "pz"],
            "projection" => "signed ensemble dot(C_j,Q)/ensemble norm(Q)^2",
            "units" => Dict("orbit" => "m", "source_kick" => "rad", "k2l" => "m^-2", "s" => "m"),
            "results_csv" => element_path,
            "direction_csv" => direction_path,
            "element_direction_csv" => element_direction_path,
            "detector_csv" => detector_path,
        ))
        TOML.print(io, metadata; sorted=true)
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

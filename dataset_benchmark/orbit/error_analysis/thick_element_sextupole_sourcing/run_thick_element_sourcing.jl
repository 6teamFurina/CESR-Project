#!/usr/bin/env julia

"""Attribute a summed nonlinear detector target to complete-element sources."""

const THICK_SOURCE_HERE = @__DIR__
const QUADRATIC_ATTRIBUTION_HERE = normpath(joinpath(
    THICK_SOURCE_HERE, "..", "quadratic_x_attribution",
))
include(joinpath(QUADRATIC_ATTRIBUTION_HERE, "run_internal_exposure_attribution.jl"))

using Dates
using LinearAlgebra
using Printf
using TOML

function parse_thick_source_args(args)
    options = Dict{String,String}(
        "ring" => "latest",
        "trials" => "100",
        "seed" => "20260804",
        "base-kick-rad" => "5e-6",
        "output-plane" => "x",
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
    output_plane = lowercase(options["output-plane"])
    output_plane in ("x", "y") || error("--output-plane must be x or y")
    isempty(options["inputs"]) &&
        (options["inputs"] = default_ring_paths(; ring).inputs)
    isempty(options["output-dir"]) &&
        (options["output-dir"] = joinpath(
            THICK_SOURCE_HERE,
            output_plane == "x" ? "horizontal_results" : "vertical_results",
            ring_artifact_id(ring),
        ))
    return options
end

"""First- and second-order derivatives after substituting the closed orbit."""
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

"""Baseline local and cumulative first-order maps at every element exit."""
function build_linear_lattice_data(
    closed_orbit,
    detectors;
    model_factory=load_ring_model,
    config=nothing,
)
    state_dimension = length(closed_orbit)
    descriptor = Descriptor(state_dimension, 1)
    variables = vars(descriptor)
    # This is a phase-space-only TPSA map.  Keep every unselected lattice
    # control primitive; promoting the complete control registry with
    # zero(variables[1]) can evaluate combined-multipole expressions such as
    # SEX_14W at sqrt(0) in the wrong numerical domain.
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
        size(current) == (state_dimension, state_dimension) || error("Unexpected cumulative map size $(size(current))")
        local_maps[element_index] = current / previous
        cumulative_maps[element_index] = current
        previous = current
        length_m = Float64(GTPSA.scalar(Beamlines.deval(element.L)))
        s_m += length_m
        element_names[element_index] = String(element.name)
        element_types[element_index] = string(element.kind)
        element_s_m[element_index] = s_m - length_m / 2
        name = uppercase(String(element.name))
        haskey(detector_lookup, name) || continue
        row = detector_lookup[name]
        found[row] && error("Detector $name occurs more than once")
        detector_element_indices[row] = element_index
        detector_maps[row] = current
        found[row] = true
    end
    all(found) || error("Missing detector maps: $(detectors[.!found])")
    return (;
        local_maps, cumulative_maps, one_turn_map=previous,
        element_names, element_types, element_s_m,
        detector_element_indices, detector_maps,
    )
end

"""Periodic detector response to a six-dimensional source at element exit."""
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

"""Reconstruct one detector coordinate from arbitrary local six-dimensional sources."""
function reconstruct_sources(linear, sources, output_coordinate, mask=nothing)
    state_dimension = size(linear.one_turn_map, 1)
    size(sources) == (state_dimension, length(linear.local_maps)) || error("Source matrix size mismatch")
    isnothing(mask) || length(mask) == length(linear.local_maps) || error("Source mask size mismatch")
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
    closure = norm(state - initial, Inf)
    return detector_values, closure
end

function source_family(element_type, element_index, inventory)
    haskey(inventory, element_index) && return "normal_sextupole"
    normalized = lowercase(strip(String(element_type)))
    normalized == "sextupole" && return "other_sextupole"
    return normalized
end

"""Exact local Hessian sources and detector targets for one direction pair."""
function direction_thick_sources(
    names, detectors, closed_orbit, h_direction, v_direction, base_kick, linear,
    output_coordinate;
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
    sources = Dict(
        :hh => zeros(state_dimension, n_elements),
        :hv => zeros(state_dimension, n_elements),
        :vv => zeros(state_dimension, n_elements),
    )
    targets = Dict(:hh => zeros(length(detectors)),
                   :hv => zeros(length(detectors)),
                   :vv => zeros(length(detectors)))
    detector_lookup = Dict(name => row for (row, name) in enumerate(detectors))
    found = falses(length(detectors))

    for (element_index, element) in enumerate(setup.model.ring.line)
        track!(bunch, element)
        exit_map = copy.(vec(bunch.coords.v))
        exit_derivatives = closed_map_derivatives(exit_map, closed_derivatives)
        local_matrix = linear.local_maps[element_index]
        sources[:hh][:, element_index] .= 0.5 .* (
            view(exit_derivatives.second, :, 1, 1) -
            local_matrix * view(entrance.second, :, 1, 1)
        )
        sources[:hv][:, element_index] .= (
            view(exit_derivatives.second, :, 1, 2) -
            local_matrix * view(entrance.second, :, 1, 2)
        )
        sources[:vv][:, element_index] .= 0.5 .* (
            view(exit_derivatives.second, :, 2, 2) -
            local_matrix * view(entrance.second, :, 2, 2)
        )

        name = uppercase(String(element.name))
        if haskey(detector_lookup, name)
            detector = detector_lookup[name]
            found[detector] && error("Detector $name occurs more than once")
            targets[:hh][detector] = 0.5 * exit_derivatives.second[output_coordinate, 1, 1]
            targets[:hv][detector] = exit_derivatives.second[output_coordinate, 1, 2]
            targets[:vv][detector] = 0.5 * exit_derivatives.second[output_coordinate, 2, 2]
            found[detector] = true
        end
        entrance = exit_derivatives
    end
    all(found) || error("Missing detector targets: $(detectors[.!found])")
    return (; sources, targets,
        first_closure_residual=closed_derivatives.first_residual,
        second_closure_residual=closed_derivatives.second_residual)
end

function relative_residual(target, reconstruction)
    denominator = norm(target)
    return denominator == 0 ? NaN : norm(target - reconstruction) / denominator
end

function assert_finite_rows(label, rows)
    isempty(rows) && error("No rows were produced for $label")
    for (row_index, row) in enumerate(rows), column in propertynames(row)
        value = getproperty(row, column)
        value isa Real || continue
        isfinite(Float64(value)) || error(
            "Non-finite $label value at row $row_index, column $column: $value",
        )
    end
    return nothing
end

function main_thick_element_sourcing(args=ARGS; model_factory=nothing, config=nothing)
    options = parse_thick_source_args(args)
    ring = Symbol(lowercase(options["ring"]))
    model_factory, config = resolve_ring_model_factory(model_factory, config; ring)
    trials = parse(Int, options["trials"])
    trials >= 1 || error("--trials must be positive")
    seed = parse(Int, options["seed"])
    base_kick = parse(Float64, options["base-kick-rad"])
    base_kick > 0 || error("--base-kick-rad must be positive")
    output_plane = lowercase(options["output-plane"])
    output_dir = abspath(options["output-dir"])

    input_path = abspath(options["inputs"])
    reference = read_samples(input_path)
    names = reference.names
    validate_control_names(names, config)
    float_model = configured_model(model_factory, config; zero_value=0.0, rf_on=true)
    detectors = configured_detector_names(float_model, config)
    closed_orbit = nominal_closed_orbit(float_model.ring)
    transverse = transverse_coordinate_indices(config; state_dimension=length(closed_orbit))
    output_coordinate = output_plane == "x" ? transverse.x : output_plane == "y" ? transverse.y :
        error("--output-plane must be x or y")
    inventory = normal_sextupole_inventory(float_model.ring)
    element_indices = collect(eachindex(float_model.ring.line))
    samples = generate_mixed_samples(names, [1.0], trials, seed, base_kick)

    println("Building baseline thick-element linear maps...")
    linear = build_linear_lattice_data(
        closed_orbit, detectors;
        model_factory, config,
    )
    length(element_indices) == length(linear.local_maps) ||
        error("Element registry changed while building baseline maps")
    println("Building periodic $(length(closed_orbit))-dimensional response for $(length(element_indices)) complete elements...")
    element_response = build_element_source_response(
        linear, element_indices, output_coordinate,
    )
    element_families = [
        source_family(linear.element_types[index], index, inventory)
        for index in eachindex(linear.element_types)
    ]
    family_names = sort!(unique(element_families))
    family_element_indices = Dict(
        family => findall(==(family), element_families)
        for family in family_names
    )
    println("Source families: $(join(family_names, ", "))")

    # Hessian blocks are only an internal route to the summed second-order
    # target.  They are not exported as paper-facing block shares.
    blocks = (:hh, :hv, :vv)
    n_sources = length(element_indices)
    total_denominator = 0.0
    total_all_residual_sq = 0.0
    total_numerator = zeros(n_sources)
    total_magnitude_sq = zeros(n_sources)
    family_total_numerator = Dict(family => 0.0 for family in family_names)
    family_total_magnitude_sq = Dict(family => 0.0 for family in family_names)
    first_closure_max = 0.0
    second_closure_max = 0.0
    periodic_source_closure_max = 0.0
    element_response_check_max = 0.0
    independent_target_check_max = 0.0
    family_partition_check_max = 0.0
    family_partition_relative_check_max = 0.0
    direction_rows = NamedTuple[]
    element_direction_rows = NamedTuple[]
    family_direction_rows = NamedTuple[]

    solve_seconds = @elapsed begin
        for trial in 1:trials
            @printf("Thick-source direction %d/%d\n", trial, trials)
            result = direction_thick_sources(
                names, detectors, closed_orbit,
                view(samples.horizontal_directions, trial, :),
                view(samples.vertical_directions, trial, :),
                base_kick, linear, output_coordinate;
                model_factory, config,
            )
            if trial <= min(trials, 3)
                independent = direction_gtpsa_q(
                    names, detectors, closed_orbit,
                    view(samples.horizontal_directions, trial, :),
                    view(samples.vertical_directions, trial, :),
                    base_kick, length(closed_orbit);
                    model_factory, config,
                )
                independent_total = sum(
                    independent.q[Symbol("$(output_plane)_$(block)")]
                    for block in blocks
                )
                independent_target_check_max = max(
                    independent_target_check_max,
                    maximum(abs, sum(result.targets[block] for block in blocks) - independent_total),
                )
            end
            first_closure_max = max(first_closure_max, result.first_closure_residual)
            second_closure_max = max(second_closure_max, result.second_closure_residual)
            all_reconstruction = Dict{Symbol,Vector{Float64}}()
            element_reconstruction = Dict{Symbol,Vector{Float64}}()
            contributions = Dict{Symbol,Matrix{Float64}}()
            family_reconstruction = Dict{Tuple{String,Symbol},Vector{Float64}}()

            for block in blocks
                all_vector, periodic_closure = reconstruct_sources(
                    linear, result.sources[block], output_coordinate,
                )
                periodic_source_closure_max = max(periodic_source_closure_max, periodic_closure)
                contribution = zeros(length(detectors), n_sources)
                for source in 1:n_sources
                    element_index = element_indices[source]
                    local_source = view(result.sources[block], :, element_index)
                    contribution[:, source] .= element_response[source] * local_source
                end
                element_vector = vec(sum(contribution; dims=2))
                element_response_check_max = max(
                    element_response_check_max, maximum(abs, element_vector - all_vector),
                )
                all_reconstruction[block] = all_vector
                element_reconstruction[block] = element_vector
                contributions[block] = contribution
                for family in family_names
                    source_columns = family_element_indices[family]
                    family_vector = isempty(source_columns) ?
                        zeros(length(detectors)) :
                        vec(sum(view(contribution, :, source_columns); dims=2))
                    family_reconstruction[(family, block)] = family_vector
                end
            end

            total_target = sum(result.targets[block] for block in blocks)
            total_all = sum(all_reconstruction[block] for block in blocks)
            total_contribution = sum(contributions[block] for block in blocks)
            total_denominator += sum(abs2, total_target)
            total_all_residual_sq += sum(abs2, total_target - total_all)
            family_total_vectors = Dict{String,Vector{Float64}}()
            for family in family_names
                vector = sum(family_reconstruction[(family, block)] for block in blocks)
                family_total_vectors[family] = vector
                family_total_numerator[family] += dot(vector, total_target)
                family_total_magnitude_sq[family] += sum(abs2, vector)
                push!(family_direction_rows, (; trial, family,
                    element_count=length(family_element_indices[family]),
                    contribution_norm_m=norm(vector),
                    projection_numerator=dot(vector, total_target)))
            end
            family_total = sum(family_total_vectors[family] for family in family_names)
            family_partition_check_max = max(
                family_partition_check_max, maximum(abs, family_total - total_all),
            )
            family_partition_relative_check_max = max(
                family_partition_relative_check_max,
                norm(family_total - total_all) / max(norm(total_target), eps(Float64)),
            )
            for source in 1:n_sources
                vector = view(total_contribution, :, source)
                total_numerator[source] += dot(vector, total_target)
                total_magnitude_sq[source] += sum(abs2, vector)
                element_index = element_indices[source]
                push!(element_direction_rows, (; trial, element_order=source,
                    element_index,
                    element_name=linear.element_names[element_index],
                    element_type=linear.element_types[element_index],
                    s_m=linear.element_s_m[element_index],
                    k2l_m2=haskey(inventory, element_index) ?
                        Float64(inventory[element_index].k2l) : 0.0,
                    projection_numerator=dot(vector, total_target),
                    contribution_norm_m=norm(vector)))
            end
            push!(direction_rows, (; trial,
                q_total_norm_m=norm(total_target),
                all_element_total_relative_closure=relative_residual(total_target, total_all),
                all_element_total_signed_projection=dot(total_all, total_target) /
                    sum(abs2, total_target),
                family_partition_absolute_error_m=norm(family_total - total_all, Inf),
                family_partition_relative_closure=norm(family_total - total_all) /
                    max(norm(total_target), eps(Float64))))
        end
    end

    element_rows = NamedTuple[]
    for source in 1:n_sources
        element_index = element_indices[source]
        push!(element_rows, (; element_order=source, element_index,
            element_name=linear.element_names[element_index],
            element_type=linear.element_types[element_index],
            s_m=linear.element_s_m[element_index],
            k2l_m2=haskey(inventory, element_index) ?
                Float64(inventory[element_index].k2l) : 0.0,
            eta_total=total_numerator[source] / total_denominator,
            magnitude_total=sqrt(total_magnitude_sq[source] / total_denominator)))
    end

    family_rows = NamedTuple[]
    for family in family_names
        push!(family_rows, (; family, element_count=length(family_element_indices[family]),
            eta_total=family_total_numerator[family] / total_denominator,
            magnitude_total=sqrt(family_total_magnitude_sq[family] / total_denominator)))
    end

    summary = (; trials, output_plane, output_coordinate,
        elements=length(linear.local_maps),
        active_normal_sextupoles=length(inventory), detectors=length(detectors),
        first_fixed_point_closure_max=first_closure_max,
        second_fixed_point_closure_max=second_closure_max,
        periodic_source_closure_max,
        element_response_check_max,
        independent_target_check_max,
        family_partition_check_max,
        family_partition_relative_check_max,
        total_all_element_relative_closure=sqrt(total_all_residual_sq / total_denominator),
        total_all_element_signed_projection=sum(total_numerator) / total_denominator,
        total_target_rms_m=sqrt(total_denominator / trials),
        solve_seconds)

    length(element_rows) == length(linear.local_maps) ||
        error("Element summary does not cover every complete lattice element")
    length(direction_rows) == trials || error("Direction summary count mismatch")
    length(element_direction_rows) == n_sources * trials ||
        error("Element direction summary does not cover every element/direction pair")
    length(family_direction_rows) == length(family_names) * trials ||
        error("Family direction summary count mismatch")
    assert_finite_rows("element contribution", element_rows)
    assert_finite_rows("family contribution", family_rows)
    assert_finite_rows("direction closure", direction_rows)
    assert_finite_rows("element direction contribution", element_direction_rows)
    assert_finite_rows("family direction contribution", family_direction_rows)
    assert_finite_rows("reconstruction summary", [summary])

    mkpath(output_dir)
    element_path = write_namedtuple_csv(joinpath(output_dir, "element_contribution_summary.csv"), element_rows)
    family_path = write_namedtuple_csv(joinpath(output_dir, "family_contribution_summary.csv"), family_rows)
    direction_path = write_namedtuple_csv(joinpath(output_dir, "direction_closure.csv"), direction_rows)
    element_direction_path = write_namedtuple_csv(joinpath(output_dir, "element_direction_contributions.csv"), element_direction_rows)
    family_direction_path = write_namedtuple_csv(joinpath(output_dir, "family_direction_contributions.csv"), family_direction_rows)
    summary_path = write_namedtuple_csv(joinpath(output_dir, "reconstruction_summary.csv"), [summary])
    metadata_path = joinpath(output_dir, "metadata.toml")
    open(metadata_path, "w") do io
        metadata = merge(ring_metadata(config; ring), Dict(
            "format" => "cesr-thick-element-hessian-sourcing-v3",
            "date" => string(Dates.today()), "trials" => trials, "seed" => seed,
            "base_kick_rad" => base_kick,
            "output_plane" => output_plane,
            "output_coordinate" => output_coordinate,
            "target" => "summed second-order nonlinear detector vector Q_$(output_plane) at $(length(detectors)) detectors",
            "target_definition" => "Q = Q_hh + Q_hv + Q_vv; Hessian blocks are internal computation only",
            "nonlinear_order" => "second_order",
            "local_source" => "exact complete-element Hessian source g_j = S_exit - A_j*S_entrance",
            "source_coordinates" => "$(length(closed_orbit))-dimensional canonical state at element exit",
            "source_scope" => "every complete element in the selected ring line",
            "source_families" => family_names,
            "element_count" => length(element_indices),
            "element_types" => sort!(unique(linear.element_types)),
            "detector_count" => length(detectors),
            "detector_names" => detectors,
            "input_csv" => input_path,
            "control_count" => length(names),
            "control_names" => names,
            "state_dimension" => length(closed_orbit),
            "projection" => "ensemble dot(C_j,Q)/ensemble norm(Q)^2",
            "closure" => "vector norm(target - sum(all element source vectors))/vector norm(target)",
        ))
        TOML.print(io, metadata; sorted=true)
    end
    @printf("Solve time: %.3f s\n", solve_seconds)
    println("Elements:           $element_path")
    println("Families:           $family_path")
    println("Directions:         $direction_path")
    println("Element directions: $element_direction_path")
    println("Family directions:  $family_direction_path")
    println("Summary:            $summary_path")
    println("Metadata:           $metadata_path")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main_thick_element_sourcing())
end

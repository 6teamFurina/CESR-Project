#!/usr/bin/env julia

"""Attribute detector Q_x to exact local Hessian sources of thick elements."""

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
        "trials" => "100",
        "seed" => "20260804",
        "base-kick-rad" => "5e-6",
        "output-dir" => joinpath(THICK_SOURCE_HERE, "results"),
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

"""First- and second-order derivatives after substituting the closed orbit."""
function closed_map_derivatives(map, closed_derivatives)
    jacobian, hessians = map_derivatives(map)
    first = jacobian * closed_derivatives.lifted
    second = zeros(6, 2, 2)
    for coordinate in 1:6, left in 1:2, right in left:2
        value = dot(
            view(jacobian, coordinate, 1:6),
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
    bunch = Bunch(v=reshape(copy.(input_map), 1, 6))
    SciBmad.BTBL.check_bl_bunch!(bunch, ring, false)
    for element in ring.line
        track!(bunch, element)
    end
    return copy.(vec(bunch.coords.v))
end

"""Baseline local and cumulative first-order maps at every element exit."""
function build_linear_lattice_data(closed_orbit, detectors)
    descriptor = Descriptor(6, 1)
    variables = vars(descriptor)
    model = load_cesr_model(zero_value=zero(variables[1]), rf_on=true)
    input_map = [closed_orbit[index] + copy(variables[index]) for index in 1:6]
    bunch = Bunch(v=reshape(input_map, 1, 6))
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

    previous = Matrix{Float64}(I, 6, 6)
    s_m = 0.0
    for (element_index, element) in enumerate(model.ring.line)
        track!(bunch, element)
        current = Matrix(GTPSA.jacobian(copy.(vec(bunch.coords.v)); include_params=true))
        size(current) == (6, 6) || error("Unexpected cumulative map size $(size(current))")
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

"""Periodic x-detector response to a six-dimensional source at element exit."""
function build_element_source_response(linear, element_indices)
    fixed_point = I - linear.one_turn_map
    responses = Matrix{Float64}[]
    for element_index in element_indices
        exit_map = linear.cumulative_maps[element_index]
        downstream_to_end = linear.one_turn_map / exit_map
        closed_start = fixed_point \ downstream_to_end
        response = zeros(length(linear.detector_maps), 6)
        for detector in eachindex(linear.detector_maps)
            detector_map = linear.detector_maps[detector]
            total = detector_map * closed_start
            if linear.detector_element_indices[detector] >= element_index
                total += detector_map / exit_map
            end
            response[detector, :] .= view(total, 1, :)
        end
        push!(responses, response)
    end
    return responses
end

"""Reconstruct detector x from arbitrary local six-dimensional sources."""
function reconstruct_sources(linear, sources, mask=nothing)
    size(sources) == (6, length(linear.local_maps)) || error("Source matrix size mismatch")
    isnothing(mask) || length(mask) == length(linear.local_maps) || error("Source mask size mismatch")
    end_source = zeros(6)
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
        iszero(detector) || (detector_values[detector] = state[1])
    end
    closure = norm(state - initial, Inf)
    return detector_values, closure
end

function source_family(element_type, element_index, inventory)
    haskey(inventory, element_index) && return "normal_sextupole"
    element_type == "Sextupole" && return "other_sextupole"
    return lowercase(element_type)
end

"""Exact local Hessian sources and detector targets for one direction pair."""
function direction_thick_sources(
    names, detectors, closed_orbit, h_direction, v_direction, base_kick, linear,
)
    setup = direction_parameterized_model(names, h_direction, v_direction, base_kick)
    input_map = [closed_orbit[index] + copy(setup.variables[index]) for index in 1:6]
    one_turn_map = one_turn_direction_map(setup.model.ring, input_map)
    closed_derivatives = implicit_closed_orbit_derivatives(one_turn_map)

    bunch = Bunch(v=reshape(copy.(input_map), 1, 6))
    SciBmad.BTBL.check_bl_bunch!(bunch, setup.model.ring, false)
    entrance = closed_map_derivatives(input_map, closed_derivatives)
    n_elements = length(setup.model.ring.line)
    sources = Dict(
        :hh => zeros(6, n_elements),
        :hv => zeros(6, n_elements),
        :vv => zeros(6, n_elements),
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
            targets[:hh][detector] = 0.5 * exit_derivatives.second[1, 1, 1]
            targets[:hv][detector] = exit_derivatives.second[1, 1, 2]
            targets[:vv][detector] = 0.5 * exit_derivatives.second[1, 2, 2]
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

function main_thick_element_sourcing(args=ARGS)
    options = parse_thick_source_args(args)
    trials = parse(Int, options["trials"])
    trials >= 1 || error("--trials must be positive")
    seed = parse(Int, options["seed"])
    base_kick = parse(Float64, options["base-kick-rad"])
    base_kick > 0 || error("--base-kick-rad must be positive")
    output_dir = abspath(options["output-dir"])

    reference = read_samples(joinpath(
        MIXED_CALCULATION_DIR, "inputs", "cesr_corrector_samples_1000.csv",
    ))
    names = reference.names
    float_model = load_cesr_model(zero_value=0.0, rf_on=true)
    detectors = detector_names(float_model.ring)
    closed_orbit = nominal_closed_orbit(float_model.ring)
    inventory = normal_sextupole_inventory(float_model.ring)
    sextupole_indices = sort!(collect(keys(inventory)))
    samples = generate_mixed_samples(names, [1.0], trials, seed, base_kick)

    println("Building baseline thick-element linear maps...")
    linear = build_linear_lattice_data(closed_orbit, detectors)
    println("Building periodic six-dimensional response for $(length(sextupole_indices)) sextupoles...")
    sextupole_response = build_element_source_response(linear, sextupole_indices)
    element_families = [
        source_family(linear.element_types[index], index, inventory)
        for index in eachindex(linear.element_types)
    ]
    family_names = sort!(unique(element_families))
    family_masks = Dict(
        family => BitVector(item == family for item in element_families)
        for family in family_names
    )
    println("Source families: $(join(family_names, ", "))")

    blocks = (:hh, :hv, :vv)
    n_sources = length(sextupole_indices)
    denominator = Dict(block => 0.0 for block in blocks)
    all_residual_sq = Dict(block => 0.0 for block in blocks)
    sext_residual_sq = Dict(block => 0.0 for block in blocks)
    numerator = Dict(block => zeros(n_sources) for block in blocks)
    magnitude_sq = Dict(block => zeros(n_sources) for block in blocks)
    total_denominator = 0.0
    total_all_residual_sq = 0.0
    total_sext_residual_sq = 0.0
    total_numerator = zeros(n_sources)
    total_magnitude_sq = zeros(n_sources)
    family_numerator = Dict((family, block) => 0.0 for family in family_names for block in blocks)
    family_magnitude_sq = Dict((family, block) => 0.0 for family in family_names for block in blocks)
    family_total_numerator = Dict(family => 0.0 for family in family_names)
    family_total_magnitude_sq = Dict(family => 0.0 for family in family_names)
    first_closure_max = 0.0
    second_closure_max = 0.0
    periodic_source_closure_max = 0.0
    sext_response_check_max = 0.0
    independent_target_check_max = 0.0
    family_partition_check_max = 0.0
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
                base_kick, linear,
            )
            if trial <= min(trials, 3)
                independent = direction_gtpsa_q(
                    names, detectors, closed_orbit,
                    view(samples.horizontal_directions, trial, :),
                    view(samples.vertical_directions, trial, :),
                    base_kick,
                )
                for block in blocks
                    independent_target_check_max = max(
                        independent_target_check_max,
                        maximum(abs, result.targets[block] - independent.q[Symbol("x_$(block)")]),
                    )
                end
            end
            first_closure_max = max(first_closure_max, result.first_closure_residual)
            second_closure_max = max(second_closure_max, result.second_closure_residual)
            all_reconstruction = Dict{Symbol,Vector{Float64}}()
            sext_reconstruction = Dict{Symbol,Vector{Float64}}()
            contributions = Dict{Symbol,Matrix{Float64}}()
            family_reconstruction = Dict{Tuple{String,Symbol},Vector{Float64}}()

            for block in blocks
                target = result.targets[block]
                all_vector, periodic_closure = reconstruct_sources(
                    linear, result.sources[block],
                )
                periodic_source_closure_max = max(periodic_source_closure_max, periodic_closure)
                contribution = zeros(length(detectors), n_sources)
                sext_sources = zeros(size(result.sources[block]))
                for source in 1:n_sources
                    element_index = sextupole_indices[source]
                    local_source = view(result.sources[block], :, element_index)
                    contribution[:, source] .= sextupole_response[source] * local_source
                    sext_sources[:, element_index] .= local_source
                    numerator[block][source] += dot(view(contribution, :, source), target)
                    magnitude_sq[block][source] += sum(abs2, view(contribution, :, source))
                end
                sext_vector = vec(sum(contribution; dims=2))
                sext_check, _ = reconstruct_sources(linear, sext_sources)
                sext_response_check_max = max(
                    sext_response_check_max, maximum(abs, sext_vector - sext_check),
                )
                all_reconstruction[block] = all_vector
                sext_reconstruction[block] = sext_vector
                contributions[block] = contribution
                denominator[block] += sum(abs2, target)
                all_residual_sq[block] += sum(abs2, target - all_vector)
                sext_residual_sq[block] += sum(abs2, target - sext_vector)
                for family in family_names
                    family_vector, _ = reconstruct_sources(
                        linear, result.sources[block], family_masks[family],
                    )
                    family_reconstruction[(family, block)] = family_vector
                    family_numerator[(family, block)] += dot(family_vector, target)
                    family_magnitude_sq[(family, block)] += sum(abs2, family_vector)
                end
            end

            total_target = sum(result.targets[block] for block in blocks)
            total_all = sum(all_reconstruction[block] for block in blocks)
            total_sext = sum(sext_reconstruction[block] for block in blocks)
            total_contribution = sum(contributions[block] for block in blocks)
            total_denominator += sum(abs2, total_target)
            total_all_residual_sq += sum(abs2, total_target - total_all)
            total_sext_residual_sq += sum(abs2, total_target - total_sext)
            family_total_vectors = Dict{String,Vector{Float64}}()
            for family in family_names
                vector = sum(family_reconstruction[(family, block)] for block in blocks)
                family_total_vectors[family] = vector
                family_total_numerator[family] += dot(vector, total_target)
                family_total_magnitude_sq[family] += sum(abs2, vector)
                push!(family_direction_rows, (; trial, family,
                    element_count=count(family_masks[family]),
                    contribution_norm_m=norm(vector),
                    projection_numerator=dot(vector, total_target)))
            end
            family_total = sum(family_total_vectors[family] for family in family_names)
            family_partition_check_max = max(
                family_partition_check_max, maximum(abs, family_total - total_all),
            )
            for source in 1:n_sources
                vector = view(total_contribution, :, source)
                total_numerator[source] += dot(vector, total_target)
                total_magnitude_sq[source] += sum(abs2, vector)
                element_index = sextupole_indices[source]
                push!(element_direction_rows, (; trial, element_order=source,
                    element_index,
                    element_name=linear.element_names[element_index],
                    s_m=linear.element_s_m[element_index],
                    k2l_m2=Float64(inventory[element_index].k2l),
                    projection_numerator=dot(vector, total_target),
                    contribution_norm_m=norm(vector)))
            end
            push!(direction_rows, (; trial,
                q_total_norm_m=norm(total_target),
                all_element_total_relative_closure=relative_residual(total_target, total_all),
                sextupole_total_relative_closure=relative_residual(total_target, total_sext),
                sextupole_total_signed_projection=dot(total_sext, total_target) / sum(abs2, total_target),
                hh_all_relative_closure=relative_residual(result.targets[:hh], all_reconstruction[:hh]),
                hv_all_relative_closure=relative_residual(result.targets[:hv], all_reconstruction[:hv]),
                vv_all_relative_closure=relative_residual(result.targets[:vv], all_reconstruction[:vv]),
                hh_sextupole_relative_closure=relative_residual(result.targets[:hh], sext_reconstruction[:hh]),
                hv_sextupole_relative_closure=relative_residual(result.targets[:hv], sext_reconstruction[:hv]),
                vv_sextupole_relative_closure=relative_residual(result.targets[:vv], sext_reconstruction[:vv])))
        end
    end

    element_rows = NamedTuple[]
    for source in 1:n_sources
        element_index = sextupole_indices[source]
        push!(element_rows, (; element_order=source, element_index,
            element_name=linear.element_names[element_index],
            s_m=linear.element_s_m[element_index],
            k2l_m2=Float64(inventory[element_index].k2l),
            eta_hh=numerator[:hh][source] / denominator[:hh],
            eta_hv=numerator[:hv][source] / denominator[:hv],
            eta_vv=numerator[:vv][source] / denominator[:vv],
            eta_total=total_numerator[source] / total_denominator,
            magnitude_hh=sqrt(magnitude_sq[:hh][source] / denominator[:hh]),
            magnitude_hv=sqrt(magnitude_sq[:hv][source] / denominator[:hv]),
            magnitude_vv=sqrt(magnitude_sq[:vv][source] / denominator[:vv]),
            magnitude_total=sqrt(total_magnitude_sq[source] / total_denominator)))
    end

    family_rows = NamedTuple[]
    for family in family_names
        push!(family_rows, (; family, element_count=count(family_masks[family]),
            eta_hh=family_numerator[(family, :hh)] / denominator[:hh],
            eta_hv=family_numerator[(family, :hv)] / denominator[:hv],
            eta_vv=family_numerator[(family, :vv)] / denominator[:vv],
            eta_total=family_total_numerator[family] / total_denominator,
            magnitude_hh=sqrt(family_magnitude_sq[(family, :hh)] / denominator[:hh]),
            magnitude_hv=sqrt(family_magnitude_sq[(family, :hv)] / denominator[:hv]),
            magnitude_vv=sqrt(family_magnitude_sq[(family, :vv)] / denominator[:vv]),
            magnitude_total=sqrt(family_total_magnitude_sq[family] / total_denominator)))
    end

    summary = (; trials, elements=length(linear.local_maps),
        active_normal_sextupoles=n_sources, detectors=length(detectors),
        first_fixed_point_closure_max=first_closure_max,
        second_fixed_point_closure_max=second_closure_max,
        periodic_source_closure_max,
        sextupole_response_check_max=sext_response_check_max,
        independent_target_check_max,
        family_partition_check_max,
        hh_all_element_relative_closure=sqrt(all_residual_sq[:hh] / denominator[:hh]),
        hv_all_element_relative_closure=sqrt(all_residual_sq[:hv] / denominator[:hv]),
        vv_all_element_relative_closure=sqrt(all_residual_sq[:vv] / denominator[:vv]),
        total_all_element_relative_closure=sqrt(total_all_residual_sq / total_denominator),
        hh_sextupole_relative_closure=sqrt(sext_residual_sq[:hh] / denominator[:hh]),
        hv_sextupole_relative_closure=sqrt(sext_residual_sq[:hv] / denominator[:hv]),
        vv_sextupole_relative_closure=sqrt(sext_residual_sq[:vv] / denominator[:vv]),
        total_sextupole_relative_closure=sqrt(total_sext_residual_sq / total_denominator),
        hh_sextupole_signed_projection=sum(numerator[:hh]) / denominator[:hh],
        hv_sextupole_signed_projection=sum(numerator[:hv]) / denominator[:hv],
        vv_sextupole_signed_projection=sum(numerator[:vv]) / denominator[:vv],
        total_sextupole_signed_projection=sum(total_numerator) / total_denominator,
        solve_seconds)

    mkpath(output_dir)
    element_path = write_namedtuple_csv(joinpath(output_dir, "thick_sextupole_contribution_summary.csv"), element_rows)
    family_path = write_namedtuple_csv(joinpath(output_dir, "family_contribution_summary.csv"), family_rows)
    direction_path = write_namedtuple_csv(joinpath(output_dir, "direction_closure.csv"), direction_rows)
    element_direction_path = write_namedtuple_csv(joinpath(output_dir, "thick_sextupole_direction_contributions.csv"), element_direction_rows)
    family_direction_path = write_namedtuple_csv(joinpath(output_dir, "family_direction_contributions.csv"), family_direction_rows)
    summary_path = write_namedtuple_csv(joinpath(output_dir, "reconstruction_summary.csv"), [summary])
    metadata_path = joinpath(output_dir, "metadata.toml")
    open(metadata_path, "w") do io
        TOML.print(io, Dict(
            "format" => "cesr-thick-element-hessian-sourcing-v1",
            "date" => string(Dates.today()), "trials" => trials, "seed" => seed,
            "base_kick_rad" => base_kick,
            "target" => "Q_x = Q_hh,x + Q_hv,x + Q_vv,x at 99 detectors",
            "local_source" => "exact complete-element Hessian source g_j = S_exit - A_j*S_entrance",
            "source_coordinates" => "six-dimensional canonical state at element exit",
            "projection" => "ensemble dot(C_j,Q)/ensemble norm(Q)^2",
        ); sorted=true)
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

#!/usr/bin/env julia

"""
Compute direction-contracted second-order CESR closed-orbit responses with
GTPSA corrector parameters.

For every fixed horizontal/vertical direction pair, two dimensionless GTPSA
parameters `a` and `b` drive the physical corrector controls as

    delta_k = base_kick * (a * h_direction + b * v_direction).

A second-order one-turn map in the six initial coordinates and the two
corrector parameters is differentiated implicitly through the closed-orbit
fixed-point equation.  The resulting detector derivatives give

    Q_hh(rho) = rho^2 / 2 * d2u/da2,
    Q_hv(rho) = rho^2     * d2u/dadb,
    Q_vv(rho) = rho^2 / 2 * d2u/db2.

Thus the calculation has no corrector finite-difference step and contains no
fourth- or higher-order contamination.
"""

const GTPSA_MIXED_HERE = @__DIR__
include(joinpath(GTPSA_MIXED_HERE, "run_mixed_term_experiment.jl"))

using GTPSA
using LinearAlgebra
using Printf

function parse_gtpsa_mixed_args(args)
    options = Dict{String,String}(
        "ring" => "latest",
        "rhos" => "0.1,0.14,0.2,0.28,0.4,0.57,0.8,1.13",
        "trials" => "100",
        "seed" => "20260804",
        "base-kick-rad" => "5e-6",
        "output-dir" => "",
        "inputs" => "",
        "finite-difference-directions" => "",
    )
    for argument in args
        startswith(argument, "--") ||
            error("Arguments must have --name=value form: $argument")
        fields = split(argument[3:end], '='; limit=2)
        length(fields) == 2 || error("Missing value in argument: $argument")
        haskey(options, fields[1]) || error("Unknown option: --$(fields[1])")
        options[fields[1]] = fields[2]
    end
    ring = Symbol(lowercase(options["ring"]))
    ring in (:latest, :latest_cesr, :repaired_latest, :legacy, :legacy_cesr, :historical) ||
        error("--ring must be latest or legacy")
    isempty(options["inputs"]) &&
        (options["inputs"] = default_ring_paths(; ring).inputs)
    artifact = ring_artifact_id(ring)
    isempty(options["output-dir"]) &&
        (options["output-dir"] = joinpath(GTPSA_MIXED_HERE, "gtpsa_results", artifact))
    isempty(options["finite-difference-directions"]) &&
        (options["finite-difference-directions"] = joinpath(
            GTPSA_MIXED_HERE, "results", artifact, "mixed_term_directions.csv",
        ))
    return options
end

function nominal_closed_orbit(ring)
    solution = find_closed_orbit(
        ring;
        coasting_beam=false,
        batch=Val{false}(),
        warn=true,
    )
    bad = findall(solution.sol.retcode .!= SciBmad.BatchSolve.RETCODE_SUCCESS)
    isempty(bad) || error("Nominal closed-orbit solve failed: $(solution.sol.retcode)")
    return vec(Float64.(solution.v0[1, :]))
end

function direction_parameterized_model(
    names,
    h_direction,
    v_direction,
    base_kick,
    state_dimension;
    model_factory=load_ring_model,
    config=nothing,
)
    descriptor = Descriptor(state_dimension, 2, 2, 2)
    variables = vars(descriptor)
    parameters = params(descriptor)
    # Keep unselected skew/group controls primitive.  Only the requested H/V
    # steering controls below carry the two GTPSA parameters; promoting the
    # complete registry can trigger the known combined-multipole sqrt(0)
    # domain failure around SEX_14W.
    model = configured_model(model_factory, config; zero_value=0.0, rf_on=true)
    for index in eachindex(names)
        value = base_kick * (
            h_direction[index] * parameters[1] +
            v_direction[index] * parameters[2]
        )
        model.controls[names[index]] = value
    end
    return (; descriptor, variables, parameters, model)
end

function track_direction_maps(ring, input_map, detectors)
    state_dimension = length(input_map)
    bunch = Bunch(v=reshape(input_map, 1, state_dimension))
    SciBmad.BTBL.check_bl_bunch!(bunch, ring, false)
    detector_index = Dict(name => index for (index, name) in enumerate(detectors))
    maps = Vector{typeof(copy.(input_map))}(undef, length(detectors))
    found = falses(length(detectors))
    for element in ring.line
        track!(bunch, element)
        name = uppercase(String(element.name))
        haskey(detector_index, name) || continue
        index = detector_index[name]
        found[index] && error("Detector $name occurs more than once")
        maps[index] = copy.(vec(bunch.coords.v))
        found[index] = true
    end
    all(found) || error("Missing detector maps: $(detectors[.!found])")
    return maps, copy.(vec(bunch.coords.v))
end

function map_derivatives(map, state_dimension)
    jacobian = Matrix(GTPSA.jacobian(map; include_params=true))
    expected_columns = state_dimension + 2
    size(jacobian) == (state_dimension, expected_columns) ||
        error("Expected a $state_dimension x $expected_columns direction-map Jacobian, got $(size(jacobian))")
    hessians = [Matrix(GTPSA.hessian(map[index]; include_params=true)) for index in 1:state_dimension]
    all(size(hessian) == (expected_columns, expected_columns) for hessian in hessians) ||
        error("Expected $expected_columns x $expected_columns direction-map Hessians")
    return jacobian, hessians
end

"""Implicit first and second derivatives of the closed-orbit fixed point."""
function implicit_closed_orbit_derivatives(one_turn_map, state_dimension)
    jacobian, hessians = map_derivatives(one_turn_map, state_dimension)
    A = jacobian[:, 1:state_dimension]
    B = jacobian[:, state_dimension + 1:state_dimension + 2]
    fixed_point_matrix = Matrix{Float64}(I, state_dimension, state_dimension) - A
    first = fixed_point_matrix \ B
    lifted = vcat(first, Matrix{Float64}(I, 2, 2))
    second = zeros(state_dimension, 2, 2)
    sources = zeros(state_dimension, 2, 2)
    for left in 1:2, right in left:2
        source = [
            dot(view(lifted, :, left), hessians[coordinate] * view(lifted, :, right))
        for coordinate in 1:state_dimension
        ]
        value = fixed_point_matrix \ source
        sources[:, left, right] .= source
        sources[:, right, left] .= source
        second[:, left, right] .= value
        second[:, right, left] .= value
    end
    first_residual = fixed_point_matrix * first - B
    second_residual = maximum(
        norm(fixed_point_matrix * view(second, :, left, right) -
             view(sources, :, left, right), Inf)
        for left in 1:2, right in 1:2
    )
    return (; first, second, lifted,
        first_residual=norm(first_residual, Inf), second_residual)
end

function detector_second_derivatives(detector_maps, closed_derivatives; config=nothing)
    count = length(detector_maps)
    horizontal = zeros(count, 2, 2)
    vertical = zeros(count, 2, 2)
    state_dimension = size(closed_derivatives.lifted, 1) - 2
    transverse = transverse_coordinate_indices(config; state_dimension)
    for detector in 1:count
        jacobian, hessians = map_derivatives(detector_maps[detector], state_dimension)
        for (coordinate, output) in ((transverse.x, horizontal), (transverse.y, vertical))
            initial_gradient = view(jacobian, coordinate, 1:state_dimension)
            hessian = hessians[coordinate]
            for left in 1:2, right in left:2
                value = dot(
                    initial_gradient,
                    view(closed_derivatives.second, :, left, right),
                ) + dot(
                    view(closed_derivatives.lifted, :, left),
                    hessian * view(closed_derivatives.lifted, :, right),
                )
                output[detector, left, right] = value
                output[detector, right, left] = value
            end
        end
    end
    return horizontal, vertical
end

function direction_gtpsa_q(
    names,
    detectors,
    closed_orbit,
    h_direction,
    v_direction,
    base_kick,
    state_dimension;
    model_factory=load_ring_model,
    config=nothing,
)
    setup = direction_parameterized_model(
        names, h_direction, v_direction, base_kick, state_dimension;
        model_factory, config,
    )
    input_map = [closed_orbit[index] + copy(setup.variables[index]) for index in 1:state_dimension]
    detector_maps, one_turn_map = track_direction_maps(
        setup.model.ring, input_map, detectors,
    )
    closed_derivatives = implicit_closed_orbit_derivatives(one_turn_map, state_dimension)
    horizontal, vertical = detector_second_derivatives(
        detector_maps, closed_derivatives; config,
    )
    q = Dict{Symbol,Vector{Float64}}(
        :x_hh => horizontal[:, 1, 1] / 2,
        :x_hv => horizontal[:, 1, 2],
        :x_vv => horizontal[:, 2, 2] / 2,
        :y_hh => vertical[:, 1, 1] / 2,
        :y_hv => vertical[:, 1, 2],
        :y_vv => vertical[:, 2, 2] / 2,
    )
    return (; q,
        first_closure_residual=closed_derivatives.first_residual,
        second_closure_residual=closed_derivatives.second_residual)
end

function q_direction_row(trial, rho, base_kick, result)
    scale = rho^2
    metrics = Dict{Symbol,Float64}()
    for plane in (:x, :y)
        squared_norms = Dict{Symbol,Float64}()
        for block in (:hh, :hv, :vv)
            vector = scale .* result.q[Symbol("$(plane)_$(block)")]
            rmse = vector_rmse(vector)
            metrics[Symbol("$(plane)_q$(block)_rmse_m")] = rmse
            squared_norms[block] = sum(abs2, vector)
        end
        denominator = sum(Base.values(squared_norms))
        for block in (:hh, :hv, :vv)
            metrics[Symbol("$(plane)_f_$(block)")] = squared_norms[block] / denominator
        end
    end
    return (;
        rho,
        trial,
        active_rms_rad=rho * base_kick,
        first_closure_residual=result.first_closure_residual,
        second_closure_residual=result.second_closure_residual,
        (key => metrics[key] for key in sort!(collect(keys(metrics)); by=String))...,
    )
end

function percentile_linear(values, probability)
    sorted = sort(collect(values))
    isempty(sorted) && return NaN
    position = 1 + probability * (length(sorted) - 1)
    lower = floor(Int, position)
    upper = ceil(Int, position)
    lower == upper && return sorted[lower]
    fraction = position - lower
    return (1 - fraction) * sorted[lower] + fraction * sorted[upper]
end

function summarize_gtpsa_rows(rows, rhos)
    metric_names = propertynames(first(rows))[6:end]
    output = NamedTuple[]
    for rho in rhos
        selected = [row for row in rows if row.rho == rho]
        statistics = Pair{Symbol,Float64}[]
        for metric in metric_names
            samples = [getproperty(row, metric) for row in selected]
            push!(statistics, Symbol("mean_$(metric)") => mean(samples))
            push!(statistics, Symbol("p10_$(metric)") => percentile_linear(samples, 0.1))
            push!(statistics, Symbol("median_$(metric)") => percentile_linear(samples, 0.5))
            push!(statistics, Symbol("p90_$(metric)") => percentile_linear(samples, 0.9))
        end
        push!(output, (; rho, direction_pairs=length(selected), statistics...))
    end
    return output
end

function read_numeric_csv(path)
    lines = readlines(path)
    isempty(lines) && error("Empty CSV: $path")
    header = Symbol.(split(lines[1], ','))
    rows = NamedTuple[]
    for line in lines[2:end]
        fields = split(line, ',')
        length(fields) == length(header) || error("Malformed CSV row in $path")
        parsed = map(fields) do field
            field == "true" && return true
            field == "false" && return false
            return parse(Float64, field)
        end
        push!(rows, NamedTuple{Tuple(header)}(Tuple(parsed)))
    end
    return rows
end

function compare_finite_difference(gtpsa_rows, finite_difference_path)
    isfile(finite_difference_path) || return NamedTuple[]
    finite_rows = read_numeric_csv(finite_difference_path)
    lookup = Dict(
        (Int(round(row.trial)), row.rho) => row
        for row in finite_rows if row.all_states_converged
    )
    metrics = [Symbol("$(plane)_q$(block)_rmse_m")
               for plane in (:x, :y) for block in (:hh, :hv, :vv)]
    output = NamedTuple[]
    for row in gtpsa_rows
        key = (row.trial, row.rho)
        haskey(lookup, key) || continue
        finite = lookup[key]
        pairs = Pair{Symbol,Float64}[]
        for metric in metrics
            exact = getproperty(row, metric)
            candidate = getproperty(finite, metric)
            relative = iszero(exact) ? abs(candidate - exact) : abs(candidate - exact) / exact
            push!(pairs, Symbol("fd_$(metric)") => candidate)
            push!(pairs, Symbol("gtpsa_$(metric)") => exact)
            push!(pairs, Symbol("relative_difference_$(metric)") => relative)
        end
        push!(output, (; rho=row.rho, trial=row.trial, pairs...))
    end
    return output
end

function main_gtpsa_mixed(args=ARGS; model_factory=nothing, config=nothing)
    options = parse_gtpsa_mixed_args(args)
    ring = Symbol(lowercase(options["ring"]))
    model_factory, config = resolve_ring_model_factory(model_factory, config; ring)
    rhos = parse.(Float64, split(options["rhos"], ','))
    !isempty(rhos) && all(isfinite, rhos) && all(>(0), rhos) ||
        error("All radii must be finite and positive")
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
    layout = observable_layout(detectors; config)
    closed_orbit = nominal_closed_orbit(float_model.ring)
    state_dimension = length(closed_orbit)
    samples = generate_mixed_samples(names, [1.0], trials, seed, base_kick; config, model=float_model)

    rows = NamedTuple[]
    solve_seconds = @elapsed begin
        for trial in 1:trials
            @printf("GTPSA direction %d/%d\n", trial, trials)
            result = direction_gtpsa_q(
                names,
                detectors,
                closed_orbit,
                view(samples.horizontal_directions, trial, :),
                view(samples.vertical_directions, trial, :),
                base_kick,
                state_dimension;
                model_factory,
                config,
            )
            for rho in rhos
                push!(rows, q_direction_row(trial, rho, base_kick, result))
            end
        end
    end

    summary = summarize_gtpsa_rows(rows, rhos)
    comparison = compare_finite_difference(
        rows, abspath(options["finite-difference-directions"]),
    )
    mkpath(output_dir)
    direction_path = write_namedtuple_csv(joinpath(output_dir, "gtpsa_direction_q.csv"), rows)
    summary_path = write_namedtuple_csv(joinpath(output_dir, "gtpsa_summary.csv"), summary)
    comparison_path = isempty(comparison) ? "" : write_namedtuple_csv(
        joinpath(output_dir, "gtpsa_finite_difference_comparison.csv"), comparison,
    )
    metadata = Dict(
        "format" => "cesr-direction-gtpsa-mixed-term-v1",
        "date" => string(Dates.today()),
        "definition" => "two second-order GTPSA corrector parameters plus implicit differentiation of the RF-on closed-orbit fixed point",
        "result_role" => "adopted final quadratic response; four-sign finite differences are an independent validation",
        "reporting_statistic" => "median [P10, P90] across fixed directions",
        "descriptor" => "Descriptor($(state_dimension), 2, 2, 2)",
        "rhos" => rhos,
        "trials" => trials,
        "seed" => seed,
        "base_kick_rad" => base_kick,
        "input_csv" => input_path,
        "control_count" => length(names),
        "control_names" => names,
        "horizontal_control_count" => length(samples.horizontal_indices),
        "vertical_control_count" => length(samples.vertical_indices),
        "state_dimension" => state_dimension,
        "detector_count" => length(detectors),
        "observable_count" => length(layout.labels),
        "observable_labels" => layout.labels,
        "solve_seconds" => solve_seconds,
        "direction_csv" => direction_path,
        "summary_csv" => summary_path,
        "finite_difference_comparison_csv" => comparison_path,
        "finite_difference_source_csv" => abspath(options["finite-difference-directions"]),
        "maximum_first_order_closure_residual" => maximum(row.first_closure_residual for row in rows),
        "maximum_second_order_closure_residual" => maximum(row.second_closure_residual for row in rows),
    )
    merge!(metadata, ring_metadata(config; ring))
    metadata_path = joinpath(output_dir, "gtpsa_metadata.toml")
    open(metadata_path, "w") do io
        TOML.print(io, metadata; sorted=true)
    end
    @printf("GTPSA solve time: %.3f s\n", solve_seconds)
    println("Direction Q: $direction_path")
    println("Summary:     $summary_path")
    !isempty(comparison_path) && println("FD comparison: $comparison_path")
    println("Metadata:    $metadata_path")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main_gtpsa_mixed())
end

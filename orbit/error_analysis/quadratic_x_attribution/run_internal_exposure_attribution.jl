#!/usr/bin/env julia

"""Diagnose the X-plane Q_vv/Q_hh imbalance from internal sextupole exposure."""

const ATTRIBUTION_HERE = @__DIR__
const MIXED_TERM_HERE = normpath(joinpath(ATTRIBUTION_HERE, "..", "mixed_terms"))
include(joinpath(MIXED_TERM_HERE, "run_mixed_term_gtpsa.jl"))

using Dates
using LinearAlgebra
using Printf
using Statistics
using TOML

function parse_exposure_args(args)
    options = Dict{String,String}(
        "ring" => "latest",
        "trials" => "100",
        "seed" => "20260804",
        "base-kick-rad" => "5e-6",
        "rho" => "1.13",
        "gtpsa-direction-csv" => "",
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
    isempty(options["inputs"]) &&
        (options["inputs"] = default_ring_paths(; ring).inputs)
    artifact = ring_artifact_id(ring)
    isempty(options["gtpsa-direction-csv"]) &&
        (options["gtpsa-direction-csv"] = joinpath(
            MIXED_TERM_HERE, "gtpsa_results", artifact, "gtpsa_direction_q.csv",
        ))
    isempty(options["output-dir"]) &&
        (options["output-dir"] = joinpath(ATTRIBUTION_HERE, "element_results", artifact))
    return options
end

function normal_sextupole_strength(element)
    kn2 = Float64(GTPSA.scalar(Beamlines.deval(element.Kn2)))
    kn2l = Float64(GTPSA.scalar(Beamlines.deval(element.Kn2L)))
    length_m = Float64(GTPSA.scalar(Beamlines.deval(element.L)))
    !iszero(kn2) && return kn2 * length_m
    return kn2l
end

function normal_sextupole_inventory(ring)
    inventory = Dict{Int,NamedTuple}()
    s_start_m = 0.0
    for (element_index, element) in enumerate(ring.line)
        length_m = Float64(GTPSA.scalar(Beamlines.deval(element.L)))
        k2l = normal_sextupole_strength(element)
        if !iszero(k2l)
            inventory[element_index] = (;
                name=String(element.name), k2l,
                s_m=s_start_m + length_m / 2,
            )
        end
        s_start_m += length_m
    end
    isempty(inventory) && error("No active normal sextupoles were found")
    return inventory
end

function all_corrector_first_order_model(
    names,
    base_kick,
    state_dimension;
    model_factory=load_ring_model,
    config=nothing,
)
    descriptor = Descriptor(state_dimension, 1, length(names), 1)
    variables = vars(descriptor)
    parameters = params(descriptor)
    model = configured_model(model_factory, config; zero_value=zero(parameters[1]), rf_on=true)
    for index in eachindex(names)
        model.controls[names[index]] = base_kick * parameters[index]
    end
    return (; variables, model)
end

function track_with_sextupole_maps(ring, input_map, inventory)
    state_dimension = length(input_map)
    state_dimension = length(input_map)
    bunch = Bunch(v=reshape(input_map, 1, state_dimension))
    SciBmad.BTBL.check_bl_bunch!(bunch, ring, false)
    sextupole_maps = NamedTuple[]
    for (element_index, element) in enumerate(ring.line)
        active = get(inventory, element_index, nothing)
        entrance = isnothing(active) ? nothing : copy.(vec(bunch.coords.v))
        track!(bunch, element)
        exit_map = copy.(vec(bunch.coords.v))
        if !isnothing(active)
            midpoint = [(entrance[index] + exit_map[index]) / 2 for index in 1:state_dimension]
            push!(sextupole_maps, (;
                element_index, name=active.name, s_m=active.s_m,
                k2l=active.k2l, map=midpoint,
            ))
        end
    end
    isempty(sextupole_maps) && error("No active normal sextupoles were found")
    return copy.(vec(bunch.coords.v)), sextupole_maps
end

function build_internal_response(
    names,
    closed_orbit,
    base_kick,
    inventory;
    model_factory=load_ring_model,
    config=nothing,
)
    state_dimension = length(closed_orbit)
    setup = all_corrector_first_order_model(
        names, base_kick, state_dimension;
        model_factory, config,
    )
    input_map = [closed_orbit[index] + copy(setup.variables[index]) for index in 1:state_dimension]
    one_turn_map, sextupole_maps = track_with_sextupole_maps(
        setup.model.ring, input_map, inventory,
    )
    jacobian = Matrix(GTPSA.jacobian(one_turn_map; include_params=true))
    expected_columns = state_dimension + length(names)
    size(jacobian) == (state_dimension, expected_columns) || error("Unexpected all-corrector Jacobian size")
    A, B = jacobian[:, 1:state_dimension], jacobian[:, state_dimension + 1:end]
    first = (I - A) \ B
    lifted = vcat(first, Matrix{Float64}(I, length(names), length(names)))
    first_residual = norm((I - A) * first - B, Inf)
    transverse = transverse_coordinate_indices(config; state_dimension)
    response_x = zeros(length(sextupole_maps), length(names))
    response_y = similar(response_x)
    k2l = zeros(length(sextupole_maps))
    element_indices = zeros(Int, length(sextupole_maps))
    element_names = Vector{String}(undef, length(sextupole_maps))
    s_m = zeros(length(sextupole_maps))
    for (row, item) in enumerate(sextupole_maps)
        internal_jacobian = Matrix(GTPSA.jacobian(item.map; include_params=true))
        size(internal_jacobian) == (state_dimension, expected_columns) ||
            error("Unexpected internal all-corrector Jacobian size")
        response_x[row, :] .= vec(view(internal_jacobian, transverse.x, :)' * lifted)
        response_y[row, :] .= vec(view(internal_jacobian, transverse.y, :)' * lifted)
        k2l[row] = item.k2l
        element_indices[row] = item.element_index
        element_names[row] = item.name
        s_m[row] = item.s_m
    end
    return (;
        response_x, response_y, k2l, element_indices, element_names, s_m,
        first_residual,
    )
end

function element_exposure_rows(response, h_direction, v_direction, trial)
    x_h = response.response_x * h_direction
    y_v = response.response_y * v_direction
    rows = NamedTuple[]
    for row in eachindex(response.k2l)
        k2l = response.k2l[row]
        orbit_sq_h = x_h[row]^2
        orbit_sq_v = y_v[row]^2
        source_h = abs(k2l) * orbit_sq_h
        source_v = abs(k2l) * orbit_sq_v
        push!(rows, (;
            trial,
            element_order=row,
            element_index=response.element_indices[row],
            element_name=response.element_names[row],
            s_m=response.s_m[row],
            k2l_m2=k2l,
            x_h_m=x_h[row],
            y_v_m=y_v[row],
            orbit_sq_h_m2=orbit_sq_h,
            orbit_sq_v_m2=orbit_sq_v,
            source_exposure_h=source_h,
            source_exposure_v=source_v,
            source_excess_v_minus_h=source_v - source_h,
            signed_local_kick_h_rad=-0.5 * k2l * orbit_sq_h,
            signed_local_kick_v_rad=0.5 * k2l * orbit_sq_v,
        ))
    end
    return rows
end

function direction_attribution(response, h_direction, v_direction, exact)
    x_h = response.response_x * h_direction
    x_v = response.response_x * v_direction
    y_h = response.response_y * h_direction
    y_v = response.response_y * v_direction
    k2l = response.k2l

    orbit_h = sqrt(mean(x_h .^ 4))
    orbit_v = sqrt(mean(y_v .^ 4))
    source_h = sum(abs.(k2l) .* x_h .^ 2)
    source_v = sum(abs.(k2l) .* y_v .^ 2)
    source_l2_h = norm(k2l .* x_h .^ 2)
    source_l2_v = norm(k2l .* y_v .^ 2)
    # Detector RMSEs come from the adopted second-order GTPSA result. A common
    # sqrt(N_detector) factor cancels from all ratios and component shares.
    norm_hh = exact.x_qhh_rmse_m
    norm_hv = exact.x_qhv_rmse_m
    norm_vv = exact.x_qvv_rmse_m
    denominator = norm_hh^2 + norm_hv^2 + norm_vv^2

    # Symmetric direction rescaling makes both source exposures equal to their
    # geometric mean. Since s_h*s_v=1, the mixed block stays unchanged.
    equal_hh = sqrt(source_v / source_h) * norm_hh
    equal_vv = sqrt(source_h / source_v) * norm_vv
    equal_denominator = equal_hh^2 + norm_hv^2 + equal_vv^2

    return (;
        q_x_hh_rmse_m=norm_hh,
        q_x_hv_rmse_m=norm_hv,
        q_x_vv_rmse_m=norm_vv,
        f_hh=norm_hh^2 / denominator,
        f_hv=norm_hv^2 / denominator,
        f_vv=norm_vv^2 / denominator,
        qnorm_ratio_vv_hh=norm_vv / norm_hh,
        orbit_exposure_h=orbit_h,
        orbit_exposure_v=orbit_v,
        orbit_exposure_ratio_v_h=orbit_v / orbit_h,
        source_exposure_h=source_h,
        source_exposure_v=source_v,
        source_exposure_ratio_v_h=source_v / source_h,
        source_l2_h,
        source_l2_v,
        source_l2_ratio_v_h=source_l2_v / source_l2_h,
        propagation_efficiency_h=norm_hh / source_h,
        propagation_efficiency_v=norm_vv / source_v,
        propagation_efficiency_ratio_v_h=(norm_vv / source_v) / (norm_hh / source_h),
        equal_source_f_hh=equal_hh^2 / equal_denominator,
        equal_source_f_hv=norm_hv^2 / equal_denominator,
        equal_source_f_vv=equal_vv^2 / equal_denominator,
        cross_plane_xv_rms_m=sqrt(mean(abs2, x_v)),
        cross_plane_yh_rms_m=sqrt(mean(abs2, y_h)),
        active_normal_sextupoles=length(k2l),
        first_closure_residual=response.first_residual,
    )
end

function summarize_attribution(rows)
    statistics = Pair{Symbol,Float64}[]
    for metric in filter(!=(:trial), propertynames(first(rows)))
        samples = Float64[getproperty(row, metric) for row in rows]
        push!(statistics, Symbol("p10_$(metric)") => percentile_linear(samples, 0.1))
        push!(statistics, Symbol("median_$(metric)") => percentile_linear(samples, 0.5))
        push!(statistics, Symbol("p90_$(metric)") => percentile_linear(samples, 0.9))
    end
    log_q = log.([row.qnorm_ratio_vv_hh for row in rows])
    log_orbit = log.([row.orbit_exposure_ratio_v_h for row in rows])
    log_source = log.([row.source_exposure_ratio_v_h for row in rows])
    push!(statistics, :cor_log_q_ratio_log_orbit_exposure_ratio => cor(log_q, log_orbit))
    push!(statistics, :cor_log_q_ratio_log_source_exposure_ratio => cor(log_q, log_source))
    return (; direction_pairs=length(rows), statistics...)
end

function main_internal_exposure(args=ARGS; model_factory=nothing, config=nothing)
    options = parse_exposure_args(args)
    ring = Symbol(lowercase(options["ring"]))
    model_factory, config = resolve_ring_model_factory(model_factory, config; ring)
    trials = parse(Int, options["trials"])
    trials >= 3 || error("--trials must be at least 3")
    seed = parse(Int, options["seed"])
    base_kick = parse(Float64, options["base-kick-rad"])
    base_kick > 0 || error("--base-kick-rad must be positive")
    rho = parse(Float64, options["rho"])
    rho > 0 || error("--rho must be positive")
    output_dir = abspath(options["output-dir"])
    gtpsa_direction_path = abspath(options["gtpsa-direction-csv"])
    exact_rows = read_numeric_csv(gtpsa_direction_path)
    exact_lookup = Dict(
        Int(round(row.trial)) => row for row in exact_rows if isapprox(row.rho, rho; atol=1e-12)
    )
    length(exact_lookup) >= trials || error(
        "GTPSA direction CSV has only $(length(exact_lookup)) directions at rho=$rho",
    )
    input_path = abspath(options["inputs"])
    input_reference = read_samples(input_path)
    names = input_reference.names
    validate_control_names(names, config)
    float_model = configured_model(model_factory, config; zero_value=0.0, rf_on=true)
    closed_orbit = nominal_closed_orbit(float_model.ring)
    inventory = normal_sextupole_inventory(float_model.ring)
    samples = generate_mixed_samples(names, [1.0], trials, seed, base_kick)
    rows = NamedTuple[]
    element_rows = NamedTuple[]
    solve_seconds = @elapsed begin
        response = build_internal_response(
            names, closed_orbit, base_kick, inventory;
            model_factory, config,
        )
        for trial in 1:trials
            @printf("Internal-exposure direction %d/%d\n", trial, trials)
            result = direction_attribution(
                response,
                view(samples.horizontal_directions, trial, :),
                view(samples.vertical_directions, trial, :),
                exact_lookup[trial],
            )
            push!(rows, (; trial, result...))
            append!(element_rows, element_exposure_rows(
                response,
                view(samples.horizontal_directions, trial, :),
                view(samples.vertical_directions, trial, :),
                trial,
            ))
        end
    end
    summary = summarize_attribution(rows)
    mkpath(output_dir)
    direction_path = write_namedtuple_csv(joinpath(output_dir, "direction_attribution.csv"), rows)
    summary_path = write_namedtuple_csv(joinpath(output_dir, "attribution_summary.csv"), [summary])
    element_path = write_namedtuple_csv(
        joinpath(output_dir, "element_exposure_directions.csv"), element_rows,
    )
    metadata = Dict(
        "format" => "cesr-x-quadratic-internal-exposure-v2",
        "date" => string(Dates.today()),
        "role" => "causal diagnostic for the X-plane Q_vv versus Q_hh imbalance",
        "method" => "one all-corrector first-order GTPSA response at sextupoles joined to adopted second-order GTPSA detector Q",
        "normalization" => "original result uses equal corrector RMS; equal-source-exposure shares are diagnostic only",
        "internal_sampling" => "average of entrance and exit TPSA maps for every element with nonzero normal Kn2",
        "source_exposure" => "sum(abs(Kn2L) * first_order_orbit^2) over active normal sextupoles",
        "trials" => trials, "seed" => seed, "base_kick_rad" => base_kick, "rho" => rho,
        "input_csv" => input_path,
        "control_count" => length(names),
        "control_names" => names,
        "state_dimension" => length(closed_orbit),
        "active_normal_sextupoles" => length(inventory),
        "gtpsa_direction_csv" => gtpsa_direction_path,
        "solve_seconds" => solve_seconds,
        "direction_csv" => direction_path, "summary_csv" => summary_path,
        "element_direction_csv" => element_path,
    )
    merge!(metadata, ring_metadata(config; ring))
    metadata_path = joinpath(output_dir, "metadata.toml")
    open(metadata_path, "w") do io
        TOML.print(io, metadata; sorted=true)
    end
    @printf("Solve time: %.3f s\n", solve_seconds)
    println("Direction attribution: $direction_path")
    println("Element directions:     $element_path")
    println("Summary:               $summary_path")
    println("Metadata:              $metadata_path")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main_internal_exposure())
end

#!/usr/bin/env julia

"""
Measure the validity range of the first-order CESR detector-orbit response.

For every error scenario, Gaussian directions are generated once and normalized
to unit RMS over the active controls. Scaling a direction by `rho * base_kick`
therefore gives every trial exactly the requested active-control RMS radius.
The same directions are reused at every rho value.

The direct first-order prediction at all configured detector observables is compared with the
converged nonlinear SciBmad RF-on closed orbit. The nonlinear reference uses
the maintained response-initialized, frozen-nominal-Jacobian solver with a
full-AD Newton fallback and an independent final closure check.
"""

const ERROR_ANALYSIS_HERE = @__DIR__
const CALCULATION_DIR = normpath(joinpath(ERROR_ANALYSIS_HERE, "..", "Orbit_Calculation"))
include(joinpath(CALCULATION_DIR, "benchmark_scibmad.jl"))
include(joinpath(ERROR_ANALYSIS_HERE, "ring_analysis_config.jl"))
using .RingErrorAnalysisConfig

using Dates
using Random

const RHO_SCENARIOS = ("all", "horizontal", "vertical")

function parse_rho_args(args)
    options = Dict{String,String}(
        "ring" => "latest",
        "rhos" => "0,0.1,0.14,0.2,0.28",
        "trials" => "600",
        "seed" => "20260803",
        "base-kick-rad" => "5e-6",
        # Keep new runs separate from the committed historical CESR output.
        "output-dir" => "",
        # Empty values are resolved after parsing so --ring=legacy selects
        # the archived paths explicitly while the default uses latest_cesr.
        "inputs" => "",
        "detector-response" => "",
        "closed-orbit-response" => "",
        "reltol" => "1e-12",
        "abstol" => "1e-13",
        "maxiter" => "100",
        "warmup-samples" => "2",
        "response-method" => "gtpsa",
        "recompute-response" => "false",
        "response-step-rad" => "1e-7",
        "response-controls-per-batch" => "8",
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
    options["response-method"] = canonical_response_method(options["response-method"])
    lowercase(options["recompute-response"]) in ("true", "false") ||
        error("--recompute-response must be true or false")
    paths = default_ring_paths(; ring, response_method=options["response-method"])
    isempty(options["output-dir"]) &&
        (options["output-dir"] = joinpath(
            ERROR_ANALYSIS_HERE,
            "response_rho_sweep_600",
            ring_artifact_id(ring),
        ))
    isempty(options["inputs"]) && (options["inputs"] = paths.inputs)
    isempty(options["detector-response"]) &&
        (options["detector-response"] = paths.detector_response)
    isempty(options["closed-orbit-response"]) &&
        (options["closed-orbit-response"] = paths.closed_orbit_response)
    return options
end

function parse_rhos(text::AbstractString)
    values = parse.(Float64, split(text, ','))
    isempty(values) && error("--rhos must contain at least one value")
    all(isfinite, values) || error("Every rho must be finite")
    all(>=(0), values) || error("Every rho must be nonnegative")
    issorted(values) || error("--rhos must be sorted in ascending order")
    length(unique(values)) == length(values) || error("--rhos contains duplicates")
    first(values) == 0.0 || error("--rhos must start with the zero baseline")
    return values
end

function unquote_csv_field(field::AbstractString)
    value = strip(field)
    if length(value) >= 2 && first(value) == '"' && last(value) == '"'
        value = value[2:end-1]
    end
    return replace(value, "\"\"" => "\"")
end

function read_detector_response(path::AbstractString, control_names, detectors; layout=nothing)
    lines = readlines(path)
    isnothing(layout) && (layout = observable_layout(detectors))
    n_observables = length(layout.labels)
    length(lines) == n_observables + 1 ||
        error("Detector response must have one header and $n_observables rows: $path")
    header = unquote_csv_field.(split(lines[1], ','))
    first(header) == "observable" ||
        error("Detector response first column must be observable")
    String.(header[2:end]) == control_names ||
        error("Detector-response control names/order do not match the orbit input")

    expected_labels = layout.labels
    response = Matrix{Float64}(undef, n_observables, length(control_names))
    labels = Vector{String}(undef, n_observables)
    for row in 1:n_observables
        fields = unquote_csv_field.(split(lines[row + 1], ','))
        length(fields) == length(header) ||
            error("Detector-response row $row has the wrong width")
        labels[row] = fields[1]
        response[row, :] .= parse.(Float64, fields[2:end])
    end
    labels == expected_labels ||
        error("Detector-response observable labels/order do not match SciBmad")
    all(isfinite, response) || error("Detector response contains non-finite values")
    return response
end

function active_control_indices(names, scenario::String; config=nothing, model=nothing)
    scenario == "all" && return collect(eachindex(names))
    scenario in ("horizontal", "vertical") || error("Unknown rho scenario: $scenario")
    return control_group_indices(names, scenario; config, model)
end

function gaussian_unit_rms_directions(
    generator::AbstractRNG,
    trials::Int,
    n_controls::Int,
    active_indices,
)
    directions = zeros(trials, n_controls)
    for trial in 1:trials
        active = randn(generator, length(active_indices))
        active_rms = sqrt(sum(abs2, active) / length(active))
        iszero(active_rms) && error("Generated an all-zero Gaussian direction")
        directions[trial, active_indices] .= active ./ active_rms
    end
    return directions
end

function generate_rho_samples(names, rhos, trials, seed, base_kick; config=nothing, model=nothing)
    n_controls = length(names)
    positive_rhos = filter(>(0), rhos)
    n_samples = 1 + length(RHO_SCENARIOS) * length(positive_rhos) * trials
    values = Matrix{Float64}(undef, n_samples, n_controls)
    scenarios = Vector{String}(undef, n_samples)
    sample_rhos = Vector{Float64}(undef, n_samples)
    trial_ids = Vector{Int}(undef, n_samples)
    active_rms = Vector{Float64}(undef, n_samples)
    global_rms = Vector{Float64}(undef, n_samples)
    maximum_kick = Vector{Float64}(undef, n_samples)

    row = 1
    values[row, :] .= 0.0
    scenarios[row] = "baseline"
    sample_rhos[row] = 0.0
    trial_ids[row] = 0
    active_rms[row] = 0.0
    global_rms[row] = 0.0
    maximum_kick[row] = 0.0
    for (scenario_index, scenario) in enumerate(RHO_SCENARIOS)
        active_indices = active_control_indices(names, scenario; config, model)
        generator = MersenneTwister(seed + scenario_index - 1)
        directions = gaussian_unit_rms_directions(
            generator,
            trials,
            n_controls,
            active_indices,
        )
        for rho in positive_rhos
            scale = rho * base_kick
            for trial in 1:trials
                row += 1
                values[row, :] .= scale .* view(directions, trial, :)
                scenarios[row] = scenario
                sample_rhos[row] = rho
                trial_ids[row] = trial
                active_rms[row] = sqrt(
                    sum(abs2, view(values, row, active_indices)) /
                    length(active_indices),
                )
                global_rms[row] = sqrt(sum(abs2, view(values, row, :)) / n_controls)
                maximum_kick[row] = maximum(abs, view(values, row, :))
            end
        end
    end
    row == n_samples || error("Internal rho-sample count mismatch")
    return (;
        sample_ids=collect(0:n_samples-1),
        values,
        scenarios,
        rhos=sample_rhos,
        trial_ids,
        active_rms,
        global_rms,
        maximum_kick,
    )
end

plane_rmse(error_row, columns) = isempty(columns) ? NaN :
    sqrt(sum(abs2, view(error_row, columns)) / length(columns))

function calculate_trial_metrics(samples, result, approximate; layout)
    difference = approximate .- result.observables
    n_samples = size(difference, 1)
    x_rmse = Vector{Float64}(undef, n_samples)
    y_rmse = similar(x_rmse)
    x_max_abs = similar(x_rmse)
    y_max_abs = similar(x_rmse)
    x_indices = plane_indices(layout, :x)
    y_indices = plane_indices(layout, :y)
    isempty(x_indices) && error("Observable layout has no horizontal/x plane")
    isempty(y_indices) && error("Observable layout has no vertical/y plane")
    for row in 1:n_samples
        x_rmse[row] = plane_rmse(view(difference, row, :), x_indices)
        y_rmse[row] = plane_rmse(view(difference, row, :), y_indices)
        x_max_abs[row] = maximum(abs, view(difference, row, x_indices))
        y_max_abs[row] = maximum(abs, view(difference, row, y_indices))
    end
    return (; difference, x_rmse, y_rmse, x_max_abs, y_max_abs)
end

function csv_number(value)
    value isa AbstractFloat && return @sprintf("%.17g", value)
    return string(value)
end

function write_trial_metrics(path, samples, result, metrics)
    mkpath(dirname(path))
    open(path, "w") do io
        println(
            io,
            "sample_id,scenario,rho,trial_id,active_rms_rad,global_rms_rad," *
            "maximum_abs_kick_rad,converged,newton_iterations,closure_norm," *
            "x_rmse_m,x_max_abs_m,y_rmse_m,y_max_abs_m",
        )
        for row in eachindex(samples.sample_ids)
            fields = (
                samples.sample_ids[row],
                samples.scenarios[row],
                samples.rhos[row],
                samples.trial_ids[row],
                samples.active_rms[row],
                samples.global_rms[row],
                samples.maximum_kick[row],
                result.converged[row],
                result.iterations[row],
                result.closure_norms[row],
                metrics.x_rmse[row],
                metrics.x_max_abs[row],
                metrics.y_rmse[row],
                metrics.y_max_abs[row],
            )
            println(io, join(csv_number.(fields), ','))
        end
    end
    return path
end

function write_summary(path, samples, result, metrics, requested_rhos)
    mkpath(dirname(path))
    open(path, "w") do io
        println(
            io,
            "scenario,rho,trials,converged_trials,mean_x_rmse_m,max_trial_x_rmse_m," *
            "mean_x_max_abs_m,global_x_max_abs_m,mean_y_rmse_m,max_trial_y_rmse_m," *
            "mean_y_max_abs_m,global_y_max_abs_m,mean_newton_iterations," *
            "max_closure_norm",
        )
        for scenario in RHO_SCENARIOS, rho in requested_rhos
            selected = if rho == 0.0
                findall(==("baseline"), samples.scenarios)
            else
                findall(
                    row -> samples.scenarios[row] == scenario && samples.rhos[row] == rho,
                    eachindex(samples.sample_ids),
                )
            end
            converged = [row for row in selected if result.converged[row]]
            x_rmse = metrics.x_rmse[converged]
            y_rmse = metrics.y_rmse[converged]
            x_max = metrics.x_max_abs[converged]
            y_max = metrics.y_max_abs[converged]
            iterations = Float64.(result.iterations[converged])
            closures = result.closure_norms[converged]
            fields = (
                scenario,
                rho,
                length(selected),
                length(converged),
                isempty(x_rmse) ? NaN : mean(x_rmse),
                isempty(x_rmse) ? NaN : maximum(x_rmse),
                isempty(x_max) ? NaN : mean(x_max),
                isempty(x_max) ? NaN : maximum(x_max),
                isempty(y_rmse) ? NaN : mean(y_rmse),
                isempty(y_rmse) ? NaN : maximum(y_rmse),
                isempty(y_max) ? NaN : mean(y_max),
                isempty(y_max) ? NaN : maximum(y_max),
                isempty(iterations) ? NaN : mean(iterations),
                isempty(closures) ? NaN : maximum(closures),
            )
            println(io, join(csv_number.(fields), ','))
        end
    end
    return path
end

function main_rho_sweep(args=ARGS; model_factory=nothing, config=nothing)
    options = parse_rho_args(args)
    ring = Symbol(lowercase(options["ring"]))
    model_factory, config = resolve_ring_model_factory(model_factory, config; ring)
    rhos = parse_rhos(options["rhos"])
    trials = parse(Int, options["trials"])
    trials >= 2 || error("--trials must be at least 2")
    seed = parse(Int, options["seed"])
    base_kick = parse(Float64, options["base-kick-rad"])
    isfinite(base_kick) && base_kick > 0 ||
        error("--base-kick-rad must be finite and positive")
    reltol = parse(Float64, options["reltol"])
    abstol = parse(Float64, options["abstol"])
    maxiter = parse(Int, options["maxiter"])
    requested_response_method = canonical_response_method(options["response-method"])
    recompute_response = lowercase(options["recompute-response"]) == "true"
    response_step_rad = parse(Float64, options["response-step-rad"])
    response_controls_per_batch =
        parse(Int, options["response-controls-per-batch"])
    isfinite(response_step_rad) && response_step_rad > 0 ||
        error("--response-step-rad must be finite and positive")
    response_controls_per_batch >= 1 ||
        error("--response-controls-per-batch must be at least 1")
    warmup_samples = parse(Int, options["warmup-samples"])
    warmup_samples >= 2 || error("--warmup-samples must be at least 2")
    output_dir = abspath(options["output-dir"])
    detector_response_path = abspath(options["detector-response"])
    closed_orbit_response_path = abspath(options["closed-orbit-response"])

    input_path = abspath(options["inputs"])
    input_reference = read_samples(input_path)
    names = input_reference.names
    validate_control_names(names, config)
    nominal_model = configured_model(model_factory, config; zero_value=0.0, rf_on=true)
    detectors = configured_detector_names(nominal_model, config)
    layout = observable_layout(detectors; config)
    detector_response = read_or_compute_detector_response(
        read_detector_response,
        detector_response_path,
        names,
        detectors;
        layout,
        model_factory,
        config,
        closed_orbit_path=closed_orbit_response_path,
        response_method=requested_response_method,
        recompute_response,
        response_step_rad,
        controls_per_batch=response_controls_per_batch,
        reltol,
        abstol,
        maxiter,
    )
    response_metadata = detector_response_cache_metadata(detector_response_path)
    response_method = if isnothing(response_metadata)
        "cache-sidecar-missing-or-legacy"
    else
        raw_method = String(get(response_metadata, "response_method", ""))
        isempty(raw_method) ? "cache-sidecar-missing-or-legacy" :
            try
                canonical_response_method(raw_method)
            catch
                raw_method
            end
    end
    actual_response_step_rad = isnothing(response_metadata) ? 0.0 :
        Float64(get(response_metadata, "response_step_rad", 0.0))
    actual_response_controls_per_batch = isnothing(response_metadata) ? 0 :
        Int(get(response_metadata, "controls_per_batch", 0))
    actual_response_reltol = isnothing(response_metadata) ? reltol :
        Float64(get(response_metadata, "reltol", reltol))
    actual_response_abstol = isnothing(response_metadata) ? abstol :
        Float64(get(response_metadata, "abstol", abstol))
    actual_response_maxiter = isnothing(response_metadata) ? maxiter :
        Int(get(response_metadata, "maxiter", maxiter))
    actual_response_pair_id = isnothing(response_metadata) ? "" :
        String(get(response_metadata, "response_pair_id", ""))
    actual_response_scibmad_version = isnothing(response_metadata) ? "" :
        String(get(response_metadata, "scibmad_version", ""))
    if ring_artifact_id(ring) == "latest_cesr"
        response_method == requested_response_method || error(
            "Latest-ring response cache method '$response_method' does not match " *
            "the requested '$requested_response_method' method",
        )
        isempty(actual_response_pair_id) && error(
            "Latest-ring response cache is missing response_pair_id provenance",
        )
    end
    size(detector_response) == (length(layout.labels), length(names)) || error(
        "Detector response shape $(size(detector_response)) does not match " *
        "the runtime layout $(length(layout.labels)) x $(length(names))",
    )

    samples = generate_rho_samples(names, rhos, trials, seed, base_kick; config, model=nominal_model)
    total_samples = size(samples.values, 1)
    BLAS.set_num_threads(max(1, Threads.nthreads()))
    @printf(
        "Orbit response rho sweep: %d samples = one shared baseline + %d scenarios x %d positive rho x %d trials\n",
        total_samples,
        length(RHO_SCENARIOS),
        count(>(0), rhos),
        trials,
    )
    @printf(
        "Julia threads: %d; BLAS threads: %d; base kick: %.3e rad\n",
        Threads.nthreads(),
        BLAS.get_num_threads(),
        base_kick,
    )

    warmup_count = min(warmup_samples, total_samples)
    warmup_seconds = @elapsed configured_simulate_batch(
        simulate_batch,
        names,
        samples.values[1:warmup_count, :];
        config,
        model_factory,
        initial_guess_mode="response-linear",
        jacobian_mode="frozen-nominal",
        response_matrix_cache=closed_orbit_response_path,
        recompute_response=false,
        response_method=requested_response_method,
        response_step_rad,
        response_controls_per_batch,
        reltol,
        abstol,
        maxiter,
    )
    @printf("Warmup (%d samples): %.3f s\n", warmup_count, warmup_seconds)

    exact_timed = @timed configured_simulate_batch(
        simulate_batch,
        names,
        samples.values;
        config,
        model_factory,
        initial_guess_mode="response-linear",
        jacobian_mode="frozen-nominal",
        response_matrix_cache=closed_orbit_response_path,
        recompute_response=false,
        response_method=requested_response_method,
        response_step_rad,
        response_controls_per_batch,
        reltol,
        abstol,
        maxiter,
    )
    result = exact_timed.value
    baseline_index = findfirst(
        row -> samples.scenarios[row] == "baseline" && samples.rhos[row] == 0.0,
        eachindex(samples.sample_ids),
    )
    isnothing(baseline_index) && error("The rho=0 all-control baseline is missing")
    result.converged[baseline_index] || error("The nominal baseline did not converge")
    nominal_observables = vec(result.observables[baseline_index, :])
    response_timed = @timed begin
        approximate = repeat(reshape(nominal_observables, 1, :), total_samples, 1)
        mul!(approximate, samples.values, transpose(detector_response), 1.0, 1.0)
        approximate
    end
    metrics = calculate_trial_metrics(samples, result, response_timed.value; layout)

    trial_path = write_trial_metrics(
        joinpath(output_dir, "rho_sweep_trial_errors.csv"),
        samples,
        result,
        metrics,
    )
    summary_path = write_summary(
        joinpath(output_dir, "rho_sweep_summary.csv"),
        samples,
        result,
        metrics,
        rhos,
    )
    metadata_path = joinpath(output_dir, "rho_sweep_metadata.toml")
    metadata = Dict(
        "format" => "cesr-orbit-response-rho-sweep-v1",
        "date" => string(Dates.today()),
        "rho_definition" => "active-control RMS kick divided by base_kick_rad",
        "direction_distribution" => "Gaussian direction normalized to exact unit RMS over active controls",
        "direction_reuse" => "the same directions are reused at every rho within each scenario",
        "scenarios" => collect(RHO_SCENARIOS),
        "rho_values" => rhos,
        "trials_per_rho_scenario" => trials,
        "rho_zero_is_one_shared_baseline" => true,
        "seed" => seed,
        "scenario_seeds" => Dict(
            scenario => seed + index - 1
            for (index, scenario) in enumerate(RHO_SCENARIOS)
        ),
        "base_kick_rad" => base_kick,
        "total_samples" => total_samples,
        "input_csv" => input_path,
        "control_count" => length(names),
        "control_names" => names,
        "horizontal_control_count" => length(active_control_indices(names, "horizontal"; config, model=nominal_model)),
        "vertical_control_count" => length(active_control_indices(names, "vertical"; config, model=nominal_model)),
        "observable_count" => length(layout.labels),
        "detector_count" => length(detectors),
        "observable_labels" => layout.labels,
        "detector_count_per_plane" => length(plane_indices(layout, :x)),
        "detector_response_csv" => detector_response_path,
        "closed_orbit_response_csv" => closed_orbit_response_path,
        "response_method" => response_method,
        "response_step_rad" => actual_response_step_rad,
        "response_controls_per_batch" =>
            actual_response_controls_per_batch,
        "response_reltol" => actual_response_reltol,
        "response_abstol" => actual_response_abstol,
        "response_maxiter" => actual_response_maxiter,
        "response_pair_id" => actual_response_pair_id,
        "response_scibmad_version" => actual_response_scibmad_version,
        "exact_reference" => "nonlinear RF-on SciBmad closed orbit with frozen nominal Jacobian, closure check, and full-AD fallback",
        "reltol" => reltol,
        "abstol" => abstol,
        "maxiter" => maxiter,
        "converged_count" => count(result.converged),
        "failed_count" => count(.!result.converged),
        "fallback_count" => result.fallback_count,
        "fallback_success_count" => result.fallback_success_count,
        "maximum_closure_norm" => maximum(result.closure_norms[result.converged]),
        "warmup_samples" => warmup_count,
        "warmup_seconds" => warmup_seconds,
        "exact_physics_seconds" => exact_timed.time,
        "response_evaluation_seconds" => response_timed.time,
        "julia_version" => string(VERSION),
        "julia_threads" => Threads.nthreads(),
        "blas_threads" => BLAS.get_num_threads(),
        "trial_errors_csv" => trial_path,
        "summary_csv" => summary_path,
    )
    merge!(metadata, ring_metadata(config; ring))
    mkpath(output_dir)
    open(metadata_path, "w") do io
        TOML.print(io, metadata; sorted=true)
    end
    @printf(
        "Exact nonlinear reference: %.3f s; linear response evaluation: %.6f s\n",
        exact_timed.time,
        response_timed.time,
    )
    @printf(
        "Converged %d/%d; fallback %d; maximum closure %.3e\n",
        count(result.converged),
        total_samples,
        result.fallback_count,
        maximum(result.closure_norms[result.converged]),
    )
    println("Trial errors: $trial_path")
    println("Summary:      $summary_path")
    println("Metadata:     $metadata_path")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main_rho_sweep())
end

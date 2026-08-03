#!/usr/bin/env julia

"""
Measure the validity range of the first-order CESR detector-orbit response.

For every error scenario, Gaussian directions are generated once and normalized
to unit RMS over the active controls. Scaling a direction by `rho * base_kick`
therefore gives every trial exactly the requested active-control RMS radius.
The same directions are reused at every rho value.

The direct first-order prediction at all 99 detectors is compared with the
converged nonlinear SciBmad RF-on closed orbit. The nonlinear reference uses
the maintained response-initialized, frozen-nominal-Jacobian solver with a
full-AD Newton fallback and an independent final closure check.
"""

include(joinpath(@__DIR__, "benchmark_scibmad.jl"))

using Dates
using Random

const RHO_SCENARIOS = ("all", "horizontal", "vertical")

function parse_rho_args(args)
    options = Dict{String,String}(
        "rhos" => "0,0.1,0.14,0.2,0.28",
        "trials" => "600",
        "seed" => "20260803",
        "base-kick-rad" => "5e-6",
        "output-dir" => joinpath(HERE, "results", "response_rho_sweep"),
        "detector-response" => joinpath(
            PROJECT_ROOT,
            "bmad_comparison",
            "bmad_control_response_rf_on",
            "scibmad_control_response_rf_on.csv",
        ),
        "closed-orbit-response" => joinpath(
            HERE,
            "reference",
            "closed_orbit_response_6x119.csv",
        ),
        "reltol" => "1e-12",
        "abstol" => "1e-13",
        "maxiter" => "100",
        "warmup-samples" => "2",
    )
    for argument in args
        startswith(argument, "--") ||
            error("Arguments must have --name=value form: $argument")
        fields = split(argument[3:end], '='; limit=2)
        length(fields) == 2 || error("Missing value in argument: $argument")
        haskey(options, fields[1]) || error("Unknown option: --$(fields[1])")
        options[fields[1]] = fields[2]
    end
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

function read_detector_response(path::AbstractString, control_names, detectors)
    lines = readlines(path)
    length(lines) == 199 ||
        error("Detector response must have one header and 198 rows: $path")
    header = unquote_csv_field.(split(lines[1], ','))
    first(header) == "observable" ||
        error("Detector response first column must be observable")
    String.(header[2:end]) == control_names ||
        error("Detector-response control names/order do not match the orbit input")

    expected_labels = vcat(
        [name * ":x" for name in detectors],
        [name * ":y" for name in detectors],
    )
    response = Matrix{Float64}(undef, 198, length(control_names))
    labels = Vector{String}(undef, 198)
    for row in 1:198
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

function active_control_indices(names, scenario::String)
    scenario == "all" && return collect(eachindex(names))
    prefix = scenario == "horizontal" ? 'H' : scenario == "vertical" ? 'V' :
        error("Unknown rho scenario: $scenario")
    indices = findall(name -> startswith(name, string(prefix)), names)
    expected = scenario == "horizontal" ? 58 : 61
    length(indices) == expected ||
        error("Expected $expected $scenario controls, found $(length(indices))")
    return indices
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

function generate_rho_samples(names, rhos, trials, seed, base_kick)
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
        active_indices = active_control_indices(names, scenario)
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

plane_rmse(error_row, columns) = sqrt(sum(abs2, view(error_row, columns)) / length(columns))

function calculate_trial_metrics(samples, result, approximate)
    difference = approximate .- result.observables
    n_samples = size(difference, 1)
    x_rmse = Vector{Float64}(undef, n_samples)
    y_rmse = similar(x_rmse)
    x_max_abs = similar(x_rmse)
    y_max_abs = similar(x_rmse)
    for row in 1:n_samples
        x_rmse[row] = plane_rmse(view(difference, row, :), 1:99)
        y_rmse[row] = plane_rmse(view(difference, row, :), 100:198)
        x_max_abs[row] = maximum(abs, view(difference, row, 1:99))
        y_max_abs[row] = maximum(abs, view(difference, row, 100:198))
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

function main_rho_sweep(args=ARGS)
    options = parse_rho_args(args)
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
    warmup_samples = parse(Int, options["warmup-samples"])
    warmup_samples >= 2 || error("--warmup-samples must be at least 2")
    output_dir = abspath(options["output-dir"])
    detector_response_path = abspath(options["detector-response"])
    closed_orbit_response_path = abspath(options["closed-orbit-response"])

    control_reference = read_response_matrix(
        closed_orbit_response_path,
        read_samples(joinpath(HERE, "inputs", "cesr_corrector_samples_1000.csv")).names,
    )
    size(control_reference) == (6, 119) || error("Invalid control reference")
    names = read_samples(joinpath(HERE, "inputs", "cesr_corrector_samples_1000.csv")).names
    nominal_model = load_cesr_model(zero_value=0.0, rf_on=true)
    detectors = detector_names(nominal_model.ring)
    detector_response = read_detector_response(
        detector_response_path,
        names,
        detectors,
    )

    samples = generate_rho_samples(names, rhos, trials, seed, base_kick)
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
    warmup_seconds = @elapsed simulate_batch(
        names,
        samples.values[1:warmup_count, :];
        initial_guess_mode="response-linear",
        jacobian_mode="frozen-nominal",
        response_matrix_cache=closed_orbit_response_path,
        recompute_response=false,
        reltol,
        abstol,
        maxiter,
    )
    @printf("Warmup (%d samples): %.3f s\n", warmup_count, warmup_seconds)

    exact_timed = @timed simulate_batch(
        names,
        samples.values;
        initial_guess_mode="response-linear",
        jacobian_mode="frozen-nominal",
        response_matrix_cache=closed_orbit_response_path,
        recompute_response=false,
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
    metrics = calculate_trial_metrics(samples, result, response_timed.value)

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
        "detector_count_per_plane" => 99,
        "detector_response_csv" => detector_response_path,
        "closed_orbit_response_csv" => closed_orbit_response_path,
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

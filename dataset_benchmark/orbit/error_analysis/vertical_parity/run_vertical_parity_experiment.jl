#!/usr/bin/env julia

"""
Separate even and odd nonlinear closed-orbit response for vertical correctors.

For every fixed random vertical-corrector direction d and positive radius rho,
the experiment evaluates both +rho*d and -rho*d.  For detector plane u,

    even     = (u(+) + u(-))/2 - u(0)
    odd_nl   = (u(+) - u(-))/2 - J_u * (+delta_k)

so even starts at second order while odd_nl starts at third order when the
linear response matrix J is consistent with the nonlinear reference.
"""

const PARITY_HERE = @__DIR__
const PARITY_ERROR_ANALYSIS_DIR = normpath(joinpath(PARITY_HERE, ".."))
const PARITY_ORBIT_ROOT = normpath(joinpath(PARITY_ERROR_ANALYSIS_DIR, ".."))
const PARITY_CALCULATION_DIR = joinpath(PARITY_ORBIT_ROOT, "Orbit_Calculation")
include(joinpath(PARITY_ERROR_ANALYSIS_DIR, "run_response_rho_sweep.jl"))

using Dates
using Random
using Statistics
using TOML

function parse_parity_args(args)
    options = Dict{String,String}(
        "rhos" => "0.05,0.075,0.1,0.14,0.2,0.28,0.4,0.57,0.8,1.13,1.6,1.8,2.0,2.26,2.5,2.8,3.2,3.6,4.0,4.53,5.0,5.6,6.0,6.4,7.5,8.8,10.05",
        "trials" => "100",
        "seed" => "20260804",
        "base-kick-rad" => "5e-6",
        "output-dir" => joinpath(PARITY_HERE, "results"),
        "detector-response" => joinpath(
            PROJECT_ROOT,
            "bmad_comparison",
            "bmad_control_response_rf_on",
            "scibmad_control_response_rf_on.csv",
        ),
        "closed-orbit-response" => joinpath(
            PARITY_ORBIT_ROOT,
            "reference",
            "closed_orbit_response_6x119.csv",
        ),
        "reltol" => "1e-12",
        "abstol" => "1e-13",
        "maxiter" => "100",
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

function generate_vertical_pairs(names, rhos, trials, seed, base_kick)
    active_indices = active_control_indices(names, "vertical")
    directions = gaussian_unit_rms_directions(
        MersenneTwister(seed),
        trials,
        length(names),
        active_indices,
    )
    n_samples = 1 + 2 * length(rhos) * trials
    values = zeros(n_samples, length(names))
    plus_rows = Matrix{Int}(undef, length(rhos), trials)
    minus_rows = similar(plus_rows)
    row = 1
    for (rho_index, rho) in enumerate(rhos), trial in 1:trials
        delta = rho * base_kick .* view(directions, trial, :)
        row += 1
        plus_rows[rho_index, trial] = row
        values[row, :] .= delta
        row += 1
        minus_rows[rho_index, trial] = row
        values[row, :] .= -delta
    end
    row == n_samples || error("Internal paired-sample count mismatch")
    return (; values, directions, active_indices, plus_rows, minus_rows)
end

vector_rmse(values) = sqrt(sum(abs2, values) / length(values))

function paired_metrics(result, response, samples, rhos, base_kick)
    baseline = vec(result.observables[1, :])
    rows = NamedTuple[]
    for (rho_index, rho) in enumerate(rhos), trial in axes(samples.directions, 1)
        plus_row = samples.plus_rows[rho_index, trial]
        minus_row = samples.minus_rows[rho_index, trial]
        both_converged = result.converged[plus_row] && result.converged[minus_row]
        if !both_converged
            push!(rows, (;
                rho,
                trial,
                active_rms_rad=rho * base_kick,
                plus_converged=result.converged[plus_row],
                minus_converged=result.converged[minus_row],
                max_closure_norm=max(
                    result.closure_norms[plus_row],
                    result.closure_norms[minus_row],
                ),
                x_even_rmse_m=NaN,
                x_odd_nl_rmse_m=NaN,
                x_mean_full_residual_rmse_m=NaN,
                y_even_rmse_m=NaN,
                y_odd_nl_rmse_m=NaN,
                y_mean_full_residual_rmse_m=NaN,
            ))
            continue
        end

        delta = rho * base_kick .* view(samples.directions, trial, :)
        linear = response * delta
        plus_delta = vec(result.observables[plus_row, :]) .- baseline
        minus_delta = vec(result.observables[minus_row, :]) .- baseline
        even = (plus_delta .+ minus_delta) ./ 2
        odd_nl = (plus_delta .- minus_delta) ./ 2 .- linear
        plus_residual = plus_delta .- linear
        minus_residual = minus_delta .+ linear

        xcols = 1:99
        ycols = 100:198
        push!(rows, (;
            rho,
            trial,
            active_rms_rad=rho * base_kick,
            plus_converged=true,
            minus_converged=true,
            max_closure_norm=max(
                result.closure_norms[plus_row],
                result.closure_norms[minus_row],
            ),
            x_even_rmse_m=vector_rmse(view(even, xcols)),
            x_odd_nl_rmse_m=vector_rmse(view(odd_nl, xcols)),
            x_mean_full_residual_rmse_m=(
                vector_rmse(view(plus_residual, xcols)) +
                vector_rmse(view(minus_residual, xcols))
            ) / 2,
            y_even_rmse_m=vector_rmse(view(even, ycols)),
            y_odd_nl_rmse_m=vector_rmse(view(odd_nl, ycols)),
            y_mean_full_residual_rmse_m=(
                vector_rmse(view(plus_residual, ycols)) +
                vector_rmse(view(minus_residual, ycols))
            ) / 2,
        ))
    end
    return rows
end

function write_namedtuple_csv(path, rows)
    isempty(rows) && error("Cannot write an empty CSV")
    mkpath(dirname(path))
    columns = propertynames(first(rows))
    open(path, "w") do io
        println(io, join(columns, ','))
        for row in rows
            println(io, join((csv_number(getproperty(row, column)) for column in columns), ','))
        end
    end
    return path
end

function summarize_pairs(rows, rhos)
    summary = NamedTuple[]
    metric_names = (
        :x_even_rmse_m,
        :x_odd_nl_rmse_m,
        :x_mean_full_residual_rmse_m,
        :y_even_rmse_m,
        :y_odd_nl_rmse_m,
        :y_mean_full_residual_rmse_m,
    )
    for rho in rhos
        selected = [row for row in rows if row.rho == rho]
        complete = [row for row in selected if row.plus_converged && row.minus_converged]
        values = Dict{Symbol,Float64}()
        for metric in metric_names
            values[Symbol("mean_", metric)] = isempty(complete) ? NaN :
                mean(getproperty(row, metric) for row in complete)
            values[Symbol("max_", metric)] = isempty(complete) ? NaN :
                maximum(getproperty(row, metric) for row in complete)
        end
        push!(summary, (;
            rho,
            requested_pairs=length(selected),
            converged_pairs=length(complete),
            maximum_closure_norm=isempty(complete) ? NaN :
                maximum(row.max_closure_norm for row in complete),
            (; values...)...,
        ))
    end
    return summary
end

function main_parity(args=ARGS)
    options = parse_parity_args(args)
    rhos = parse.(Float64, split(options["rhos"], ','))
    all(>(0), rhos) || error("All parity radii must be positive")
    issorted(rhos) || error("Parity radii must be sorted")
    trials = parse(Int, options["trials"])
    trials >= 2 || error("--trials must be at least 2")
    seed = parse(Int, options["seed"])
    base_kick = parse(Float64, options["base-kick-rad"])
    reltol = parse(Float64, options["reltol"])
    abstol = parse(Float64, options["abstol"])
    maxiter = parse(Int, options["maxiter"])
    output_dir = abspath(options["output-dir"])

    input_reference = read_samples(
        joinpath(PARITY_CALCULATION_DIR, "inputs", "cesr_corrector_samples_1000.csv"),
    )
    names = input_reference.names
    model = load_cesr_model(zero_value=0.0, rf_on=true)
    detectors = detector_names(model.ring)
    response = read_detector_response(
        abspath(options["detector-response"]),
        names,
        detectors,
    )
    samples = generate_vertical_pairs(names, rhos, trials, seed, base_kick)

    @printf(
        "Vertical parity experiment: %d nonlinear states = baseline + 2 signs x %d radii x %d directions\n",
        size(samples.values, 1),
        length(rhos),
        trials,
    )
    timed = @timed simulate_batch(
        names,
        samples.values;
        initial_guess_mode="response-linear",
        jacobian_mode="frozen-nominal",
        response_matrix_cache=abspath(options["closed-orbit-response"]),
        recompute_response=false,
        reltol,
        abstol,
        maxiter,
    )
    result = timed.value
    result.converged[1] || error("Nominal baseline did not converge")
    pairs = paired_metrics(result, response, samples, rhos, base_kick)
    summary = summarize_pairs(pairs, rhos)
    pair_path = write_namedtuple_csv(joinpath(output_dir, "vertical_parity_pairs.csv"), pairs)
    summary_path = write_namedtuple_csv(joinpath(output_dir, "vertical_parity_summary.csv"), summary)

    metadata = Dict(
        "format" => "cesr-vertical-corrector-parity-v1",
        "date" => string(Dates.today()),
        "definition" => "paired + and - vertical-corrector directions; even=(u(+)+u(-))/2-u(0); odd_nl=(u(+)-u(-))/2-J*delta_k",
        "rhos" => rhos,
        "trials_per_rho" => trials,
        "seed" => seed,
        "base_kick_rad" => base_kick,
        "total_nonlinear_states" => size(samples.values, 1),
        "converged_states" => count(result.converged),
        "failed_states" => count(.!result.converged),
        "fallback_count" => result.fallback_count,
        "fallback_success_count" => result.fallback_success_count,
        "maximum_closure_norm" => maximum(result.closure_norms[result.converged]),
        "solve_seconds" => timed.time,
        "reltol" => reltol,
        "abstol" => abstol,
        "maxiter" => maxiter,
        "pair_csv" => pair_path,
        "summary_csv" => summary_path,
    )
    mkpath(output_dir)
    metadata_path = joinpath(output_dir, "vertical_parity_metadata.toml")
    open(metadata_path, "w") do io
        TOML.print(io, metadata; sorted=true)
    end

    @printf(
        "Converged %d/%d states; fallback %d; maximum closure %.3e; solve %.2f s\n",
        count(result.converged),
        length(result.converged),
        result.fallback_count,
        maximum(result.closure_norms[result.converged]),
        timed.time,
    )
    println("Paired metrics: $pair_path")
    println("Summary:        $summary_path")
    println("Metadata:       $metadata_path")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main_parity())
end

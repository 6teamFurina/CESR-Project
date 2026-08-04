#!/usr/bin/env julia

"""
Decompose the all-corrector nonlinear orbit residual into pure-horizontal,
pure-vertical, and mixed horizontal--vertical even-sign components.

For each fixed pair of unit-RMS directions h and v and each positive rho, this
experiment evaluates the pure +/-h and +/-v states and all four joint sign
states. With all orbits expressed relative to the shared zero-input orbit,

    Q_hh = (u(+h,0) + u(-h,0))/2
    Q_vv = (u(0,+v) + u(0,-v))/2
    Q_hv = (u(+h,+v) - u(+h,-v) - u(-h,+v) + u(-h,-v))/4.

The signed second-order reconstruction is Q_hh + Q_vv + sh*sv*Q_hv.
"""

const MIXED_HERE = @__DIR__
const MIXED_ERROR_ANALYSIS_DIR = normpath(joinpath(MIXED_HERE, ".."))
const MIXED_ORBIT_ROOT = normpath(joinpath(MIXED_ERROR_ANALYSIS_DIR, ".."))
const MIXED_CALCULATION_DIR = joinpath(MIXED_ORBIT_ROOT, "Orbit_Calculation")
include(joinpath(MIXED_ERROR_ANALYSIS_DIR, "run_response_rho_sweep.jl"))

using Dates
using Random
using Statistics
using TOML

function parse_mixed_args(args)
    options = Dict{String,String}(
        "rhos" => "0.1,0.14,0.2,0.28,0.4,0.57,0.8,1.13",
        "trials" => "100",
        "seed" => "20260804",
        "base-kick-rad" => "5e-6",
        "output-dir" => joinpath(MIXED_HERE, "results"),
        "detector-response" => joinpath(
            PROJECT_ROOT, "bmad_comparison", "bmad_control_response_rf_on",
            "scibmad_control_response_rf_on.csv",
        ),
        "closed-orbit-response" => joinpath(
            MIXED_ORBIT_ROOT, "reference", "closed_orbit_response_6x119.csv",
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

function generate_mixed_samples(names, rhos, trials, seed, base_kick)
    horizontal_indices = active_control_indices(names, "horizontal")
    vertical_indices = active_control_indices(names, "vertical")
    horizontal_directions = gaussian_unit_rms_directions(
        MersenneTwister(seed), trials, length(names), horizontal_indices,
    )
    vertical_directions = gaussian_unit_rms_directions(
        MersenneTwister(seed + 1), trials, length(names), vertical_indices,
    )

    n_samples = 1 + 8 * length(rhos) * trials
    values = zeros(n_samples, length(names))
    hp = Matrix{Int}(undef, length(rhos), trials)
    hm = similar(hp)
    vp = similar(hp)
    vm = similar(hp)
    pp = similar(hp)
    pm = similar(hp)
    mp = similar(hp)
    mm = similar(hp)
    row_arrays = (; hp, hm, vp, vm, pp, pm, mp, mm)
    row = 1
    for (rho_index, rho) in enumerate(rhos), trial in 1:trials
        h = rho * base_kick .* view(horizontal_directions, trial, :)
        v = rho * base_kick .* view(vertical_directions, trial, :)
        for (name, candidate) in (
            (:hp, h), (:hm, -h), (:vp, v), (:vm, -v),
            (:pp, h + v), (:pm, h - v), (:mp, -h + v), (:mm, -h - v),
        )
            row += 1
            getproperty(row_arrays, name)[rho_index, trial] = row
            values[row, :] .= candidate
        end
    end
    row == n_samples || error("Internal mixed-term sample count mismatch")
    return (;
        values,
        horizontal_directions,
        vertical_directions,
        horizontal_indices,
        vertical_indices,
        rows=row_arrays,
    )
end

vector_rmse(values) = sqrt(sum(abs2, values) / length(values))

function aggregate_signed_rmse(vectors, columns)
    return sqrt(
        sum(sum(abs2, view(vector, columns)) for vector in vectors) /
        (length(vectors) * length(columns)),
    )
end

function plane_metrics(qhh, qvv, qhv, residuals, columns)
    pure = qhh + qvv
    reconstructions = [pure + sign * qhv for sign in (1, -1, -1, 1)]
    pure_remainders = [residual - pure for residual in residuals]
    full_remainders = [
        residual - reconstruction
        for (residual, reconstruction) in zip(residuals, reconstructions)
    ]
    qhh_rmse = vector_rmse(view(qhh, columns))
    qvv_rmse = vector_rmse(view(qvv, columns))
    qhv_rmse = vector_rmse(view(qhv, columns))
    pure_rmse = vector_rmse(view(pure, columns))
    second_order_rmse = sqrt(pure_rmse^2 + qhv_rmse^2)
    exact_rmse = aggregate_signed_rmse(residuals, columns)
    pure_remainder_rmse = aggregate_signed_rmse(pure_remainders, columns)
    full_remainder_rmse = aggregate_signed_rmse(full_remainders, columns)
    mixed_energy_share = qhv_rmse^2 / (pure_rmse^2 + qhv_rmse^2)
    reconstruction_improvement = 1 - full_remainder_rmse^2 / pure_remainder_rmse^2
    return (;
        qhh_rmse_m=qhh_rmse,
        qvv_rmse_m=qvv_rmse,
        qhv_rmse_m=qhv_rmse,
        pure_rmse_m=pure_rmse,
        second_order_rmse_m=second_order_rmse,
        exact_residual_rmse_m=exact_rmse,
        pure_only_remainder_rmse_m=pure_remainder_rmse,
        full_reconstruction_remainder_rmse_m=full_remainder_rmse,
        mixed_to_pure_ratio=qhv_rmse / pure_rmse,
        mixed_energy_share,
        reconstruction_improvement,
        relative_reconstruction_remainder=full_remainder_rmse / exact_rmse,
    )
end

function mixed_metrics(result, response, samples, rhos, base_kick)
    baseline = vec(result.observables[1, :])
    rows = NamedTuple[]
    joint_names = (:pp, :pm, :mp, :mm)
    joint_signs = ((1, 1), (1, -1), (-1, 1), (-1, -1))
    for (rho_index, rho) in enumerate(rhos), trial in axes(samples.horizontal_directions, 1)
        state_rows = [getproperty(samples.rows, name)[rho_index, trial]
                      for name in (:hp, :hm, :vp, :vm, :pp, :pm, :mp, :mm)]
        all_converged = all(result.converged[state_rows])
        max_closure = maximum(result.closure_norms[state_rows])
        if !all_converged
            nan_metrics = (;
                (Symbol("$(plane)_$(metric)") => NaN
                 for plane in (:x, :y)
                 for metric in (
                     :qhh_rmse_m, :qvv_rmse_m, :qhv_rmse_m, :pure_rmse_m,
                     :second_order_rmse_m, :exact_residual_rmse_m,
                     :pure_only_remainder_rmse_m,
                     :full_reconstruction_remainder_rmse_m,
                     :mixed_to_pure_ratio, :mixed_energy_share,
                     :reconstruction_improvement, :relative_reconstruction_remainder,
                 ))...,
            )
            push!(rows, (; rho, trial, active_rms_rad=rho * base_kick,
                all_states_converged=false, max_closure_norm=max_closure,
                nan_metrics...))
            continue
        end

        orbit(name) = vec(result.observables[
            getproperty(samples.rows, name)[rho_index, trial], :
        ])
        qhh = (orbit(:hp) + orbit(:hm)) / 2 - baseline
        qvv = (orbit(:vp) + orbit(:vm)) / 2 - baseline
        qhv = (orbit(:pp) - orbit(:pm) - orbit(:mp) + orbit(:mm)) / 4
        h = rho * base_kick .* view(samples.horizontal_directions, trial, :)
        v = rho * base_kick .* view(samples.vertical_directions, trial, :)
        linear_h = response * h
        linear_v = response * v
        residuals = [
            orbit(name) - baseline - sh * linear_h - sv * linear_v
            for (name, (sh, sv)) in zip(joint_names, joint_signs)
        ]
        x_metrics = plane_metrics(qhh, qvv, qhv, residuals, 1:99)
        y_metrics = plane_metrics(qhh, qvv, qhv, residuals, 100:198)
        prefixed = (;
            (Symbol("x_$(name)") => value for (name, value) in pairs(x_metrics))...,
            (Symbol("y_$(name)") => value for (name, value) in pairs(y_metrics))...,
        )
        push!(rows, (; rho, trial, active_rms_rad=rho * base_kick,
            all_states_converged=true, max_closure_norm=max_closure, prefixed...))
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

function summarize_mixed(rows, rhos)
    summary = NamedTuple[]
    metric_names = propertynames(first(rows))[6:end]
    for rho in rhos
        selected = [row for row in rows if row.rho == rho]
        complete = [row for row in selected if row.all_states_converged]
        means = (;
            (Symbol("mean_$(metric)") => (isempty(complete) ? NaN :
                mean(getproperty(row, metric) for row in complete))
             for metric in metric_names)...,
        )
        push!(summary, (;
            rho,
            requested_direction_pairs=length(selected),
            converged_direction_pairs=length(complete),
            maximum_closure_norm=isempty(complete) ? NaN :
                maximum(row.max_closure_norm for row in complete),
            means...,
        ))
    end
    return summary
end

function main_mixed(args=ARGS)
    options = parse_mixed_args(args)
    rhos = parse.(Float64, split(options["rhos"], ','))
    !isempty(rhos) && all(isfinite, rhos) && all(>(0), rhos) ||
        error("All mixed-term radii must be finite and positive")
    issorted(rhos) && length(unique(rhos)) == length(rhos) ||
        error("Mixed-term radii must be sorted and unique")
    trials = parse(Int, options["trials"])
    trials >= 2 || error("--trials must be at least 2")
    seed = parse(Int, options["seed"])
    base_kick = parse(Float64, options["base-kick-rad"])
    isfinite(base_kick) && base_kick > 0 ||
        error("--base-kick-rad must be finite and positive")
    reltol = parse(Float64, options["reltol"])
    abstol = parse(Float64, options["abstol"])
    maxiter = parse(Int, options["maxiter"])
    output_dir = abspath(options["output-dir"])

    input_reference = read_samples(joinpath(
        MIXED_CALCULATION_DIR, "inputs", "cesr_corrector_samples_1000.csv",
    ))
    names = input_reference.names
    model = load_cesr_model(zero_value=0.0, rf_on=true)
    detectors = detector_names(model.ring)
    response = read_detector_response(
        abspath(options["detector-response"]), names, detectors,
    )
    samples = generate_mixed_samples(names, rhos, trials, seed, base_kick)

    @printf(
        "Mixed H/V experiment: %d nonlinear states = baseline + 8 states x %d radii x %d direction pairs\n",
        size(samples.values, 1), length(rhos), trials,
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
    direction_rows = mixed_metrics(result, response, samples, rhos, base_kick)
    summary = summarize_mixed(direction_rows, rhos)
    direction_path = write_namedtuple_csv(
        joinpath(output_dir, "mixed_term_directions.csv"), direction_rows,
    )
    summary_path = write_namedtuple_csv(
        joinpath(output_dir, "mixed_term_summary.csv"), summary,
    )

    metadata = Dict(
        "format" => "cesr-all-corrector-mixed-term-v1",
        "date" => string(Dates.today()),
        "definition" => "four-sign finite-difference decomposition into Q_hh, Q_vv, and Q_hv with signed vector-level reconstruction",
        "rhos" => rhos,
        "trials_per_rho" => trials,
        "seed" => seed,
        "horizontal_direction_seed" => seed,
        "vertical_direction_seed" => seed + 1,
        "direction_distribution" => "independent Gaussian H and V directions, each normalized to exact unit RMS in its active family and reused at every rho",
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
        "detector_response_csv" => abspath(options["detector-response"]),
        "closed_orbit_response_csv" => abspath(options["closed-orbit-response"]),
        "direction_csv" => direction_path,
        "summary_csv" => summary_path,
    )
    mkpath(output_dir)
    metadata_path = joinpath(output_dir, "mixed_term_metadata.toml")
    open(metadata_path, "w") do io
        TOML.print(io, metadata; sorted=true)
    end

    @printf(
        "Converged %d/%d states; fallback %d; maximum closure %.3e; solve %.2f s\n",
        count(result.converged), length(result.converged), result.fallback_count,
        maximum(result.closure_norms[result.converged]), timed.time,
    )
    println("Direction metrics: $direction_path")
    println("Summary:           $summary_path")
    println("Metadata:          $metadata_path")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main_mixed())
end

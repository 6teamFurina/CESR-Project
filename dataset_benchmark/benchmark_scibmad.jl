#!/usr/bin/env julia

"""
Generate exact CESR closed-orbit samples with SciBmad batch parameters.

The timed physics region includes assigning all sampled control arrays, solving
all closed orbits, and tracking the solved orbits to all DET_* markers. Model
loading, Julia compilation warmup, and CSV writing are reported separately.
"""

using Beamlines
using LinearAlgebra
using Printf
using SciBmad
using Statistics
using TOML

const HERE = @__DIR__
const PROJECT_ROOT = normpath(joinpath(HERE, ".."))
include(joinpath(PROJECT_ROOT, "cesr_model.jl"))

function parse_args(args)
    options = Dict{String,String}(
        "inputs" => joinpath(HERE, "inputs", "cesr_corrector_samples_1000.csv"),
        "output" => joinpath(HERE, "results", "formal_1000", "scibmad", "scibmad_rf_on_samples.csv"),
        "metadata" => joinpath(HERE, "results", "formal_1000", "scibmad", "scibmad_rf_on_metadata.toml"),
        "mode" => "rf_on",
        "reltol" => "1e-13",
        "abstol" => "1e-13",
        "maxiter" => "100",
        "warmup-samples" => "2",
        "initial-guess" => "zero",
        "jacobian-mode" => "full",
        "fallback-full-newton" => "true",
    )
    for argument in args
        startswith(argument, "--") ||
            error("Arguments must have --name=value form: $argument")
        fields = split(argument[3:end], "="; limit=2)
        length(fields) == 2 || error("Missing value in argument: $argument")
        haskey(options, fields[1]) || error("Unknown option: --$(fields[1])")
        options[fields[1]] = fields[2]
    end
    options["mode"] == "rf_on" ||
        error("The first benchmark release supports --mode=rf_on only")
    options["initial-guess"] in ("zero", "nominal-z0") ||
        error("--initial-guess must be zero or nominal-z0")
    options["jacobian-mode"] in ("full", "frozen-nominal") ||
        error("--jacobian-mode must be full or frozen-nominal")
    if options["jacobian-mode"] == "frozen-nominal"
        options["initial-guess"] == "nominal-z0" ||
            error("--jacobian-mode=frozen-nominal requires --initial-guess=nominal-z0")
    end
    lowercase(options["fallback-full-newton"]) in ("true", "false") ||
        error("--fallback-full-newton must be true or false")
    return options
end

function read_samples(path::AbstractString)
    lines = readlines(path)
    length(lines) >= 2 || error("Sample CSV has no data: $path")
    header = split(lines[1], ',')
    first(header) == "sample_id" || error("First CSV column must be sample_id")
    names = String.(header[2:end])
    length(names) == 119 || error("Expected 119 controls, found $(length(names))")
    length(unique(names)) == length(names) || error("Control names are not unique")

    sample_ids = Vector{Int}(undef, length(lines) - 1)
    values = Matrix{Float64}(undef, length(lines) - 1, length(names))
    for (row, line) in enumerate(lines[2:end])
        fields = split(line, ',')
        length(fields) == length(header) ||
            error("CSV row $row has $(length(fields)) fields; expected $(length(header))")
        sample_ids[row] = parse(Int, fields[1])
        for column in eachindex(names)
            values[row, column] = parse(Float64, fields[column + 1])
        end
    end
    return (; sample_ids, names, values)
end

detector_names(ring) = [
    uppercase(String(element.name))
    for element in ring.line
    if startswith(uppercase(String(element.name)), "DET_")
]

function prepare_batch_model(
    names::Vector{String},
    values::Matrix{Float64},
)
    n_samples, n_controls = size(values)
    n_controls == length(names) || error("Control matrix width does not match labels")

    model = load_cesr_model(zero_value=BatchParam(0.0), rf_on=true)
    for (column, name) in enumerate(names)
        model.controls[name] = BatchParam(view(values, :, column))
    end
    return model
end

function solve_and_track(
    model,
    n_samples::Int;
    initial_v0::Union{Nothing,AbstractMatrix}=nothing,
    reltol::Float64,
    abstol::Float64,
    maxiter::Int,
)
    detectors = detector_names(model.ring)
    length(detectors) == 99 || error("Expected 99 detectors, found $(length(detectors))")
    v0 = isnothing(initial_v0) ? zeros(n_samples, 6) : copy(initial_v0)
    size(v0) == (n_samples, 6) ||
        error("Initial closed-orbit guess must have size ($n_samples, 6)")
    solve_seconds = @elapsed begin
        solution = find_closed_orbit(
            model.ring;
            v0,
            coasting_beam=false,
            batch=Val{true}(),
            reltol,
            abstol,
            maxiter,
            warn=false,
        )
    end

    converged = Array(solution.sol.retcode .== SciBmad.BatchSolve.RETCODE_SUCCESS)
    iterations = vec(Array(solution.sol.iters))
    horizontal = similar(solution.v0, n_samples, length(detectors))
    vertical = similar(horizontal)
    closure_norms = zeros(n_samples)
    track_seconds = @elapsed begin
        bunch = Bunch(v=copy(solution.v0))
        SciBmad.BTBL.check_bl_bunch!(bunch, model.ring, false)
        detector_index = 0
        for element in model.ring.line
            track!(bunch, element)
            name = uppercase(String(element.name))
            startswith(name, "DET_") || continue
            detector_index += 1
            name == detectors[detector_index] ||
                error("Detector order changed at $name")
            horizontal[:, detector_index] .= bunch.coords.v[:, 1]
            vertical[:, detector_index] .= bunch.coords.v[:, 3]
        end
        detector_index == length(detectors) || error("Not all detectors were tracked")
        closure_norms .= sqrt.(
            vec(sum(abs2, Array(bunch.coords.v) .- solution.v0; dims=2)),
        )
    end

    return (;
        observables=Array(hcat(horizontal, vertical)),
        detectors,
        converged,
        iterations,
        solution,
        final_v0=copy(solution.v0),
        solve_seconds,
        track_seconds,
        factorization_seconds=0.0,
        closure_norms,
        fallback_count=0,
        fallback_success_count=0,
        fallback_seconds=0.0,
        fallback_iterations_max=0,
    )
end

function frozen_solve_and_track(
    model,
    n_samples::Int,
    frozen_jacobian::AbstractMatrix;
    initial_v0::AbstractMatrix,
    reltol::Float64,
    abstol::Float64,
    maxiter::Int,
)
    detectors = detector_names(model.ring)
    length(detectors) == 99 || error("Expected 99 detectors, found $(length(detectors))")
    size(initial_v0) == (n_samples, 6) ||
        error("Initial closed-orbit guess must have size ($n_samples, 6)")
    size(frozen_jacobian) == (6, 6) ||
        error("Frozen closed-orbit Jacobian must have size (6, 6)")

    v = copy(initial_v0)
    residual = similar(v)
    v_cache = similar(v)
    rhs = zeros(eltype(v), 6, n_samples)
    step = similar(v)
    active = trues(n_samples)
    converged = falses(n_samples)
    iterations = fill(maxiter, n_samples)
    set_kernel! = SciBmad.set_v!(SciBmad.KA.get_backend(v))
    sub_kernel! = SciBmad.sub_v!(SciBmad.KA.get_backend(v))
    factorization_seconds = @elapsed factorization = lu(Matrix(frozen_jacobian))

    solve_seconds = @elapsed begin
        for iteration in 1:maxiter
            SciBmad._co_res!(
                residual,
                v,
                model.ring,
                set_kernel!,
                sub_kernel!,
                v_cache,
            )

            for sample in 1:n_samples
                active[sample] || continue
                if sum(abs2, view(residual, sample, :)) < abstol^2
                    active[sample] = false
                    converged[sample] = true
                    iterations[sample] = iteration - 1
                end
            end
            all(.!active) && break

            fill!(rhs, 0)
            for sample in 1:n_samples
                active[sample] || continue
                rhs[:, sample] .= -view(residual, sample, :)
            end
            ldiv!(factorization, rhs)
            step .= transpose(rhs)

            for sample in 1:n_samples
                active[sample] || continue
                sample_step = view(step, sample, :)
                if any(value -> !isfinite(value), sample_step)
                    active[sample] = false
                    iterations[sample] = iteration - 1
                    continue
                end
                view(v, sample, :) .+= sample_step
                if sum(abs2, sample_step) <
                   reltol^2 * sum(abs2, view(v, sample, :))
                    active[sample] = false
                    converged[sample] = true
                    iterations[sample] = iteration
                end
            end
            all(.!active) && break
        end
    end

    horizontal = similar(v, n_samples, length(detectors))
    vertical = similar(horizontal)
    closure_norms = zeros(n_samples)
    track_seconds = @elapsed begin
        bunch = Bunch(v=copy(v))
        SciBmad.BTBL.check_bl_bunch!(bunch, model.ring, false)
        detector_index = 0
        for element in model.ring.line
            track!(bunch, element)
            name = uppercase(String(element.name))
            startswith(name, "DET_") || continue
            detector_index += 1
            name == detectors[detector_index] ||
                error("Detector order changed at $name")
            horizontal[:, detector_index] .= bunch.coords.v[:, 1]
            vertical[:, detector_index] .= bunch.coords.v[:, 3]
        end
        detector_index == length(detectors) || error("Not all detectors were tracked")
        closure_norms .= sqrt.(vec(sum(abs2, Array(bunch.coords.v) .- v; dims=2)))
    end

    return (;
        observables=Array(hcat(horizontal, vertical)),
        detectors,
        converged,
        iterations,
        solution=nothing,
        final_v0=copy(v),
        solve_seconds,
        track_seconds,
        factorization_seconds,
        closure_norms,
        fallback_count=0,
        fallback_success_count=0,
        fallback_seconds=0.0,
        fallback_iterations_max=0,
    )
end

function apply_full_newton_fallback(
    result,
    names::Vector{String},
    values::Matrix{Float64};
    reltol::Float64,
    abstol::Float64,
    maxiter::Int,
)
    fallback_indices = findall(
        .!result.converged .|
        .!isfinite.(result.closure_norms) .|
        (result.closure_norms .> abstol),
    )
    isempty(fallback_indices) && return result

    # BatchParam requires at least two entries. Duplicate a single failed lane,
    # then keep only the first result after the full-AD solve.
    solve_indices = length(fallback_indices) == 1 ?
        [fallback_indices[1], fallback_indices[1]] : fallback_indices
    fallback_result = nothing
    fallback_seconds = @elapsed begin
        fallback_model = prepare_batch_model(
            names,
            Matrix(values[solve_indices, :]),
        )
        fallback_result = solve_and_track(
            fallback_model,
            length(solve_indices);
            initial_v0=Matrix(result.final_v0[solve_indices, :]),
            reltol,
            abstol,
            maxiter,
        )
    end

    observables = copy(result.observables)
    converged = copy(result.converged)
    iterations = copy(result.iterations)
    final_v0 = copy(result.final_v0)
    closure_norms = copy(result.closure_norms)
    success_count = 0
    for (fallback_row, original_row) in enumerate(fallback_indices)
        observables[original_row, :] .=
            view(fallback_result.observables, fallback_row, :)
        final_v0[original_row, :] .=
            view(fallback_result.final_v0, fallback_row, :)
        closure_norms[original_row] =
            fallback_result.closure_norms[fallback_row]
        good =
            fallback_result.converged[fallback_row] &&
            isfinite(closure_norms[original_row]) &&
            closure_norms[original_row] <= abstol
        converged[original_row] = good
        success_count += good
        iterations[original_row] +=
            fallback_result.iterations[fallback_row]
    end

    return merge(
        result,
        (;
            observables,
            converged,
            iterations,
            final_v0,
            closure_norms,
            fallback_count=length(fallback_indices),
            fallback_success_count=success_count,
            fallback_seconds,
            fallback_iterations_max=maximum(fallback_result.iterations),
        ),
    )
end

function prepare_initial_guess(
    names::Vector{String},
    n_samples::Int,
    mode::String;
    reltol::Float64,
    abstol::Float64,
    maxiter::Int,
)
    if mode == "zero"
        return (;
            v0=zeros(n_samples, 6),
            nominal_orbit=zeros(6),
            nominal_jacobian=zeros(6, 6),
            nominal_model_setup_seconds=0.0,
            nominal_solve_seconds=0.0,
            nominal_iterations=0,
        )
    end

    mode == "nominal-z0" || error("Unsupported initial-guess mode: $mode")
    setup_seconds = @elapsed nominal_model = load_cesr_model(
        zero_value=0.0,
        rf_on=true,
    )
    solve_seconds = @elapsed nominal_solution = find_closed_orbit(
        nominal_model.ring;
        v0=zeros(1, 6),
        coasting_beam=false,
        batch=Val{false}(),
        reltol,
        abstol,
        maxiter,
        warn=false,
    )
    converged =
        nominal_solution.sol.retcode == SciBmad.BatchSolve.RETCODE_SUCCESS
    converged || error("Nominal closed orbit did not converge")
    nominal_orbit = vec(copy(nominal_solution.v0))
    nominal_jacobian = Matrix(nominal_solution.sol.jac)
    nominal_iterations = Int(nominal_solution.sol.iters)
    return (;
        v0=repeat(reshape(nominal_orbit, 1, 6), n_samples, 1),
        nominal_orbit,
        nominal_jacobian,
        nominal_model_setup_seconds=setup_seconds,
        nominal_solve_seconds=solve_seconds,
        nominal_iterations,
    )
end

function simulate_batch(
    names::Vector{String},
    values::Matrix{Float64};
    initial_guess_mode::String,
    jacobian_mode::String,
    fallback_full_newton::Bool,
    reltol::Float64,
    abstol::Float64,
    maxiter::Int,
)
    guess = prepare_initial_guess(
        names,
        size(values, 1),
        initial_guess_mode;
        reltol,
        abstol,
        maxiter,
    )
    model_setup_seconds = @elapsed model = prepare_batch_model(names, values)
    result = if jacobian_mode == "full"
        solve_and_track(
            model,
            size(values, 1);
            initial_v0=guess.v0,
            reltol,
            abstol,
            maxiter,
        )
    elseif jacobian_mode == "frozen-nominal"
        frozen_solve_and_track(
            model,
            size(values, 1),
            guess.nominal_jacobian;
            initial_v0=guess.v0,
            reltol,
            abstol,
            maxiter,
        )
    else
        error("Unsupported Jacobian mode: $jacobian_mode")
    end
    if jacobian_mode == "frozen-nominal" && fallback_full_newton
        result = apply_full_newton_fallback(
            result,
            names,
            values;
            reltol,
            abstol,
            maxiter,
        )
    end
    return merge(result, guess, (; model_setup_seconds))
end

function write_outputs(path, sample_ids, result)
    mkpath(dirname(path))
    labels = vcat(
        [name * ":x" for name in result.detectors],
        [name * ":y" for name in result.detectors],
    )
    open(path, "w") do io
        println(io, join(vcat("sample_id", "converged", labels), ','))
        for row in eachindex(sample_ids)
            fields = Any[sample_ids[row], result.converged[row]]
            append!(fields, result.observables[row, :])
            println(io, join(fields, ','))
        end
    end
    return labels
end

function main(args=ARGS)
    options = parse_args(args)
    inputs = abspath(options["inputs"])
    output = abspath(options["output"])
    metadata_path = abspath(options["metadata"])
    samples = read_samples(inputs)
    reltol = parse(Float64, options["reltol"])
    abstol = parse(Float64, options["abstol"])
    maxiter = parse(Int, options["maxiter"])
    initial_guess_mode = options["initial-guess"]
    jacobian_mode = options["jacobian-mode"]
    fallback_full_newton =
        lowercase(options["fallback-full-newton"]) == "true"
    warmup_samples = min(parse(Int, options["warmup-samples"]), size(samples.values, 1))
    warmup_samples >= 1 || error("--warmup-samples must be positive")

    @printf(
        "SciBmad CESR batch benchmark: %d samples x %d controls\n",
        size(samples.values)...,
    )
    warmup_elapsed = @elapsed simulate_batch(
        samples.names,
        samples.values[1:warmup_samples, :];
        initial_guess_mode,
        jacobian_mode,
        fallback_full_newton,
        reltol,
        abstol,
        maxiter,
    )
    @printf("Warmup/compilation batch (%d samples): %.3f s\n", warmup_samples, warmup_elapsed)

    guess = prepare_initial_guess(
        samples.names,
        size(samples.values, 1),
        initial_guess_mode;
        reltol,
        abstol,
        maxiter,
    )
    if initial_guess_mode == "nominal-z0"
        @printf(
            "Nominal z0: setup %.3f s + solve %.3f s, %d Newton iterations\n",
            guess.nominal_model_setup_seconds,
            guess.nominal_solve_seconds,
            guess.nominal_iterations,
        )
        @printf(
            "Nominal z0 [x, px, y, py, z, pz] = [%s]\n",
            join((@sprintf("%.16e", value) for value in guess.nominal_orbit), ", "),
        )
    end

    model_timed = @timed prepare_batch_model(samples.names, samples.values)
    model = model_timed.value
    timed = @timed begin
        current_result = if jacobian_mode == "full"
            solve_and_track(
                model,
                size(samples.values, 1);
                initial_v0=guess.v0,
                reltol,
                abstol,
                maxiter,
            )
        else
            frozen_solve_and_track(
                model,
                size(samples.values, 1),
                guess.nominal_jacobian;
                initial_v0=guess.v0,
                reltol,
                abstol,
                maxiter,
            )
        end
        if jacobian_mode == "frozen-nominal" && fallback_full_newton
            current_result = apply_full_newton_fallback(
                current_result,
                samples.names,
                samples.values;
                reltol,
                abstol,
                maxiter,
            )
        end
        current_result
    end
    result = timed.value
    physics_seconds = timed.time
    @printf(
        "Model setup: %.3f s\nPhysics: %.3f s (solve %.3f + track %.3f), %.3f samples/s, converged %d/%d\n",
        model_timed.time,
        physics_seconds,
        result.solve_seconds,
        result.track_seconds,
        size(samples.values, 1) / physics_seconds,
        count(result.converged),
        length(result.converged),
    )
    @printf(
        "Newton iterations: min %d, median %.1f, mean %.3f, max %d\n",
        minimum(result.iterations),
        median(result.iterations),
        mean(result.iterations),
        maximum(result.iterations),
    )
    if jacobian_mode == "frozen-nominal"
        @printf(
            "Frozen Jacobian: factorization %.6f s, closure norm median %.3e, max %.3e\n",
            result.factorization_seconds,
            median(result.closure_norms),
            maximum(result.closure_norms),
        )
        @printf(
            "Full-AD fallback: %d attempted, %d succeeded, %.3f s\n",
            result.fallback_count,
            result.fallback_success_count,
            result.fallback_seconds,
        )
    end

    write_seconds = @elapsed labels = write_outputs(output, samples.sample_ids, result)
    mkpath(dirname(metadata_path))
    metadata = Dict(
        "format" => "cesr-dataset-benchmark-v1",
        "engine" => "SciBmad",
        "device" => "cpu",
        "mode" => options["mode"],
        "input_csv" => inputs,
        "output_csv" => output,
        "sample_count" => size(samples.values, 1),
        "control_count" => size(samples.values, 2),
        "observable_count" => length(labels),
        "detector_count" => length(result.detectors),
        "converged_count" => count(result.converged),
        "failed_count" => count(.!result.converged),
        "warmup_sample_count" => warmup_samples,
        "warmup_seconds" => warmup_elapsed,
        "initial_guess_mode" => initial_guess_mode,
        "jacobian_mode" => jacobian_mode,
        "fallback_full_newton_enabled" => fallback_full_newton,
        "nominal_model_setup_seconds" => guess.nominal_model_setup_seconds,
        "nominal_closed_orbit_seconds" => guess.nominal_solve_seconds,
        "nominal_closed_orbit_iterations" => guess.nominal_iterations,
        "nominal_closed_orbit" => guess.nominal_orbit,
        "nominal_jacobian_condition_number" => (
            jacobian_mode == "frozen-nominal" ?
            cond(guess.nominal_jacobian) : 0.0
        ),
        "model_setup_seconds" => model_timed.time,
        "model_setup_allocated_bytes" => model_timed.bytes,
        "physics_seconds" => physics_seconds,
        "closed_orbit_seconds" => result.solve_seconds,
        "newton_iterations_min" => minimum(result.iterations),
        "newton_iterations_median" => median(result.iterations),
        "newton_iterations_mean" => mean(result.iterations),
        "newton_iterations_max" => maximum(result.iterations),
        "detector_tracking_seconds" => result.track_seconds,
        "samples_per_second" => size(samples.values, 1) / physics_seconds,
        "write_seconds" => write_seconds,
        "allocated_bytes" => timed.bytes,
        "gc_seconds" => timed.gctime,
        "julia_version" => string(VERSION),
        "julia_threads" => Threads.nthreads(),
        "reltol" => reltol,
        "abstol" => abstol,
        "maxiter" => maxiter,
        "execution_model" => (
            jacobian_mode == "full" ?
            "one BatchParam array per control; full AD Jacobian each Newton iteration" :
            "one BatchParam array per control; one nominal 6x6 Jacobian reused for all samples and iterations"
        ),
        "timed_region" => "closed-orbit solve + detector tracking",
    )
    if jacobian_mode == "frozen-nominal"
        metadata["frozen_jacobian_factorization_seconds"] =
            result.factorization_seconds
        metadata["final_closure_norm_median"] = median(result.closure_norms)
        metadata["final_closure_norm_max"] = maximum(result.closure_norms)
        metadata["fallback_count"] = result.fallback_count
        metadata["fallback_success_count"] = result.fallback_success_count
        metadata["fallback_seconds"] = result.fallback_seconds
        metadata["fallback_iterations_max"] =
            result.fallback_iterations_max
    end
    open(metadata_path, "w") do io
        TOML.print(io, metadata; sorted=true)
    end
    println("Output:   $output")
    println("Metadata: $metadata_path")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end

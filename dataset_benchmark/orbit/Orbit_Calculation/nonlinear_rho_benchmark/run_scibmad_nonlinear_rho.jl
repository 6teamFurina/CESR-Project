#!/usr/bin/env julia

const HERE_NL = @__DIR__
const CALCULATION_DIR_NL = normpath(joinpath(HERE_NL, ".."))
include(joinpath(CALCULATION_DIR_NL, "benchmark_scibmad.jl"))

using Dates

const INPUT_PATH_NL = joinpath(HERE_NL, "shared_input", "nonlinear_rho_correctors.csv")
const MANIFEST_PATH_NL = joinpath(HERE_NL, "shared_input", "sample_manifest.csv")
const RESULT_DIR_NL = joinpath(HERE_NL, "results", "scibmad")
const RESPONSE_PATH_NL = joinpath(ORBIT_ROOT, "reference", "closed_orbit_response_6x119.csv")
const RELTOL_NL = 1.0e-8
const ABSTOL_NL = 1.0e-10
const MAXITER_NL = 100

function read_manifest(path)
    lines = readlines(path)
    sample_id = Int[]
    scenario = String[]
    rho = Float64[]
    trial_id = Int[]
    for line in lines[2:end]
        fields = split(line, ',')
        push!(sample_id, parse(Int, fields[1]))
        push!(scenario, fields[2])
        push!(rho, parse(Float64, fields[3]))
        push!(trial_id, parse(Int, fields[4]))
    end
    return (; sample_id, scenario, rho, trial_id)
end

function ordered_groups(manifest)
    groups = Tuple{String,Float64}[]
    for row in eachindex(manifest.sample_id)
        manifest.scenario[row] == "baseline" && continue
        key = (manifest.scenario[row], manifest.rho[row])
        key in groups || push!(groups, key)
    end
    return groups
end

function write_outputs_nl(path, samples, observables, converged, detectors)
    labels = vcat(["$name:x" for name in detectors], ["$name:y" for name in detectors])
    open(path, "w") do io
        println(io, join(vcat("sample_id", "converged", labels), ','))
        for row in eachindex(samples.sample_ids)
            fields = Any[samples.sample_ids[row], converged[row]]
            append!(fields, observables[row, :])
            println(
                io,
                join(
                    (value isa AbstractFloat ? @sprintf("%.17g", value) : string(value) for value in fields),
                    ',',
                ),
            )
        end
    end
end

function main_nl()
    samples = read_samples(INPUT_PATH_NL)
    manifest = read_manifest(MANIFEST_PATH_NL)
    samples.sample_ids == manifest.sample_id || error("Input and manifest sample IDs differ")
    mkpath(RESULT_DIR_NL)

    warmup_sample_count = 600
    warmup_seconds = @elapsed simulate_batch(
        samples.names,
        samples.values[2:(1 + warmup_sample_count), :];
        initial_guess_mode="response-linear",
        jacobian_mode="frozen-nominal",
        response_matrix_cache=RESPONSE_PATH_NL,
        recompute_response=false,
        reltol=RELTOL_NL,
        abstol=ABSTOL_NL,
        maxiter=MAXITER_NL,
    )

    guess_timed = @timed prepare_initial_guess(
        samples.names,
        samples.values,
        "response-linear";
        response_matrix_cache=RESPONSE_PATH_NL,
        recompute_response=false,
        reltol=RELTOL_NL,
        abstol=ABSTOL_NL,
        maxiter=MAXITER_NL,
    )
    guess = guess_timed.value
    n_samples = length(samples.sample_ids)
    observables = Matrix{Float64}(undef, n_samples, 198)
    converged = falses(n_samples)
    iterations = fill(MAXITER_NL, n_samples)
    closure_norms = fill(Inf, n_samples)
    detectors = String[]
    group_rows = NamedTuple[]

    baseline_result = simulate_batch(
        samples.names,
        vcat(samples.values[1:1, :], samples.values[1:1, :]);
        initial_guess_mode="response-linear",
        jacobian_mode="frozen-nominal",
        response_matrix_cache=RESPONSE_PATH_NL,
        recompute_response=false,
        reltol=RELTOL_NL,
        abstol=ABSTOL_NL,
        maxiter=MAXITER_NL,
    )
    observables[1, :] .= baseline_result.observables[1, :]
    converged[1] = baseline_result.converged[1]
    iterations[1] = baseline_result.iterations[1]
    closure_norms[1] = baseline_result.closure_norms[1]
    append!(detectors, baseline_result.detectors)

    for (scenario, rho) in ordered_groups(manifest)
        selected = findall(
            row -> manifest.scenario[row] == scenario && manifest.rho[row] == rho,
            eachindex(manifest.sample_id),
        )
        run_indices = selected
        values = Matrix(samples.values[run_indices, :])
        model_timed = @timed prepare_batch_model(samples.names, values)
        exact_timed = @timed begin
            current = frozen_solve_and_track(
                model_timed.value,
                length(run_indices),
                guess.nominal_jacobian;
                initial_v0=Matrix(guess.v0[run_indices, :]),
                reltol=RELTOL_NL,
                abstol=ABSTOL_NL,
                maxiter=MAXITER_NL,
            )
            apply_full_newton_fallback(
                current,
                samples.names,
                values;
                reltol=RELTOL_NL,
                abstol=ABSTOL_NL,
                maxiter=MAXITER_NL,
            )
        end
        result = exact_timed.value
        observables[run_indices, :] .= result.observables
        converged[run_indices] .= result.converged
        iterations[run_indices] .= result.iterations
        closure_norms[run_indices] .= result.closure_norms
        selected_local = eachindex(run_indices)
        push!(group_rows, (;
            scenario,
            rho,
            samples=length(selected),
            converged=count(result.converged[selected_local]),
            model_setup_seconds=model_timed.time,
            physics_seconds=exact_timed.time,
            solve_seconds=result.solve_seconds,
            track_seconds=result.track_seconds,
            samples_per_physics_second=length(selected) / exact_timed.time,
            samples_per_setup_plus_physics_second=length(selected) / (model_timed.time + exact_timed.time),
            mean_iterations=mean(result.iterations[selected_local]),
            max_iterations=maximum(result.iterations[selected_local]),
            max_closure_norm=maximum(result.closure_norms[selected_local]),
            fallback_count=result.fallback_count,
            fallback_success_count=result.fallback_success_count,
            allocated_bytes=exact_timed.bytes,
        ))
        @printf(
            "SciBmad %-10s rho=%4.2f: %d/%d, %.3f s, %.2f samples/s, mean iter %.2f\n",
            scenario, rho, group_rows[end].converged, length(selected), exact_timed.time,
            group_rows[end].samples_per_physics_second, group_rows[end].mean_iterations,
        )
    end

    output_path = joinpath(RESULT_DIR_NL, "scibmad_samples.csv")
    diagnostics_path = joinpath(RESULT_DIR_NL, "scibmad_sample_diagnostics.csv")
    timing_path = joinpath(RESULT_DIR_NL, "scibmad_group_timings.csv")
    write_outputs_nl(output_path, samples, observables, converged, detectors)
    open(diagnostics_path, "w") do io
        println(io, "sample_id,scenario,rho,trial_id,converged,newton_iterations,closure_norm")
        for row in eachindex(samples.sample_ids)
            @printf(
                io, "%d,%s,%.17g,%d,%s,%d,%.17g\n",
                samples.sample_ids[row], manifest.scenario[row], manifest.rho[row],
                manifest.trial_id[row], converged[row], iterations[row], closure_norms[row],
            )
        end
    end
    open(timing_path, "w") do io
        println(io, join(string.(keys(first(group_rows))), ','))
        for row in group_rows
            println(io, join(string.(values(row)), ','))
        end
    end
    metadata = Dict(
        "format" => "cesr-nonlinear-rho-scibmad-v1",
        "date" => string(Dates.today()),
        "engine" => "SciBmad",
        "input_csv" => INPUT_PATH_NL,
        "sample_count" => n_samples,
        "converged_count" => count(converged),
        "failed_count" => count(.!converged),
        "reltol" => RELTOL_NL,
        "abstol" => ABSTOL_NL,
        "maxiter" => MAXITER_NL,
        "initial_guess" => "nominal first-order closed-orbit response",
        "jacobian" => "frozen nominal 6x6 Jacobian with shared LU per group",
        "full_ad_fallback" => true,
        "warmup_seconds" => warmup_seconds,
        "warmup_sample_count" => warmup_sample_count,
        "shared_initial_guess_setup_seconds" => guess_timed.time,
        "group_physics_seconds_sum" => sum(row.physics_seconds for row in group_rows),
        "group_model_setup_seconds_sum" => sum(row.model_setup_seconds for row in group_rows),
        "output_csv" => output_path,
        "diagnostics_csv" => diagnostics_path,
        "group_timings_csv" => timing_path,
        "julia_version" => string(VERSION),
        "julia_threads" => Threads.nthreads(),
        "blas_threads" => BLAS.get_num_threads(),
    )
    open(joinpath(RESULT_DIR_NL, "scibmad_metadata.toml"), "w") do io
        TOML.print(io, metadata; sorted=true)
    end
end

main_nl()

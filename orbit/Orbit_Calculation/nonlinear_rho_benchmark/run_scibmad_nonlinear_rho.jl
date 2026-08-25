#!/usr/bin/env julia

const HERE_NL = @__DIR__
const CALCULATION_DIR_NL = normpath(joinpath(HERE_NL, ".."))
include(joinpath(CALCULATION_DIR_NL, "benchmark_scibmad.jl"))

using Dates

function selected_nl_ring()
    ring = :latest
    for argument in ARGS
        startswith(argument, "--ring=") || continue
        ring = Symbol(lowercase(split(argument, "="; limit=2)[2]))
    end
    return canonical_ring_id(ring)
end

function selected_cpu_multithreading()
    for argument in ARGS
        startswith(argument, "--cpu-multithreading=") || continue
        value = lowercase(split(argument, "="; limit=2)[2])
        value in ("true", "false") ||
            error("--cpu-multithreading must be true or false")
        return value == "true"
    end
    return false
end

function selected_reuse_batch_model()
    for argument in ARGS
        startswith(argument, "--reuse-batch-model=") || continue
        value = lowercase(split(argument, "="; limit=2)[2])
        value in ("true", "false") ||
            error("--reuse-batch-model must be true or false")
        return value == "true"
    end
    return false
end

const ARTIFACT_RING_NL = selected_nl_ring()
const CPU_MULTITHREADING_NL = selected_cpu_multithreading()
const REUSE_BATCH_MODEL_NL = selected_reuse_batch_model()
const RING_NL = ARTIFACT_RING_NL == :latest_cesr ? :latest : :legacy
const INPUT_PATH_NL = joinpath(
    HERE_NL,
    "shared_input",
    String(ARTIFACT_RING_NL),
    "nonlinear_rho_correctors.csv",
)
const MANIFEST_PATH_NL = joinpath(
    HERE_NL,
    "shared_input",
    String(ARTIFACT_RING_NL),
    "sample_manifest.csv",
)
const RESULT_VARIANT_NL = string(
    CPU_MULTITHREADING_NL ? "scibmad_threads$(Threads.nthreads())" : "scibmad",
    REUSE_BATCH_MODEL_NL ? "_reuse" : "",
)
const RESULT_DIR_NL = joinpath(
    HERE_NL,
    "results",
    String(ARTIFACT_RING_NL),
    RESULT_VARIANT_NL,
)
const RESPONSE_DIR_NL = joinpath(
    ORBIT_ROOT,
    "reference",
    String(ARTIFACT_RING_NL),
    ARTIFACT_RING_NL == :latest_cesr ? "gtpsa" : "",
)
const RESPONSE_PATH_NL = joinpath(
    RESPONSE_DIR_NL,
    "closed_orbit_response.csv",
)
const RESPONSE_METHOD_NL = ARTIFACT_RING_NL == :latest_cesr ? "gtpsa" : "central-difference"
const MODEL_FACTORY_NL = (; kwargs...) -> load_ring_model(; ring=RING_NL, kwargs...)
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
    CPU_MULTITHREADING_NL && BLAS.set_num_threads(1)
    samples = read_samples(INPUT_PATH_NL)
    manifest = read_manifest(MANIFEST_PATH_NL)
    samples.sample_ids == manifest.sample_id || error("Input and manifest sample IDs differ")
    mkpath(RESULT_DIR_NL)

    warmup_sample_count = 600
    warmup_sample_count = min(warmup_sample_count, max(1, size(samples.values, 1) - 1))
    warmup_seconds = @elapsed simulate_batch(
        samples.names,
        samples.values[2:(1 + warmup_sample_count), :];
        initial_guess_mode="response-linear",
        jacobian_mode="frozen-nominal",
        response_matrix_cache=RESPONSE_PATH_NL,
        recompute_response=false,
        response_method=RESPONSE_METHOD_NL,
        reltol=RELTOL_NL,
        abstol=ABSTOL_NL,
        maxiter=MAXITER_NL,
        model_factory=MODEL_FACTORY_NL,
        use_cpu_multithreading=CPU_MULTITHREADING_NL,
    )

    guess_timed = @timed prepare_initial_guess(
        samples.names,
        samples.values,
        "response-linear";
        response_matrix_cache=RESPONSE_PATH_NL,
        recompute_response=false,
        response_method=RESPONSE_METHOD_NL,
        reltol=RELTOL_NL,
        abstol=ABSTOL_NL,
        maxiter=MAXITER_NL,
        model_factory=MODEL_FACTORY_NL,
    )
    guess = guess_timed.value
    n_samples = length(samples.sample_ids)
    group_rows = NamedTuple[]
    groups = ordered_groups(manifest)

    baseline_result = simulate_batch(
        samples.names,
        vcat(samples.values[1:1, :], samples.values[1:1, :]);
        initial_guess_mode="response-linear",
        jacobian_mode="frozen-nominal",
        response_matrix_cache=RESPONSE_PATH_NL,
        recompute_response=false,
        response_method=RESPONSE_METHOD_NL,
        reltol=RELTOL_NL,
        abstol=ABSTOL_NL,
        maxiter=MAXITER_NL,
        model_factory=MODEL_FACTORY_NL,
        use_cpu_multithreading=CPU_MULTITHREADING_NL,
    )
    observables = Matrix{Float64}(
        undef,
        n_samples,
        size(baseline_result.observables, 2),
    )
    converged = falses(n_samples)
    iterations = fill(MAXITER_NL, n_samples)
    closure_norms = fill(Inf, n_samples)
    detectors = copy(baseline_result.detectors)
    observables[1, :] .= baseline_result.observables[1, :]
    converged[1] = baseline_result.converged[1]
    iterations[1] = baseline_result.iterations[1]
    closure_norms[1] = baseline_result.closure_norms[1]

    shared_model = nothing
    shared_model_setup_seconds = 0.0
    if REUSE_BATCH_MODEL_NL
        first_scenario, first_rho = first(groups)
        first_indices = findall(
            row -> manifest.scenario[row] == first_scenario &&
                manifest.rho[row] == first_rho,
            eachindex(manifest.sample_id),
        )
        shared_model_timed = @timed prepare_batch_model(
            samples.names,
            Matrix(samples.values[first_indices, :]);
            model_factory=MODEL_FACTORY_NL,
        )
        shared_model = shared_model_timed.value
        shared_model_setup_seconds = shared_model_timed.time
    end

    for (group_index, (scenario, rho)) in enumerate(groups)
        selected = findall(
            row -> manifest.scenario[row] == scenario && manifest.rho[row] == rho,
            eachindex(manifest.sample_id),
        )
        run_indices = selected
        values = Matrix(samples.values[run_indices, :])
        model_setup_seconds = 0.0
        control_update_seconds = 0.0
        model = if REUSE_BATCH_MODEL_NL
            if group_index > 1
                control_update_seconds = @elapsed update_batch_model_controls!(
                    shared_model,
                    samples.names,
                    values,
                )
            else
                model_setup_seconds = shared_model_setup_seconds
            end
            shared_model
        else
            model_timed = @timed prepare_batch_model(
                samples.names,
                values;
                model_factory=MODEL_FACTORY_NL,
            )
            model_setup_seconds = model_timed.time
            model_timed.value
        end
        exact_timed = @timed begin
            current = frozen_solve_and_track(
                model,
                length(run_indices),
                guess.nominal_jacobian;
                initial_v0=Matrix(guess.v0[run_indices, :]),
                reltol=RELTOL_NL,
                abstol=ABSTOL_NL,
                maxiter=MAXITER_NL,
                use_cpu_multithreading=CPU_MULTITHREADING_NL,
            )
            apply_full_newton_fallback(
                current,
                samples.names,
                values;
                reltol=RELTOL_NL,
                abstol=ABSTOL_NL,
                maxiter=MAXITER_NL,
                model_factory=MODEL_FACTORY_NL,
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
            model_setup_seconds,
            control_update_seconds,
            physics_seconds=exact_timed.time,
            solve_seconds=result.solve_seconds,
            track_seconds=result.track_seconds,
            samples_per_physics_second=length(selected) / exact_timed.time,
            samples_per_setup_plus_physics_second=length(selected) /
                (model_setup_seconds + control_update_seconds + exact_timed.time),
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
        "format" => "ring-nonlinear-rho-scibmad-v2",
        "ring" => String(ARTIFACT_RING_NL),
        "date" => string(Dates.today()),
        "engine" => "SciBmad",
        "input_csv" => INPUT_PATH_NL,
        "lattice_path" => String(
            ARTIFACT_RING_NL == :latest_cesr ?
            ORBIT_ADAPTER_LATEST_LATTICE : ORBIT_ADAPTER_LEGACY_MODEL,
        ),
        "response_method" => RESPONSE_METHOD_NL,
        "response_matrix_cache" => RESPONSE_PATH_NL,
        "sample_count" => n_samples,
        "converged_count" => count(converged),
        "failed_count" => count(.!converged),
        "reltol" => RELTOL_NL,
        "abstol" => ABSTOL_NL,
        "maxiter" => MAXITER_NL,
        "initial_guess" => "nominal first-order closed-orbit response",
        "jacobian" => "frozen nominal phase-space Jacobian with shared LU per group",
        "full_ad_fallback" => true,
        "cpu_multithreaded_tracking" => CPU_MULTITHREADING_NL,
        "batch_model_reused" => REUSE_BATCH_MODEL_NL,
        "warmup_seconds" => warmup_seconds,
        "warmup_sample_count" => warmup_sample_count,
        "shared_initial_guess_setup_seconds" => guess_timed.time,
        "group_physics_seconds_sum" => sum(row.physics_seconds for row in group_rows),
        "group_model_setup_seconds_sum" => sum(row.model_setup_seconds for row in group_rows),
        "group_control_update_seconds_sum" =>
            sum(row.control_update_seconds for row in group_rows),
        "setup_plus_physics_seconds" => sum(
            row.model_setup_seconds + row.control_update_seconds + row.physics_seconds
            for row in group_rows
        ),
        "all_runtime_setup_plus_physics_seconds" => guess_timed.time + sum(
            row.model_setup_seconds + row.control_update_seconds + row.physics_seconds
            for row in group_rows
        ),
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

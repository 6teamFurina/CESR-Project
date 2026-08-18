#!/usr/bin/env julia

"""Run exact RF-off chromatic optics on the 9,001 nonlinear-rho inputs.

The calculation is checkpointed by the original (scenario, rho) cells.  A
failure in the coasting closed-orbit or Twiss stage is recorded per sample and
does not discard already completed cells.
"""

include(joinpath(@__DIR__, "..", "benchmark_scibmad_chromatic_optics.jl"))
using Dates

const INPUT_PATH = normpath(joinpath(
    @__DIR__, "..", "..", "orbit", "Orbit_Calculation",
    "nonlinear_rho_benchmark", "shared_input", "nonlinear_rho_correctors.csv",
))
const MANIFEST_PATH = normpath(joinpath(
    @__DIR__, "..", "..", "orbit", "Orbit_Calculation",
    "nonlinear_rho_benchmark", "shared_input", "sample_manifest.csv",
))
const RESPONSE_PATH = normpath(joinpath(
    @__DIR__, "..", "..", "orbit", "reference", "closed_orbit_response_6x119.csv",
))
const DEFAULT_OUTPUT_DIR = joinpath(@__DIR__, "results", "scibmad")

function parse_runner_args(args)
    options = Dict(
        "output-dir" => DEFAULT_OUTPUT_DIR,
        "max-groups" => "16",
        "reltol" => "1e-8",
        "abstol" => "1e-10",
        "maxiter" => "100",
    )
    for argument in args
        startswith(argument, "--") || error("Expected --name=value: $argument")
        fields = split(argument[3:end], '='; limit=2)
        length(fields) == 2 || error("Missing value in $argument")
        haskey(options, fields[1]) || error("Unknown option --$(fields[1])")
        options[fields[1]] = fields[2]
    end
    return options
end

function read_nonlinear_manifest(path)
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

function manifest_groups(manifest)
    groups = NamedTuple[]
    baseline = findall(==("baseline"), manifest.scenario)
    push!(groups, (; scenario="baseline", rho=0.0, rows=baseline))
    keys_seen = Tuple{String,Float64}[]
    for row in eachindex(manifest.sample_id)
        manifest.scenario[row] == "baseline" && continue
        key = (manifest.scenario[row], manifest.rho[row])
        key in keys_seen && continue
        push!(keys_seen, key)
        rows = findall(
            index -> manifest.scenario[index] == key[1] && manifest.rho[index] == key[2],
            eachindex(manifest.sample_id),
        )
        push!(groups, (; scenario=key[1], rho=key[2], rows))
    end
    return groups
end

rho_label(rho) = replace(@sprintf("%.2f", rho), "." => "p")
function group_label(index, scenario, rho)
    scenario == "baseline" && return @sprintf("%02d_baseline", index - 1)
    return @sprintf("%02d_%s_rho_%s", index - 1, scenario, rho_label(rho))
end

sanitize_error(exception) = replace(
    sprint(showerror, exception), '\n' => ' ', '\r' => ' ', ',' => ';',
)

function solve_rf_on_group(names, values, initial_v0, nominal_jacobian; reltol, abstol, maxiter)
    duplicated = size(values, 1) == 1
    run_values = duplicated ? vcat(values, values) : values
    run_v0 = duplicated ? vcat(initial_v0, initial_v0) : initial_v0
    model = prepare_batch_model(names, run_values)
    result = frozen_solve_and_track(
        model,
        size(run_values, 1),
        nominal_jacobian;
        initial_v0=run_v0,
        reltol,
        abstol,
        maxiter,
    )
    result = apply_full_newton_fallback(
        result,
        names,
        run_values;
        reltol,
        abstol,
        maxiter,
    )
    return result
end

function write_status_header(path)
    open(path, "w") do io
        println(
            io,
            "sample_id,scenario,rho,trial_id,rf_on_converged,rf_on_iterations," *
            "rf_on_closure_norm,coasting_converged,coasting_closure_norm," *
            "twiss_converged,coasting_seconds,twiss_seconds,error_stage,error_message",
        )
    end
end

function main(args=ARGS)
    options = parse_runner_args(args)
    output_dir = abspath(options["output-dir"])
    max_groups = parse(Int, options["max-groups"])
    reltol = parse(Float64, options["reltol"])
    abstol = parse(Float64, options["abstol"])
    maxiter = parse(Int, options["maxiter"])
    mkpath(output_dir)

    samples = read_samples(INPUT_PATH)
    manifest = read_nonlinear_manifest(MANIFEST_PATH)
    samples.sample_ids == manifest.sample_id || error("Input and manifest sample IDs differ")
    groups = manifest_groups(manifest)[1:min(max_groups, 16)]

    guess_timed = @timed prepare_initial_guess(
        samples.names,
        samples.values,
        "response-linear";
        response_matrix_cache=RESPONSE_PATH,
        recompute_response=false,
        reltol,
        abstol,
        maxiter,
    )
    guess = guess_timed.value
    optics_setup_timed = @timed prepare_chromatic_optics_model()
    optics_setup = optics_setup_timed.value
    reusable_workspace = nothing

    status_path = joinpath(output_dir, "scibmad_sample_status.csv")
    group_path = joinpath(output_dir, "scibmad_group_summary.csv")
    write_status_header(status_path)
    open(group_path, "w") do io
        println(
            io,
            "group_index,scenario,rho,samples,rf_on_converged,coasting_converged," *
            "twiss_converged,rf_on_seconds,coasting_seconds,twiss_seconds," *
            "maximum_rf_on_closure_norm,maximum_coasting_closure_norm",
        )
    end

    total_rf_on = 0
    total_coasting = 0
    total_twiss = 0
    total_rf_seconds = 0.0
    total_coasting_seconds = 0.0
    total_twiss_seconds = 0.0

    for (group_index, group) in enumerate(groups)
        rows = group.rows
        values = Matrix(samples.values[rows, :])
        initial_v0 = Matrix(guess.v0[rows, :])
        rf_timed = @timed solve_rf_on_group(
            samples.names,
            values,
            initial_v0,
            guess.nominal_jacobian;
            reltol,
            abstol,
            maxiter,
        )
        rf = rf_timed.value
        group_count = length(rows)
        rf_converged = rf.converged[1:group_count]
        rf_iterations = rf.iterations[1:group_count]
        rf_closure = rf.closure_norms[1:group_count]
        rf_orbits = Matrix(rf.final_v0[1:group_count, :])

        successful_ids = Int[]
        successful_optics = Any[]
        successful_orbits = Vector{Vector{Float64}}()
        successful_residuals = Float64[]
        successful_seconds = Float64[]
        status_rows = String[]
        coasting_group_seconds = 0.0
        twiss_group_seconds = 0.0
        coasting_count = 0
        twiss_count = 0

        for local_row in 1:group_count
            global_row = rows[local_row]
            sample_id = samples.sample_ids[global_row]
            coasting_ok = false
            twiss_ok = false
            coasting_residual = Inf
            coasting_seconds = 0.0
            twiss_seconds = 0.0
            error_stage = ""
            error_message = ""
            coasting_orbit = zeros(6)

            if rf_converged[local_row]
                try
                    coasting = solve_coasting_closed_orbits!(
                        optics_setup.model,
                        samples.names,
                        values[local_row:local_row, :],
                        rf_orbits[local_row:local_row, :];
                        reltol,
                        abstol,
                        maxiter,
                    )
                    coasting_orbit .= coasting.orbits[1, :]
                    coasting_residual = coasting.residuals[1]
                    coasting_seconds = coasting.sample_seconds[1]
                    coasting_ok = isfinite(coasting_residual)
                    if isnothing(reusable_workspace)
                        reusable_workspace = prepare_reusable_twiss_workspace(
                            optics_setup.model.ring,
                            optics_setup.detectors,
                            optics_setup.descriptor,
                            coasting_orbit,
                        )
                    end
                    calculation = calculate_one_chromatic_optics(
                        optics_setup.model,
                        optics_setup.detectors,
                        samples.names,
                        view(values, local_row, :),
                        coasting_orbit,
                        optics_setup.descriptor,
                        reusable_workspace,
                    )
                    twiss_seconds = calculation.twiss_seconds
                    twiss_ok = true
                    push!(successful_ids, sample_id)
                    push!(successful_optics, calculation.optics)
                    push!(successful_orbits, copy(coasting_orbit))
                    push!(successful_residuals, coasting_residual)
                    push!(successful_seconds, calculation.control_update_seconds + calculation.twiss_seconds)
                catch exception
                    error_stage = coasting_ok ? "twiss" : "coasting"
                    error_message = sanitize_error(exception)
                    if error_stage == "twiss"
                        reusable_workspace = nothing
                    end
                end
            else
                error_stage = "rf_on"
                error_message = "RF-on closed orbit did not converge"
            end

            coasting_group_seconds += coasting_seconds
            twiss_group_seconds += twiss_seconds
            coasting_count += coasting_ok
            twiss_count += twiss_ok
            push!(
                status_rows,
                join((
                    sample_id,
                    manifest.scenario[global_row],
                    @sprintf("%.17g", manifest.rho[global_row]),
                    manifest.trial_id[global_row],
                    rf_converged[local_row],
                    rf_iterations[local_row],
                    @sprintf("%.17g", rf_closure[local_row]),
                    coasting_ok,
                    @sprintf("%.17g", coasting_residual),
                    twiss_ok,
                    @sprintf("%.17g", coasting_seconds),
                    @sprintf("%.17g", twiss_seconds),
                    error_stage,
                    error_message,
                ), ','),
            )
        end

        chunk_dir = joinpath(
            output_dir,
            "chunks",
            group_label(group_index, group.scenario, group.rho),
        )
        mkpath(chunk_dir)
        if !isempty(successful_ids)
            write_chromatic_detector_optics(
                joinpath(chunk_dir, "scibmad_detector_chromatic_twiss.csv"),
                successful_ids,
                successful_optics,
            )
            write_chromatic_ring_optics(
                joinpath(chunk_dir, "scibmad_ring_chromatic_twiss.csv"),
                successful_ids,
                successful_optics,
                successful_seconds,
            )
            orbit_matrix = reduce(vcat, permutedims.(successful_orbits))
            write_start_orbits(
                joinpath(chunk_dir, "scibmad_start_closed_orbits.csv"),
                successful_ids,
                orbit_matrix,
                successful_residuals,
            )
        end
        open(status_path, "a") do io
            for row in status_rows
                println(io, row)
            end
        end
        open(group_path, "a") do io
            @printf(
                io,
                "%d,%s,%.17g,%d,%d,%d,%d,%.17g,%.17g,%.17g,%.17g,%.17g\n",
                group_index - 1,
                group.scenario,
                group.rho,
                group_count,
                count(rf_converged),
                coasting_count,
                twiss_count,
                rf_timed.time,
                coasting_group_seconds,
                twiss_group_seconds,
                maximum(rf_closure),
                isempty(successful_residuals) ? Inf : maximum(successful_residuals),
            )
        end

        total_rf_on += count(rf_converged)
        total_coasting += coasting_count
        total_twiss += twiss_count
        total_rf_seconds += rf_timed.time
        total_coasting_seconds += coasting_group_seconds
        total_twiss_seconds += twiss_group_seconds
        @printf(
            "SciBmad optics %-10s rho=%4.2f: RF %d/%d, coasting %d/%d, Twiss %d/%d\n",
            group.scenario,
            group.rho,
            count(rf_converged),
            group_count,
            coasting_count,
            group_count,
            twiss_count,
            group_count,
        )
        flush(stdout)
    end

    metadata = Dict(
        "format" => "cesr-nonlinear-rho-optics-scibmad-v1",
        "date" => string(Dates.today()),
        "engine" => "SciBmad",
        "input_csv" => INPUT_PATH,
        "manifest_csv" => MANIFEST_PATH,
        "group_count" => length(groups),
        "sample_count" => sum(length(group.rows) for group in groups),
        "rf_on_converged_count" => total_rf_on,
        "coasting_converged_count" => total_coasting,
        "twiss_converged_count" => total_twiss,
        "shared_initial_guess_seconds" => guess_timed.time,
        "optics_model_setup_seconds" => optics_setup_timed.time,
        "rf_on_seconds" => total_rf_seconds,
        "coasting_seconds" => total_coasting_seconds,
        "twiss_seconds" => total_twiss_seconds,
        "twiss_mode" => "reusable exact pointwise workspace",
        "rf_mode" => "off (coasting) for optics",
        "reltol" => reltol,
        "abstol" => abstol,
        "maxiter" => maxiter,
        "julia_version" => string(VERSION),
        "julia_threads" => Threads.nthreads(),
        "status_csv" => status_path,
        "group_summary_csv" => group_path,
    )
    open(joinpath(output_dir, "scibmad_metadata.toml"), "w") do io
        TOML.print(io, metadata; sorted=true)
    end
    return total_twiss == sum(length(group.rows) for group in groups) ? 0 : 2
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end

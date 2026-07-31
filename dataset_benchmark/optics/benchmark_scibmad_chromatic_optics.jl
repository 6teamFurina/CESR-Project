#!/usr/bin/env julia

"""
Benchmark second-order GTPSA periodic optics at the 99 CESR DET_* elements.

`Descriptor(6, 2)` retains quadratic terms in the six phase-space variables.
Consequently, the first-order coefficient with index 6 gives the derivative
with respect to relative momentum deviation delta.  It is not a derivative
with respect to any of the 119 corrector settings.

The CESR model, detector list, and GTPSA descriptor are constructed once and
reused for all samples.  The validated RF-on batch closed orbits seed a short
RF-off/coasting solve. Those solved orbits are supplied to `twiss`, so `twiss`
does not solve them again.
"""

include(joinpath(@__DIR__, "optics_common.jl"))
include(joinpath(PROJECT_DIR, "scibmad_coasting_forwarddiff_patch.jl"))
include(joinpath(@__DIR__, "twiss_reuse.jl"))
using .SciBmadCoastingForwardDiffPatch

const DELTA_INDEX = 6

constant_term(value) = Float64(GTPSA.scalar(value))
function delta_derivative(value)
    try
        return Float64(value[DELTA_INDEX])
    catch exception
        exception isa BoundsError || rethrow()
        return 0.0
    end
end

function prepare_chromatic_optics_model()
    # With RF off, delta is a fixed energy parameter. With RF on it is a
    # synchrotron coordinate, and its linear coefficient is not chromaticity.
    model = load_cesr_model(zero_value=0.0, rf_on=false)
    detectors = detector_elements(model.ring)
    length(detectors) == 99 || error("Expected 99 detectors, found $(length(detectors))")
    descriptor = Descriptor(6, 2)
    return (; model, detectors, descriptor)
end

const CHROMATIC_TWISS_COLUMNS = (
    :phi_1,
    :beta_1,
    :alpha_1,
    :phi_2,
    :beta_2,
    :alpha_2,
    :phi_3,
    :gamma_c,
    :c11,
    :c12,
    :c21,
    :c22,
)

const CHROMATIC_ORBIT_COLUMNS = (
    :orbit_x,
    :orbit_px,
    :orbit_y,
    :orbit_py,
    :orbit_z,
    :orbit_pz,
)

function write_chromatic_detector_optics(path, sample_ids, optics_results)
    mkpath(dirname(path))
    open(path, "w") do io
        header = String["sample_id", "s", "beamline_index", "name"]
        for column in CHROMATIC_TWISS_COLUMNS
            push!(header, String(column))
            push!(header, "d$(column)_ddelta")
        end
        for column in CHROMATIC_ORBIT_COLUMNS
            push!(header, String(column))
            push!(header, "d$(column)_ddelta")
        end
        println(io, join(header, ','))

        for (sample_row, optics) in enumerate(optics_results)
            table = optics.table
            for detector_row in eachindex(table.name)
                fields = Any[
                    sample_ids[sample_row],
                    table.s[detector_row],
                    table.beamline_index[detector_row],
                    table.name[detector_row],
                ]
                for column in CHROMATIC_TWISS_COLUMNS
                    value = getproperty(table, column)[detector_row]
                    push!(fields, constant_term(value))
                    push!(fields, delta_derivative(value))
                end
                for column in CHROMATIC_ORBIT_COLUMNS
                    value = getproperty(table, column)[detector_row]
                    push!(fields, constant_term(value))
                    push!(fields, delta_derivative(value))
                end
                println(io, join(fields, ','))
            end
        end
    end
    return path
end

function write_chromatic_ring_optics(path, sample_ids, optics_results, sample_seconds)
    mkpath(dirname(path))
    open(path, "w") do io
        println(
            io,
            "sample_id,Q1_signed,Q2_signed,Qx_fractional,Qy_fractional," *
            "slip_tps_constant,xi_1,xi_2,slip_factor,twiss_seconds",
        )
        for row in eachindex(optics_results)
            tune_tps = optics_results[row].tunes
            tunes = constant_term.(tune_tps)
            chromaticities = delta_derivative.(tune_tps)
            comparison_tunes = (mod(tunes[1], 1), mod(tunes[2], 1))
            println(
                io,
                join(
                    (
                        sample_ids[row],
                        tunes[1],
                        tunes[2],
                        comparison_tunes...,
                        tunes[3],
                        chromaticities[1],
                        chromaticities[2],
                        chromaticities[3],
                        sample_seconds[row],
                    ),
                    ',',
                ),
            )
        end
    end
    return path
end

function solve_coasting_closed_orbits!(
    model,
    names,
    values,
    initial_guesses;
    reltol,
    abstol,
    maxiter,
)
    sample_count = size(values, 1)
    orbits = zeros(sample_count, 6)
    residuals = zeros(sample_count)
    sample_seconds = zeros(sample_count)
    for row in 1:sample_count
        apply_sample!(model, names, view(values, row, :))
        guess = copy(view(initial_guesses, row, :))
        guess[5:6] .= 0.0
        solution = nothing
        sample_seconds[row] = @elapsed begin
            solution = find_closed_orbit_coasting_forwarddiff(
                model.ring;
                v0=guess,
                coasting_beam=true,
                z=0.0,
                pz=0.0,
                reltol,
                abstol,
                maxiter,
            )
        end
        orbits[row, :] .= solution.orbit
        residuals[row] = solution.residual
    end
    return (; orbits, residuals, sample_seconds)
end

function calculate_one_chromatic_optics(
    model,
    detectors,
    names,
    values,
    closed_orbit,
    descriptor,
    reusable_workspace=nothing,
)
    control_update_seconds = @elapsed apply_sample!(model, names, values)
    orbit_matrix = reshape(copy(closed_orbit), 1, 6)
    twiss_seconds = @elapsed begin
        if isnothing(reusable_workspace)
            optics = twiss(
                model.ring;
                GTPSA_descriptor=descriptor,
                at=detectors,
                v0_and_coast=(orbit_matrix, true),
                spin=false,
                RDTs=false,
                normalizing_map=false,
            )
        else
            optics = twiss!(reusable_workspace, closed_orbit)
        end
    end
    length(optics.table.name) == 99 ||
        error("Twiss returned $(length(optics.table.name)) detector rows")
    return (; optics, control_update_seconds, twiss_seconds)
end

function main_chromatic_optics(args=ARGS)
    options = parse_optics_args(args)
    if !any(startswith(argument, "--output-dir=") for argument in args)
        options["output-dir"] = joinpath(OPTICS_DIR, "results", "chromatic_test_10")
    end
    inputs = abspath(options["inputs"])
    output_dir = abspath(options["output-dir"])
    response_matrix_cache = abspath(options["response-matrix-cache"])
    sample_count = parse(Int, options["sample-count"])
    reltol = parse(Float64, options["reltol"])
    abstol = parse(Float64, options["abstol"])
    maxiter = parse(Int, options["maxiter"])
    run_warmup = lowercase(options["warmup"]) == "true"
    twiss_mode = options["twiss-mode"]

    all_samples = read_samples(inputs)
    1 <= sample_count <= size(all_samples.values, 1) || error(
        "--sample-count must be between 1 and $(size(all_samples.values, 1))",
    )
    sample_ids = all_samples.sample_ids[1:sample_count]
    values = Matrix(all_samples.values[1:sample_count, :])

    @printf(
        "SciBmad Descriptor(6,2) chromatic-optics benchmark (%s Twiss): %d samples x %d controls x 99 detectors\n",
        twiss_mode,
        sample_count,
        length(all_samples.names),
    )

    closed_orbit_warmup_seconds = 0.0
    if run_warmup
        warmup_rows = min(2, sample_count)
        warmup_values = warmup_rows == 1 ? repeat(values[1:1, :], 2, 1) : values[1:warmup_rows, :]
        closed_orbit_warmup_seconds = @elapsed solve_input_closed_orbits(
            all_samples.names,
            Matrix(warmup_values);
            response_matrix_cache,
            reltol,
            abstol,
            maxiter,
        )
        @printf(
            "Closed-orbit warmup/compilation (%d lanes): %.3f s\n",
            size(warmup_values, 1),
            closed_orbit_warmup_seconds,
        )
    end

    closed_orbit_timed = @timed solve_input_closed_orbits(
        all_samples.names,
        values;
        response_matrix_cache,
        reltol,
        abstol,
        maxiter,
    )
    closed_orbit_data = closed_orbit_timed.value
    closed_orbits = Matrix(closed_orbit_data.result.final_v0)
    @printf(
        "Closed orbit: %.3f s wall (Newton %.3f s), closure max %.3e, fallback %d\n",
        closed_orbit_timed.time,
        closed_orbit_data.result.solve_seconds,
        maximum(closed_orbit_data.result.closure_norms),
        closed_orbit_data.result.fallback_count,
    )

    optics_setup_timed = @timed prepare_chromatic_optics_model()
    optics_setup = optics_setup_timed.value
    @printf(
        "Persistent optics model + detector list + Descriptor(6,2): %.3f s\n",
        optics_setup_timed.time,
    )

    coasting_warmup_seconds = 0.0
    if run_warmup
        coasting_warmup_seconds = @elapsed solve_coasting_closed_orbits!(
            optics_setup.model,
            all_samples.names,
            values[1:1, :],
            closed_orbits[1:1, :];
            reltol,
            abstol,
            maxiter,
        )
        @printf("Coasting closed-orbit warmup/compilation: %.3f s\n", coasting_warmup_seconds)
    end
    coasting_timed = @timed solve_coasting_closed_orbits!(
        optics_setup.model,
        all_samples.names,
        values,
        closed_orbits;
        reltol,
        abstol,
        maxiter,
    )
    coasting_data = coasting_timed.value
    closed_orbits = coasting_data.orbits
    @printf(
        "Coasting closed orbits: %.3f s, residual max %.3e\n",
        coasting_timed.time,
        maximum(coasting_data.residuals),
    )

    reusable_workspace = nothing
    reusable_workspace_setup_seconds = 0.0
    if twiss_mode == "reuse"
        reusable_workspace_setup_seconds = @elapsed begin
            reusable_workspace = prepare_reusable_twiss_workspace(
                optics_setup.model.ring,
                optics_setup.detectors,
                optics_setup.descriptor,
                view(closed_orbits, 1, :),
            )
        end
        @printf(
            "Reusable Twiss DAMap/TPS workspace: %.3f s\n",
            reusable_workspace_setup_seconds,
        )
    end

    warmup_seconds = 0.0
    if run_warmup
        warmup_seconds = @elapsed calculate_one_chromatic_optics(
            optics_setup.model,
            optics_setup.detectors,
            all_samples.names,
            view(values, 1, :),
            view(closed_orbits, 1, :),
            optics_setup.descriptor,
            reusable_workspace,
        )
        @printf("Second-order Twiss warmup/compilation: %.3f s\n", warmup_seconds)
    end

    optics_results = Vector{Any}(undef, sample_count)
    sample_seconds = zeros(sample_count)
    control_update_seconds = zeros(sample_count)
    twiss_kernel_seconds = zeros(sample_count)
    twiss_timed = @timed begin
        for row in 1:sample_count
            calculation = nothing
            sample_seconds[row] = @elapsed begin
                calculation = calculate_one_chromatic_optics(
                    optics_setup.model,
                    optics_setup.detectors,
                    all_samples.names,
                    view(values, row, :),
                    view(closed_orbits, row, :),
                    optics_setup.descriptor,
                    reusable_workspace,
                )
            end
            optics_results[row] = calculation.optics
            control_update_seconds[row] = calculation.control_update_seconds
            twiss_kernel_seconds[row] = calculation.twiss_seconds
            @printf(
                "Second-order Twiss %d/%d (sample_id=%d): %.3f s (controls %.4f + twiss %.3f)\n",
                row,
                sample_count,
                sample_ids[row],
                sample_seconds[row],
                control_update_seconds[row],
                twiss_kernel_seconds[row],
            )
        end
    end

    detector_path = joinpath(output_dir, "scibmad_detector_chromatic_twiss.csv")
    ring_path = joinpath(output_dir, "scibmad_ring_chromatic_twiss.csv")
    orbit_path = joinpath(output_dir, "scibmad_start_closed_orbits.csv")
    metadata_path = joinpath(output_dir, "scibmad_chromatic_optics_metadata.toml")
    write_seconds = @elapsed begin
        write_chromatic_detector_optics(detector_path, sample_ids, optics_results)
        write_chromatic_ring_optics(ring_path, sample_ids, optics_results, sample_seconds)
        write_start_orbits(
            orbit_path,
            sample_ids,
            closed_orbits,
            coasting_data.residuals,
        )
    end

    metadata = Dict(
        "format" => "cesr-chromatic-optics-benchmark-v1",
        "engine" => "SciBmad",
        "method" => "second-order GTPSA one-turn map and periodic normal modes at DET_* beginnings",
        "twiss_mode" => twiss_mode,
        "twiss_reuses_damap_tps_workspace" => twiss_mode == "reuse",
        "twiss_workspace_setup_seconds" => reusable_workspace_setup_seconds,
        "rf_mode" => "off (coasting), required for delta to represent a fixed momentum offset",
        "derivative_coordinate" => "relative momentum deviation delta (phase-space variable 6)",
        "corrector_derivatives_included" => false,
        "closed_orbit_source" => "batched response-linear initial guess + frozen nominal Jacobian + full-AD fallback",
        "twiss_closed_orbit_resolved_internally" => false,
        "input_csv" => inputs,
        "output_directory" => output_dir,
        "sample_count" => sample_count,
        "control_count" => length(all_samples.names),
        "detector_count" => 99,
        "detector_row_count" => 99 * sample_count,
        "closed_orbit_wall_seconds" => closed_orbit_timed.time,
        "closed_orbit_warmup_seconds" => closed_orbit_warmup_seconds,
        "closed_orbit_newton_seconds" => closed_orbit_data.result.solve_seconds,
        "closed_orbit_fallback_seconds" => closed_orbit_data.result.fallback_seconds,
        "closed_orbit_fallback_count" => closed_orbit_data.result.fallback_count,
        "closed_orbit_maximum_residual" => maximum(closed_orbit_data.result.closure_norms),
        "coasting_closed_orbit_warmup_seconds" => coasting_warmup_seconds,
        "coasting_closed_orbit_wall_seconds" => coasting_timed.time,
        "coasting_closed_orbit_seconds_per_sample" => coasting_data.sample_seconds,
        "coasting_closed_orbit_maximum_residual" => maximum(coasting_data.residuals),
        "twiss_warmup_seconds" => warmup_seconds,
        "optics_model_setup_seconds" => optics_setup_timed.time,
        "optics_model_setup_allocated_bytes" => optics_setup_timed.bytes,
        "persistent_optics_model" => true,
        "persistent_detector_list" => true,
        "persistent_GTPSA_descriptor" => true,
        "twiss_physics_seconds" => twiss_timed.time,
        "twiss_seconds_per_sample" => sample_seconds,
        "control_update_seconds_per_sample" => control_update_seconds,
        "control_update_seconds_total" => sum(control_update_seconds),
        "twiss_kernel_seconds_per_sample" => twiss_kernel_seconds,
        "twiss_kernel_seconds_total" => sum(twiss_kernel_seconds),
        "twiss_seconds_per_sample_mean" => mean(sample_seconds),
        "twiss_seconds_per_sample_median" => median(sample_seconds),
        "twiss_samples_per_second" => sample_count / twiss_timed.time,
        "write_seconds" => write_seconds,
        "allocated_bytes" => twiss_timed.bytes,
        "gc_seconds" => twiss_timed.gctime,
        "julia_version" => string(VERSION),
        "julia_threads" => Threads.nthreads(),
        "reltol" => reltol,
        "abstol" => abstol,
        "maxiter" => maxiter,
        "GTPSA_descriptor" => "Descriptor(6, 2)",
        "saved_at" => "beginning of 99 DET_* elements",
    )
    mkpath(output_dir)
    open(metadata_path, "w") do io
        TOML.print(io, metadata; sorted=true)
    end

    @printf(
        "Second-order Twiss physics: %.3f s, %.3f samples/s (mean %.3f s/sample)\n",
        twiss_timed.time,
        sample_count / twiss_timed.time,
        mean(sample_seconds),
    )
    println("Detector chromatic optics: $detector_path")
    println("Ring chromatic optics:     $ring_path")
    println("Input orbits:              $orbit_path")
    println("Metadata:                  $metadata_path")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main_chromatic_optics())
end

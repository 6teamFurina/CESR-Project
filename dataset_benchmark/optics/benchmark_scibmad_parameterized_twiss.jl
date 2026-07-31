#!/usr/bin/env julia

"""
Compute the full CESR chromatic-optics dataset from one parameterized Twiss
calculation.

The descriptor keeps the six phase-space variables to configurable order
(third order by default) and adds all 119 corrector controls as first-order
GTPSA parameters. `twiss` is called once at the zero-corrector machine. Its TPS
outputs are then evaluated at all requested control samples. This is a local
first-order surrogate in corrector space. Third phase-space order is required
to retain the mixed `delta * corrector` terms needed for corrector-dependent
chromaticity.
"""

include(joinpath(@__DIR__, "benchmark_scibmad_chromatic_optics.jl"))

function parse_parameterized_args(args)
    options = Dict{String,String}(
        "inputs" => joinpath(DATASET_DIR, "inputs", "cesr_corrector_samples_1000.csv"),
        "sample-count" => "10",
        "output-dir" => joinpath(OPTICS_DIR, "results", "methods_10", "scibmad_parameterized"),
        "warmup" => "true",
        "phase-space-order" => "3",
        "reltol" => "1e-8",
        "abstol" => "1e-10",
        "maxiter" => "100",
        "response-matrix-cache" => joinpath(DATASET_DIR, "reference", "closed_orbit_response_6x119.csv"),
    )
    for argument in args
        startswith(argument, "--") ||
            error("Arguments must have --name=value form: $argument")
        fields = split(argument[3:end], "="; limit=2)
        length(fields) == 2 || error("Missing value in argument: $argument")
        haskey(options, fields[1]) || error("Unknown option: --$(fields[1])")
        options[fields[1]] = fields[2]
    end
    lowercase(options["warmup"]) in ("true", "false") ||
        error("--warmup must be true or false")
    parse(Int, options["phase-space-order"]) in (2, 3) ||
        error("--phase-space-order must be 2 or 3")
    return options
end

function nominal_coasting_orbit(names; response_matrix_cache, reltol, abstol, maxiter)
    zero_controls = zeros(1, length(names))
    initial = prepare_initial_guess(
        names,
        zero_controls,
        "nominal-z0";
        response_matrix_cache,
        recompute_response=false,
        reltol,
        abstol,
        maxiter,
    )
    scalar_setup = prepare_chromatic_optics_model()
    coast = solve_coasting_closed_orbits!(
        scalar_setup.model,
        names,
        zero_controls,
        initial.v0;
        reltol,
        abstol,
        maxiter,
    )
    return (; orbit=vec(coast.orbits), residual=only(coast.residuals))
end

function prepare_parameterized_optics_model(names, phase_space_order)
    control_count = length(names)
    descriptor = Descriptor(6, phase_space_order, control_count, 1)
    corrector_parameters = params(descriptor)
    model = load_cesr_model(
        zero_value=zero(corrector_parameters[1]),
        rf_on=false,
    )
    for (index, name) in enumerate(names)
        model.controls[name] = corrector_parameters[index]
    end
    detectors = detector_elements(model.ring)
    length(detectors) == 99 || error("Expected 99 detectors, found $(length(detectors))")
    return (; model, detectors, descriptor)
end

function calculate_parameterized_twiss(setup, nominal_orbit)
    return twiss(
        setup.model.ring;
        GTPSA_descriptor=setup.descriptor,
        at=setup.detectors,
        v0_and_coast=(reshape(copy(nominal_orbit), 1, 6), true),
        spin=false,
        RDTs=false,
        normalizing_map=false,
    )
end

function collect_detector_series(optics)
    columns = (CHROMATIC_TWISS_COLUMNS..., CHROMATIC_ORBIT_COLUMNS...)
    series = [
        getproperty(optics.table, column)[detector]
        for detector in eachindex(optics.table.name)
        for column in columns
    ]
    return (; columns, series, delta_series=deriv.(series, DELTA_INDEX))
end

function evaluate_parameterized_outputs(optics, values)
    sample_count, control_count = size(values)
    detector = collect_detector_series(optics)
    ring_series = collect(optics.tunes)
    ring_delta_series = deriv.(ring_series, DELTA_INDEX)

    # Parameter order is exactly one, so at zero phase-space amplitude every
    # requested quantity is affine in the corrector controls.  Extract the
    # coefficient matrices once and evaluate every sample with four dense
    # matrix products instead of thousands of scalar TPS evaluations.
    control_columns = 7:(6 + control_count)
    detector_base = constant_term.(detector.series)
    detector_delta_base = constant_term.(detector.delta_series)
    detector_control = Matrix(
        GTPSA.jacobian(detector.series; include_params=true)[:, control_columns],
    )
    detector_delta_control = Matrix(
        GTPSA.jacobian(detector.delta_series; include_params=true)[:, control_columns],
    )
    ring_base = constant_term.(ring_series)
    ring_delta_base = constant_term.(ring_delta_series)
    ring_control = Matrix(
        GTPSA.jacobian(ring_series; include_params=true)[:, control_columns],
    )
    ring_delta_control = Matrix(
        GTPSA.jacobian(ring_delta_series; include_params=true)[:, control_columns],
    )

    detector_values = values * transpose(detector_control)
    detector_values .+= transpose(detector_base)
    detector_delta = values * transpose(detector_delta_control)
    detector_delta .+= transpose(detector_delta_base)
    ring_values = values * transpose(ring_control)
    ring_values .+= transpose(ring_base)
    # With RF off, the constant term of the longitudinal one-turn coefficient
    # is identically zero. Enforce that invariant after affine evaluation so a
    # tiny parameter-expansion artifact is not reported as a physical value.
    ring_values[:, 3] .= 0.0
    ring_delta = values * transpose(ring_delta_control)
    ring_delta .+= transpose(ring_delta_base)
    sample_seconds = fill(0.0, sample_count)
    return (;
        detector,
        detector_values,
        detector_delta,
        ring_values,
        ring_delta,
        sample_seconds,
    )
end

function write_parameterized_detector(path, sample_ids, optics, evaluated)
    mkpath(dirname(path))
    open(path, "w") do io
        header = String["sample_id", "s", "beamline_index", "name"]
        for column in evaluated.detector.columns
            push!(header, String(column))
            push!(header, "d$(column)_ddelta")
        end
        println(io, join(header, ','))

        column_count = length(evaluated.detector.columns)
        for sample_row in eachindex(sample_ids)
            for detector_row in eachindex(optics.table.name)
                fields = Any[
                    sample_ids[sample_row],
                    optics.table.s[detector_row],
                    optics.table.beamline_index[detector_row],
                    optics.table.name[detector_row],
                ]
                offset = (detector_row - 1) * column_count
                for column in 1:column_count
                    index = offset + column
                    push!(fields, evaluated.detector_values[sample_row, index])
                    push!(fields, evaluated.detector_delta[sample_row, index])
                end
                println(io, join(fields, ','))
            end
        end
    end
    return path
end

function write_parameterized_ring(
    path,
    sample_ids,
    evaluated,
    parameterized_twiss_seconds,
)
    mkpath(dirname(path))
    shared_seconds = parameterized_twiss_seconds / length(sample_ids)
    open(path, "w") do io
        println(
            io,
            "sample_id,Q1_signed,Q2_signed,Qx_fractional,Qy_fractional," *
            "slip_tps_constant,xi_1,xi_2,slip_factor,twiss_seconds",
        )
        for row in eachindex(sample_ids)
            tunes = view(evaluated.ring_values, row, :)
            chromaticities = view(evaluated.ring_delta, row, :)
            println(
                io,
                join(
                    (
                        sample_ids[row],
                        tunes[1],
                        tunes[2],
                        mod(tunes[1], 1),
                        mod(tunes[2], 1),
                        tunes[3],
                        chromaticities[1],
                        chromaticities[2],
                        chromaticities[3],
                        shared_seconds + evaluated.sample_seconds[row],
                    ),
                    ',',
                ),
            )
        end
    end
    return path
end

function main_parameterized_twiss(args=ARGS)
    options = parse_parameterized_args(args)
    inputs = abspath(options["inputs"])
    output_dir = abspath(options["output-dir"])
    sample_count = parse(Int, options["sample-count"])
    run_warmup = lowercase(options["warmup"]) == "true"
    reltol = parse(Float64, options["reltol"])
    abstol = parse(Float64, options["abstol"])
    maxiter = parse(Int, options["maxiter"])
    response_matrix_cache = abspath(options["response-matrix-cache"])
    phase_space_order = parse(Int, options["phase-space-order"])

    samples = read_samples(inputs)
    1 <= sample_count <= size(samples.values, 1) || error(
        "--sample-count must be between 1 and $(size(samples.values, 1))",
    )
    sample_ids = samples.sample_ids[1:sample_count]
    values = Matrix(samples.values[1:sample_count, :])

    @printf(
        "SciBmad single parameterized Twiss: %d samples x %d controls x 99 detectors\n",
        sample_count,
        length(samples.names),
    )
    nominal_timed = @timed nominal_coasting_orbit(
        samples.names;
        response_matrix_cache,
        reltol,
        abstol,
        maxiter,
    )
    nominal = nominal_timed.value
    @printf(
        "Nominal coasting orbit setup: %.3f s, residual %.3e\n",
        nominal_timed.time,
        nominal.residual,
    )

    setup_timed = @timed prepare_parameterized_optics_model(
        samples.names,
        phase_space_order,
    )
    setup = setup_timed.value
    @printf(
        "Parameterized model + Descriptor(6,%d,119,1): %.3f s\n",
        phase_space_order,
        setup_timed.time,
    )

    warmup_seconds = 0.0
    evaluation_warmup_seconds = 0.0
    if run_warmup
        warmup_optics = nothing
        warmup_seconds = @elapsed begin
            warmup_optics = calculate_parameterized_twiss(setup, nominal.orbit)
        end
        evaluation_warmup_seconds = @elapsed evaluate_parameterized_outputs(
            warmup_optics,
            values[1:1, :],
        )
        @printf("Parameterized Twiss warmup/compilation: %.3f s\n", warmup_seconds)
        @printf("Affine evaluation warmup/compilation: %.3f s\n", evaluation_warmup_seconds)
    end

    twiss_timed = @timed calculate_parameterized_twiss(setup, nominal.orbit)
    optics = twiss_timed.value
    length(optics.table.name) == 99 ||
        error("Twiss returned $(length(optics.table.name)) detector rows")
    @printf("One parameterized Twiss physics call: %.3f s\n", twiss_timed.time)

    evaluation_timed = @timed evaluate_parameterized_outputs(optics, values)
    evaluated = evaluation_timed.value
    physics_seconds = twiss_timed.time + evaluation_timed.time
    @printf(
        "TPS evaluation: %.3f s; total physics %.3f s; %.3f samples/s\n",
        evaluation_timed.time,
        physics_seconds,
        sample_count / physics_seconds,
    )

    detector_path = joinpath(output_dir, "scibmad_detector_chromatic_twiss.csv")
    ring_path = joinpath(output_dir, "scibmad_ring_chromatic_twiss.csv")
    metadata_path = joinpath(output_dir, "scibmad_chromatic_optics_metadata.toml")
    write_seconds = @elapsed begin
        write_parameterized_detector(detector_path, sample_ids, optics, evaluated)
        write_parameterized_ring(
            ring_path,
            sample_ids,
            evaluated,
            twiss_timed.time,
        )
    end

    metadata = Dict(
        "format" => "cesr-chromatic-optics-benchmark-v1",
        "engine" => "SciBmad",
        "method" => "one parameterized Twiss map, first order in 119 corrector parameters",
        "twiss_mode" => "single-parameterized",
        "parameterized_twiss_is_local_surrogate" => true,
        "phase_space_order" => phase_space_order,
        "corrector_parameter_order" => 1,
        "GTPSA_descriptor" => "Descriptor(6, $phase_space_order, 119, 1)",
        "input_csv" => inputs,
        "output_directory" => output_dir,
        "sample_count" => sample_count,
        "control_count" => length(samples.names),
        "detector_count" => 99,
        "detector_row_count" => 99 * sample_count,
        "nominal_coasting_setup_seconds" => nominal_timed.time,
        "nominal_coasting_residual" => nominal.residual,
        "parameterized_model_setup_seconds" => setup_timed.time,
        "parameterized_model_setup_allocated_bytes" => setup_timed.bytes,
        "twiss_warmup_seconds" => warmup_seconds,
        "parameter_evaluation_warmup_seconds" => evaluation_warmup_seconds,
        "parameterized_twiss_seconds" => twiss_timed.time,
        "parameterized_twiss_allocated_bytes" => twiss_timed.bytes,
        "parameterized_twiss_gc_seconds" => twiss_timed.gctime,
        "parameter_evaluation_seconds" => evaluation_timed.time,
        "parameter_evaluation_method" => "affine coefficient extraction plus dense matrix multiplication",
        "parameter_evaluation_seconds_per_sample" => evaluated.sample_seconds,
        "twiss_physics_seconds" => physics_seconds,
        "twiss_samples_per_second" => sample_count / physics_seconds,
        "write_seconds" => write_seconds,
        "julia_version" => string(VERSION),
        "julia_threads" => Threads.nthreads(),
        "saved_at" => "beginning of 99 DET_* elements",
    )
    mkpath(output_dir)
    open(metadata_path, "w") do io
        TOML.print(io, metadata; sorted=true)
    end
    println("Detector chromatic optics: $detector_path")
    println("Ring chromatic optics:     $ring_path")
    println("Metadata:                  $metadata_path")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main_parameterized_twiss())
end

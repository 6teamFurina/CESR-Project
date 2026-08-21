#!/usr/bin/env julia

"""
Compute all 119 RF-off CESR corrector derivatives with one GTPSA Twiss call.

The periodic closed-orbit derivative is obtained first from the implicit
fixed-point equation X = (I-A)\\B.  That parameterized periodic orbit is then
passed to one Twiss calculation, so detector orbit and optics derivatives both
include the response of the initial closed orbit.
"""

include(joinpath(@__DIR__, "..", "benchmark_scibmad_parameterized_twiss.jl"))

function parse_jacobian_args(args)
    options = Dict{String,String}(
        "inputs" => joinpath(ORBIT_CALCULATION_DIR, "inputs", "cesr_corrector_samples_1000.csv"),
        "output-dir" => joinpath(@__DIR__, "results", "scibmad"),
        "phase-space-order" => "2",
        "warmup" => "true",
    )
    for argument in args
        startswith(argument, "--") || error("Arguments must have --name=value form: $argument")
        fields = split(argument[3:end], '='; limit=2)
        length(fields) == 2 || error("Missing value in argument: $argument")
        haskey(options, fields[1]) || error("Unknown option: --$(fields[1])")
        options[fields[1]] = fields[2]
    end
    parse(Int, options["phase-space-order"]) in (2, 3) ||
        error("--phase-space-order must be 2 or 3")
    lowercase(options["warmup"]) in ("true", "false") ||
        error("--warmup must be true or false")
    return options
end

function track_map(ring, coordinates)
    bunch = Bunch(v=reshape(coordinates, 1, 6))
    SciBmad.BTBL.check_bl_bunch!(bunch, ring, false)
    track!(bunch, ring)
    return vec(bunch.coords.v)
end

function implicit_coasting_response(model, descriptor, closed_orbit, control_count)
    variables = vars(descriptor)
    input_map = [closed_orbit[index] + copy(variables[index]) for index in 1:6]
    output_map = track_map(model.ring, input_map)
    full_jacobian = Matrix(GTPSA.jacobian(output_map; include_params=true))
    size(full_jacobian) == (6, 6 + control_count) ||
        error("Unexpected one-turn Jacobian size $(size(full_jacobian))")
    A = full_jacobian[:, 1:6]
    B = full_jacobian[:, 7:end]
    response = zeros(6, control_count)
    response[1:4, :] .= (I - A[1:4, 1:4]) \ B[1:4, :]
    residual = (I - A[1:4, 1:4]) * response[1:4, :] - B[1:4, :]
    return (; response, residual)
end

function calculate_jacobian_twiss(setup, closed_orbit)
    # The Twiss normal form supplies its own parameter-dependent fixed-point
    # translation. The Float64 orbit is only the nominal expansion reference;
    # inserting the implicit response here as well would double count it.
    return calculate_parameterized_twiss(setup, closed_orbit)
end

function detector_series(optics)
    columns = (CHROMATIC_TWISS_COLUMNS..., CHROMATIC_ORBIT_COLUMNS...)
    series = [
        getproperty(optics.table, column)[detector]
        for detector in eachindex(optics.table.name)
        for column in columns
    ]
    labels = [
        "$(lowercase(String(optics.table.name[detector]))):$(String(column))"
        for detector in eachindex(optics.table.name)
        for column in columns
    ]
    return (; columns, series, labels)
end

function write_labeled_matrix(path, row_labels, column_labels, values)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, join(vcat("observable", column_labels), ','))
        for row in eachindex(row_labels)
            println(io, join(vcat(row_labels[row], values[row, :]), ','))
        end
    end
    return path
end

function main_scibmad_jacobian(args=ARGS)
    options = parse_jacobian_args(args)
    inputs = abspath(options["inputs"])
    output_dir = abspath(options["output-dir"])
    phase_space_order = parse(Int, options["phase-space-order"])
    run_warmup = lowercase(options["warmup"]) == "true"
    samples = read_samples(inputs)
    names = samples.names
    control_count = length(names)

    baseline_timed = @timed begin
        scalar_model = load_cesr_model(rf_on=false)
        solution = find_closed_orbit_coasting_forwarddiff(
            scalar_model.ring;
            coasting_beam=true,
        )
        (orbit=solution.orbit, residual=solution.residual)
    end
    baseline = baseline_timed.value

    setup_timed = @timed prepare_parameterized_optics_model(names, phase_space_order)
    setup = setup_timed.value
    response_warmup_seconds = 0.0
    if run_warmup
        response_warmup_seconds = @elapsed implicit_coasting_response(
            setup.model,
            setup.descriptor,
            baseline.orbit,
            control_count,
        )
        @printf("SciBmad implicit closed-orbit warmup: %.3f s\n", response_warmup_seconds)
    end
    response_timed = @timed implicit_coasting_response(
        setup.model,
        setup.descriptor,
        baseline.orbit,
        control_count,
    )
    response = response_timed.value
    warmup_seconds = 0.0
    extraction_warmup_seconds = 0.0
    if run_warmup
        warmup_optics = nothing
        warmup_seconds = @elapsed warmup_optics = calculate_jacobian_twiss(
            setup,
            baseline.orbit,
        )
        extraction_warmup_seconds = @elapsed begin
            warmup_detector = detector_series(warmup_optics)
            Matrix(GTPSA.jacobian(warmup_detector.series; include_params=true))[:, 7:end]
        end
        @printf("SciBmad parameterized Twiss warmup: %.3f s\n", warmup_seconds)
    end

    twiss_timed = @timed calculate_jacobian_twiss(
        setup,
        baseline.orbit,
    )
    optics = twiss_timed.value
    extraction_timed = @timed begin
        detector = detector_series(optics)
        detector_jacobian = Matrix(
            GTPSA.jacobian(detector.series; include_params=true)[:, 7:end],
        )
        tune_jacobian = Matrix(
            GTPSA.jacobian(collect(optics.tunes); include_params=true)[:, 7:end],
        )
        (; detector, detector_jacobian, tune_jacobian)
    end
    extracted = extraction_timed.value
    physics_seconds = response_timed.time + twiss_timed.time + extraction_timed.time

    detector_path = joinpath(output_dir, "scibmad_detector_optics_jacobian.csv")
    orbit_path = joinpath(output_dir, "scibmad_closed_orbit_jacobian_6x119.csv")
    ring_path = joinpath(output_dir, "scibmad_ring_tune_jacobian.csv")
    metadata_path = joinpath(output_dir, "scibmad_corrector_jacobian_metadata.toml")
    write_timed = @timed begin
        write_labeled_matrix(detector_path, extracted.detector.labels, names, extracted.detector_jacobian)
        write_labeled_matrix(
            orbit_path,
            ["start:x", "start:px", "start:y", "start:py", "start:z", "start:pz"],
            names,
            response.response,
        )
        write_labeled_matrix(ring_path, ["ring:Q1", "ring:Q2", "ring:slip"], names, extracted.tune_jacobian)
    end

    metadata = Dict(
        "format" => "cesr-corrector-optics-jacobian-v1",
        "engine" => "SciBmad/GTPSA",
        "method" => "one parameterized Twiss plus implicit periodic closed-orbit differentiation",
        "rf_mode" => "off (4D coasting)",
        "GTPSA_descriptor" => "Descriptor(6, $phase_space_order, $control_count, 1)",
        "phase_space_order" => phase_space_order,
        "corrector_parameter_order" => 1,
        "control_count" => control_count,
        "detector_count" => length(optics.table.name),
        "detector_quantity_count" => length(extracted.detector.columns),
        "detector_jacobian_shape" => collect(size(extracted.detector_jacobian)),
        "closed_orbit_jacobian_shape" => collect(size(response.response)),
        "baseline_closed_orbit_seconds" => baseline_timed.time,
        "baseline_closed_orbit_residual" => baseline.residual,
        "parameterized_model_setup_seconds" => setup_timed.time,
        "implicit_response_warmup_seconds" => response_warmup_seconds,
        "implicit_response_seconds" => response_timed.time,
        "implicit_response_max_abs_residual" => maximum(abs, response.residual),
        "warmup_enabled" => run_warmup,
        "twiss_warmup_seconds" => warmup_seconds,
        "extraction_warmup_seconds" => extraction_warmup_seconds,
        "parameterized_twiss_seconds" => twiss_timed.time,
        "jacobian_extraction_seconds" => extraction_timed.time,
        "physics_seconds" => physics_seconds,
        "write_seconds" => write_timed.time,
        "timed_region" => "implicit 6x119 periodic-orbit response + one parameterized Twiss + coefficient extraction",
        "input_csv" => inputs,
        "julia_version" => string(VERSION),
        "julia_threads" => Threads.nthreads(),
    )
    mkpath(output_dir)
    open(metadata_path, "w") do io
        TOML.print(io, metadata; sorted=true)
    end

    @printf("SciBmad corrector Jacobian physics: %.3f s\n", physics_seconds)
    @printf("  implicit closed orbit: %.3f s\n", response_timed.time)
    @printf("  one parameterized Twiss: %.3f s\n", twiss_timed.time)
    @printf("  coefficient extraction: %.3f s\n", extraction_timed.time)
    println("Detector Jacobian: $detector_path")
    println("Closed-orbit Jacobian: $orbit_path")
    println("Metadata: $metadata_path")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main_scibmad_jacobian())
end

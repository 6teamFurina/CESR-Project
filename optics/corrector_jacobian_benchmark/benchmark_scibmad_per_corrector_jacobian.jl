#!/usr/bin/env julia

"""
Compute the 119-column RF-off CESR corrector-to-optics Jacobian using the
SciBmad-recommended narrow-parameter strategy: one `Descriptor(6, 2, 1, 1)`
Twiss calculation per corrector.

The descriptor and typed lattice are constructed once. For each column, only
the selected corrector is assigned the single GTPSA parameter; the preceding
corrector is reset to a typed zero. The independent columns are assembled into
the same labeled matrices as the existing wide-descriptor and Bmad benchmarks.
"""

include(joinpath(@__DIR__, "benchmark_scibmad_corrector_jacobian.jl"))

function parse_per_corrector_args(args)
    options = Dict{String,String}(
        "inputs" => joinpath(
            ORBIT_CALCULATION_DIR,
            "inputs",
            "cesr_corrector_samples_1000.csv",
        ),
        "output-dir" => joinpath(@__DIR__, "results", "scibmad_per_corrector"),
        "phase-space-order" => "2",
        "control-count" => "119",
        "warmup" => "true",
    )
    for argument in args
        startswith(argument, "--") ||
            error("Arguments must have --name=value form: $argument")
        fields = split(argument[3:end], '='; limit=2)
        length(fields) == 2 || error("Missing value in argument: $argument")
        haskey(options, fields[1]) || error("Unknown option: --$(fields[1])")
        options[fields[1]] = fields[2]
    end
    parse(Int, options["phase-space-order"]) == 2 ||
        error("This ordinary-optics benchmark requires --phase-space-order=2")
    parse(Int, options["control-count"]) in 1:119 ||
        error("--control-count must be between 1 and 119")
    lowercase(options["warmup"]) in ("true", "false") ||
        error("--warmup must be true or false")
    return options
end

function prepare_single_parameter_model(phase_space_order)
    descriptor = Descriptor(6, phase_space_order, 1, 1)
    parameter = only(params(descriptor))
    typed_zero = zero(parameter)
    model = load_cesr_model(zero_value=typed_zero, rf_on=false)
    detectors = detector_elements(model.ring)
    length(detectors) == 99 || error("Expected 99 detectors, found $(length(detectors))")
    setup = (; model, detectors, descriptor)
    return (; setup, parameter, typed_zero)
end

function activate_single_corrector!(single, current_name, previous_name=nothing)
    if !isnothing(previous_name)
        single.setup.model.controls[previous_name] = single.typed_zero
    end
    single.setup.model.controls[current_name] = single.parameter
    return single
end

function extract_single_column(optics)
    detector = detector_series(optics)
    detector_full = GTPSA.jacobian(detector.series; include_params=true)
    tune_full = GTPSA.jacobian(collect(optics.tunes); include_params=true)
    size(detector_full, 2) == 7 ||
        error("Expected six phase-space columns plus one parameter column")
    size(tune_full, 2) == 7 ||
        error("Unexpected tune Jacobian width $(size(tune_full, 2))")
    return (;
        detector,
        detector_column=Float64.(detector_full[:, 7]),
        tune_column=Float64.(tune_full[:, 7]),
    )
end

function finite_column(values, label)
    all(isfinite, values) || error("Non-finite values in $label")
    return values
end

function distribution_summary(values)
    return Dict(
        "minimum" => minimum(values),
        "median" => median(values),
        "maximum" => maximum(values),
        "mean" => mean(values),
    )
end

function main_scibmad_per_corrector(args=ARGS)
    options = parse_per_corrector_args(args)
    inputs = abspath(options["inputs"])
    output_dir = abspath(options["output-dir"])
    phase_space_order = parse(Int, options["phase-space-order"])
    requested_count = parse(Int, options["control-count"])
    run_warmup = lowercase(options["warmup"]) == "true"

    samples = read_samples(inputs)
    names = samples.names[1:requested_count]
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

    setup_timed = @timed prepare_single_parameter_model(phase_space_order)
    single = setup_timed.value

    response_warmup_seconds = 0.0
    twiss_warmup_seconds = 0.0
    extraction_warmup_seconds = 0.0
    if run_warmup
        activate_single_corrector!(single, first(names))
        response_warmup_seconds = @elapsed implicit_coasting_response(
            single.setup.model,
            single.setup.descriptor,
            baseline.orbit,
            1,
        )
        warmup_optics = nothing
        twiss_warmup_seconds = @elapsed warmup_optics = calculate_jacobian_twiss(
            single.setup,
            baseline.orbit,
        )
        extraction_warmup_seconds = @elapsed extract_single_column(warmup_optics)
        single.setup.model.controls[first(names)] = single.typed_zero
        @printf(
            "P=1 warmup: response %.3f s, Twiss %.3f s, extraction %.3f s\n",
            response_warmup_seconds,
            twiss_warmup_seconds,
            extraction_warmup_seconds,
        )
        flush(stdout)
        GC.gc()
    end

    detector_jacobian = Matrix{Float64}(undef, 0, control_count)
    orbit_jacobian = Matrix{Float64}(undef, 6, control_count)
    tune_jacobian = Matrix{Float64}(undef, 0, control_count)
    detector_labels = String[]
    detector_columns = ()
    activation_seconds = zeros(control_count)
    response_seconds = zeros(control_count)
    twiss_seconds = zeros(control_count)
    extraction_seconds = zeros(control_count)
    response_allocated_bytes = zeros(Int, control_count)
    twiss_allocated_bytes = zeros(Int, control_count)
    extraction_allocated_bytes = zeros(Int, control_count)
    response_residuals = zeros(control_count)

    loop_timed = @timed begin
        previous_name = nothing
        for (column, name) in enumerate(names)
            activation_seconds[column] = @elapsed activate_single_corrector!(
                single,
                name,
                previous_name,
            )

            response_timed = @timed implicit_coasting_response(
                single.setup.model,
                single.setup.descriptor,
                baseline.orbit,
                1,
            )
            response = response_timed.value
            response_seconds[column] = response_timed.time
            response_allocated_bytes[column] = response_timed.bytes
            response_residuals[column] = maximum(abs, response.residual)
            orbit_jacobian[:, column] .= finite_column(
                vec(response.response),
                "$name periodic-orbit response",
            )

            twiss_timed = @timed calculate_jacobian_twiss(
                single.setup,
                baseline.orbit,
            )
            twiss_seconds[column] = twiss_timed.time
            twiss_allocated_bytes[column] = twiss_timed.bytes

            extraction_timed = @timed extract_single_column(twiss_timed.value)
            extracted = extraction_timed.value
            extraction_seconds[column] = extraction_timed.time
            extraction_allocated_bytes[column] = extraction_timed.bytes

            if column == 1
                detector_labels = extracted.detector.labels
                detector_columns = extracted.detector.columns
                detector_jacobian = Matrix{Float64}(
                    undef,
                    length(extracted.detector_column),
                    control_count,
                )
                tune_jacobian = Matrix{Float64}(
                    undef,
                    length(extracted.tune_column),
                    control_count,
                )
            else
                extracted.detector.labels == detector_labels ||
                    error("Detector labels changed for $name")
            end
            detector_jacobian[:, column] .= finite_column(
                extracted.detector_column,
                "$name detector Jacobian",
            )
            tune_jacobian[:, column] .= finite_column(
                extracted.tune_column,
                "$name tune Jacobian",
            )

            previous_name = name
            if column == 1 || column % 10 == 0 || column == control_count
                @printf(
                    "P=1 corrector %d/%d %-6s: response %.3f s, Twiss %.3f s\n",
                    column,
                    control_count,
                    name,
                    response_seconds[column],
                    twiss_seconds[column],
                )
                flush(stdout)
            end
        end
        if !isnothing(previous_name)
            single.setup.model.controls[previous_name] = single.typed_zero
        end
    end

    optics_only_seconds =
        sum(activation_seconds) + sum(twiss_seconds) + sum(extraction_seconds)
    physics_with_response_seconds = optics_only_seconds + sum(response_seconds)

    detector_path = joinpath(
        output_dir,
        "scibmad_per_corrector_detector_optics_jacobian.csv",
    )
    orbit_path = joinpath(
        output_dir,
        "scibmad_per_corrector_closed_orbit_jacobian.csv",
    )
    ring_path = joinpath(
        output_dir,
        "scibmad_per_corrector_ring_tune_jacobian.csv",
    )
    metadata_path = joinpath(
        output_dir,
        "scibmad_per_corrector_jacobian_metadata.toml",
    )

    write_timed = @timed begin
        write_labeled_matrix(detector_path, detector_labels, names, detector_jacobian)
        write_labeled_matrix(
            orbit_path,
            ["start:x", "start:px", "start:y", "start:py", "start:z", "start:pz"],
            names,
            orbit_jacobian,
        )
        write_labeled_matrix(
            ring_path,
            ["ring:Q1", "ring:Q2", "ring:slip"],
            names,
            tune_jacobian,
        )
    end

    metadata = Dict(
        "format" => "cesr-per-corrector-optics-jacobian-v1",
        "engine" => "SciBmad/GTPSA",
        "method" => "one Descriptor(6,2,1,1) parameterized Twiss per corrector",
        "rf_mode" => "off (4D coasting)",
        "GTPSA_descriptor" => "Descriptor(6, $phase_space_order, 1, 1)",
        "phase_space_order" => phase_space_order,
        "corrector_parameter_order" => 1,
        "control_count" => control_count,
        "controls" => names,
        "detector_count" => div(length(detector_labels), length(detector_columns)),
        "detector_quantity_count" => length(detector_columns),
        "detector_jacobian_shape" => collect(size(detector_jacobian)),
        "closed_orbit_jacobian_shape" => collect(size(orbit_jacobian)),
        "ring_tune_jacobian_shape" => collect(size(tune_jacobian)),
        "baseline_closed_orbit_seconds" => baseline_timed.time,
        "baseline_closed_orbit_residual" => baseline.residual,
        "model_setup_seconds" => setup_timed.time,
        "warmup_enabled" => run_warmup,
        "implicit_response_warmup_seconds" => response_warmup_seconds,
        "twiss_warmup_seconds" => twiss_warmup_seconds,
        "extraction_warmup_seconds" => extraction_warmup_seconds,
        "parameter_activation_seconds" => sum(activation_seconds),
        "implicit_response_seconds" => sum(response_seconds),
        "parameterized_twiss_seconds" => sum(twiss_seconds),
        "jacobian_extraction_seconds" => sum(extraction_seconds),
        "optics_only_seconds" => optics_only_seconds,
        "physics_seconds_with_separate_start_orbit_response" => physics_with_response_seconds,
        "loop_wall_seconds" => loop_timed.time,
        "write_seconds" => write_timed.time,
        "twiss_seconds_per_corrector" => twiss_seconds,
        "implicit_response_seconds_per_corrector" => response_seconds,
        "extraction_seconds_per_corrector" => extraction_seconds,
        "twiss_seconds_distribution" => distribution_summary(twiss_seconds),
        "implicit_response_seconds_distribution" => distribution_summary(response_seconds),
        "twiss_allocated_bytes" => sum(twiss_allocated_bytes),
        "implicit_response_allocated_bytes" => sum(response_allocated_bytes),
        "extraction_allocated_bytes" => sum(extraction_allocated_bytes),
        "implicit_response_max_abs_residual" => maximum(response_residuals),
        "timed_region" => "119 serial P=1 parameter activations, implicit responses, Twiss calls, and coefficient extractions",
        "input_csv" => inputs,
        "julia_version" => string(VERSION),
        "julia_threads" => Threads.nthreads(),
    )
    mkpath(output_dir)
    open(metadata_path, "w") do io
        TOML.print(io, metadata; sorted=true)
    end

    @printf("SciBmad per-corrector P=1 optics-only: %.3f s\n", optics_only_seconds)
    @printf(
        "SciBmad per-corrector P=1 with separate orbit response: %.3f s\n",
        physics_with_response_seconds,
    )
    @printf("  parameter activation: %.3f s\n", sum(activation_seconds))
    @printf("  implicit responses: %.3f s\n", sum(response_seconds))
    @printf("  parameterized Twiss: %.3f s\n", sum(twiss_seconds))
    @printf("  coefficient extraction: %.3f s\n", sum(extraction_seconds))
    println("Detector Jacobian: $detector_path")
    println("Metadata: $metadata_path")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main_scibmad_per_corrector())
end

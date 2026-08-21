#!/usr/bin/env julia

"""Benchmark SciBmad optics Jacobians for correctors plus Kn1/Kn2 strengths."""

include(joinpath(@__DIR__, "benchmark_scibmad_corrector_jacobian.jl"))

function parse_extended_args(args)
    options = Dict{String,String}(
        "inputs" => joinpath(ORBIT_CALCULATION_DIR, "inputs", "cesr_corrector_samples_1000.csv"),
        "case" => "both",
        "output-root" => joinpath(@__DIR__, "results", "extended"),
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
    options["case"] in ("correctors_quads", "correctors_quads_sextupoles", "both") ||
        error("--case must be correctors_quads, correctors_quads_sextupoles, or both")
    parse(Int, options["phase-space-order"]) == 2 ||
        error("This ordinary-optics benchmark requires --phase-space-order=2")
    lowercase(options["warmup"]) in ("true", "false") ||
        error("--warmup must be true or false")
    return options
end

base_bmad_name(element) = first(split(uppercase(String(element.name)), '!'))

function active_strength_bases(ring, attribute::Symbol)
    result = Set{String}()
    for element in ring.line
        strength = Float64(GTPSA.scalar(Beamlines.deval(getproperty(element, attribute))))
        iszero(strength) || push!(result, base_bmad_name(element))
    end
    return sort!(collect(result))
end

function add_strength_parameter!(ring, base_name, attribute, parameter)
    matches = 0
    for element in ring.line
        base_bmad_name(element) == base_name || continue
        baseline = Beamlines.deval(getproperty(element, attribute))
        iszero(Float64(GTPSA.scalar(baseline))) && continue
        setproperty!(element, attribute, baseline + parameter)
        matches += 1
    end
    matches > 0 || error("No active $attribute slices found for $base_name")
    return matches
end

function prepare_extended_model(corrector_names, case, phase_space_order)
    scalar = load_cesr_model(rf_on=false)
    quadrupoles = active_strength_bases(scalar.ring, :Kn1)
    sextupoles = active_strength_bases(scalar.ring, :Kn2)
    length(quadrupoles) == 106 || error("Expected 106 active quadrupoles, found $(length(quadrupoles))")
    length(sextupoles) == 76 || error("Expected 76 active sextupoles, found $(length(sextupoles))")
    selected_sextupoles = case == "correctors_quads" ? String[] : sextupoles
    parameter_count = length(corrector_names) + length(quadrupoles) + length(selected_sextupoles)
    descriptor = Descriptor(6, phase_space_order, parameter_count, 1)
    parameters = params(descriptor)
    model = load_cesr_model(zero_value=zero(parameters[1]), rf_on=false)

    labels = String[]
    families = String[]
    parameter_index = 0
    for name in corrector_names
        parameter_index += 1
        model.controls[name] = parameters[parameter_index]
        push!(labels, "COR:$name")
        push!(families, "corrector")
    end
    for name in quadrupoles
        parameter_index += 1
        add_strength_parameter!(model.ring, name, :Kn1, parameters[parameter_index])
        push!(labels, "K1:$name")
        push!(families, "quadrupole")
    end
    for name in selected_sextupoles
        parameter_index += 1
        add_strength_parameter!(model.ring, name, :Kn2, parameters[parameter_index])
        push!(labels, "K2:$name")
        push!(families, "sextupole")
    end
    parameter_index == parameter_count || error("Extended parameter count mismatch")
    detectors = detector_elements(model.ring)
    length(detectors) == 99 || error("Expected 99 detectors, found $(length(detectors))")
    return (; model, detectors, descriptor, labels, families, quadrupoles, sextupoles=selected_sextupoles)
end

function run_extended_case(
    corrector_names,
    baseline,
    case,
    output_dir,
    phase_space_order,
    run_warmup,
    baseline_seconds,
)
    mkpath(output_dir)
    progress_path = joinpath(output_dir, "scibmad_extended_jacobian_progress.toml")
    progress = Dict{String,Any}(
        "case" => case,
        "status" => "setup",
    )
    function checkpoint!(status; fields...)
        progress["status"] = status
        for (key, value) in pairs(fields)
            progress[String(key)] = value
        end
        open(progress_path, "w") do io
            TOML.print(io, progress; sorted=true)
        end
        println("Checkpoint $case: $status")
        flush(stdout)
    end

    setup_timed = @timed prepare_extended_model(corrector_names, case, phase_space_order)
    setup = setup_timed.value
    parameter_count = length(setup.labels)
    @printf(
        "SciBmad extended Jacobian %s: P=%d (119 correctors + %d quads + %d sextupoles)\n",
        case,
        parameter_count,
        length(setup.quadrupoles),
        length(setup.sextupoles),
    )
    flush(stdout)
    checkpoint!("setup_complete"; parameter_count, setup_seconds=setup_timed.time)

    response_warmup_seconds = 0.0
    twiss_warmup_seconds = 0.0
    extraction_warmup_seconds = 0.0
    if run_warmup
        warmup_response = nothing
        response_warmup_seconds = @elapsed warmup_response = implicit_coasting_response(
                setup.model,
                setup.descriptor,
                baseline.orbit,
                parameter_count,
            )
        checkpoint!("response_warmup_complete"; response_warmup_seconds)
        warmup_optics = nothing
        twiss_warmup_seconds = @elapsed warmup_optics = calculate_jacobian_twiss(
            setup,
            baseline.orbit,
        )
        checkpoint!("twiss_warmup_complete"; response_warmup_seconds, twiss_warmup_seconds)
        extraction_warmup_seconds = @elapsed begin
            warmup_detector = detector_series(warmup_optics)
            Matrix(GTPSA.jacobian(warmup_detector.series; include_params=true))[:, 7:end]
        end
        # Keep the fixed-point solve alive through the warmup block so its
        # allocations are not optimized away by future refactors.
        maximum(abs, warmup_response.residual) < 1e-10 || error("Warmup response failed")
        checkpoint!(
            "extraction_warmup_complete";
            response_warmup_seconds,
            twiss_warmup_seconds,
            extraction_warmup_seconds,
        )
    end

    response_timed = @timed implicit_coasting_response(
        setup.model,
        setup.descriptor,
        baseline.orbit,
        parameter_count,
    )
    response = response_timed.value
    checkpoint!("stable_response_complete"; implicit_response_seconds=response_timed.time)
    twiss_timed = @timed calculate_jacobian_twiss(setup, baseline.orbit)
    optics = twiss_timed.value
    checkpoint!(
        "stable_twiss_complete";
        implicit_response_seconds=response_timed.time,
        parameterized_twiss_seconds=twiss_timed.time,
    )
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
    checkpoint!(
        "stable_extraction_complete";
        implicit_response_seconds=response_timed.time,
        parameterized_twiss_seconds=twiss_timed.time,
        jacobian_extraction_seconds=extraction_timed.time,
    )
    physics_seconds = response_timed.time + twiss_timed.time + extraction_timed.time
    optics_only_seconds = twiss_timed.time + extraction_timed.time

    detector_path = joinpath(output_dir, "scibmad_detector_optics_jacobian.csv")
    orbit_path = joinpath(output_dir, "scibmad_closed_orbit_jacobian.csv")
    ring_path = joinpath(output_dir, "scibmad_ring_tune_jacobian.csv")
    metadata_path = joinpath(output_dir, "scibmad_extended_jacobian_metadata.toml")
    write_timed = @timed begin
        write_labeled_matrix(
            detector_path,
            extracted.detector.labels,
            setup.labels,
            extracted.detector_jacobian,
        )
        write_labeled_matrix(
            orbit_path,
            ["start:x", "start:px", "start:y", "start:py", "start:z", "start:pz"],
            setup.labels,
            response.response,
        )
        write_labeled_matrix(
            ring_path,
            ["ring:Q1", "ring:Q2", "ring:slip"],
            setup.labels,
            extracted.tune_jacobian,
        )
    end

    family_counts = Dict(
        family => count(==(family), setup.families)
        for family in ("corrector", "quadrupole", "sextupole")
    )
    metadata = Dict(
        "format" => "cesr-extended-optics-jacobian-v1",
        "engine" => "SciBmad/GTPSA",
        "case" => case,
        "method" => "one first-order-parameter Twiss plus implicit periodic closed-orbit response",
        "rf_mode" => "off (4D coasting)",
        "GTPSA_descriptor" => "Descriptor(6, $phase_space_order, $parameter_count, 1)",
        "phase_space_order" => phase_space_order,
        "parameter_count" => parameter_count,
        "family_counts" => family_counts,
        "parameters" => setup.labels,
        "detector_count" => 99,
        "detector_quantity_count" => length(extracted.detector.columns),
        "detector_jacobian_shape" => collect(size(extracted.detector_jacobian)),
        "baseline_closed_orbit_seconds" => baseline_seconds,
        "baseline_closed_orbit_residual" => baseline.residual,
        "parameterized_model_setup_seconds" => setup_timed.time,
        "parameterized_model_setup_allocated_bytes" => setup_timed.bytes,
        "implicit_response_warmup_seconds" => response_warmup_seconds,
        "twiss_warmup_seconds" => twiss_warmup_seconds,
        "extraction_warmup_seconds" => extraction_warmup_seconds,
        "implicit_response_seconds" => response_timed.time,
        "implicit_response_allocated_bytes" => response_timed.bytes,
        "implicit_response_max_abs_residual" => maximum(abs, response.residual),
        "parameterized_twiss_seconds" => twiss_timed.time,
        "parameterized_twiss_allocated_bytes" => twiss_timed.bytes,
        "jacobian_extraction_seconds" => extraction_timed.time,
        "optics_only_seconds" => optics_only_seconds,
        "physics_seconds_with_separate_start_orbit_response" => physics_seconds,
        "write_seconds" => write_timed.time,
        "julia_version" => string(VERSION),
        "julia_threads" => Threads.nthreads(),
    )
    mkpath(output_dir)
    open(metadata_path, "w") do io
        TOML.print(io, metadata; sorted=true)
    end
    checkpoint!("complete"; metadata_path)
    @printf(
        "SciBmad P=%d: optics %.3f s; separate orbit response %.3f s; total %.3f s\n",
        parameter_count,
        optics_only_seconds,
        response_timed.time,
        physics_seconds,
    )
    return metadata
end

function main_extended(args=ARGS)
    options = parse_extended_args(args)
    inputs = abspath(options["inputs"])
    output_root = abspath(options["output-root"])
    phase_space_order = parse(Int, options["phase-space-order"])
    run_warmup = lowercase(options["warmup"]) == "true"
    samples = read_samples(inputs)

    baseline_timed = @timed begin
        scalar_model = load_cesr_model(rf_on=false)
        solution = find_closed_orbit_coasting_forwarddiff(
            scalar_model.ring;
            coasting_beam=true,
        )
        (orbit=solution.orbit, residual=solution.residual)
    end
    baseline = baseline_timed.value
    cases = options["case"] == "both" ?
        ("correctors_quads", "correctors_quads_sextupoles") :
        (options["case"],)
    for case in cases
        run_extended_case(
            samples.names,
            baseline,
            case,
            joinpath(output_root, case, "scibmad"),
            phase_space_order,
            run_warmup,
            baseline_timed.time,
        )
    end
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main_extended())
end

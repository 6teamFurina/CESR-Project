#!/usr/bin/env julia

"""Four-corner K2/wiggler ablation for the vertical detector-y cubic vector."""

include(joinpath(@__DIR__, "run_sextupole_cascade_experiment.jl"))

function parse_corner_args(args)
    options = Dict{String,String}(
        "ring" => "latest",
        "rhos" => "0.4,0.57,0.8,1.13,1.6,2.26",
        "trials" => "100",
        "seed" => "20260804",
        "base-kick-rad" => "5e-6",
        "output-dir" => "",
        "inputs" => "",
        "reltol" => "1e-12",
        "abstol" => "1e-13",
        "maxiter" => "100",
    )
    for argument in args
        startswith(argument, "--") || error("Arguments must have --name=value form: $argument")
        fields = split(argument[3:end], '='; limit=2)
        length(fields) == 2 || error("Missing value in argument: $argument")
        haskey(options, fields[1]) || error("Unknown option: --$(fields[1])")
        options[fields[1]] = fields[2]
    end
    ring = Symbol(lowercase(options["ring"]))
    ring in (:latest, :latest_cesr, :repaired_latest, :legacy, :legacy_cesr, :historical) ||
        error("--ring must be latest or legacy")
    isempty(options["inputs"]) &&
        (options["inputs"] = default_ring_paths(; ring).inputs)
    isempty(options["output-dir"]) &&
        (options["output-dir"] = joinpath(
            @__DIR__, "wiggler_corner_results", ring_artifact_id(ring),
        ))
    return options
end

function write_corner_vectors(path, variants, detectors)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "lambda2,wiggler_scale,trial,detector,c3_y_m,c5_y_m")
        for variant in variants, trial in axes(variant.c3_vectors, 1), detector in eachindex(detectors)
            println(io, join((
                csv_number(variant.lambda2), csv_number(variant.wiggler_scale),
                trial, detectors[detector],
                csv_number(variant.c3_vectors[trial, detector]),
                csv_number(variant.c5_vectors[trial, detector]),
            ), ','))
        end
    end
    return path
end

function main_corner(args=ARGS; model_factory=nothing, config=nothing)
    options = parse_corner_args(args)
    ring = Symbol(lowercase(options["ring"]))
    model_factory, config = resolve_ring_model_factory(model_factory, config; ring)
    rhos = parse_sorted_unique(options["rhos"], "--rhos")
    all(>(0), rhos) || error("Every rho must be positive")
    trials = parse(Int, options["trials"])
    trials >= 2 || error("--trials must be at least 2")
    seed = parse(Int, options["seed"])
    base_kick = parse(Float64, options["base-kick-rad"])
    reltol = parse(Float64, options["reltol"])
    abstol = parse(Float64, options["abstol"])
    maxiter = parse(Int, options["maxiter"])
    output_dir = abspath(options["output-dir"])
    mkpath(output_dir)

    input_path = abspath(options["inputs"])
    input_reference = read_samples(input_path)
    names = input_reference.names
    validate_control_names(names, config)
    base_model = configured_model(model_factory, config; zero_value=0.0, rf_on=true)
    samples = generate_vertical_pairs(names, rhos, trials, seed, base_kick; config, model=base_model)
    corners = ((0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (1.0, 1.0))
    variants = NamedTuple[]
    detectors = String[]
    layout = nothing
    maximum_closure = 0.0
    total_fallback = 0
    total_seconds = 0.0

    @printf("K2/wiggler corner experiment: 4 variants x %d states\n", size(samples.values, 1))
    for (lambda2, wiggler_scale) in corners
        factory = ScaledCESRFactory(
            lambda2,
            wiggler_scale;
            base_factory=model_factory,
            config,
        )
        @printf("K2=%.0f, wiggler=%.0f: response map...\n", lambda2, wiggler_scale)
        linear = nominal_and_detector_response(factory, names; reltol, abstol, maxiter)
        if isempty(detectors)
            detectors = linear.detectors
            layout = linear.layout
        else
            detectors == linear.detectors || error("Detector order changed")
        end
        cache = joinpath(
            output_dir,
            "closed_response_k2_$(lambda_token(lambda2))_wig_$(lambda_token(wiggler_scale)).csv",
        )
        timed = @timed configured_simulate_batch(
            simulate_batch,
            names, samples.values;
            config,
            initial_guess_mode="response-linear",
            jacobian_mode="frozen-nominal",
            response_matrix_cache=cache,
            recompute_response=true,
            reltol, abstol, maxiter,
            model_factory=factory,
        )
        result = timed.value
        all(result.converged) || error("A corner variant has nonconverged states")
        extracted = extract_variant(
            lambda2, result, linear.detector_response, samples, rhos, base_kick, linear.layout,
        )
        push!(variants, (; lambda2, wiggler_scale, extracted...))
        maximum_closure = max(maximum_closure, maximum(result.closure_norms))
        total_fallback += result.fallback_count
        total_seconds += timed.time
        @printf(
            "K2=%.0f, wiggler=%.0f complete: fallback=%d closure=%.3e\n",
            lambda2, wiggler_scale, result.fallback_count, maximum(result.closure_norms),
        )
    end

    vector_path = write_corner_vectors(
        joinpath(output_dir, "wiggler_corner_vectors.csv"), variants, detectors,
    )
    pair_rows = NamedTuple[]
    coefficient_rows = NamedTuple[]
    for variant in variants
        append!(pair_rows, [merge((wiggler_scale=variant.wiggler_scale,), row) for row in variant.pair_rows])
        append!(coefficient_rows, [merge((wiggler_scale=variant.wiggler_scale,), row) for row in variant.coefficient_rows])
    end
    pair_path = write_namedtuple_csv(joinpath(output_dir, "wiggler_corner_pairs.csv"), pair_rows)
    coefficient_path = write_namedtuple_csv(
        joinpath(output_dir, "wiggler_corner_coefficients.csv"), coefficient_rows,
    )
    metadata = Dict(
        "format" => "cesr-k2-wiggler-cubic-corners-v1",
        "date" => string(Dates.today()),
        "rhos" => rhos,
        "trials_per_corner" => trials,
        "seed" => seed,
        "base_kick_rad" => base_kick,
        "input_csv" => input_path,
        "control_count" => length(names),
        "control_names" => names,
        "detector_count" => length(detectors),
        "observable_count" => isnothing(layout) ? 0 : length(layout.labels),
        "total_nonlinear_states" => length(corners) * size(samples.values, 1),
        "total_solve_seconds" => total_seconds,
        "total_fallback_count" => total_fallback,
        "maximum_closure_norm" => maximum_closure,
        "vectors_csv" => vector_path,
        "pairs_csv" => pair_path,
        "coefficients_csv" => coefficient_path,
    )
    merge!(metadata, ring_metadata(config; ring))
    metadata_path = joinpath(output_dir, "wiggler_corner_metadata.toml")
    open(metadata_path, "w") do io
        TOML.print(io, metadata; sorted=true)
    end
    println("Vectors:  $vector_path")
    println("Metadata: $metadata_path")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main_corner())
end

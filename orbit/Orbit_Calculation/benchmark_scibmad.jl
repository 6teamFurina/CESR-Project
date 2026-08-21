#!/usr/bin/env julia

"""
Generate exact CESR closed-orbit samples with SciBmad batch parameters.

The timed physics region includes assigning all sampled control arrays, solving
all closed orbits, and tracking the solved orbits to the configured detector
registry. Model loading, Julia compilation warmup, and CSV writing are
reported separately.
"""

using Beamlines
using GTPSA
using LinearAlgebra
using Printf
using SciBmad
using Statistics
using TOML

const HERE = @__DIR__
const ORBIT_ROOT = normpath(joinpath(HERE, ".."))
const PROJECT_ROOT = normpath(joinpath(HERE, "..", ".."))
# SciBmad's Bunch state is the six-component accelerator phase-space vector
# (x, px, y, py, z, pz).  Control and detector dimensions are ring metadata;
# this is the one physical interface dimension retained by the runner.
const ORBIT_PHASE_SPACE_DIMENSION = 6
const ORBIT_PHASE_SPACE_LABELS = ("x", "px", "y", "py", "z", "pz")
const DEFAULT_RESPONSE_STEP_RAD = 1.0e-7
const DEFAULT_RESPONSE_CONTROLS_PER_BATCH = 8
const DEFAULT_RESPONSE_METHOD = "gtpsa"

function canonical_response_method(method::AbstractString)
    value = lowercase(replace(strip(String(method)), '_' => '-'))
    value in ("central-difference", "batch-central-finite-difference") &&
        return "central-difference"
    value in ("gtpsa", "gtpsa-implicit", "gtpsa-implicit-closed-orbit") &&
        return "gtpsa"
    error("Unsupported response method '$method'; use gtpsa or central-difference")
end

if !isdefined(@__MODULE__, :load_ring_model)
    include(joinpath(HERE, "ring_model_adapter.jl"))
end

function parse_args(args)
    options = Dict{String,String}(
        # Resolve these after parsing --ring so an explicit legacy run cannot
        # write into the maintained latest-ring tree.
        "inputs" => "",
        "output" => "",
        "metadata" => "",
        "ring" => "latest",
        "mode" => "rf_on",
        "reltol" => "1e-8",
        "abstol" => "1e-10",
        "maxiter" => "100",
        "warmup-samples" => "2",
        "initial-guess" => "response-linear",
        "jacobian-mode" => "frozen-nominal",
        "response-matrix-cache" => "",
        "recompute-response" => "false",
        "response-method" => DEFAULT_RESPONSE_METHOD,
        "response-step-rad" => string(DEFAULT_RESPONSE_STEP_RAD),
        "response-controls-per-batch" =>
            string(DEFAULT_RESPONSE_CONTROLS_PER_BATCH),
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
    options["ring"] in ("latest", "latest_cesr", "repaired_latest", "legacy", "legacy_cesr", "historical") ||
        error("--ring must be latest or explicitly legacy")
    ring = Symbol(options["ring"])
    artifact_ring = String(canonical_ring_id(ring))
    if isempty(options["inputs"])
        options["inputs"] = artifact_ring == "latest_cesr" ?
            joinpath(HERE, "inputs", "latest_cesr", "corrector_samples.csv") :
            joinpath(HERE, "inputs", "cesr_corrector_samples_1000.csv")
    end
    if isempty(options["output"])
        options["output"] = joinpath(
            HERE,
            "results",
            artifact_ring,
            "formal_1000",
            "scibmad_response_initial_frozen_fallback_bmad_tolerance",
            "scibmad_rf_on_samples.csv",
        )
    end
    if isempty(options["metadata"])
        options["metadata"] = joinpath(
            HERE,
            "results",
            artifact_ring,
            "formal_1000",
            "scibmad_response_initial_frozen_fallback_bmad_tolerance",
            "scibmad_rf_on_metadata.toml",
        )
    end
    options["initial-guess"] in ("zero", "nominal-z0", "response-linear") ||
        error("--initial-guess must be zero, nominal-z0, or response-linear")
    options["jacobian-mode"] in ("full", "frozen-nominal") ||
        error("--jacobian-mode must be full or frozen-nominal")
    if options["jacobian-mode"] == "frozen-nominal"
        options["initial-guess"] in ("nominal-z0", "response-linear") ||
            error("--jacobian-mode=frozen-nominal requires --initial-guess=nominal-z0 or response-linear")
    end
    lowercase(options["recompute-response"]) in ("true", "false") ||
        error("--recompute-response must be true or false")
    options["response-method"] = canonical_response_method(
        options["response-method"],
    )
    parse(Float64, options["response-step-rad"]) > 0 ||
        error("--response-step-rad must be positive for central-difference validation")
    parse(Int, options["response-controls-per-batch"]) >= 1 ||
        error("--response-controls-per-batch must be at least 1")
    return options
end

function read_samples(path::AbstractString)
    lines = readlines(path)
    length(lines) >= 2 || error("Sample CSV has no data: $path")
    header = split(lines[1], ',')
    first(header) == "sample_id" || error("First CSV column must be sample_id")
    names = String.(header[2:end])
    isempty(names) && error("Sample CSV must contain at least one control")
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

function detector_registry(model)
    metadata = hasproperty(model, :metadata) ? model.metadata : nothing
    if !isnothing(metadata) &&
       hasproperty(metadata, :detector_names) &&
       hasproperty(metadata, :detector_element_indices)
        return (
            String.(collect(metadata.detector_names)),
            Int.(collect(metadata.detector_element_indices)),
        )
    end
    ring = hasproperty(model, :ring) ? model.ring : model
    names = String[]
    indices = Int[]
    for (element_index, element) in enumerate(ring.line)
        name = uppercase(String(element.name))
        startswith(name, "DET_") || continue
        push!(names, name)
        push!(indices, element_index)
    end
    return names, indices
end

detector_names(model) = first(detector_registry(model))

function orbit_coordinate_count(model)
    metadata = hasproperty(model, :metadata) ? model.metadata : nothing
    if !isnothing(metadata) && hasproperty(metadata, :phase_space_dimension)
        return Int(metadata.phase_space_dimension)
    end
    return ORBIT_PHASE_SPACE_DIMENSION
end

function orbit_coordinate_labels(model, count::Int)
    metadata = hasproperty(model, :metadata) ? model.metadata : nothing
    labels = if !isnothing(metadata) && hasproperty(metadata, :coordinate_labels)
        String.(collect(metadata.coordinate_labels))
    else
        collect(ORBIT_PHASE_SPACE_LABELS)
    end
    length(labels) >= count || error("Model provides fewer coordinate labels than its orbit dimension")
    return labels[1:count]
end

function orbit_transverse_coordinate_indices(model, coordinate_count::Int)
    metadata = hasproperty(model, :metadata) ? model.metadata : nothing
    indices = if !isnothing(metadata) && hasproperty(metadata, :transverse_coordinate_indices)
        metadata.transverse_coordinate_indices
    else
        (x=1, y=3)
    end
    x_index = Int(getproperty(indices, :x))
    y_index = Int(getproperty(indices, :y))
    all(index -> 1 <= index <= coordinate_count, (x_index, y_index)) ||
        error("Transverse coordinate indices ($x_index, $y_index) are outside the $coordinate_count-component orbit state")
    return (; x=x_index, y=y_index)
end

function prepare_batch_model(
    names::Vector{String},
    values::Matrix{Float64},
    ;
    model_factory=load_ring_model,
)
    n_samples, n_controls = size(values)
    n_controls == length(names) || error("Control matrix width does not match labels")

    model = model_factory(zero_value=BatchParam(0.0), rf_on=true)
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
    detectors, detector_indices = detector_registry(model)
    coordinate_count = orbit_coordinate_count(model)
    transverse_indices = orbit_transverse_coordinate_indices(model, coordinate_count)
    v0 = isnothing(initial_v0) ? zeros(n_samples, coordinate_count) : copy(initial_v0)
    size(v0) == (n_samples, coordinate_count) ||
        error("Initial closed-orbit guess must have size ($n_samples, $coordinate_count)")
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
        for (element_index, element) in enumerate(model.ring.line)
            track!(bunch, element)
            name = uppercase(String(element.name))
            element_index in detector_indices || continue
            detector_index += 1
            name == detectors[detector_index] ||
                error("Detector order changed at $name")
            horizontal[:, detector_index] .= bunch.coords.v[:, transverse_indices.x]
            vertical[:, detector_index] .= bunch.coords.v[:, transverse_indices.y]
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

"""Compute first-order closed-orbit and detector responses with GTPSA.

The one-turn map is differentiated with respect to the six phase-space
coordinates and the selected steering controls, then the closed-orbit
derivative is obtained from the implicit fixed-point equation
`(I - A) * dz = B * dk`.  The model is initialized with primitive Float64
zeros for the complete control registry; only `names` are assigned GTPSA
parameters.  This selected-control construction is important for the latest
ring because unused Overlay/Group globals should not be promoted to GTPSA
values.  The method is the default first-order response backend.

Signed central differences remain available through the explicit
`response_method="central-difference"` option for validation or for a caller
whose requested observable cannot yet propagate GTPSA values.

For the repaired CESR export this GTPSA backend intentionally accepts the
default normal H/V steering subset.  Skew and Group controls are retained in
the registry, but are rejected here with a diagnostic rather than being
silently promoted into an all-control GTPSA map.
"""
function gtpsa_first_order_responses(
    names::Vector{String};
    nominal_orbit::Union{Nothing,AbstractVector}=nothing,
    rf_on::Bool=true,
    reltol::Float64=1.0e-8,
    abstol::Float64=1.0e-10,
    maxiter::Int=100,
    model_factory=load_ring_model,
)
    n_controls = length(names)
    n_controls > 0 || error("Cannot compute a response with no controls")
    length(unique(names)) == n_controls ||
        error("GTPSA response control names are not unique")

    nominal_model_setup_seconds = 0.0
    nominal_solve_seconds = 0.0
    nominal = nothing
    coordinate_count = ORBIT_PHASE_SPACE_DIMENSION
    if isnothing(nominal_orbit)
        nominal_model_setup_seconds = @elapsed nominal_model = model_factory(
            zero_value=0.0,
            rf_on=rf_on,
        )
        coordinate_count = orbit_coordinate_count(nominal_model)
        coordinate_count == ORBIT_PHASE_SPACE_DIMENSION || error(
            "SciBmad GTPSA response requires the $ORBIT_PHASE_SPACE_DIMENSION-component " *
            "Bunch state; model reports $coordinate_count",
        )
        nominal_solve_seconds = @elapsed nominal_solution = find_closed_orbit(
            nominal_model.ring;
            v0=zeros(1, coordinate_count),
            coasting_beam=false,
            batch=Val{false}(),
            reltol,
            abstol,
            maxiter,
            warn=false,
        )
        nominal_solution.sol.retcode == SciBmad.BatchSolve.RETCODE_SUCCESS ||
            error("Nominal closed orbit did not converge while computing the GTPSA response")
        nominal = vec(Float64.(copy(nominal_solution.v0)))
    else
        nominal = vec(Float64.(collect(nominal_orbit)))
        coordinate_count = length(nominal)
    end
    coordinate_count == ORBIT_PHASE_SPACE_DIMENSION || error(
        "SciBmad GTPSA response requires the $ORBIT_PHASE_SPACE_DIMENSION-component " *
        "Bunch state; model reports $coordinate_count",
    )
    all(isfinite, nominal) || error("Nominal orbit contains a non-finite value")

    descriptor = Descriptor(coordinate_count, 1, n_controls, 1)
    variables = vars(descriptor)
    parameters = params(descriptor)
    response_model = nothing
    response_model_setup_seconds = @elapsed begin
        response_model = parameterized_ring_model(
            names,
            parameters;
            model_factory,
            zero_value=0.0,
            rf_on,
        )
    end
    if hasproperty(response_model.metadata, :control_plane)
        unsupported = [
            name for name in names if
            get(response_model.metadata.control_plane, name, :other)
                ∉ (:horizontal, :vertical)
        ]
        isempty(unsupported) || error(
            "GTPSA response currently supports only normal H/V steering controls; " *
            "incompatible selected controls: $(join(unsupported, ", ")). " *
            "Use --response-method=central-difference for this explicit subset.",
        )
    end

    input_map = [
        nominal[index] + copy(variables[index])
        for index in 1:coordinate_count
    ]
    map_bunch = Bunch(v=reshape(input_map, 1, coordinate_count))
    SciBmad.BTBL.check_bl_bunch!(map_bunch, response_model.ring, false)
    map_seconds = @elapsed track!(map_bunch, response_model.ring)
    output_map = vec(map_bunch.coords.v)
    full_jacobian = Matrix(GTPSA.jacobian(output_map; include_params=true))
    size(full_jacobian) == (coordinate_count, coordinate_count + n_controls) ||
        error("Unexpected GTPSA one-turn Jacobian size: $(size(full_jacobian)); expected " *
              "($coordinate_count, $(coordinate_count + n_controls))")
    all(isfinite, full_jacobian) ||
        error("GTPSA one-turn Jacobian contains a non-finite value")

    A = full_jacobian[:, 1:coordinate_count]
    B = full_jacobian[:, coordinate_count + 1:end]
    closed_orbit_response = zeros(Float64, coordinate_count, n_controls)
    if rf_on
        closed_orbit_response .= (I - A) \ B
        closure_residual = (I - A) * closed_orbit_response - B
    else
        active = 1:min(4, coordinate_count)
        closed_orbit_response[active, :] .=
            (I - A[active, active]) \ B[active, :]
        closure_residual = (I - A[active, active]) *
            closed_orbit_response[active, :] - B[active, :]
    end
    all(isfinite, closed_orbit_response) ||
        error("GTPSA closed-orbit response contains a non-finite value")

    detectors, detector_indices = detector_registry(response_model)
    transverse_indices = orbit_transverse_coordinate_indices(
        response_model,
        coordinate_count,
    )
    detector_response = zeros(Float64, 2 * length(detectors), n_controls)
    detector_bunch = Bunch(v=reshape([
        nominal[index] + sum(
            closed_orbit_response[index, control] * parameters[control]
            for control in 1:n_controls
        )
        for index in 1:coordinate_count
    ], 1, coordinate_count))
    SciBmad.BTBL.check_bl_bunch!(detector_bunch, response_model.ring, false)
    detector_index = 0
    detector_track_seconds = @elapsed begin
        for (element_index, element) in enumerate(response_model.ring.line)
            track!(detector_bunch, element)
            element_index in detector_indices || continue
            detector_index += 1
            name = uppercase(String(element.name))
            name == detectors[detector_index] || error(
                "Detector order changed at $name; expected $(detectors[detector_index])",
            )
            jacobian = Matrix(GTPSA.jacobian(
                vec(detector_bunch.coords.v);
                include_params=true,
            ))
            size(jacobian) == (coordinate_count, coordinate_count + n_controls) ||
                error("Unexpected detector GTPSA Jacobian size at $name: $(size(jacobian))")
            detector_response[detector_index, :] .=
                jacobian[transverse_indices.x, coordinate_count + 1:end]
            detector_response[length(detectors) + detector_index, :] .=
                jacobian[transverse_indices.y, coordinate_count + 1:end]
        end
    end
    detector_index == length(detectors) || error(
        "Not all configured detectors were tracked ($(detector_index)/$(length(detectors)))",
    )
    all(isfinite, detector_response) ||
        error("GTPSA detector response contains a non-finite value")
    observable_labels = if hasproperty(response_model.metadata, :observable_labels)
        String.(collect(response_model.metadata.observable_labels))
    else
        vcat(
            [name * ":x" for name in detectors],
            [name * ":y" for name in detectors],
        )
    end
    length(observable_labels) == size(detector_response, 1) || error(
        "Detector response labels do not match response rows",
    )
    return (;
        descriptor,
        variables,
        parameters,
        model=response_model,
        A,
        B,
        output_map,
        closed_orbit_response,
        detector_response,
        detectors,
        observable_labels,
        closure_residual,
        closure_norm_max=maximum(abs, closure_residual),
        response_method="gtpsa",
        response_step_rad=0.0,
        controls_per_batch=0,
        chunk_count=1,
        model_setup_seconds=response_model_setup_seconds,
        map_seconds,
        track_seconds=map_seconds + detector_track_seconds,
        total_seconds=map_seconds + detector_track_seconds +
            response_model_setup_seconds,
        nominal_model_setup_seconds,
        nominal_solve_seconds,
        iteration_max=0,
    )
end

"""Compute closed-orbit and detector responses with explicit central differences.

This is an independent validation/fallback backend.  It is not the default:
the production first-order response path uses `gtpsa_first_order_responses`.
"""
function central_finite_difference_responses(
    names::Vector{String};
    nominal_orbit::Union{Nothing,AbstractVector}=nothing,
    response_step_rad::Float64=DEFAULT_RESPONSE_STEP_RAD,
    controls_per_batch::Int=DEFAULT_RESPONSE_CONTROLS_PER_BATCH,
    reltol::Float64,
    abstol::Float64,
    maxiter::Int,
    model_factory=load_ring_model,
)
    n_controls = length(names)
    n_controls > 0 || error("Cannot compute a response with no controls")
    isfinite(response_step_rad) && response_step_rad > 0 ||
        error("The response finite-difference step must be finite and positive")
    controls_per_batch >= 1 ||
        error("The response controls-per-batch value must be at least 1")

    nominal_model_setup_seconds = 0.0
    nominal_solve_seconds = 0.0
    if isnothing(nominal_orbit)
        nominal_model_setup_seconds = @elapsed nominal_model = model_factory(
            zero_value=0.0,
            rf_on=true,
        )
        coordinate_count = orbit_coordinate_count(nominal_model)
        nominal_solve_seconds = @elapsed nominal_solution = find_closed_orbit(
            nominal_model.ring;
            v0=zeros(1, coordinate_count),
            coasting_beam=false,
            batch=Val{false}(),
            reltol,
            abstol,
            maxiter,
            warn=false,
        )
        nominal_solution.sol.retcode == SciBmad.BatchSolve.RETCODE_SUCCESS ||
            error("Nominal closed orbit did not converge while computing the response")
        nominal = vec(Float64.(copy(nominal_solution.v0)))
    else
        nominal = vec(Float64.(collect(nominal_orbit)))
        coordinate_count = length(nominal)
    end
    coordinate_count == ORBIT_PHASE_SPACE_DIMENSION ||
        error("SciBmad response generation requires the $ORBIT_PHASE_SPACE_DIMENSION-component Bunch state; model reports $coordinate_count")

    closed_orbit_response = zeros(coordinate_count, n_controls)
    detector_response = Matrix{Float64}(undef, 0, n_controls)
    response_detectors = String[]
    model_setup_seconds = 0.0
    solve_seconds = 0.0
    track_seconds = 0.0
    closure_norm_max = 0.0
    iteration_max = 0
    chunk_count = 0
    total_seconds = @elapsed begin
        for first_control in 1:controls_per_batch:n_controls
            last_control = min(first_control + controls_per_batch - 1, n_controls)
            chunk_controls = first_control:last_control
            lanes = 2 * length(chunk_controls)
            perturbations = zeros(lanes, n_controls)
            for (local_control, global_control) in enumerate(chunk_controls)
                perturbations[2 * local_control - 1, global_control] = response_step_rad
                perturbations[2 * local_control, global_control] = -response_step_rad
            end

            setup = @elapsed model = prepare_batch_model(
                names,
                perturbations;
                model_factory,
            )
            initial_v0 = repeat(reshape(nominal, 1, coordinate_count), lanes, 1)
            result = solve_and_track(
                model,
                lanes;
                initial_v0,
                reltol,
                abstol,
                maxiter,
            )
            all(result.converged) || error(
                "Only $(count(result.converged))/$lanes finite-difference lanes converged for controls $first_control:$last_control",
            )
            all(isfinite, result.closure_norms) ||
                error("A finite-difference response lane has a non-finite closure norm")
            maximum(result.closure_norms) <= abstol || error(
                "Finite-difference response closure norm $(maximum(result.closure_norms)) exceeds abstol=$abstol",
            )

            if isempty(response_detectors)
                response_detectors = copy(result.detectors)
                detector_response = zeros(
                    size(result.observables, 2),
                    n_controls,
                )
            else
                result.detectors == response_detectors ||
                    error("Detector registry changed between response chunks")
                size(result.observables, 2) == size(detector_response, 1) ||
                    error("Observable count changed between response chunks")
            end

            inverse_span = inv(2 * response_step_rad)
            for (local_control, global_control) in enumerate(chunk_controls)
                plus_lane = 2 * local_control - 1
                minus_lane = 2 * local_control
                closed_orbit_response[:, global_control] .=
                    (view(result.final_v0, plus_lane, :) .-
                     view(result.final_v0, minus_lane, :)) .* inverse_span
                detector_response[:, global_control] .=
                    (view(result.observables, plus_lane, :) .-
                     view(result.observables, minus_lane, :)) .* inverse_span
            end
            model_setup_seconds += setup
            solve_seconds += result.solve_seconds
            track_seconds += result.track_seconds
            closure_norm_max = max(
                closure_norm_max,
                maximum(result.closure_norms),
            )
            iteration_max = max(iteration_max, maximum(result.iterations))
            chunk_count += 1
        end
    end

    all(isfinite, closed_orbit_response) ||
        error("Closed-orbit finite-difference response contains a non-finite value")
    all(isfinite, detector_response) ||
        error("Detector finite-difference response contains a non-finite value")
    observable_labels = vcat(
        [name * ":x" for name in response_detectors],
        [name * ":y" for name in response_detectors],
    )
    return (;
        closed_orbit_response,
        detector_response,
        detectors=response_detectors,
        observable_labels,
        response_method="central-difference",
        response_step_rad,
        controls_per_batch,
        chunk_count,
        model_setup_seconds,
        solve_seconds,
        track_seconds,
        total_seconds,
        closure_norm_max,
        iteration_max,
        nominal_model_setup_seconds,
        nominal_solve_seconds,
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
    detectors, detector_indices = detector_registry(model)
    coordinate_count = orbit_coordinate_count(model)
    transverse_indices = orbit_transverse_coordinate_indices(model, coordinate_count)
    size(initial_v0) == (n_samples, coordinate_count) ||
        error("Initial closed-orbit guess must have size ($n_samples, $coordinate_count)")
    size(frozen_jacobian) == (coordinate_count, coordinate_count) ||
        error("Frozen closed-orbit Jacobian must have size ($coordinate_count, $coordinate_count)")

    v = copy(initial_v0)
    residual = similar(v)
    v_cache = similar(v)
    rhs = zeros(eltype(v), coordinate_count, n_samples)
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
        for (element_index, element) in enumerate(model.ring.line)
            track!(bunch, element)
            name = uppercase(String(element.name))
            element_index in detector_indices || continue
            detector_index += 1
            name == detectors[detector_index] ||
                error("Detector order changed at $name")
            horizontal[:, detector_index] .= bunch.coords.v[:, transverse_indices.x]
            vertical[:, detector_index] .= bunch.coords.v[:, transverse_indices.y]
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
    model_factory=load_ring_model,
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
            ;
            model_factory,
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
    values::Matrix{Float64},
    mode::String;
    response_matrix_cache::AbstractString,
    recompute_response::Bool,
    reltol::Float64,
    abstol::Float64,
    maxiter::Int,
    response_step_rad::Float64=DEFAULT_RESPONSE_STEP_RAD,
    response_controls_per_batch::Int=DEFAULT_RESPONSE_CONTROLS_PER_BATCH,
    response_method::String=DEFAULT_RESPONSE_METHOD,
    model_factory=load_ring_model,
)
    n_samples, n_controls = size(values)
    n_controls == length(names) || error("Control matrix width does not match labels")
    requested_response_method = canonical_response_method(response_method)
    coordinate_count = ORBIT_PHASE_SPACE_DIMENSION
    if mode == "zero"
        return (;
            v0=zeros(n_samples, coordinate_count),
            nominal_orbit=zeros(coordinate_count),
            nominal_jacobian=zeros(coordinate_count, coordinate_count),
            nominal_model_setup_seconds=0.0,
            nominal_solve_seconds=0.0,
            nominal_iterations=0,
            response_matrix=zeros(coordinate_count, n_controls),
            response_model_setup_seconds=0.0,
            response_map_seconds=0.0,
            response_load_seconds=0.0,
            response_cache_write_seconds=0.0,
            response_closure_residual_max=0.0,
            detector_response=zeros(0, n_controls),
            detector_response_labels=String[],
            response_detectors=String[],
            response_source="not-used",
            response_method="not-used",
            response_step_rad=0.0,
            response_controls_per_batch=0,
            response_chunk_count=0,
        )
    end

    mode in ("nominal-z0", "response-linear") ||
        error("Unsupported initial-guess mode: $mode")
    setup_seconds = @elapsed nominal_model = model_factory(
        zero_value=0.0,
        rf_on=true,
    )
    coordinate_count = orbit_coordinate_count(nominal_model)
    coordinate_count == ORBIT_PHASE_SPACE_DIMENSION ||
        error("SciBmad closed-orbit runner requires the $ORBIT_PHASE_SPACE_DIMENSION-component Bunch state; model reports $coordinate_count")
    solve_seconds = @elapsed nominal_solution = find_closed_orbit(
        nominal_model.ring;
        v0=zeros(1, coordinate_count),
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

    response_matrix = zeros(length(nominal_orbit), n_controls)
    response_model_setup_seconds = 0.0
    response_map_seconds = 0.0
    response_load_seconds = 0.0
    response_cache_write_seconds = 0.0
    response_closure_residual_max = 0.0
    detector_response = Matrix{Float64}(undef, 0, n_controls)
    detector_response_labels = String[]
    response_detectors = String[]
    response_source = "not-used"
    actual_response_method = "not-used"
    effective_response_step_rad = 0.0
    effective_response_controls_per_batch = 0
    response_chunk_count = 0
    if mode == "response-linear"
        if isfile(response_matrix_cache) && !recompute_response
            response_load_seconds = @elapsed begin
                response_matrix .= read_response_matrix(
                    response_matrix_cache,
                    names,
                )
            end
            response_source = "loaded"
            cache_metadata = read_response_cache_metadata(
                response_matrix_cache,
                names,
                response_matrix;
                model=nominal_model,
                requested_response_method=requested_response_method,
                requested_response_step_rad=response_step_rad,
                requested_response_controls_per_batch=response_controls_per_batch,
                required_reltol=reltol,
                required_abstol=abstol,
                required_maxiter=maxiter,
            )
            if isnothing(cache_metadata)
                selected_ring_id = hasproperty(nominal_model, :metadata) &&
                    hasproperty(nominal_model.metadata, :ring_id) ?
                    lowercase(String(nominal_model.metadata.ring_id)) : ""
                startswith(selected_ring_id, "legacy") || error(
                    "Latest/custom response cache is missing its provenance sidecar; use --recompute-response=true: $(response_cache_metadata_path(response_matrix_cache))",
                )
                actual_response_method = "cached-labeled-csv-unknown-generator"
            else
                actual_response_method = String(get(
                    cache_metadata,
                    "response_method",
                    "cached-labeled-csv-unknown-generator",
                ))
                effective_response_step_rad = Float64(get(
                    cache_metadata,
                    "response_step_rad",
                    0.0,
                ))
                effective_response_controls_per_batch = Int(get(
                    cache_metadata,
                    "controls_per_batch",
                    0,
                ))
                response_chunk_count = Int(get(
                    cache_metadata,
                    "chunk_count",
                    0,
                ))
                response_closure_residual_max = Float64(get(
                    cache_metadata,
                    "closure_norm_max",
                    0.0,
                ))
            end
        else
            response = if requested_response_method == "gtpsa"
                gtpsa_first_order_responses(
                    names;
                    nominal_orbit,
                    rf_on=true,
                    reltol,
                    abstol,
                    maxiter,
                    model_factory,
                )
            else
                central_finite_difference_responses(
                    names;
                    nominal_orbit,
                    response_step_rad,
                    controls_per_batch=response_controls_per_batch,
                    reltol,
                    abstol,
                    maxiter,
                    model_factory,
                )
            end
            response_matrix .= response.closed_orbit_response
            detector_response = response.detector_response
            detector_response_labels = response.observable_labels
            response_detectors = response.detectors
            response_model_setup_seconds = response.model_setup_seconds
            response_map_seconds = response.total_seconds
            response_closure_residual_max = response.closure_norm_max
            effective_response_step_rad = response.response_step_rad
            effective_response_controls_per_batch = response.controls_per_batch
            response_chunk_count = response.chunk_count
            response_cache_write_seconds = @elapsed begin
                write_response_matrix(
                    response_matrix_cache,
                    names,
                    response_matrix,
                )
                write_response_cache_metadata(
                    response_matrix_cache,
                    names,
                    response_matrix;
                    method=String(response.response_method),
                    response_step_rad=effective_response_step_rad,
                    controls_per_batch=effective_response_controls_per_batch,
                    chunk_count=response_chunk_count,
                    closure_norm_max=response_closure_residual_max,
                    reltol,
                    abstol,
                    maxiter,
                    model=nominal_model,
                )
            end
            response_source = "computed"
            actual_response_method = String(response.response_method)
        end
    end

    v0 = repeat(reshape(nominal_orbit, 1, coordinate_count), n_samples, 1)
    if mode == "response-linear"
        v0 .+= values * transpose(response_matrix)
    end
    return (;
        v0,
        nominal_orbit,
        nominal_jacobian,
        nominal_model_setup_seconds=setup_seconds,
        nominal_solve_seconds=solve_seconds,
        nominal_iterations,
        response_matrix,
        response_model_setup_seconds,
        response_map_seconds,
        response_load_seconds,
        response_cache_write_seconds,
        response_closure_residual_max,
        detector_response,
        detector_response_labels,
        response_detectors,
        response_source,
        response_method=actual_response_method,
        response_step_rad=effective_response_step_rad,
        response_controls_per_batch=effective_response_controls_per_batch,
        response_chunk_count,
    )
end

function simulate_batch(
    names::Vector{String},
    values::Matrix{Float64};
    initial_guess_mode::String,
    jacobian_mode::String,
    response_matrix_cache::AbstractString,
    recompute_response::Bool,
    reltol::Float64,
    abstol::Float64,
    maxiter::Int,
    response_step_rad::Float64=DEFAULT_RESPONSE_STEP_RAD,
    response_controls_per_batch::Int=DEFAULT_RESPONSE_CONTROLS_PER_BATCH,
    response_method::String=DEFAULT_RESPONSE_METHOD,
    model_factory=load_ring_model,
)
    guess = prepare_initial_guess(
        names,
        values,
        initial_guess_mode;
        response_matrix_cache,
        recompute_response,
        reltol,
        abstol,
        maxiter,
        response_step_rad,
        response_controls_per_batch,
        response_method,
        model_factory,
    )
    model_setup_seconds = @elapsed model = prepare_batch_model(
        names,
        values;
        model_factory,
    )
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
    if jacobian_mode == "frozen-nominal"
        result = apply_full_newton_fallback(
            result,
            names,
            values;
            reltol,
            abstol,
            maxiter,
            model_factory,
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

function write_response_matrix(path, names, response_matrix)
    size(response_matrix, 2) == length(names) ||
        error("Closed-orbit response matrix has an unexpected size")
    mkpath(dirname(path))
    coordinate_labels = collect(ORBIT_PHASE_SPACE_LABELS)
    size(response_matrix, 1) <= length(coordinate_labels) ||
        error("No coordinate labels are available for response rows")
    open(path, "w") do io
        println(io, join(vcat("coordinate", names), ','))
        for row in axes(response_matrix, 1)
            fields = Any[coordinate_labels[row]]
            append!(fields, response_matrix[row, :])
            println(io, join(fields, ','))
        end
    end
    return path
end

"""Write a dynamically labeled observable-by-control response matrix."""
function write_labeled_response_matrix(
    path::AbstractString,
    row_labels,
    names,
    response_matrix;
    row_key::AbstractString="observable",
)
    labels = String.(collect(row_labels))
    control_names = String.(collect(names))
    size(response_matrix) == (length(labels), length(control_names)) || error(
        "Labeled response shape $(size(response_matrix)) does not match " *
        "$(length(labels)) rows x $(length(control_names)) controls",
    )
    length(unique(labels)) == length(labels) ||
        error("Observable response labels are not unique")
    length(unique(control_names)) == length(control_names) ||
        error("Observable response control names are not unique")
    all(isfinite, response_matrix) ||
        error("Observable response contains a non-finite value")
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, join(vcat(row_key, control_names), ','))
        for row in axes(response_matrix, 1)
            println(io, join(vcat(labels[row], response_matrix[row, :]), ','))
        end
    end
    return path
end

response_cache_metadata_path(path::AbstractString) = path * ".metadata.toml"

"""Write provenance beside a reusable closed-orbit response cache."""
function write_response_cache_metadata(
    path::AbstractString,
    names,
    response_matrix;
    method::AbstractString,
    response_step_rad::Float64,
    controls_per_batch::Int,
    chunk_count::Int,
    closure_norm_max::Float64,
    pair_id::AbstractString="",
    reltol::Union{Nothing,Float64}=nothing,
    abstol::Union{Nothing,Float64}=nothing,
    maxiter::Union{Nothing,Int}=nothing,
    model=nothing,
)
    metadata = Dict{String,Any}(
        "format" => "scibmad-closed-orbit-response-cache-v1",
        "engine" => "SciBmad",
        "rf_on" => true,
        "response_method" => String(method),
        "response_parameterization" =>
            (startswith(String(method), "gtpsa") ? "selected-controls-only" : "batch-controls"),
        "unparameterized_controls_are_primitive_zero" =>
            startswith(String(method), "gtpsa"),
        "response_step_rad" => response_step_rad,
        "controls_per_batch" => controls_per_batch,
        "chunk_count" => chunk_count,
        "closure_norm_max" => closure_norm_max,
        "control_count" => length(names),
        "control_names" => String.(names),
        "response_shape" => [size(response_matrix)...],
        "coordinate_labels" => collect(ORBIT_PHASE_SPACE_LABELS)[
            1:size(response_matrix, 1)
        ],
        "scibmad_version" => string(Base.pkgversion(SciBmad)),
    )
    isempty(pair_id) || (metadata["response_pair_id"] = String(pair_id))
    isnothing(reltol) || (metadata["reltol"] = reltol)
    isnothing(abstol) || (metadata["abstol"] = abstol)
    isnothing(maxiter) || (metadata["maxiter"] = maxiter)
    if !isnothing(model) && hasproperty(model, :metadata)
        model_metadata = model.metadata
        hasproperty(model_metadata, :ring_id) &&
            (metadata["ring_id"] = String(model_metadata.ring_id))
        hasproperty(model_metadata, :lattice_path) &&
            (metadata["lattice_path"] = String(model_metadata.lattice_path))
        hasproperty(model_metadata, :branch) &&
            (metadata["branch"] = Int(model_metadata.branch))
        if hasproperty(model_metadata, :rf_voltage) &&
           !isnothing(model_metadata.rf_voltage)
            metadata["rf_voltage"] = Float64(model_metadata.rf_voltage)
        end
    end
    metadata_path = response_cache_metadata_path(path)
    mkpath(dirname(metadata_path))
    open(metadata_path, "w") do io
        TOML.print(io, metadata; sorted=true)
    end
    return metadata_path
end

"""Read and validate optional provenance for a reusable response cache."""
function read_response_cache_metadata(
    path::AbstractString,
    names,
    response_matrix;
    model=nothing,
    requested_response_method::Union{Nothing,AbstractString}=nothing,
    requested_response_step_rad::Union{Nothing,Float64}=nothing,
    requested_response_controls_per_batch::Union{Nothing,Int}=nothing,
    required_reltol::Union{Nothing,Float64}=nothing,
    required_abstol::Union{Nothing,Float64}=nothing,
    required_maxiter::Union{Nothing,Int}=nothing,
)
    metadata_path = response_cache_metadata_path(path)
    isfile(metadata_path) || return nothing
    metadata = TOML.parsefile(metadata_path)
    get(metadata, "format", "") == "scibmad-closed-orbit-response-cache-v1" ||
        error("Unsupported response-cache metadata format: $metadata_path")
    model_ring_id = if !isnothing(model) &&
                       hasproperty(model, :metadata) &&
                       hasproperty(model.metadata, :ring_id) &&
                       !isnothing(model.metadata.ring_id)
        lowercase(String(model.metadata.ring_id))
    else
        ""
    end
    cached_ring_id = lowercase(String(get(metadata, "ring_id", "")))
    legacy_cache = startswith(model_ring_id, "legacy") ||
        startswith(cached_ring_id, "legacy")
    if !legacy_cache
        for key in (
            "response_method",
            "response_step_rad",
            "controls_per_batch",
            "chunk_count",
            "closure_norm_max",
            "reltol",
            "abstol",
            "maxiter",
            "scibmad_version",
        )
            haskey(metadata, key) || error(
                "Response cache is missing required provenance field '$key': $metadata_path",
            )
        end
    end
    String.(get(metadata, "control_names", String[])) == String.(names) ||
        error("Response-cache metadata control names/order do not match: $metadata_path")
    Int.(get(metadata, "response_shape", Int[])) == collect(size(response_matrix)) ||
        error("Response-cache metadata shape does not match the CSV: $metadata_path")
    get(metadata, "engine", "") == "SciBmad" ||
        error("Primary response cache is not marked as SciBmad: $metadata_path")
    Bool(get(metadata, "rf_on", false)) ||
        error("Response cache is not marked RF-on: $metadata_path")
    cached_version = String(get(metadata, "scibmad_version", ""))
    if !legacy_cache
        cached_version == string(Base.pkgversion(SciBmad)) ||
            error("Response cache SciBmad version does not match the runtime: $metadata_path")
    end
    raw_method = String(get(metadata, "response_method", ""))
    method = if legacy_cache && isempty(raw_method)
        "legacy"
    else
        canonical_response_method(raw_method)
    end
    if !isnothing(requested_response_method) && !legacy_cache
        requested_method = canonical_response_method(requested_response_method)
        method == requested_method || error(
            "Response cache method '$method' does not match requested method " *
            "'$(requested_response_method)'; use --recompute-response=true: $metadata_path",
        )
    end
    if !isnothing(requested_response_step_rad) &&
       method == "central-difference"
        cached_step = Float64(get(metadata, "response_step_rad", 0.0))
        isapprox(
            cached_step,
            requested_response_step_rad;
            rtol=8 * eps(Float64),
            atol=0.0,
        ) || error(
            "Response cache step $cached_step does not match requested step $requested_response_step_rad; use --recompute-response=true: $metadata_path",
        )
    end
    if !legacy_cache
        cached_controls = Int(get(metadata, "controls_per_batch", -1))
        cached_chunks = Int(get(metadata, "chunk_count", 0))
        if method == "central-difference"
            cached_controls >= 1 || error(
                "Central-difference response cache has invalid controls_per_batch: $metadata_path",
            )
            cached_chunks >= 1 || error(
                "Central-difference response cache has invalid chunk_count: $metadata_path",
            )
            if !isnothing(requested_response_controls_per_batch)
                cached_controls == requested_response_controls_per_batch || error(
                    "Response cache controls_per_batch $cached_controls does not match requested " *
                    "$requested_response_controls_per_batch; use --recompute-response=true: $metadata_path",
                )
            end
        else
            cached_controls == 0 || error(
                "GTPSA response cache must record controls_per_batch=0: $metadata_path",
            )
            cached_chunks == 1 || error(
                "GTPSA response cache must record chunk_count=1: $metadata_path",
            )
        end
    end
    if !isnothing(required_reltol) && haskey(metadata, "reltol")
        Float64(metadata["reltol"]) <= required_reltol || error(
            "Response cache reltol is looser than the requested solve: $metadata_path",
        )
    end
    if !isnothing(required_abstol) && haskey(metadata, "abstol")
        Float64(metadata["abstol"]) <= required_abstol || error(
            "Response cache abstol is looser than the requested solve: $metadata_path",
        )
    end
    if !isnothing(required_maxiter) && haskey(metadata, "maxiter")
        Int(metadata["maxiter"]) >= required_maxiter || error(
            "Response cache maxiter is lower than the requested solve: $metadata_path",
        )
    end
    if !isnothing(model) && hasproperty(model, :metadata)
        model_metadata = model.metadata
        if hasproperty(model_metadata, :ring_id)
            haskey(metadata, "ring_id") || legacy_cache || error(
                "Response cache is missing ring_id provenance: $metadata_path",
            )
            if haskey(metadata, "ring_id")
                String(metadata["ring_id"]) == String(model_metadata.ring_id) ||
                    error("Response cache ring_id does not match the selected model: $metadata_path")
            end
        end
        if hasproperty(model_metadata, :lattice_path) &&
           !isnothing(model_metadata.lattice_path) &&
           !isempty(String(model_metadata.lattice_path))
            haskey(metadata, "lattice_path") || legacy_cache || error(
                "Response cache is missing lattice-path provenance: $metadata_path",
            )
            if haskey(metadata, "lattice_path")
                normpath(String(metadata["lattice_path"])) ==
                    normpath(String(model_metadata.lattice_path)) ||
                    error("Response cache lattice path does not match the selected model: $metadata_path")
            end
        end
        if hasproperty(model_metadata, :branch)
            haskey(metadata, "branch") || legacy_cache || error(
                "Response cache is missing branch provenance: $metadata_path",
            )
            if haskey(metadata, "branch")
                Int(metadata["branch"]) == Int(model_metadata.branch) ||
                    error("Response cache branch does not match the selected model: $metadata_path")
            end
        end
        if hasproperty(model_metadata, :rf_voltage) &&
           !isnothing(model_metadata.rf_voltage)
            haskey(metadata, "rf_voltage") || legacy_cache || error(
                "Response cache is missing RF-voltage provenance: $metadata_path",
            )
            if haskey(metadata, "rf_voltage")
                Float64(metadata["rf_voltage"]) == Float64(model_metadata.rf_voltage) ||
                    error("Response cache RF voltage does not match the selected model: $metadata_path")
            end
        end
    end
    return metadata
end

function default_response_cache_path(
    artifact_ring::Symbol,
    input_path::AbstractString;
    response_method::AbstractString=DEFAULT_RESPONSE_METHOD,
)
    method = canonical_response_method(response_method)
    input_name = lowercase(basename(input_path))
    if artifact_ring == :latest_cesr && input_name == "corrector_samples.csv"
        reference_dir = joinpath(ORBIT_ROOT, "reference", "latest_cesr")
        return joinpath(
            method == "gtpsa" ? joinpath(reference_dir, "gtpsa") : reference_dir,
            "closed_orbit_response.csv",
        )
    elseif artifact_ring == :legacy && input_name == "cesr_corrector_samples_1000.csv"
        # Preserve the historical cache for explicit legacy reproduction.
        return joinpath(ORBIT_ROOT, "reference", "closed_orbit_response_6x119.csv")
    end
    stem = splitext(basename(input_path))[1]
    safe_stem = replace(stem, r"[^A-Za-z0-9_.-]" => "_")
    isempty(safe_stem) && (safe_stem = "controls")
    return joinpath(
        ORBIT_ROOT,
        "reference",
        String(artifact_ring),
        "inputs",
        safe_stem,
        method,
        "closed_orbit_response.csv",
    )
end

function read_response_matrix(path, names)
    lines = readlines(path)
    length(lines) >= 2 ||
        error("Cached closed-orbit response must have a header and at least one row: $path")
    header = split(lines[1], ',')
    first(header) == "coordinate" ||
        error("Cached response first column must be coordinate: $path")
    String.(header[2:end]) == names ||
        error("Cached response control names/order do not match the input CSV: $path")
    coordinate_labels = collect(ORBIT_PHASE_SPACE_LABELS)
    n_coordinates = length(lines) - 1
    n_coordinates >= 1 || error("Cached response has no coordinate rows: $path")
    response_matrix = Matrix{Float64}(undef, n_coordinates, length(names))
    for row in 1:n_coordinates
        fields = split(lines[row + 1], ',')
        length(fields) == length(header) ||
            error("Cached response row $row has the wrong width: $path")
        row <= length(coordinate_labels) ||
            error("Cached response has an unsupported coordinate row $row: $path")
        fields[1] == coordinate_labels[row] ||
            error("Cached response coordinate order is invalid at row $row: $path")
        for column in eachindex(names)
            response_matrix[row, column] = parse(Float64, fields[column + 1])
        end
    end
    all(isfinite, response_matrix) ||
        error("Cached response contains a non-finite value: $path")
    return response_matrix
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
    response_method = options["response-method"]
    ring = Symbol(options["ring"])
    artifact_ring = canonical_ring_id(ring)
    model_factory = (; kwargs...) -> load_ring_model(; ring, kwargs...)
    response_matrix_cache = isempty(options["response-matrix-cache"]) ?
        default_response_cache_path(
            artifact_ring,
            inputs;
            response_method,
        ) :
        abspath(options["response-matrix-cache"])
    recompute_response =
        lowercase(options["recompute-response"]) == "true"
    response_step_rad = parse(Float64, options["response-step-rad"])
    response_controls_per_batch =
        parse(Int, options["response-controls-per-batch"])
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
        response_matrix_cache,
        recompute_response,
        reltol,
        abstol,
        maxiter,
        response_step_rad,
        response_controls_per_batch,
        response_method,
        model_factory,
    )
    @printf("Warmup/compilation batch (%d samples): %.3f s\n", warmup_samples, warmup_elapsed)

    guess = prepare_initial_guess(
        samples.names,
        samples.values,
        initial_guess_mode;
        response_matrix_cache,
        recompute_response,
        reltol,
        abstol,
        maxiter,
        response_step_rad,
        response_controls_per_batch,
        response_method,
        model_factory,
    )
    if initial_guess_mode != "zero"
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
    if initial_guess_mode == "response-linear"
        if guess.response_source == "loaded"
            @printf(
                "Closed-orbit response %dx%d: loaded cache in %.6f s\n",
                size(guess.response_matrix)...,
                guess.response_load_seconds,
            )
        else
            @printf(
                "Closed-orbit response %dx%d: %.3f s via %s (chunks=%d, h=%.3e rad, controls/chunk=%d), closure max %.3e\n",
                size(guess.response_matrix)...,
                guess.response_map_seconds,
                guess.response_method,
                guess.response_chunk_count,
                guess.response_step_rad,
                guess.response_controls_per_batch,
                guess.response_closure_residual_max,
            )
        end
        println("Response cache: $response_matrix_cache")
    end

    model_timed = @timed prepare_batch_model(
        samples.names,
        samples.values;
        model_factory,
    )
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
        if jacobian_mode == "frozen-nominal"
            current_result = apply_full_newton_fallback(
                current_result,
                samples.names,
                samples.values;
                reltol,
                abstol,
                maxiter,
                model_factory,
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
    response_matrix_path = ""
    detector_response_path = ""
    if initial_guess_mode == "response-linear"
        response_matrix_path = joinpath(
            dirname(metadata_path),
            "closed_orbit_response_$(size(guess.response_matrix, 1))x$(size(guess.response_matrix, 2)).csv",
        )
        write_response_matrix(
            response_matrix_path,
            samples.names,
            guess.response_matrix,
        )
        if !isempty(guess.detector_response_labels)
            detector_response_path = joinpath(
                dirname(metadata_path),
                "detector_response_$(size(guess.detector_response, 1))x$(size(guess.detector_response, 2)).csv",
            )
            write_labeled_response_matrix(
                detector_response_path,
                guess.detector_response_labels,
                samples.names,
                guess.detector_response,
            )
        end
    end
    metadata = Dict(
        "format" => "ring-dataset-benchmark-v2",
        "engine" => "SciBmad",
        "device" => "cpu",
        "mode" => options["mode"],
        "ring" => String(ring),
        "artifact_ring" => String(artifact_ring),
        "ring_id" => String(model.metadata.ring_id),
        "lattice_path" => model.metadata.lattice_path,
        "scibmad_version" => model.metadata.scibmad_version,
        "input_csv" => inputs,
        "output_csv" => output,
        "sample_count" => size(samples.values, 1),
        "control_count" => size(samples.values, 2),
        "observable_count" => length(labels),
        "detector_count" => length(result.detectors),
        "control_names" => samples.names,
        "all_control_names" => model.metadata.all_control_names,
        "steering_control_names" => model.metadata.steering_control_names,
        "control_groups" => Dict(
            String(group) => names for (group, names) in model.metadata.control_groups
        ),
        "control_plane" => Dict(
            name => String(plane) for (name, plane) in model.metadata.control_plane
        ),
        "detector_names" => result.detectors,
        "detector_order" => result.detectors,
        "detector_element_indices" => model.metadata.detector_element_indices,
        "observable_labels" => labels,
        "coordinate_labels" => collect(model.metadata.coordinate_labels),
        "transverse_coordinate_indices" => [
            model.metadata.transverse_coordinate_indices.x,
            model.metadata.transverse_coordinate_indices.y,
        ],
        "converged_count" => count(result.converged),
        "failed_count" => count(.!result.converged),
        "warmup_sample_count" => warmup_samples,
        "warmup_seconds" => warmup_elapsed,
        "initial_guess_mode" => initial_guess_mode,
        "jacobian_mode" => jacobian_mode,
        "fallback_full_newton_enabled" =>
            jacobian_mode == "frozen-nominal",
        "nominal_model_setup_seconds" => guess.nominal_model_setup_seconds,
        "nominal_closed_orbit_seconds" => guess.nominal_solve_seconds,
        "nominal_closed_orbit_iterations" => guess.nominal_iterations,
        "nominal_closed_orbit" => guess.nominal_orbit,
        "response_matrix_path" => response_matrix_path,
        "response_matrix_cache" => response_matrix_cache,
        "response_matrix_source" => guess.response_source,
        "response_matrix_method" => guess.response_method,
        "response_method_requested" => response_method,
        "response_parameterized_control_count" => length(samples.names),
        "response_unparameterized_controls_remain_primitive_zero" => true,
        "response_matrix_shape" => [size(guess.response_matrix)...],
        "detector_response_path" => detector_response_path,
        "detector_response_shape" => [size(guess.detector_response)...],
        "detector_response_labels" => guess.detector_response_labels,
        "response_step_rad" => guess.response_step_rad,
        "response_controls_per_batch" =>
            guess.response_controls_per_batch,
        "response_chunk_count" => guess.response_chunk_count,
        "response_model_setup_seconds" => guess.response_model_setup_seconds,
        "response_map_seconds" => guess.response_map_seconds,
        "response_load_seconds" => guess.response_load_seconds,
        "response_cache_write_seconds" =>
            guess.response_cache_write_seconds,
        "response_closure_residual_max" =>
            guess.response_closure_residual_max,
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
        "final_closure_norm_median" => median(result.closure_norms),
        "final_closure_norm_max" => maximum(result.closure_norms),
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
            "one BatchParam array per control; one nominal " *
            "$(size(guess.nominal_jacobian, 1))x$(size(guess.nominal_jacobian, 2)) " *
            "Jacobian reused for all samples and iterations"
        ),
        "timed_region" => "closed-orbit solve + detector tracking",
    )
    if jacobian_mode == "frozen-nominal"
        metadata["frozen_jacobian_factorization_seconds"] =
            result.factorization_seconds
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
    isempty(response_matrix_path) ||
        println("Response: $response_matrix_path")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end

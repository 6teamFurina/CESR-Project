"""Shared configuration and layout helpers for orbit error-analysis studies.

The original studies were written for one CESR export and encoded its control
and detector counts in several independent places.  This file keeps the
historical prefix-based behaviour as a fallback, while allowing a caller to
provide a ring configuration (or a model factory that already captures one)
with named control groups and observable-plane metadata.
"""

module RingErrorAnalysisConfig

using Beamlines
using GTPSA
using LinearAlgebra
using SciBmad
using TOML
using UUIDs

export ObservableLayout,
    config_value,
    configured_model,
    configured_simulate_batch,
    canonical_response_method,
    configured_detector_names,
    validate_control_names,
    control_group_indices,
    detector_response_from_ring,
    detector_response_cache_metadata,
    read_or_compute_detector_response,
    default_ring_config,
    default_ring_paths,
    ring_artifact_id,
    observable_layout,
    plane_indices,
    observable_labels,
    ring_metadata,
    resolve_ring_configuration,
    resolve_ring_model_factory,
    coordinate_index,
    transverse_coordinate_indices,
    transverse_momentum_indices

const _RING_ADAPTER_PATH = normpath(joinpath(
    @__DIR__,
    "..",
    "Orbit_Calculation",
    "ring_model_adapter.jl",
))

# The core adapter lives next to the shared closed-orbit runner.  Include it in
# Main so its model factory can be passed to the existing runners without
# copying the latest-lattice inventory into every error-analysis script.  The
# guard also makes nested includes (parity -> rho sweep -> helper) harmless.
if !isdefined(Main, :load_ring_model) && isfile(_RING_ADAPTER_PATH)
    Base.include(Main, _RING_ADAPTER_PATH)
end

const _MISSING = gensym(:missing)

"""Normalize response-method labels used by old and new cache sidecars."""
function _canonical_response_method(method)
    value = lowercase(replace(strip(String(method)), '_' => '-'))
    value in ("central-difference", "batch-central-finite-difference") &&
        return "central-difference"
    value in ("gtpsa", "gtpsa-implicit", "gtpsa-implicit-closed-orbit") &&
        return "gtpsa"
    error("Unsupported response method '$method' in response-cache metadata")
end

canonical_response_method(method) = _canonical_response_method(method)

"""Read a key from a NamedTuple, Dict-like object, or struct."""
function config_value(config, key::Symbol, default=nothing)
    isnothing(config) && return default
    if config isa NamedTuple
        return haskey(config, key) ? getproperty(config, key) : default
    elseif config isa AbstractDict
        if haskey(config, key)
            return config[key]
        elseif haskey(config, String(key))
            return config[String(key)]
        end
        return default
    elseif hasproperty(config, key)
        return getproperty(config, key)
    end
    return default
end

"""Return the authoritative ring config supplied by the core adapter."""
function default_ring_config(; ring::Symbol=:latest)
    if isdefined(Main, :ring_analysis_config)
        return Main.ring_analysis_config(; ring)
    end
    return nothing
end

"""Resolve an omitted config to the maintained latest ring (or explicit legacy)."""
function resolve_ring_configuration(config=nothing; ring::Symbol=:latest)
    return isnothing(config) ? default_ring_config(; ring) : config
end

"""Resolve the model factory from ring metadata, preserving explicit overrides."""
function resolve_ring_model_factory(model_factory=nothing, config=nothing; ring::Symbol=:latest)
    if !isnothing(model_factory)
        # A caller-supplied factory without a config is intentionally opaque;
        # do not attach the default CESR registry to an unrelated ring.
        resolved_config = isnothing(config) ? nothing :
            resolve_ring_configuration(config; ring)
        return model_factory, resolved_config
    end
    resolved_config = resolve_ring_configuration(config; ring)
    configured_factory = config_value(resolved_config, :model_factory, nothing)
    if !isnothing(configured_factory)
        return configured_factory, resolved_config
    end
    # This fallback is only reached when the core adapter is unavailable.  It
    # retains the historical loader for old checkouts while normal checkouts
    # always resolve the latest config above.
    fallback = isdefined(Main, :load_cesr_model) ? Main.load_cesr_model : nothing
    isnothing(fallback) && error("No ring model factory is available")
    return fallback, resolved_config
end

"""Return ring-scoped input/cache defaults without encoding dimensions."""
function default_ring_paths(
    ;
    ring::Symbol=:latest,
    config=nothing,
    response_method::AbstractString="gtpsa",
)
    resolved_ring = Symbol(config_value(
        config,
        :ring_id,
        config_value(config, :ring, ring),
    ))
    method = _canonical_response_method(response_method)
    artifact = ring_artifact_id(resolved_ring)
    if resolved_ring in (:legacy, :legacy_cesr, :historical) || artifact == "legacy"
        orbit_root = normpath(joinpath(@__DIR__, ".."))
        calculation_dir = joinpath(orbit_root, "Orbit_Calculation")
        project_root = normpath(joinpath(orbit_root, ".."))
        if method == "gtpsa"
            response_dir = joinpath(orbit_root, "reference", "legacy", "gtpsa")
            return (
                inputs=joinpath(calculation_dir, "inputs", "cesr_corrector_samples_1000.csv"),
                detector_response=joinpath(response_dir, "detector_response.csv"),
                closed_orbit_response=joinpath(response_dir, "closed_orbit_response.csv"),
            )
        end
        return (
            inputs=joinpath(calculation_dir, "inputs", "cesr_corrector_samples_1000.csv"),
            detector_response=joinpath(
                project_root,
                "older_ring_version",
                "bmad_comparison",
                "bmad_control_response_rf_on",
                "scibmad_control_response_rf_on.csv",
            ),
            closed_orbit_response=joinpath(
                orbit_root,
                "reference",
                "closed_orbit_response_6x119.csv",
            ),
        )
    end

    orbit_root = normpath(joinpath(@__DIR__, ".."))
    calculation_dir = joinpath(orbit_root, "Orbit_Calculation")
    reference_dir = joinpath(orbit_root, "reference", artifact)
    response_dir = method == "gtpsa" ? joinpath(reference_dir, "gtpsa") : reference_dir
    input_dir = joinpath(calculation_dir, "inputs", artifact)
    return (
        inputs=joinpath(input_dir, "corrector_samples.csv"),
        detector_response=joinpath(response_dir, "detector_response.csv"),
        closed_orbit_response=joinpath(response_dir, "closed_orbit_response.csv"),
    )
end

"""Map a CLI ring alias to the directory name used for generated artifacts."""
function ring_artifact_id(ring::Symbol)
    ring in (:latest, :latest_cesr, :repaired_latest) && return "latest_cesr"
    ring in (:legacy, :legacy_cesr, :historical) && return "legacy"
    artifact = replace(String(ring), r"[^A-Za-z0-9_.-]" => "_")
    isempty(artifact) && error("Ring identifier cannot be empty")
    return artifact
end

ring_artifact_id(ring) = ring_artifact_id(Symbol(ring))

"""Return provenance fields common to every latest/legacy result metadata file."""
function ring_metadata(config=nothing; ring::Symbol=:latest)
    raw_ring = Symbol(config_value(
        config,
        :ring_id,
        config_value(config, :ring, ring),
    ))
    ring_id = ring_artifact_id(raw_ring)
    scibmad_version = try
        string(Base.pkgversion(SciBmad))
    catch
        "unknown"
    end
    rf_value = config_value(config, :rf_on, true)
    isnothing(rf_value) && (rf_value = true)
    metadata = Dict{String,Any}(
        "ring_id" => ring_id,
        "lattice_path" => String(config_value(config, :lattice_path, "")),
        "branch" => config_value(config, :branch, 0),
        "rf_on" => Bool(rf_value),
        "scibmad_version" => scibmad_version,
        "control_names_from_config" => config_value(config, :control_names, String[]),
        "observable_labels_from_config" => config_value(config, :observable_labels, String[]),
    )
    rf_voltage = config_value(config, :rf_voltage, nothing)
    isnothing(rf_voltage) || (metadata["rf_voltage"] = Float64(rf_voltage))
    return metadata
end

"""Return the detector registry from ring metadata, with legacy fallback.

The core benchmark's historical `detector_names` helper discovers `DET_*`
elements.  New ring configs carry an authoritative ordered detector registry
so custom detector naming/prefixes do not silently change observable order.
"""
function configured_detector_names(model=nothing, config=nothing)
    for source in (config, model, config_value(model, :metadata, nothing))
        configured = config_value(source, :detector_names, nothing)
        isnothing(configured) && (configured = config_value(source, :detector_order, nothing))
        if !isnothing(configured)
            names = uppercase.(String.(collect(configured)))
            isempty(names) && error("Configured detector registry is empty")
            length(unique(names)) == length(names) || error("Configured detector names are not unique")
            return names
        end
    end
    model_ring = config_value(model, :ring, nothing)
    if !isnothing(model_ring) && isdefined(Main, :detector_names)
        names = uppercase.(String.(Main.detector_names(model_ring)))
        isempty(names) && error("No detectors were found in the selected ring")
        return names
    end
    error("No detector registry is available for the selected ring")
end

"""Validate an input subset/order against the selected ring's named registry."""
function validate_control_names(names, config=nothing)
    unique_names = String.(names)
    length(unique(unique_names)) == length(unique_names) ||
        error("Control names are not unique")
    configured = config_value(config, :all_control_names, nothing)
    isnothing(configured) && (configured = config_value(config, :control_names, nothing))
    isnothing(configured) && return unique_names
    allowed = Set(String.(collect(configured)))
    unknown = [name for name in unique_names if !(name in allowed)]
    isempty(unknown) || error(
        "Input controls are not present in the selected ring registry: $(join(unknown, ", "))",
    )
    return unique_names
end

"""Resolve one named phase-space coordinate from ring metadata.

The maintained adapter records coordinate labels and the transverse position
indices.  A custom ring may provide either of those fields; only the standard
SciBmad ordering remains as a compatibility fallback for archived configs.
"""
function coordinate_index(
    config,
    label::Symbol;
    state_dimension::Union{Nothing,Int}=nothing,
)
    configured = config_value(config, :coordinate_indices, nothing)
    isnothing(configured) && (configured = config_value(config, :transverse_indices, nothing))
    coordinate = config_value(configured, label, nothing)
    labels = config_value(config, :coordinate_labels, nothing)
    if isnothing(coordinate) && !isnothing(labels)
        lowered = lowercase.(String.(collect(labels)))
        coordinate = findfirst(==(lowercase(String(label))), lowered)
    end
    # The fallback is intentionally centralized here so every study has the
    # same explicit historical compatibility behavior.
    fallback = Dict(:x => 1, :px => 2, :y => 3, :py => 4, :z => 5, :pz => 6)
    coordinate = isnothing(coordinate) ? get(fallback, label, nothing) : Int(coordinate)
    isnothing(coordinate) && error("No phase-space coordinate metadata for '$label'")
    if !isnothing(state_dimension)
        1 <= coordinate <= state_dimension ||
            error("Configured '$label' coordinate index $coordinate is outside the $state_dimension-state orbit")
    end
    return coordinate
end

"""Resolve the configured x/y coordinates, retaining the SciBmad fallback."""
function transverse_coordinate_indices(config=nothing; state_dimension::Union{Nothing,Int}=nothing)
    return (
        x=coordinate_index(config, :x; state_dimension),
        y=coordinate_index(config, :y; state_dimension),
    )
end

"""Resolve the configured px/py coordinates for local transverse kicks."""
function transverse_momentum_indices(config=nothing; state_dimension::Union{Nothing,Int}=nothing)
    return (
        px=coordinate_index(config, :px; state_dimension),
        py=coordinate_index(config, :py; state_dimension),
    )
end

"""Write a detector-response matrix using the common labeled CSV contract."""
function _write_detector_response(path, layout, names, response)
    size(response) == (length(layout.labels), length(names)) ||
        error("Detector response has shape $(size(response)); expected $(length(layout.labels)) x $(length(names))")
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, join(vcat("observable", String.(names)), ','))
        for row in eachindex(layout.labels)
            println(io, join(vcat(layout.labels[row], response[row, :]), ','))
        end
    end
    return path
end

function _finite_difference_response_pair(
    model_factory,
    config,
    names,
    detectors,
    layout;
    rf_on=true,
    reltol=1e-12,
    abstol=1e-13,
    maxiter=100,
    response_step_rad=1.0e-7,
    controls_per_batch=8,
)
    rf_on || error("The maintained detector response is defined for RF-on closed orbits")
    control_names = String.(collect(names))
    n_controls = length(control_names)
    n_controls > 0 || error("Cannot compute a detector response with no controls")
    isdefined(Main, :central_finite_difference_responses) || error(
        "The shared SciBmad runner must be loaded before computing detector responses",
    )
    response_factory = (; kwargs...) -> configured_model(
        model_factory,
        config;
        kwargs...,
    )
    pair = Main.central_finite_difference_responses(
        control_names;
        response_step_rad,
        controls_per_batch,
        reltol,
        abstol,
        maxiter,
        model_factory=response_factory,
    )
    return _select_response_pair(pair, detectors, layout)
end

"""Select the configured detector/observable rows from a shared response pair."""
function _select_response_pair(pair, detectors, layout)
    uppercase.(pair.detectors) == uppercase.(String.(collect(detectors))) ||
        error("Response detector registry does not match the configured order")
    source_index = Dict(
        label => index for (index, label) in enumerate(pair.observable_labels)
    )
    requested_indices = Int[]
    for label in layout.labels
        haskey(source_index, label) || error(
            "Configured observable '$label' is not an x/y detector observable produced by the shared SciBmad runner",
        )
        push!(requested_indices, source_index[label])
    end
    detector_response = pair.detector_response[requested_indices, :]
    return merge(pair, (; detector_response))
end

"""Compute the selected-control GTPSA response pair through the core runner."""
function _gtpsa_response_pair(
    model_factory,
    config,
    names,
    detectors,
    layout;
    rf_on=true,
    reltol=1e-12,
    abstol=1e-13,
    maxiter=100,
)
    rf_on || error("The maintained GTPSA detector response is defined for RF-on closed orbits")
    isdefined(Main, :gtpsa_first_order_responses) || error(
        "The shared SciBmad runner must be loaded before computing GTPSA detector responses",
    )
    control_names = String.(collect(names))
    isempty(control_names) && error("Cannot compute a detector response with no controls")
    response_factory = (; kwargs...) -> configured_model(
        model_factory,
        config;
        kwargs...,
    )
    pair = Main.gtpsa_first_order_responses(
        control_names;
        rf_on,
        reltol,
        abstol,
        maxiter,
        model_factory=response_factory,
    )
    return _select_response_pair(pair, detectors, layout)
end

"""Compute a labeled detector response for the selected ring.

GTPSA implicit closed-orbit differentiation is the default backend and
publishes the paired closed-orbit/detector caches.  Use
`response_method="central-difference"` explicitly for an independent
validation or fallback calculation.  The default latest CESR H/V steering
subset is supported by GTPSA; selecting known skew/group controls such as
`SK_Q14W` remains an explicit unsupported-control case.
"""
function detector_response_from_ring(
    model_factory,
    config,
    names,
    detectors,
    layout;
    response_method::AbstractString="gtpsa",
    kwargs...,
)
    method = _canonical_response_method(response_method)
    pair = if method == "gtpsa"
        _gtpsa_response_pair(
            model_factory,
            config,
            names,
            detectors,
            layout;
            kwargs...,
        )
    else
        _finite_difference_response_pair(
            model_factory,
            config,
            names,
            detectors,
            layout;
            kwargs...,
        )
    end
    return pair.detector_response
end

function _write_detector_response_metadata(
    path,
    layout,
    names,
    pair,
    config,
    ;
    pair_id::AbstractString,
    reltol::Float64,
    abstol::Float64,
    maxiter::Int,
)
    metadata = isnothing(config) ? Dict{String,Any}() : ring_metadata(config)
    merge!(metadata, Dict{String,Any}(
        "format" => "scibmad-detector-response-cache-v1",
        "engine" => "SciBmad",
        "rf_on" => true,
        "scibmad_version" => try
            string(Base.pkgversion(SciBmad))
        catch
            "unknown"
        end,
        "response_method" => String(pair.response_method),
        "response_pair_id" => String(pair_id),
        "response_step_rad" => pair.response_step_rad,
        "controls_per_batch" => pair.controls_per_batch,
        "chunk_count" => pair.chunk_count,
        "closure_norm_max" => pair.closure_norm_max,
        "control_count" => length(names),
        "control_names" => String.(names),
        "detector_count" => length(layout.detectors),
        "detector_names" => layout.detectors,
        "observable_count" => length(layout.labels),
        "observable_labels" => layout.labels,
        "response_shape" => [size(pair.detector_response)...],
        "reltol" => reltol,
        "abstol" => abstol,
        "maxiter" => maxiter,
    ))
    metadata_path = String(path) * ".metadata.toml"
    open(metadata_path, "w") do io
        TOML.print(io, metadata; sorted=true)
    end
    return metadata_path
end

function _validate_detector_response_metadata(
    path,
    layout,
    names,
    config,
    ;
    requested_response_method=nothing,
    requested_response_step_rad::Union{Nothing,Float64}=nothing,
    requested_response_controls_per_batch::Union{Nothing,Int}=nothing,
    required_reltol::Union{Nothing,Float64}=nothing,
    required_abstol::Union{Nothing,Float64}=nothing,
    required_maxiter::Union{Nothing,Int}=nothing,
)
    metadata_path = String(path) * ".metadata.toml"
    isfile(metadata_path) || return nothing
    metadata = TOML.parsefile(metadata_path)
    get(metadata, "format", "") == "scibmad-detector-response-cache-v1" ||
        error("Unsupported detector-response metadata format: $metadata_path")
    selected_ring = config_value(
        config,
        :ring_id,
        config_value(config, :ring, ""),
    )
    selected_ring_id = isnothing(selected_ring) ? "" : lowercase(String(selected_ring))
    cached_ring_id = lowercase(String(get(metadata, "ring_id", "")))
    legacy_cache = startswith(selected_ring_id, "legacy") ||
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
                "Detector-response cache is missing required provenance field '$key': $metadata_path",
            )
        end
    end
    get(metadata, "engine", "") == "SciBmad" ||
        error("Primary detector response is not marked as SciBmad: $metadata_path")
    Bool(get(metadata, "rf_on", false)) ||
        error("Detector response is not marked RF-on: $metadata_path")
    cached_version = String(get(metadata, "scibmad_version", ""))
    if !legacy_cache
        cached_version == string(Base.pkgversion(SciBmad)) ||
            error("Detector-response SciBmad version does not match the runtime: $metadata_path")
    end
    String.(get(metadata, "control_names", String[])) == String.(names) ||
        error("Detector-response metadata control names/order do not match: $metadata_path")
    String.(get(metadata, "detector_names", String[])) == layout.detectors ||
        error("Detector-response metadata detector names/order do not match: $metadata_path")
    String.(get(metadata, "observable_labels", String[])) == layout.labels ||
        error("Detector-response metadata observable labels/order do not match: $metadata_path")
    Int.(get(metadata, "response_shape", Int[])) ==
        [length(layout.labels), length(names)] ||
        error("Detector-response metadata shape does not match the requested layout: $metadata_path")
    raw_method = String(get(metadata, "response_method", ""))
    method = if legacy_cache && isempty(raw_method)
        "legacy"
    else
        _canonical_response_method(raw_method)
    end
    if !legacy_cache
        if !isnothing(requested_response_method)
            requested_method = _canonical_response_method(requested_response_method)
            method == requested_method || error(
                "Detector-response method '$raw_method' does not match requested " *
                "'$requested_response_method': $metadata_path",
            )
        end
        method in ("central-difference", "gtpsa") || error(
            "Unsupported detector-response method '$raw_method': $metadata_path",
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
            "Detector-response step $cached_step does not match requested step $requested_response_step_rad: $metadata_path",
        )
    end
    if !legacy_cache
        cached_controls = Int(get(metadata, "controls_per_batch", -1))
        cached_chunks = Int(get(metadata, "chunk_count", 0))
        if method == "central-difference"
            cached_controls >= 1 || error(
                "Central-difference detector response has invalid controls_per_batch: $metadata_path",
            )
            cached_chunks >= 1 || error(
                "Central-difference detector response has invalid chunk_count: $metadata_path",
            )
            if !isnothing(requested_response_controls_per_batch)
                cached_controls == requested_response_controls_per_batch || error(
                    "Detector-response controls_per_batch $cached_controls does not match requested " *
                    "$requested_response_controls_per_batch: $metadata_path",
                )
            end
        else
            cached_controls == 0 || error(
                "GTPSA detector response must record controls_per_batch=0: $metadata_path",
            )
            cached_chunks == 1 || error(
                "GTPSA detector response must record chunk_count=1: $metadata_path",
            )
        end
    end
    if !isnothing(required_reltol) && haskey(metadata, "reltol")
        Float64(metadata["reltol"]) <= required_reltol ||
            error("Detector-response reltol is looser than requested: $metadata_path")
    end
    if !isnothing(required_abstol) && haskey(metadata, "abstol")
        Float64(metadata["abstol"]) <= required_abstol ||
            error("Detector-response abstol is looser than requested: $metadata_path")
    end
    if !isnothing(required_maxiter) && haskey(metadata, "maxiter")
        Int(metadata["maxiter"]) >= required_maxiter ||
            error("Detector-response maxiter is lower than requested: $metadata_path")
    end

    if !isnothing(config)
        expected = ring_metadata(config)
        expected_ring = String(get(expected, "ring_id", ""))
        cached_ring = String(get(metadata, "ring_id", ""))
        if !isempty(expected_ring)
            !isempty(cached_ring) || legacy_cache || error(
                "Detector-response cache is missing ring_id provenance: $metadata_path",
            )
            isempty(cached_ring) || expected_ring == cached_ring || error(
                "Detector-response ring_id does not match the selected config: $metadata_path",
            )
        end
        expected_lattice = String(get(expected, "lattice_path", ""))
        cached_lattice = String(get(metadata, "lattice_path", ""))
        if !isempty(expected_lattice)
            !isempty(cached_lattice) || legacy_cache || error(
                "Detector-response cache is missing lattice-path provenance: $metadata_path",
            )
            isempty(cached_lattice) ||
                normpath(expected_lattice) == normpath(cached_lattice) || error(
                    "Detector-response lattice path does not match the selected config: $metadata_path",
                )
        end
        if haskey(expected, "branch")
            haskey(metadata, "branch") || legacy_cache || error(
                "Detector-response cache is missing branch provenance: $metadata_path",
            )
            if haskey(metadata, "branch")
                Int(expected["branch"]) == Int(metadata["branch"]) ||
                    error("Detector-response branch does not match the selected config: $metadata_path")
            end
        end
        expected_rf_voltage = config_value(config, :rf_voltage, nothing)
        if !isnothing(expected_rf_voltage)
            haskey(metadata, "rf_voltage") || legacy_cache || error(
                "Detector-response cache is missing RF-voltage provenance: $metadata_path",
            )
            if haskey(metadata, "rf_voltage")
                Float64(expected_rf_voltage) == Float64(metadata["rf_voltage"]) ||
                    error("Detector-response RF voltage does not match the selected config: $metadata_path")
            end
        end
    end
    return metadata
end

function _validate_response_pair(
    detector_metadata,
    closed_metadata,
    detector_path,
    closed_path,
)
    if isnothing(detector_metadata) || isnothing(closed_metadata)
        isnothing(detector_metadata) && isnothing(closed_metadata) && return nothing
        error(
            "Only one response cache has provenance metadata; recompute the pair: $detector_path, $closed_path",
        )
    end
    for key in (
        "response_method",
        "response_step_rad",
        "controls_per_batch",
        "chunk_count",
        "closure_norm_max",
        "reltol",
        "abstol",
        "maxiter",
        "ring_id",
        "lattice_path",
        "branch",
        "rf_on",
        "rf_voltage",
        "scibmad_version",
    )
        if haskey(detector_metadata, key) || haskey(closed_metadata, key)
            haskey(detector_metadata, key) && haskey(closed_metadata, key) ||
                error("Response cache pair differs in metadata field '$key'")
            if key == "response_method"
                detector_method = String(detector_metadata[key])
                closed_method = String(closed_metadata[key])
                detector_method = try
                    _canonical_response_method(detector_method)
                catch
                    detector_method
                end
                closed_method = try
                    _canonical_response_method(closed_method)
                catch
                    closed_method
                end
                detector_method == closed_method || error(
                    "Response cache pair has mismatched metadata field '$key'",
                )
            else
                detector_metadata[key] == closed_metadata[key] ||
                    error("Response cache pair has mismatched metadata field '$key'")
            end
        end
    end
    detector_pair_id = String(get(detector_metadata, "response_pair_id", ""))
    closed_pair_id = String(get(closed_metadata, "response_pair_id", ""))
    detector_ring_id = lowercase(String(get(detector_metadata, "ring_id", "")))
    closed_ring_id = lowercase(String(get(closed_metadata, "ring_id", "")))
    identified_nonlegacy = (!isempty(detector_ring_id) &&
                            !startswith(detector_ring_id, "legacy")) ||
        (!isempty(closed_ring_id) && !startswith(closed_ring_id, "legacy"))
    if identified_nonlegacy
        !isempty(detector_pair_id) && detector_pair_id == closed_pair_id ||
            error("Closed-orbit and detector responses do not share a pair id")
    elseif !isempty(detector_pair_id) || !isempty(closed_pair_id)
        !isempty(detector_pair_id) && detector_pair_id == closed_pair_id ||
            error("Closed-orbit and detector responses do not share a pair id")
    end
    return nothing
end

"""Return optional detector-response provenance without requiring a cache."""
function detector_response_cache_metadata(path)
    metadata_path = String(path) * ".metadata.toml"
    return isfile(metadata_path) ? TOML.parsefile(metadata_path) : nothing
end

"""Stage both response artifacts, then publish their shared pair id.

There is no portable two-file rename primitive, so a pending marker makes a
partially published pair unreadable.  Each CSV and sidecar is first written to
a sibling temporary path; readers either see the previous complete pair or
wait for the next call to repair the marked pair.
"""
function _publish_response_pair(
    detector_path,
    closed_orbit_path,
    layout,
    names,
    pair,
    metadata_model,
    metadata_source;
    pair_id::AbstractString,
    reltol::Float64,
    abstol::Float64,
    maxiter::Int,
)
    detector_path = String(detector_path)
    closed_path = isnothing(closed_orbit_path) ? nothing : String(closed_orbit_path)
    !isnothing(closed_path) && detector_path == closed_path && error(
        "Detector and closed-orbit response paths must be different",
    )
    mkpath(dirname(detector_path))
    isnothing(closed_path) || mkpath(dirname(closed_path))
    pending_path = isnothing(closed_path) ?
        detector_path * ".pair.pending" : closed_path * ".pair.pending"
    temporary_paths = String[]
    detector_tmp = tempname(dirname(detector_path))
    detector_metadata_tmp = detector_tmp * ".metadata.toml"
    push!(temporary_paths, detector_tmp, detector_metadata_tmp)
    closed_tmp = nothing
    closed_metadata_tmp = nothing
    if !isnothing(closed_path)
        closed_tmp = tempname(dirname(closed_path))
        closed_metadata_tmp = closed_tmp * ".metadata.toml"
        push!(temporary_paths, closed_tmp, closed_metadata_tmp)
        isdefined(Main, :write_response_matrix) || error(
            "The shared SciBmad response writer is not loaded",
        )
        isdefined(Main, :write_response_cache_metadata) || error(
            "The shared SciBmad response metadata writer is not loaded",
        )
    end
    open(pending_path, "w") do io
        println(io, pair_id)
    end
    published = false
    try
        _write_detector_response(detector_tmp, layout, names, pair.detector_response)
        _write_detector_response_metadata(
            detector_tmp,
            layout,
            names,
            pair,
            metadata_source;
            pair_id,
            reltol,
            abstol,
            maxiter,
        )
        if !isnothing(closed_path)
            Main.write_response_matrix(
                closed_tmp,
                String.(names),
                pair.closed_orbit_response,
            )
            Main.write_response_cache_metadata(
                closed_tmp,
                String.(names),
                pair.closed_orbit_response;
                method=String(pair.response_method),
                response_step_rad=pair.response_step_rad,
                controls_per_batch=pair.controls_per_batch,
                chunk_count=pair.chunk_count,
                closure_norm_max=pair.closure_norm_max,
                pair_id,
                reltol,
                abstol,
                maxiter,
                model=metadata_model,
            )
        end

        # Publish data before sidecars; the marker prevents readers from using
        # the intermediate state.  The final sidecar moves complete the pair.
        mv(detector_tmp, detector_path; force=true)
        !isnothing(closed_path) && mv(closed_tmp, closed_path; force=true)
        mv(detector_metadata_tmp, detector_path * ".metadata.toml"; force=true)
        !isnothing(closed_path) &&
            mv(closed_metadata_tmp, closed_path * ".metadata.toml"; force=true)
        published = true
    finally
        published && isfile(pending_path) && rm(pending_path; force=true)
        for temporary_path in temporary_paths
            isfile(temporary_path) && rm(temporary_path; force=true)
        end
    end
    return nothing
end

"""Read a labeled response or compute it into the ring-scoped path on demand."""
function read_or_compute_detector_response(
    reader,
    path::AbstractString,
    names,
    detectors;
    layout,
    model_factory,
    config=nothing,
    closed_orbit_path=nothing,
    rf_on=true,
    reltol::Float64=1e-12,
    abstol::Float64=1e-13,
    maxiter::Int=100,
    response_step_rad::Float64=1.0e-7,
    controls_per_batch::Int=8,
    response_method="gtpsa",
    recompute_response::Bool=false,
)
    requested_method = isnothing(response_method) ? nothing :
        _canonical_response_method(response_method)
    pending_path = isnothing(closed_orbit_path) ?
        String(path) * ".pair.pending" : String(closed_orbit_path) * ".pair.pending"
    pair_pending = isfile(pending_path)
    detector_exists = !recompute_response && isfile(path) && !pair_pending
    closed_exists = !recompute_response &&
        (isnothing(closed_orbit_path) || isfile(closed_orbit_path)) &&
        !pair_pending
    detector_metadata = if detector_exists && closed_exists
        _validate_detector_response_metadata(
            path,
            layout,
            names,
            config;
            requested_response_method=requested_method,
            requested_response_step_rad=response_step_rad,
            requested_response_controls_per_batch=controls_per_batch,
            required_reltol=reltol,
            required_abstol=abstol,
            required_maxiter=maxiter,
        )
    else
        nothing
    end
    if detector_exists && closed_exists && isnothing(detector_metadata)
        selected_ring = config_value(
            config,
            :ring_id,
            config_value(config, :ring, ""),
        )
        legacy_selection = !isnothing(selected_ring) &&
            startswith(lowercase(String(selected_ring)), "legacy")
        legacy_selection || error(
            "Latest/custom detector-response cache is missing its provenance sidecar; " *
            "recompute the response pair: $(String(path)).metadata.toml",
        )
    end
    if detector_exists && closed_exists
        if !isnothing(closed_orbit_path)
            isdefined(Main, :read_response_matrix) &&
                isdefined(Main, :read_response_cache_metadata) ||
                error("The shared SciBmad response cache reader is not loaded")
            closed_response = Main.read_response_matrix(
                closed_orbit_path,
                String.(names),
            )
            closed_metadata = Main.read_response_cache_metadata(
                closed_orbit_path,
                String.(names),
                closed_response;
                requested_response_method=requested_method,
                requested_response_step_rad=response_step_rad,
                requested_response_controls_per_batch=controls_per_batch,
                required_reltol=reltol,
                required_abstol=abstol,
                required_maxiter=maxiter,
            )
            _validate_response_pair(
                detector_metadata,
                closed_metadata,
                path,
                closed_orbit_path,
            )
        end
        return reader(path, names, detectors; layout)
    end
    pair = if requested_method == "gtpsa"
        _gtpsa_response_pair(
            model_factory,
            config,
            names,
            detectors,
            layout;
            rf_on,
            reltol,
            abstol,
            maxiter,
        )
    else
        _finite_difference_response_pair(
            model_factory,
            config,
            names,
            detectors,
            layout;
            rf_on,
            reltol,
            abstol,
            maxiter,
            response_step_rad,
            controls_per_batch,
        )
    end
    pair_id = string(uuid4())
    metadata_model = configured_model(
        model_factory,
        config;
        zero_value=0.0,
        rf_on=true,
    )
    metadata_source = config_value(metadata_model, :metadata, config)
    _publish_response_pair(
        path,
        closed_orbit_path,
        layout,
        names,
        pair,
        metadata_model,
        metadata_source;
        pair_id,
        reltol,
        abstol,
        maxiter,
    )
    return pair.detector_response
end

"""Construct a model while accepting the optional `config` hook.

The fallback is intentional: historical `load_cesr_model` does not accept a
configuration keyword, whereas newer model factories may.  A custom factory
can therefore be used unchanged by old and new scripts alike.
"""
function configured_model(factory, config=nothing; kwargs...)
    isnothing(config) && return factory(; kwargs...)
    try
        return factory(; kwargs..., config=config)
    catch exception
        unsupported_config_keyword = exception isa MethodError &&
            exception.f === Core.kwcall &&
            length(exception.args) >= 2 &&
            exception.args[2] === factory
        unsupported_config_keyword || rethrow()
        return factory(; kwargs...)
    end
end

"""Call the core batch runner with an optional ring configuration."""
function configured_simulate_batch(
    runner,
    names,
    values;
    config=nothing,
    model_factory,
    kwargs...,
)
    effective_factory = isnothing(config) ? model_factory :
        (; runtime_kwargs...) -> configured_model(
            model_factory,
            config;
            runtime_kwargs...,
        )
    return runner(names, values; kwargs..., model_factory=effective_factory)
end

"""Return the index set for a named control plane/group.

`config.control_groups` may map a group name to control names or integer
indices.  `config.control_plane` may instead map each control name to a plane.
The legacy H/V name-prefix convention is retained only as an explicit fallback
for archived CESR inputs.
"""
function control_group_indices(names, group::AbstractString; config=nothing, model=nothing)
    group_name = lowercase(String(group))
    groups = config_value(config, :control_groups, nothing)
    if isnothing(groups) && !isnothing(model)
        groups = config_value(model, :control_groups, nothing)
    end
    candidate = nothing
    if !isnothing(groups)
        for key in (Symbol(group_name), group_name, Symbol(uppercase(group_name)))
            if groups isa AbstractDict && haskey(groups, key)
                candidate = groups[key]
                break
            elseif groups isa NamedTuple && haskey(groups, Symbol(key))
                candidate = getproperty(groups, Symbol(key))
                break
            end
        end
    end
    if !isnothing(candidate)
        if all(item -> item isa Integer, candidate)
            indices = Int.(collect(candidate))
            registry = config_value(config, :all_control_names, nothing)
            isnothing(registry) && (registry = config_value(config_value(model, :metadata, nothing), :all_control_names, nothing))
            if !isnothing(registry)
                registry_names = String.(collect(registry))
                all(index -> 1 <= index <= length(registry_names), indices) ||
                    error("Control group '$group' contains an out-of-range registry index")
                wanted = Set(registry_names[indices])
                selected = findall(name -> name in wanted, names)
                isempty(selected) && error("Control group '$group' has no controls in the input subset")
                return selected
            end
            all(index -> 1 <= index <= length(names), indices) ||
                error("Control group '$group' contains an out-of-range input index")
            return indices
        end
        wanted = Set(String.(collect(candidate)))
        indices = findall(name -> name in wanted, names)
        isempty(indices) && error("Control group '$group' has no controls in the input subset")
        return indices
    end

    plane_map = config_value(config, :control_plane, nothing)
    if !isnothing(plane_map)
        indices = findall(names) do name
            value = if plane_map isa AbstractDict
                get(plane_map, name, get(plane_map, String(name), nothing))
            else
                nothing
            end
            isnothing(value) ? false : lowercase(String(value)) == group_name
        end
        !isempty(indices) && return indices
    end

    # Historical CESR fallback.  This is deliberately not a count assertion.
    prefix = group_name == "horizontal" ? "H" :
        group_name == "vertical" ? "V" : nothing
    isnothing(prefix) && error("No control metadata is available for group '$group'")
    indices = findall(name -> startswith(uppercase(String(name)), prefix), names)
    isempty(indices) && error("No controls found for group '$group'")
    return indices
end

struct ObservableLayout
    labels::Vector{String}
    detectors::Vector{String}
    planes::Dict{Symbol,Vector{Int}}
end

function _as_index_vector(value, labels, detectors)
    isnothing(value) && return Int[]
    items = collect(value)
    if all(item -> item isa Integer, items)
        indices = Int.(items)
        all(index -> 1 <= index <= length(labels), indices) ||
            error("Observable-plane index is outside the response vector")
        return indices
    end
    wanted = Set(String.(items))
    indices = findall(label -> label in wanted, labels)
    length(indices) == length(wanted) ||
        error("Observable-plane labels do not match the response labels")
    return indices
end

"""Build observable labels and named detector-plane index sets."""
function observable_layout(detectors; config=nothing, labels=nothing)
    detector_names = String.(detectors)
    configured_labels = isnothing(labels) ? config_value(config, :observable_labels, nothing) : labels
    labels_vector = if isnothing(configured_labels)
        vcat(["$(name):x" for name in detector_names], ["$(name):y" for name in detector_names])
    else
        String.(collect(configured_labels))
    end
    length(unique(labels_vector)) == length(labels_vector) ||
        error("Observable labels are not unique")

    configured_planes = config_value(config, :observable_planes, nothing)
    if isnothing(configured_planes)
        configured_planes = config_value(config, :detector_plane_indices, nothing)
    end
    planes = Dict{Symbol,Vector{Int}}()
    if !isnothing(configured_planes)
        for (key, value) in pairs(configured_planes)
            planes[Symbol(lowercase(String(key)))] =
                _as_index_vector(value, labels_vector, detector_names)
        end
    else
        for plane in (:x, :y)
            suffix = ":$(plane)"
            indices = findall(label -> endswith(lowercase(label), suffix), labels_vector)
            !isempty(indices) && (planes[plane] = indices)
        end
    end
    # Keep the standard two-plane layout usable when a custom ring supplies
    # labels without suffixes but preserves the conventional block ordering.
    if isempty(planes) && length(labels_vector) == 2 * length(detector_names)
        planes[:x] = collect(1:length(detector_names))
        planes[:y] = collect(length(detector_names) + 1:length(labels_vector))
    end
    isempty(planes) && error("No observable plane metadata could be inferred")
    return ObservableLayout(labels_vector, detector_names, planes)
end

observable_labels(layout::ObservableLayout) = layout.labels
plane_indices(layout::ObservableLayout, plane::Symbol) =
    get(layout.planes, plane, Int[])

end # module RingErrorAnalysisConfig

"""
Ring/lattice adapter used by the orbit benchmark.

The original benchmark was tied to the historical `cesr_model.jl` loader and
assumed fixed control and detector/output dimensions.  The repaired latest
lattice exports its Bmad Overlay/Group controls as module globals, so a model
owns a private Julia namespace and a Dict-like control registry.  This keeps
the deferred expressions live while allowing Float64, BatchParam, and GTPSA
values to be assigned without converting them to Float64.

The latest repaired lattice is the default.  Use `ring=:legacy` explicitly for
historical CESR reproductions.  A custom lattice can be supplied with
`lattice_path`, `ring_symbol`, and `control_specs`.
"""

using Beamlines
using SciBmad

"""Map accepted CLI/API aliases to the stable artifact directory name."""
function canonical_ring_id(ring::Symbol)
    if ring in (:latest, :latest_cesr, :repaired_latest)
        return :latest_cesr
    elseif ring in (:legacy, :legacy_cesr, :historical)
        return :legacy
    end
    error("Unknown ring '$ring'; use :latest or explicitly opt into :legacy")
end

const ORBIT_ADAPTER_HERE = @__DIR__
const ORBIT_ADAPTER_PROJECT_ROOT = normpath(joinpath(ORBIT_ADAPTER_HERE, "..", "..", ".."))
const ORBIT_ADAPTER_LATEST_DIR = joinpath(ORBIT_ADAPTER_PROJECT_ROOT, "Latest_Lattice")
const ORBIT_ADAPTER_LATEST_LATTICE = joinpath(
    ORBIT_ADAPTER_LATEST_DIR,
    "latest_cesr_scibmad_repaired.jl",
)
const ORBIT_ADAPTER_LATEST_CONTROLS = joinpath(
    ORBIT_ADAPTER_LATEST_DIR,
    "bmad_control_tracking_reference",
    "controls.csv",
)
const ORBIT_ADAPTER_LEGACY_MODEL = joinpath(ORBIT_ADAPTER_PROJECT_ROOT, "cesr_model.jl")
const ORBIT_ADAPTER_WIGGLER_SOURCE = joinpath(
    ORBIT_ADAPTER_PROJECT_ROOT,
    "wigglers",
    "wiggler.jl",
)

# The potential is a stateless numerical function, unlike the mutable control
# globals and lattice elements that must remain namespace-local.  Loading a
# fresh copy from inside every dynamically-created lattice module makes that
# function newer than already-compiled BatchSolve/BeamTracking call sites and
# triggers a Julia world-age error during ForwardDiff tracking.  Define it once
# before the benchmark methods are compiled, then let isolated lattice modules
# import this stable implementation.
if !isdefined(@__MODULE__, :WigglerModels)
    Base.include(@__MODULE__, ORBIT_ADAPTER_WIGGLER_SOURCE)
end
const ORBIT_ADAPTER_WIGGLER_MODULE = getfield(@__MODULE__, :WigglerModels)

"""Metadata for one writable lattice control."""
struct OrbitControlSpec
    name::String
    key::String
    variable::String
    symbol::Symbol
    step::Float64
    relation_count::Int
    observation_count::Int
    plane::Symbol
end

"""A Dict-like registry that writes a value into a lattice module global."""
mutable struct OrbitControlRegistry <: AbstractDict{String,Any}
    owner::Module
    specs::Vector{OrbitControlSpec}
    by_name::Dict{String,Int}
    values::Dict{String,Any}
end

function OrbitControlRegistry(
    owner::Module,
    specs::AbstractVector{<:OrbitControlSpec},
)
    ordered = collect(specs)
    names = getfield.(ordered, :name)
    length(unique(names)) == length(names) ||
        error("Control names are not unique")
    by_name = Dict(name => index for (index, name) in enumerate(names))
    return OrbitControlRegistry(owner, ordered, by_name, Dict{String,Any}())
end

Base.length(registry::OrbitControlRegistry) = length(registry.specs)
Base.keys(registry::OrbitControlRegistry) = getfield.(registry.specs, :name)
Base.values(registry::OrbitControlRegistry) =
    [get(registry.values, spec.name, nothing) for spec in registry.specs]
Base.haskey(registry::OrbitControlRegistry, name) =
    haskey(registry.by_name, String(name))
Base.getindex(registry::OrbitControlRegistry, name::AbstractString) =
    registry.values[String(name)]

function Base.iterate(registry::OrbitControlRegistry, state::Int=1)
    state > length(registry.specs) && return nothing
    spec = registry.specs[state]
    return (spec.name => get(registry.values, spec.name, nothing), state + 1)
end

function _control_spec(registry::OrbitControlRegistry, name::AbstractString)
    index = get(registry.by_name, String(name), 0)
    index > 0 || throw(KeyError(String(name)))
    return registry.specs[index]
end

"""Assign a scalar, BatchParam, or GTPSA value without narrowing its type."""
function Base.setindex!(
    registry::OrbitControlRegistry,
    value,
    name::AbstractString,
)
    spec = _control_spec(registry, name)
    # Expr can hold an arbitrary Julia object as a literal.  This is important
    # for BatchParam arrays and GTPSA parameter objects, which must remain
    # typed values in the DefExpr closure rather than Float64 approximations.
    Core.eval(registry.owner, Expr(:(=), spec.symbol, value))
    registry.values[spec.name] = value
    return registry
end

Base.setindex!(registry::OrbitControlRegistry, value, name::Symbol) =
    setindex!(registry, value, String(name))

"""Set every known control to `value` and return the registry."""
function set_all_controls!(registry::OrbitControlRegistry, value)
    for spec in registry.specs
        registry[spec.name] = value
    end
    return registry
end

function _split_csv_line(line::AbstractString)
    # The generated control metadata contains no quoted commas.  Keep this
    # small parser local so the adapter does not add a CSV package dependency.
    return split(chomp(line), ','; keepempty=true)
end

function _control_plane(
    name::AbstractString,
    key::AbstractString,
    variable::AbstractString="",
)
    variable_name = uppercase(String(variable))
    variable_name == "HKICK" && return :horizontal
    variable_name == "VKICK" && return :vertical
    uppercase(String(key)) == "GROUP" && return :group
    upper = uppercase(String(name))
    startswith(upper, "SK_") && return :skew
    startswith(upper, "V") && return :vertical
    startswith(upper, "H") && return :horizontal
    startswith(upper, "RAW_") && return :group
    return :other
end

function _adapter_config_value(config, key::Symbol, default=nothing)
    isnothing(config) && return default
    if config isa NamedTuple
        return haskey(config, key) ? getproperty(config, key) : default
    elseif config isa AbstractDict
        if haskey(config, key)
            return config[key]
        elseif haskey(config, String(key))
            return config[String(key)]
        end
    elseif hasproperty(config, key)
        return getproperty(config, key)
    end
    return default
end

"""Read ordered Overlay/Group metadata from the latest-lattice CSV."""
function read_control_specs(path::AbstractString=ORBIT_ADAPTER_LATEST_CONTROLS)
    lines = readlines(path)
    length(lines) >= 2 || error("Control metadata has no rows: $path")
    header = _split_csv_line(lines[1])
    header == [
        "lord_id",
        "lord_name",
        "lord_key",
        "variable",
        "step",
        "relation_count",
        "observation_count",
    ] || error("Unexpected control metadata header in $path")

    specs = OrbitControlSpec[]
    for (row, line) in enumerate(lines[2:end])
        fields = _split_csv_line(line)
        length(fields) == length(header) ||
            error("Control metadata row $row has $(length(fields)) fields")
        name = String(fields[2])
        key = String(fields[3])
        variable = uppercase(String(fields[4]))
        symbol = Symbol(name * "_" * lowercase(variable))
        push!(specs, OrbitControlSpec(
            name,
            key,
            variable,
            symbol,
            parse(Float64, fields[5]),
            parse(Int, fields[6]),
            parse(Int, fields[7]),
            _control_plane(name, key, variable),
        ))
    end
    return specs
end

function _coerce_control_spec(spec)
    spec isa OrbitControlSpec && return spec
    name = String(_adapter_config_value(
        spec,
        :name,
        _adapter_config_value(spec, :lord_name, ""),
    ))
    isempty(name) && error("Custom control spec has no name: $spec")
    key = String(_adapter_config_value(
        spec,
        :key,
        _adapter_config_value(spec, :lord_key, "Overlay"),
    ))
    variable = uppercase(String(_adapter_config_value(spec, :variable, "COMMAND")))
    symbol = Symbol(_adapter_config_value(
        spec,
        :symbol,
        name * "_" * lowercase(variable),
    ))
    return OrbitControlSpec(
        name,
        key,
        variable,
        symbol,
        Float64(_adapter_config_value(spec, :step, 0.0)),
        Int(_adapter_config_value(spec, :relation_count, 0)),
        Int(_adapter_config_value(spec, :observation_count, 0)),
        Symbol(lowercase(String(_adapter_config_value(
            spec,
            :plane,
            _control_plane(name, key, variable),
        )))),
    )
end

function _custom_specs(specs)
    converted = OrbitControlSpec[_coerce_control_spec(spec) for spec in specs]
    length(unique(getfield.(converted, :name))) == length(converted) ||
        error("Custom control names are not unique")
    return converted
end

"""Return labels and index groups without assuming detector/control counts."""
function orbit_ring_metadata(
    ring,
    specs::AbstractVector{<:OrbitControlSpec};
    ring_id::Symbol,
    lattice_path::AbstractString,
    detector_prefix::AbstractString="DET_",
    transverse_coordinate_indices=(x=1, y=3),
)
    prefix = uppercase(String(detector_prefix))
    detector_element_indices = Int[]
    detectors = String[]
    for (element_index, element) in enumerate(ring.line)
        name = uppercase(String(element.name))
        startswith(name, prefix) || continue
        push!(detector_element_indices, element_index)
        push!(detectors, name)
    end
    length(unique(detectors)) == length(detectors) ||
        error("Detector names are not unique")
    length(unique(detector_element_indices)) == length(detector_element_indices) ||
        error("Detector element indices are not unique")
    labels = vcat(
        ["$name:x" for name in detectors],
        ["$name:y" for name in detectors],
    )
    length(unique(labels)) == length(labels) ||
        error("Observable labels are not unique")
    names = getfield.(specs, :name)
    control_plane = Dict(spec.name => spec.plane for spec in specs)
    groups = Dict{Symbol,Vector{String}}(
        :all => copy(names),
        :horizontal => [spec.name for spec in specs if spec.plane == :horizontal],
        :vertical => [spec.name for spec in specs if spec.plane == :vertical],
        :skew => [spec.name for spec in specs if spec.plane == :skew],
        :group => [spec.name for spec in specs if spec.plane == :group],
        :overlay => [spec.name for spec in specs if uppercase(spec.key) == "OVERLAY"],
    )
    # The default orbit-input registry follows the lattice/control CSV order.
    # Plane-specific groups retain their filtered order, but concatenating H
    # then V here would silently reorder interleaved exports.
    steering_names = [
        spec.name for spec in specs if spec.plane in (:horizontal, :vertical)
    ]
    return (
        ring_id=ring_id,
        lattice_path=abspath(lattice_path),
        control_specs=collect(specs),
        control_names=names,
        all_control_names=names,
        steering_control_names=steering_names,
        control_plane=control_plane,
        control_groups=groups,
        detector_names=detectors,
        detector_order=copy(detectors),
        detector_element_indices=detector_element_indices,
        detector_prefix=prefix,
        detector_count=length(detectors),
        observable_labels=labels,
        observable_count=length(labels),
        phase_space_dimension=6,
        coordinate_labels=("x", "px", "y", "py", "z", "pz"),
        transverse_coordinate_indices=(
            x=Int(_adapter_config_value(transverse_coordinate_indices, :x, 1)),
            y=Int(_adapter_config_value(transverse_coordinate_indices, :y, 3)),
        ),
        scibmad_version=string(Base.pkgversion(SciBmad)),
        branch=0,
    )
end

"""Install the include helper and aliases needed by an isolated lattice module."""
function _new_lattice_namespace(lattice_path::AbstractString)
    owner = Module(gensym(:OrbitLattice), true, true)
    # A bare dynamically-created module has no one-argument include binding.
    # Intercept only the known stateless wiggler source so the latest-lattice
    # support file reuses the stable module above. All other includes remain
    # private to this lattice namespace.
    wiggler_source = normpath(abspath(ORBIT_ADAPTER_WIGGLER_SOURCE))
    wiggler_module = ORBIT_ADAPTER_WIGGLER_MODULE
    Core.eval(owner, :(const WigglerModels = $wiggler_module))
    Core.eval(owner, quote
        function include(path::AbstractString)
            Base.normpath(Base.abspath(path)) == $wiggler_source &&
                return $wiggler_module
            return Base.include($owner, path)
        end
    end)
    Core.eval(owner, :(using Beamlines))
    Core.eval(owner, :(using SciBmad))
    if !isdefined(Beamlines, :PhaseRef) && isdefined(Beamlines, :PhaseReference)
        Core.eval(owner, :(const PhaseRef = Beamlines.PhaseReference))
    elseif isdefined(Beamlines, :PhaseRef)
        Core.eval(owner, :(const PhaseRef = Beamlines.PhaseRef))
    end
    Base.include(owner, lattice_path)
    return owner
end

function _latest_model(
    ;
    zero_value,
    rf_on::Union{Nothing,Bool},
    rf_voltage::Real,
    lattice_path::AbstractString,
    ring_symbol::Symbol,
    control_metadata_path::AbstractString,
    control_specs=nothing,
    detector_prefix::AbstractString="DET_",
    ring_id::Symbol=:latest_cesr,
    transverse_coordinate_indices=(x=1, y=3),
)
    owner = _new_lattice_namespace(lattice_path)
    Base.invokelatest(isdefined, owner, ring_symbol) ||
        error("Latest/custom lattice does not define ring symbol $ring_symbol")
    ring = Base.invokelatest(getfield, owner, ring_symbol)
    specs = isnothing(control_specs) ?
        read_control_specs(control_metadata_path) : _custom_specs(control_specs)
    symbols = getfield.(specs, :symbol)
    length(unique(symbols)) == length(symbols) ||
        error("Control symbols are not unique")
    for spec in specs
        Base.invokelatest(isdefined, owner, spec.symbol) ||
            error("Lattice module does not define control symbol $(spec.symbol) for $(spec.name)")
    end
    controls = OrbitControlRegistry(owner, specs)
    set_all_controls!(controls, zero_value)

    if !isnothing(rf_on)
        selected_voltage = rf_on ? Float64(rf_voltage) : 0.0
        cavities = [
            element for element in ring.line
            if uppercase(String(element.kind)) in ("RFCAVITY", "CRABCAVITY")
        ]
        for cavity in cavities
            cavity.voltage = selected_voltage
        end
    end

    metadata = orbit_ring_metadata(
        ring,
        specs;
        ring_id,
        lattice_path,
        detector_prefix,
        transverse_coordinate_indices,
    )
    metadata = merge(
        metadata,
        (;
            rf_on,
            rf_voltage=isnothing(rf_on) ? nothing :
                (rf_on ? Float64(rf_voltage) : 0.0),
        ),
    )
    return (; ring, controls, metadata, namespace=owner)
end

function _legacy_model(; zero_value, rf_on, rf_voltage, ring_id::Symbol=:legacy_cesr, kwargs...)
    if !isdefined(Main, :load_cesr_model)
        Base.include(Main, ORBIT_ADAPTER_LEGACY_MODEL)
    end
    # The first call can dynamically include `cesr_model.jl` above.  Calling
    # the newly defined method through the current world age avoids the
    # MethodError seen when `ring=:legacy` is selected for the first time.
    legacy_loader = Base.invokelatest(getfield, Main, :load_cesr_model)
    model = Base.invokelatest(legacy_loader; zero_value, rf_on, rf_voltage)
    specs = OrbitControlSpec[]
    # The historical loader exposes the control names directly.  Keep the
    # registry read/write contract while retaining its historical ordering.
    for name in keys(model.controls.overlays)
        push!(specs, OrbitControlSpec(
            String(name),
            "Overlay",
            "COMMAND",
            Symbol(name),
            0.0,
            0,
            0,
            _control_plane(String(name), "Overlay"),
        ))
    end
    for name in keys(model.controls.groups)
        push!(specs, OrbitControlSpec(
            String(name),
            "Group",
            "COMMAND",
            Symbol(name),
            0.0,
            0,
            0,
            :group,
        ))
    end
    # Legacy controls are already typed scalar objects; mirror them in a
    # normal Dict-like registry only for metadata/consumer compatibility.
    metadata = orbit_ring_metadata(
        model.ring,
        specs;
        ring_id,
        lattice_path=ORBIT_ADAPTER_LEGACY_MODEL,
    )
    metadata = merge(
        metadata,
        (;
            rf_on,
            rf_voltage=isnothing(rf_on) ? nothing :
                (rf_on ? Float64(rf_voltage) : 0.0),
        ),
    )
    return merge(model, (; metadata))
end

"""
    load_ring_model(; ring=:latest, zero_value=0.0, rf_on=true, ...)

Load an isolated model.  `ring=:latest` is the maintained repaired lattice;
`ring=:legacy` is opt-in historical reproduction.  The latest loader accepts
Float64, `BatchParam`, and GTPSA values through `zero_value` and the returned
`controls` registry.  For first-order latest-ring maps, keep the complete
registry primitive and parameterize the normal H/V steering subset with
`parameterized_ring_model`; some skew/group combined-multipole controls are
not valid GTPSA inputs at their zero operating point.
"""
function load_ring_model(
    ;
    ring::Symbol=:latest,
    zero_value=0.0,
    rf_on::Union{Nothing,Bool}=nothing,
    rf_voltage::Real=1.5e6,
    lattice_path::AbstractString=ORBIT_ADAPTER_LATEST_LATTICE,
    ring_symbol::Symbol=:cesr,
    control_metadata_path::AbstractString=ORBIT_ADAPTER_LATEST_CONTROLS,
    control_specs=nothing,
    detector_prefix::AbstractString="DET_",
    ring_id::Union{Nothing,Symbol}=nothing,
    transverse_coordinate_indices=(x=1, y=3),
    kwargs...,
)
    artifact_ring = canonical_ring_id(ring)
    resolved_ring_id = isnothing(ring_id) ? artifact_ring : ring_id
    if artifact_ring == :latest_cesr
        return _latest_model(
            ;
            zero_value,
            rf_on,
            rf_voltage,
            lattice_path,
            ring_symbol,
            control_metadata_path,
            control_specs,
            detector_prefix,
            ring_id=resolved_ring_id,
            transverse_coordinate_indices,
        )
    elseif artifact_ring == :legacy
        return _legacy_model(; zero_value, rf_on, rf_voltage, ring_id=resolved_ring_id, kwargs...)
    end
    error("Unknown ring '$ring'; use :latest or explicitly opt into :legacy")
end

"""Assign only the requested controls to a typed first-order model.

`load_ring_model` deliberately initializes every control with the primitive
`zero_value` before this helper assigns the selected controls.  In particular,
callers constructing a GTPSA response model should leave `zero_value=0.0` and
assign only the steering subset to `GTPSA` parameters.  Initializing the full
registry with `zero(parameters[1])` needlessly promotes inactive Overlay/Group
controls and can send otherwise-unused element expressions outside the
primitive numerical domain.  This is a control-parameterization issue, not a
lattice or wiggler limitation.

The model factory may be `load_ring_model` itself or a closure capturing a
custom ring/lattice configuration.  The returned model owns its isolated
lattice namespace, so parameterized models do not share live control globals.
"""
function parameterize_controls!(
    model,
    names::AbstractVector,
    parameters;
)
    length(names) == length(parameters) || error(
        "The parameterized control name count ($(length(names))) does not match " *
        "the parameter count ($(length(parameters)))",
    )
    names_as_strings = String.(names)
    length(unique(names_as_strings)) == length(names_as_strings) ||
        error("Parameterized control names are not unique")
    for (index, name) in enumerate(names_as_strings)
        haskey(model.controls, name) ||
            error("Control '$name' is not present in the selected ring registry")
        model.controls[name] = parameters[index]
    end
    return model
end

"""Build an isolated ring model with only `names` represented by parameters.

This is the shared entry point for first-order GTPSA maps.  It keeps all
unselected controls at `zero_value` (Float64 by default) and supports the same
`Float64`, `BatchParam`, and GTPSA value types as `load_ring_model`.
"""
function parameterized_ring_model(
    names::AbstractVector,
    parameters;
    model_factory=load_ring_model,
    zero_value=0.0,
    rf_on::Union{Nothing,Bool}=true,
    kwargs...,
)
    model = model_factory(; zero_value, rf_on, kwargs...)
    return parameterize_controls!(model, names, parameters)
end

"""Compact config consumed by orbit error-analysis helpers."""
function ring_analysis_config(; ring::Symbol=:latest, kwargs...)
    static_kwargs = (; kwargs...)
    model = load_ring_model(; ring, zero_value=0.0, rf_on=true, static_kwargs...)
    metadata = model.metadata
    model_factory = (; runtime_kwargs...) -> begin
        merged_kwargs = merge(static_kwargs, (; runtime_kwargs...))
        load_ring_model(; ring, merged_kwargs...)
    end
    return merge(
        metadata,
        (
            ring=ring,
            model_factory,
            # Orbit studies steer only normal horizontal/vertical correctors
            # by default.  The complete control registry remains
            # available as `all_control_names` and through `control_groups`.
            control_names=metadata.steering_control_names,
            all_control_names=metadata.all_control_names,
            steering_control_names=metadata.steering_control_names,
            control_groups=metadata.control_groups,
            control_plane=metadata.control_plane,
            observable_labels=metadata.observable_labels,
            observable_planes=(
                x=collect(1:metadata.detector_count),
                y=collect(metadata.detector_count + 1:metadata.observable_count),
            ),
            coordinate_indices=metadata.transverse_coordinate_indices,
            detector_element_indices=metadata.detector_element_indices,
            detector_order=metadata.detector_order,
        ),
    )
end

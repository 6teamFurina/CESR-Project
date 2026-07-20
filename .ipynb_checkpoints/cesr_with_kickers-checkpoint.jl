module CESRWithKickers

import Beamlines

export CESRKickerSpec,
       KickTarget,
       apply_cesr_kickers!,
       cesr_kicker_spec,
       load_cesr_with_kickers

const _HERE = @__DIR__
const _SCIBMAD_LATTICE = joinpath(_HERE, "cesr_5p289gev.jl")
const _BMAD_LATTICE = joinpath(_HERE, "cesr.bmad")

# Compatibility shims for constructs emitted by the current Bmad exporter but not
# implemented by the Beamlines 0.9 version pinned by this tutorial.
Wiggler(; kwargs...) = Beamlines.Drift(; kwargs...)

function ELSeparator(; kwargs...)
    kept = Dict{Symbol,Any}(kwargs)
    pop!(kept, :En0, nothing)
    pop!(kept, :Es0, nothing)
    return Beamlines.Kicker(; kept...)
end

function RFCavity(; kwargs...)
    kept = Dict{Symbol,Any}(kwargs)
    is_on = pop!(kept, :is_on, true)
    is_on || (kept[:voltage] = 0.0)
    return Beamlines.RFCavity(; kept...)
end

function _skip_broken_group_expression(expression)
    text = repr(expression)
    if occursin("DefExpr", text) && occursin(r"RAW_[A-Za-z0-9_]+_command", text)
        return :(nothing)
    end
    return expression
end

# Group controllers are explicitly unsupported by the Bmad exporter, but it still
# emits DefExpr assignments that reference undefined RAW_*_command variables. The
# element declarations already contain Bmad's evaluated static strengths, so skip
# only those broken assignments while loading the template.
const _CESR_TEMPLATE = Base.include(
    _skip_broken_group_expression,
    @__MODULE__,
    _SCIBMAD_LATTICE,
)

struct KickTarget
    element::String
    plane::Symbol
    coefficient::Float64
end

struct DirectKick
    element::String
    kind::String
    hkick::Float64
    vkick::Float64
end

struct CESRKickerSpec
    overlays::Dict{String,Vector{KickTarget}}
    overlay_defaults::Dict{String,Float64}
    direct_kicks::Vector{DirectKick}
    groups::Dict{String,Vector{KickTarget}}
    group_defaults::Dict{String,Float64}
end

const _NUMBER_PATTERN = raw"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][-+]?\d+)?"

_canonical_name(name) = uppercase(String(name))
_plane(name) = uppercase(String(name)) == "HKICK" ? :H : :V
_parse_number(text) = parse(Float64, replace(strip(text), r"[dD]" => "e"))

function _logical_statements(path)
    statements = String[]
    current = ""
    for raw_line in eachline(path)
        line = strip(first(split(raw_line, '!'; limit=2)))
        isempty(line) && continue
        if occursin(r"^[A-Za-z][A-Za-z0-9_]*\s*:", line)
            isempty(current) || push!(statements, current)
            current = line
        elseif !isempty(current)
            current *= " " * line
        end
    end
    isempty(current) || push!(statements, current)
    return statements
end

function _coefficient(expression, variable)
    isnothing(expression) && return 1.0
    expr = replace(strip(expression), " " => "")
    uppercase(expr) == variable && return 1.0
    match_result = match(
        Regex("^(" * _NUMBER_PATTERN * ")\\*" * variable * raw"$", "i"),
        expr,
    )
    isnothing(match_result) && error("Unsupported Bmad control expression: $expression")
    return _parse_number(match_result.captures[1])
end

function _attribute_value(statement, attribute; default=0.0)
    result = match(
        Regex("\\b" * attribute * "\\s*=\\s*(" * _NUMBER_PATTERN * ")", "i"),
        statement,
    )
    return isnothing(result) ? default : _parse_number(result.captures[1])
end

function _parse_overlay_targets(body, variable)
    targets = KickTarget[]
    for item in split(body, ',')
        result = match(r"^\s*([A-Za-z0-9_]+)(?:\s*:\s*(.+))?\s*$", item)
        isnothing(result) && error("Cannot parse Bmad overlay target: $item")
        coefficient = _coefficient(result.captures[2], variable)
        push!(targets, KickTarget(_canonical_name(result.captures[1]), _plane(variable), coefficient))
    end
    return targets
end

function _parse_group_targets(body)
    targets = KickTarget[]
    for item in split(body, ',')
        result = match(
            r"^\s*([A-Za-z0-9_]+)\s*\[\s*(HKICK|VKICK)\s*\](?:\s*:\s*(.+))?\s*$"i,
            item,
        )
        isnothing(result) && continue
        coefficient = _coefficient(result.captures[3], "COMMAND")
        push!(
            targets,
            KickTarget(
                _canonical_name(result.captures[1]),
                _plane(result.captures[2]),
                coefficient,
            ),
        )
    end
    return targets
end

"""
    cesr_kicker_spec([path]) -> CESRKickerSpec

Parse the direct `HKICK`/`VKICK` values, corrector overlays, and kick-related
groups from the Bmad CESR lattice. Bmad control names are stored case-insensitively.
"""
function cesr_kicker_spec(path=_BMAD_LATTICE)
    isfile(path) || error("Bmad CESR lattice not found: $path")

    overlays = Dict{String,Vector{KickTarget}}()
    overlay_defaults = Dict{String,Float64}()
    direct_kicks = DirectKick[]
    groups = Dict{String,Vector{KickTarget}}()
    group_defaults = Dict{String,Float64}()

    for statement in _logical_statements(path)
        overlay_match = match(
            r"^([A-Za-z0-9_]+)\s*:\s*overlay\s*=\s*\{([^}]*)\}\s*,\s*var\s*=\s*\{(HKICK|VKICK)\}"i,
            statement,
        )
        if !isnothing(overlay_match)
            name = _canonical_name(overlay_match.captures[1])
            variable = uppercase(overlay_match.captures[3])
            overlays[name] = _parse_overlay_targets(overlay_match.captures[2], variable)
            overlay_defaults[name] = 0.0
            continue
        end

        group_match = match(
            r"^([A-Za-z0-9_]+)\s*:\s*group\s*=\s*\{([^}]*)\}\s*,\s*var\s*=\s*\{COMMAND\}"i,
            statement,
        )
        if !isnothing(group_match)
            targets = _parse_group_targets(group_match.captures[2])
            if !isempty(targets)
                name = _canonical_name(group_match.captures[1])
                groups[name] = targets
                group_defaults[name] = 0.0
            end
            continue
        end

        element_match = match(r"^([A-Za-z0-9_]+)\s*:\s*([A-Za-z_]+)\s*,"i, statement)
        isnothing(element_match) && continue
        hkick = _attribute_value(statement, "HKICK")
        vkick = _attribute_value(statement, "VKICK")
        (hkick == 0.0 && vkick == 0.0) && continue
        push!(
            direct_kicks,
            DirectKick(
                _canonical_name(element_match.captures[1]),
                lowercase(element_match.captures[2]),
                hkick,
                vkick,
            ),
        )
    end

    # Honor explicit controller assignments if a future CESR file contains them.
    assignment_pattern = Regex(
        "^\\s*([A-Za-z0-9_]+)\\s*\\[\\s*(HKICK|VKICK|COMMAND)\\s*\\]\\s*=\\s*(" *
        _NUMBER_PATTERN * ")",
        "i",
    )
    for raw_line in eachline(path)
        line = first(split(raw_line, '!'; limit=2))
        result = match(assignment_pattern, line)
        isnothing(result) && continue
        name = _canonical_name(result.captures[1])
        value = _parse_number(result.captures[3])
        if haskey(overlay_defaults, name)
            overlay_defaults[name] = value
        elseif haskey(group_defaults, name)
            group_defaults[name] = value
        end
    end

    return CESRKickerSpec(overlays, overlay_defaults, direct_kicks, groups, group_defaults)
end

function _element_lookup(ring)
    lookup = Dict{String,Vector{Any}}()
    for element in ring.line
        name = _canonical_name(element.name)
        elements = get!(lookup, name, Any[])
        any(other -> other === element, elements) || push!(elements, element)
    end
    return lookup
end

function _normalized_values(values)
    result = Dict{String,Float64}()
    for (name, value) in pairs(values)
        result[_canonical_name(name)] = Float64(value)
    end
    return result
end

const _APPLIED_KICKS = IdDict{Any,IdDict{Any,Tuple{Float64,Float64}}}()

function _remove_previous_kicks!(ring)
    previous = pop!(_APPLIED_KICKS, ring, nothing)
    isnothing(previous) && return
    for (element, (delta_kn0l, delta_ks0l)) in previous
        element.Kn0L -= delta_kn0l
        element.Ks0L -= delta_ks0l
    end
end

function _add_kick!(applied, element, plane, bmad_angle; charge_sign=nothing)
    # SciBmad: delta(px) = -Kn0L and delta(py) = Ks0L.
    # A magnetic Bmad HKICK/VKICK is the desired delta(px)/delta(py).
    # An ELSeparator kick is specified for a positive particle, hence charge_sign.
    if isnothing(charge_sign)
        delta_kn0l = plane == :H ? -bmad_angle : 0.0
        delta_ks0l = plane == :V ? bmad_angle : 0.0
    else
        delta_kn0l = plane == :H ? -charge_sign * bmad_angle : 0.0
        delta_ks0l = plane == :V ? charge_sign * bmad_angle : 0.0
    end

    element.Kn0L += delta_kn0l
    element.Ks0L += delta_ks0l
    old_kn0l, old_ks0l = get(applied, element, (0.0, 0.0))
    applied[element] = (old_kn0l + delta_kn0l, old_ks0l + delta_ks0l)
end

"""
    apply_cesr_kickers!(ring; kicker_values=Dict(), group_commands=Dict(),
                        spec=cesr_kicker_spec())

Apply the CESR Bmad corrector model to a SciBmad ring.

`kicker_values` overrides Bmad overlay knobs such as `"H01W"` or `"V01W"` in
radians. `group_commands` overrides kick-related group knobs such as
`"RAW_PRETZING_1"`. Names are case-insensitive.

The function is idempotent for a given ring: kicks from a previous call are removed
before the new settings are applied. It returns a report named tuple.
"""
function apply_cesr_kickers!(
    ring;
    kicker_values=Dict{String,Float64}(),
    group_commands=Dict{String,Float64}(),
    spec=cesr_kicker_spec(),
)
    _remove_previous_kicks!(ring)
    lookup = _element_lookup(ring)
    knob_values = copy(spec.overlay_defaults)
    merge!(knob_values, _normalized_values(kicker_values))
    commands = copy(spec.group_defaults)
    merge!(commands, _normalized_values(group_commands))

    missing_knobs = sort!(collect(setdiff(keys(_normalized_values(kicker_values)), keys(spec.overlays))))
    isempty(missing_knobs) || error("Unknown CESR kicker knobs: $(join(missing_knobs, ", "))")
    missing_groups = sort!(collect(setdiff(keys(_normalized_values(group_commands)), keys(spec.groups))))
    isempty(missing_groups) || error("Unknown CESR kicker groups: $(join(missing_groups, ", "))")

    direct_group_kicks = Dict{Tuple{String,Symbol},Float64}()
    for (group_name, targets) in spec.groups
        command = commands[group_name]
        for target in targets
            delta = target.coefficient * command
            if haskey(knob_values, target.element)
                knob_values[target.element] += delta
            else
                key = (target.element, target.plane)
                direct_group_kicks[key] = get(direct_group_kicks, key, 0.0) + delta
            end
        end
    end

    applied = IdDict{Any,Tuple{Float64,Float64}}()
    missing_elements = String[]
    charge_sign = sign(Beamlines.chargeof(ring.species_ref))

    # The pinned Beamlines version has no electric multipoles. ELSeparators were
    # loaded as Kicker elements, so reproduce their reference-particle kicks here.
    for direct in spec.direct_kicks
        elements = get(lookup, direct.element, Any[])
        if isempty(elements)
            push!(missing_elements, direct.element)
            continue
        end
        separator_sign = direct.kind == "elseparator" ? charge_sign : nothing
        for element in elements
            direct.hkick == 0.0 || _add_kick!(applied, element, :H, direct.hkick; charge_sign=separator_sign)
            direct.vkick == 0.0 || _add_kick!(applied, element, :V, direct.vkick; charge_sign=separator_sign)
        end
    end

    for ((element_name, plane), angle) in direct_group_kicks
        elements = get(lookup, element_name, Any[])
        if isempty(elements)
            push!(missing_elements, element_name)
            continue
        end
        direct = findfirst(item -> item.element == element_name, spec.direct_kicks)
        separator_sign = !isnothing(direct) && spec.direct_kicks[direct].kind == "elseparator" ? charge_sign : nothing
        for element in elements
            angle == 0.0 || _add_kick!(applied, element, plane, angle; charge_sign=separator_sign)
        end
    end

    for (knob_name, targets) in spec.overlays
        value = knob_values[knob_name]
        for target in targets
            elements = get(lookup, target.element, Any[])
            if isempty(elements)
                push!(missing_elements, target.element)
                continue
            end
            angle = target.coefficient * value
            for element in elements
                angle == 0.0 || _add_kick!(applied, element, target.plane, angle)
            end
        end
    end

    unique!(missing_elements)
    isempty(missing_elements) || error(
        "CESR kicker targets missing from the SciBmad ring: $(join(sort!(missing_elements), ", "))",
    )
    _APPLIED_KICKS[ring] = applied

    return (
        overlay_knobs=length(spec.overlays),
        overlay_targets=sum(length, values(spec.overlays)),
        direct_kick_elements=length(spec.direct_kicks),
        kick_groups=length(spec.groups),
        nonzero_overlay_knobs=count(!iszero, values(knob_values)),
        nonzero_group_commands=count(!iszero, values(commands)),
        elements_modified=length(applied),
    )
end

"""
    load_cesr_with_kickers(; kicker_values=Dict(), group_commands=Dict())

Return an independent copy of the exported CESR ring with Bmad direct and overlay
kicks applied. Use `kicker_values` and `group_commands` to set optimization knobs.
"""
function load_cesr_with_kickers(;
    kicker_values=Dict{String,Float64}(),
    group_commands=Dict{String,Float64}(),
    spec=cesr_kicker_spec(),
)
    # Rebuild the Beamline so every copied element receives fresh BeamlineParams.
    # A plain deepcopy copies the elements but leaves their reference-data links unset.
    elements = Beamlines.deepcopy_no_beamline.(_CESR_TEMPLATE.line)
    ring = Beamlines.Beamline(
        elements;
        p_over_q_ref=_CESR_TEMPLATE.p_over_q_ref,
        species_ref=_CESR_TEMPLATE.species_ref,
    )
    report = apply_cesr_kickers!(
        ring;
        kicker_values=kicker_values,
        group_commands=group_commands,
        spec=spec,
    )
    return ring, report
end

end # module CESRWithKickers

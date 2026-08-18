using Beamlines
using Dates
using GTPSA
using LinearAlgebra
using Printf
using Random
using SciBmad
using Statistics
using TOML

const EXACT11_HERE = @__DIR__
const AFFINITY_HERE = normpath(joinpath(EXACT11_HERE, ".."))
const PROJECT_DIR = normpath(joinpath(AFFINITY_HERE, "..", "..", "..", ".."))

include(joinpath(AFFINITY_HERE, "generate_scibmad_affinity_responses.jl"))

function parse_exact11_options(defaults, args)
    options = copy(defaults)
    for argument in args
        startswith(argument, "--") || error("Arguments must have --name=value form: $argument")
        fields = split(argument[3:end], '='; limit=2)
        length(fields) == 2 || error("Missing value in argument: $argument")
        haskey(options, fields[1]) || error("Unknown option: --$(fields[1])")
        options[fields[1]] = fields[2]
    end
    return options
end

function corrector_axis(name, label, length_m)
    length_m > 0 || return nothing
    upper_name = uppercase(name)
    upper_label = uppercase(label)
    if startswith(upper_name, "HX") || startswith(upper_name, "HKICK") ||
       endswith(upper_name, "_H") || occursin("HORZ", upper_label) ||
       occursin("CHICANE", upper_name) || occursin("DOG_LEG", upper_name)
        return :Kn0
    elseif startswith(upper_name, "VX") || startswith(upper_name, "VKICK") ||
           endswith(upper_name, "_V") || occursin("VERT", upper_label)
        return :Ks0
    end
    return nothing
end

function independent_corrector_inventory(ring)
    grouped = Dict{Tuple{String,Symbol},Vector{NamedTuple}}()
    order = Tuple{String,Symbol}[]
    for (index, element) in enumerate(ring.line)
        string(element.kind) == "Kicker" || continue
        name = String(base_name(element))
        length_m = constant_term(element.L)
        axis = corrector_axis(name, String(element.label), length_m)
        isnothing(axis) && continue
        key = (name, axis)
        haskey(grouped, key) || push!(order, key)
        push!(get!(grouped, key, NamedTuple[]), (;
            index,
            runtime_name=String(element.name),
            length_m,
            original=axis == :Kn0 ? element.Kn0 : element.Ks0,
        ))
    end
    controls = NamedTuple[]
    for key in order
        name, axis = key
        slices = grouped[key]
        push!(controls, (;
            name,
            axis,
            indices=getproperty.(slices, :index),
            runtime_names=getproperty.(slices, :runtime_name),
            total_length_m=sum(getproperty.(slices, :length_m)),
            originals=getproperty.(slices, :original),
        ))
    end
    isempty(controls) && error("No independent correctors found")
    return controls
end

function set_corrector_values!(ring, control, values)
    parameter = BatchParam(values)
    for index in control.indices
        if control.axis == :Kn0
            ring.line[index].Kn0 = parameter
        else
            ring.line[index].Ks0 = parameter
        end
    end
    return nothing
end

function set_corrector_scalar!(ring, control, value)
    for index in control.indices
        if control.axis == :Kn0
            ring.line[index].Kn0 = value
        else
            ring.line[index].Ks0 = value
        end
    end
    return nothing
end

function scalar_coordinate(value)
    array = Array(value)
    return ndims(array) == 1 ? Float64(array[1]) : Float64(array[1, 1])
end

function track_orbits_at_names(ring, closed, requested_names)
    requested = Set(uppercase.(requested_names))
    state_count = size(closed.v0, 1)
    horizontal = Dict(name => zeros(state_count) for name in requested)
    vertical = Dict(name => zeros(state_count) for name in requested)
    seen = Set{String}()
    bunch = Bunch(v=copy(closed.v0))
    SciBmad.BTBL.check_bl_bunch!(bunch, ring, false)
    for element in ring.line
        name = String(base_name(element))
        if name in requested && !(name in seen)
            coordinates = Array(bunch.coords.v)
            horizontal[name] .= coordinates[:, 1]
            vertical[name] .= coordinates[:, 3]
            push!(seen, name)
        end
        track!(bunch, element)
    end
    missing = setdiff(requested, seen)
    isempty(missing) || error("Requested orbit locations were not tracked: $(join(sort!(collect(missing)), ", "))")
    return (; horizontal, vertical)
end

function solve_batch_closed_orbit(ring, state_count; initial_v0=nothing)
    v0 = isnothing(initial_v0) ? zeros(state_count, 6) : copy(initial_v0)
    solution = find_closed_orbit(
        ring;
        v0,
        coasting_beam=false,
        batch=Val{true}(),
        warn=false,
    )
    converged = Array(solution.sol.retcode .== SciBmad.BatchSolve.RETCODE_SUCCESS)
    all(converged) || error("Batch closed-orbit solve failed for $(count(!, converged)) states")
    return solution
end

function read_five_candidate_sets(path)
    rows = read_simple_csv(path)
    result = Dict{String,Vector{String}}()
    for row in rows
        result[uppercase(row["sextupole"])] = [row["candidate_$index"] for index in 1:5]
    end
    return result
end

function read_bump_knobs(path, target)
    rows = read_simple_csv(path)
    selected = [row for row in rows if uppercase(row["target_sextupole"]) == uppercase(target)]
    isempty(selected) && error("No bump knob rows for $target in $path")
    return selected
end

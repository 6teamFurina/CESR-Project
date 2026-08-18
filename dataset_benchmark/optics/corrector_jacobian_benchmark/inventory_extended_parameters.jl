#!/usr/bin/env julia

"""Print the nominal CESR quadrupole/sextupole strength inventory."""

include(joinpath(@__DIR__, "..", "optics_common.jl"))

base_bmad_name(element) = first(split(uppercase(String(element.name)), '!'))

function strength_inventory(ring, attribute::Symbol)
    grouped = Dict{String,Vector{NamedTuple}}()
    for (index, element) in enumerate(ring.line)
        strength = Float64(GTPSA.scalar(Beamlines.deval(getproperty(element, attribute))))
        iszero(strength) && continue
        name = base_bmad_name(element)
        push!(get!(grouped, name, NamedTuple[]), (;
            index,
            runtime_name=String(element.name),
            kind=string(element.kind),
            strength,
            length_m=Float64(GTPSA.scalar(Beamlines.deval(element.L))),
        ))
    end
    return grouped
end

function print_inventory(label, inventory)
    println("$label base-name count: $(length(inventory))")
    println("$label runtime-slice count: $(sum(length, values(inventory)))")
    kinds = Dict{String,Int}()
    for entries in values(inventory), entry in entries
        kinds[entry.kind] = get(kinds, entry.kind, 0) + 1
    end
    println("$label kinds: $(sort!(collect(kinds)))")
    for name in sort!(collect(keys(inventory)))
        entries = inventory[name]
        strengths = unique(entry.strength for entry in entries)
        println(join((
            label,
            name,
            length(entries),
            join(unique(entry.kind for entry in entries), '|'),
            join(strengths, '|'),
            join((entry.runtime_name for entry in entries), '|'),
        ), ','))
    end
end

function main_inventory()
    model = load_cesr_model(rf_on=false)
    print_inventory("Kn1", strength_inventory(model.ring, :Kn1))
    print_inventory("Kn2", strength_inventory(model.ring, :Kn2))
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main_inventory())
end

using Beamlines
using Dates
using GTPSA
using LinearAlgebra
using Printf
using Random
using SciBmad
using Statistics
using TOML

const INVERSION_HERE = @__DIR__
const ALIGNMENT_STUDY_DIR = normpath(joinpath(INVERSION_HERE, ".."))
const PROJECT_DIR = normpath(joinpath(INVERSION_HERE, "..", ".."))
const RESPONSE_MAP_DIR = joinpath(ALIGNMENT_STUDY_DIR, "response_map")
const CORRECTOR_INPUT_CSV = joinpath(
    PROJECT_DIR,
    "orbit",
    "Orbit_Calculation",
    "inputs",
    "cesr_corrector_samples_1000.csv",
)

include(joinpath(PROJECT_DIR, "older_ring_version", "cesr_model.jl"))

const DETECTOR_COLUMNS = (
    :orbit_x,
    :orbit_y,
    :phi_1,
    :phi_2,
    :beta_1,
    :beta_2,
    :alpha_1,
    :alpha_2,
    :c11,
    :c12,
    :c21,
    :c22,
)

constant_term(value) = Float64(GTPSA.scalar(value))
base_name(element) = first(split(uppercase(String(element.name)), '!'))

function active_sextupole_inventory(ring)
    entries = NamedTuple[]
    s_m = 0.0
    for (index, element) in enumerate(ring.line)
        length_m = constant_term(Beamlines.deval(element.L))
        if string(element.kind) == "Sextupole"
            kn2 = constant_term(Beamlines.deval(element.Kn2))
            if !iszero(kn2)
                push!(entries, (;
                    index,
                    name=base_name(element),
                    runtime_name=String(element.name),
                    s_m=s_m + length_m / 2,
                    kn2_m3=kn2,
                    length_m,
                    x_offset_m=constant_term(Beamlines.deval(element.x_offset)),
                    y_offset_m=constant_term(Beamlines.deval(element.y_offset)),
                ))
            end
        end
        s_m += length_m
    end
    sort!(entries; by=entry -> entry.name)
    length(entries) == 76 ||
        error("Expected 76 active normal sextupoles, found $(length(entries))")
    return entries
end

function detector_elements(ring)
    elements = [
        element for element in ring.line
        if startswith(uppercase(String(element.name)), "DET_")
    ]
    length(elements) == 99 || error("Expected 99 detectors, found $(length(elements))")
    return elements
end

function corrector_names(path=CORRECTOR_INPUT_CSV)
    header = split(strip(first(eachline(path))), ',')
    first(header) == "sample_id" || error("Unexpected corrector CSV header")
    names = String.(header[2:end])
    length(names) == 119 || error("Expected 119 correctors, found $(length(names))")
    return names
end

function find_inventory_entry(inventory, name)
    normalized = uppercase(name)
    index = findfirst(entry -> entry.name == normalized, inventory)
    isnothing(index) && error("Unknown active sextupole: $name")
    return inventory[index]
end

function solve_rf_on_closed_orbit(ring)
    solution = find_closed_orbit(
        ring;
        coasting_beam=false,
        batch=Val{false}(),
        warn=true,
    )
    all(solution.sol.retcode .== SciBmad.BatchSolve.RETCODE_SUCCESS) ||
        error("RF-on closed-orbit solve failed: $(solution.sol.retcode)")
    return solution
end

function calculate_scalar_observables(model, target)
    closed = solve_rf_on_closed_orbit(model.ring)
    detectors = detector_elements(model.ring)
    target_element = model.ring.line[target.index]
    optics = twiss(
        model.ring;
        GTPSA_descriptor=Descriptor(6, 3),
        at=vcat(detectors, [target_element]),
        v0=closed.v0,
        v0_and_coast=(closed.v0, closed.coasting_beam),
        spin=false,
        RDTs=false,
        normalizing_map=false,
    )
    names = lowercase.(String.(optics.table.name))
    detector_indices = findall(name -> startswith(name, "det_"), names)
    length(detector_indices) == 99 ||
        error("Expected 99 detector rows in Twiss output, found $(length(detector_indices))")
    reference_index = findfirst(==("det_00w"), names)
    isnothing(reference_index) && error("DET_00W phase reference is missing")
    target_index = findfirst(name -> uppercase(first(split(name, '!'))) == target.name, names)
    isnothing(target_index) && error("Target $(target.name) is missing from Twiss output")

    rows = NamedTuple[]
    for index in detector_indices
        for column in DETECTOR_COLUMNS
            value = constant_term(getproperty(optics.table, column)[index])
            if column in (:phi_1, :phi_2)
                value -= constant_term(getproperty(optics.table, column)[reference_index])
            end
            push!(rows, (;
                observation_scope="detector",
                observation_name=names[index],
                observable=String(column),
                value,
            ))
        end
    end
    for tune_index in 1:2
        push!(rows, (;
            observation_scope="ring",
            observation_name="ring",
            observable="tune_$(tune_index)",
            value=constant_term(optics.tunes[tune_index]),
        ))
    end
    return (;
        rows,
        target_orbit_x_m=constant_term(optics.table.orbit_x[target_index]),
        target_orbit_y_m=constant_term(optics.table.orbit_y[target_index]),
        closed_orbit=closed,
    )
end

function csv_value(value)
    text = string(value)
    occursin(r"[,\"\n]", text) || return text
    return "\"" * replace(text, "\"" => "\"\"") * "\""
end

function write_rows(path, rows)
    isempty(rows) && error("Refusing to write empty CSV: $path")
    mkpath(dirname(path))
    columns = propertynames(first(rows))
    open(path, "w") do io
        println(io, join(string.(columns), ','))
        for row in rows
            println(io, join((csv_value(getproperty(row, column)) for column in columns), ','))
        end
    end
    return path
end

function parse_key_value_args(defaults, args)
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


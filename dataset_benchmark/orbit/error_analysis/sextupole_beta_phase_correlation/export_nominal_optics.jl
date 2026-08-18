#!/usr/bin/env julia

"""Export nominal RF-on Twiss data for the active sextupole sources and BPM detectors.

The maintained thick-element Hessian attribution inserts each complete-element
source at the element exit.  The sextupole Twiss point exported here therefore
uses the entrance of the next lattice element, which is the same longitudinal
boundary.  Detector markers are evaluated at their own entrance.
"""

using Beamlines
using Dates
using GTPSA
using Printf
using SciBmad
using TOML

const HERE = @__DIR__
const PROJECT_DIR = normpath(joinpath(HERE, "..", "..", "..", ".."))
include(joinpath(PROJECT_DIR, "cesr_model.jl"))

constant_term(value) = Float64(GTPSA.scalar(value))

function parse_options(args)
    options = Dict{String,String}(
        "output-dir" => joinpath(HERE, "results"),
    )
    for argument in args
        startswith(argument, "--") || error("Arguments must have --name=value form: $argument")
        fields = split(argument[3:end], '='; limit=2)
        length(fields) == 2 || error("Missing value in argument: $argument")
        haskey(options, fields[1]) || error("Unknown option: --$(fields[1])")
        options[fields[1]] = fields[2]
    end
    return options
end

function integrated_normal_sextupole_strength(element)
    kn2 = constant_term(Beamlines.deval(element.Kn2))
    kn2l = constant_term(Beamlines.deval(element.Kn2L))
    length_m = constant_term(Beamlines.deval(element.L))
    return iszero(kn2) ? kn2l : kn2 * length_m
end

function csv_value(value)
    text = string(value)
    occursin(r"[,\"\n]", text) || return text
    return "\"" * replace(text, "\"" => "\"\"") * "\""
end

function write_rows(path, rows)
    isempty(rows) && error("Refusing to write an empty CSV: $path")
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

function main(args=ARGS)
    options = parse_options(args)
    output_dir = abspath(options["output-dir"])

    model = load_cesr_model(zero_value=0.0, rf_on=true)
    ring = model.ring
    closed = find_closed_orbit(ring)
    descriptor = Descriptor(6, 1)
    optics = twiss(
        ring;
        GTPSA_descriptor=descriptor,
        at=:,
        v0=closed.v0,
        v0_and_coast=(closed.v0, closed.coasting_beam),
        spin=false,
        RDTs=false,
        normalizing_map=false,
    )
    table = optics.table

    # The first Twiss row for a beamline index is its entrance.  Thick elements
    # may have additional internal rows, so retain only the first occurrence.
    entrance_row = Dict{Int,Int}()
    for row in eachindex(table.beamline_index)
        index = Int(table.beamline_index[row])
        index > 0 || continue
        get!(entrance_row, index, row)
    end

    element_start_s = zeros(length(ring.line))
    element_exit_s = zeros(length(ring.line))
    s_m = 0.0
    for (index, element) in enumerate(ring.line)
        element_start_s[index] = s_m
        s_m += constant_term(Beamlines.deval(element.L))
        element_exit_s[index] = s_m
    end

    rows = NamedTuple[]
    sextupole_count = 0
    detector_count = 0
    for (index, element) in enumerate(ring.line)
        name = String(element.name)
        upper_name = uppercase(name)
        k2l = integrated_normal_sextupole_strength(element)
        if !iszero(k2l)
            index < length(ring.line) || error("Active sextupole at final lattice index")
            reference_index = index + 1
            haskey(entrance_row, reference_index) || error(
                "No Twiss entrance row for sextupole-exit reference index $reference_index",
            )
            table_row = entrance_row[reference_index]
            sextupole_count += 1
            push!(rows, (;
                point_type="sextupole_exit",
                element_index=index,
                element_name=name,
                s_m=element_exit_s[index],
                k2l_m2=k2l,
                beta_1_m=constant_term(table.beta_1[table_row]),
                beta_2_m=constant_term(table.beta_2[table_row]),
                phi_1_turn=constant_term(table.phi_1[table_row]),
                phi_2_turn=constant_term(table.phi_2[table_row]),
                twiss_reference_index=reference_index,
                twiss_reference_name=String(table.name[table_row]),
                twiss_reference_s_m=constant_term(table.s[table_row]),
            ))
        end
        if startswith(upper_name, "DET_")
            haskey(entrance_row, index) || error("No Twiss row for detector index $index")
            table_row = entrance_row[index]
            detector_count += 1
            push!(rows, (;
                point_type="detector",
                element_index=index,
                element_name=name,
                s_m=element_start_s[index],
                k2l_m2=0.0,
                beta_1_m=constant_term(table.beta_1[table_row]),
                beta_2_m=constant_term(table.beta_2[table_row]),
                phi_1_turn=constant_term(table.phi_1[table_row]),
                phi_2_turn=constant_term(table.phi_2[table_row]),
                twiss_reference_index=index,
                twiss_reference_name=String(table.name[table_row]),
                twiss_reference_s_m=constant_term(table.s[table_row]),
            ))
        end
    end
    sextupole_count == 76 || error("Expected 76 active normal sextupoles, found $sextupole_count")
    detector_count == 99 || error("Expected 99 detectors, found $detector_count")

    # SciBmad table phases are accumulated in turns.  The final rows contain
    # the full accumulated tunes, whereas optics.tunes are signed eigentunes.
    full_tune_1 = constant_term(table.phi_1[end])
    full_tune_2 = constant_term(table.phi_2[end])
    points_path = write_rows(joinpath(output_dir, "nominal_optics_points.csv"), rows)
    metadata = Dict(
        "format" => "cesr-sextupole-beta-phase-nominal-optics-v1",
        "date" => string(Dates.today()),
        "lattice_mode" => "RF-on nominal CESR",
        "phase_units" => "turn",
        "sextupole_reference" => "complete-element exit boundary",
        "detector_reference" => "detector-marker entrance",
        "active_normal_sextupoles" => sextupole_count,
        "detectors" => detector_count,
        "ring_length_m" => s_m,
        "full_tune_1_turn" => full_tune_1,
        "full_tune_2_turn" => full_tune_2,
        "signed_eigentune_1" => constant_term(optics.tunes[1]),
        "signed_eigentune_2" => constant_term(optics.tunes[2]),
        "points_csv" => points_path,
    )
    metadata_path = joinpath(output_dir, "nominal_optics_metadata.toml")
    mkpath(output_dir)
    open(metadata_path, "w") do io
        TOML.print(io, metadata; sorted=true)
    end
    @printf("Nominal optics: %s\n", points_path)
    @printf("Full accumulated tunes: Q1=%.12g, Q2=%.12g\n", full_tune_1, full_tune_2)
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end

#!/usr/bin/env julia

"""Compare every repaired SciBmad element with the Bmad branch-0 local map."""

using LinearAlgebra
using Printf
using SciBmad
using GTPSA
import Beamlines

const HERE = @__DIR__
const LATTICE_DIR = normpath(joinpath(HERE, ".."))
const REFERENCE = joinpath(LATTICE_DIR, "bmad_reference", "branch0", "local_maps.csv")
const OUTPUT = joinpath(LATTICE_DIR, "scibmad_validation", "scibmad_local_map_comparison.csv")
const SUMMARY = joinpath(LATTICE_DIR, "scibmad_validation", "LOCAL_MAP_COMPARISON.md")

include(joinpath(LATTICE_DIR, "latest_cesr_scibmad_repaired.jl"))

struct ReferenceElement
    index::Int
    name::String
    key::String
    length::Float64
    orbit_out::Vector{Float64}
    R::Matrix{Float64}
    vec0::Vector{Float64}
end

function read_reference(path)
    lines = readlines(path)
    header = split(first(lines), ',')
    expected = [
        "index", "name", "key", "length", "x", "px", "y", "py", "z", "pz",
        ["r$(row)$(column)" for row in 1:6 for column in 1:6]...,
        ["v$row" for row in 1:6]...,
    ]
    header == expected || error("Unexpected local-map CSV columns")

    result = ReferenceElement[]
    for line in @view lines[2:end]
        values = split(line, ',')
        length(values) == length(header) || error("Malformed local-map row: $line")
        R = Matrix{Float64}(undef, 6, 6)
        for row in 1:6, column in 1:6
            R[row, column] = parse(Float64, values[10 + 6(row - 1) + column])
        end
        push!(result, ReferenceElement(
            parse(Int, values[1]),
            values[2],
            values[3],
            parse(Float64, values[4]),
            parse.(Float64, values[5:10]),
            R,
            parse.(Float64, values[47:52]),
        ))
    end
    getproperty.(result, :index) == collect(0:length(result)-1) ||
        error("Bmad reference indices are not contiguous")
    return result
end

function normalized_name(name)
    return lowercase(replace(String(name), "#" => "!s"))
end

function element_kind(element)
    hasproperty(element, :kind) && return String(element.kind)
    return string(nameof(typeof(element)))
end

function local_map(element, entrance_orbit)
    # A harmon_master cavity derives its frequency from the containing ring.
    # Preserve that resolved full-ring frequency before placing a copy in the
    # one-element diagnostic beamline, whose circumference is intentionally
    # different.
    resolved_rf_frequency = if element_kind(element) == "RFCavity" && element.harmon_master
        Float64(element.rf_frequency)
    else
        nothing
    end
    copied = Beamlines.deepcopy_no_beamline(element)
    if !isnothing(resolved_rf_frequency)
        # Store an explicit rate before clearing harmon_master; otherwise the
        # setter tries to resolve the copied element's harmonic without a ring.
        copied.rf_frequency = resolved_rf_frequency
        copied.harmon_master = false
    end
    line = Beamlines.Beamline(
        [copied];
        p_over_q_ref=cesr.p_over_q_ref,
        species_ref=cesr.species_ref,
    )
    map = WigglerModels.gtpsa_transport_map(line; v0=entrance_orbit, order=1)
    R = Matrix(GTPSA.jacobian(map.v))
    output = GTPSA.scalar.(map.v)
    vec0 = output - R * entrance_orbit
    return R, vec0, output
end

function csv_quote(value)
    text = replace(string(value), '"' => "\"\"")
    return '"' * text * '"'
end

function main(args=ARGS)
    reference = read_reference(REFERENCE)
    length(reference) == length(cesr.line) + 1 || error(
        "Expected one BEGINNING plus $(length(cesr.line)) Bmad elements; " *
        "found $(length(reference)) rows",
    )
    max_elements = isempty(args) ? length(cesr.line) : parse(Int, only(args))
    n = min(max_elements, length(cesr.line))

    rows = NamedTuple[]
    cumulative_bmad = Matrix{Float64}(I, 6, 6)
    cumulative_scibmad = Matrix{Float64}(I, 6, 6)
    for i in 1:n
        bmad = reference[i + 1]
        entrance = reference[i].orbit_out
        element = cesr.line[i]
        R, vec0, output = local_map(element, entrance)
        delta_R = R - bmad.R
        cumulative_bmad = bmad.R * cumulative_bmad
        cumulative_scibmad = R * cumulative_scibmad
        worst = argmax(abs.(delta_R))
        name_match = normalized_name(element.name) == normalized_name(bmad.name)
        push!(rows, (
            index=i,
            bmad_name=bmad.name,
            scibmad_name=String(element.name),
            name_match,
            bmad_key=bmad.key,
            scibmad_kind=element_kind(element),
            bmad_length=bmad.length,
            scibmad_length=Float64(element.L),
            length_error=abs(Float64(element.L) - bmad.length),
            local_R_max=maximum(abs, delta_R),
            local_R_fro=norm(delta_R),
            local_R_relative_fro=norm(delta_R) / max(norm(bmad.R), eps()),
            local_R_row=worst[1],
            local_R_column=worst[2],
            vec0_max=maximum(abs, vec0 - bmad.vec0),
            orbit_out_max=maximum(abs, output - bmad.orbit_out),
            cumulative_R_max=maximum(abs, cumulative_scibmad - cumulative_bmad),
        ))
        (i == 1 || i % 100 == 0 || i == n) && println("Compared $i/$n elements")
    end

    columns = propertynames(first(rows))
    open(OUTPUT, "w") do io
        println(io, join(string.(columns), ','))
        for row in rows
            println(io, join(csv_quote.(getproperty.(Ref(row), columns)), ','))
        end
    end

    worst_local = sort(rows; by=row -> row.local_R_max, rev=true)[1:min(25, n)]
    worst_vec0 = sort(rows; by=row -> row.vec0_max, rev=true)[1:min(10, n)]
    max_length = maximum(row.length_error for row in rows)
    mismatched_names = count(!, getproperty.(rows, :name_match))

    open(SUMMARY, "w") do io
        println(io, "# Latest CESR local-map comparison")
        println(io)
        println(io, "- Bmad/Tao reference: `20260814-0`, branch 0, RF on.")
        println(io, "- Compared elements: `$n`.")
        println(io, "- Name mismatches after `#` to `!s` normalization: `$mismatched_names`.")
        @printf(io, "- Maximum element-length mismatch: `%.6e m`.\n", max_length)
        @printf(io, "- Maximum local-matrix entry mismatch: `%.6e`.\n", first(worst_local).local_R_max)
        @printf(io, "- Maximum affine-vector mismatch: `%.6e`.\n", first(worst_vec0).vec0_max)
        println(io)
        println(io, "## Largest local matrix discrepancies")
        println(io)
        println(io, "| Index | Bmad name | Bmad key | SciBmad kind | max abs dR | relative Frobenius | entry |")
        println(io, "|---:|---|---|---|---:|---:|---|")
        for row in worst_local
            @printf(
                io,
                "| %d | `%s` | `%s` | `%s` | %.6e | %.6e | R%d%d |\n",
                row.index, row.bmad_name, row.bmad_key, row.scibmad_kind,
                row.local_R_max, row.local_R_relative_fro,
                row.local_R_row, row.local_R_column,
            )
        end
        println(io)
        println(io, "## Largest affine-vector discrepancies")
        println(io)
        println(io, "| Index | Bmad name | Bmad key | SciBmad kind | max abs dvec0 | exit-orbit mismatch |")
        println(io, "|---:|---|---|---|---:|---:|")
        for row in worst_vec0
            @printf(
                io,
                "| %d | `%s` | `%s` | `%s` | %.6e | %.6e |\n",
                row.index, row.bmad_name, row.bmad_key, row.scibmad_kind,
                row.vec0_max, row.orbit_out_max,
            )
        end
        println(io)
        println(io, "## Split-wiggler interpretation")
        println(io)
        println(io, "The three `ID_S1A` slices must not be interpreted independently: the " *
                    "continuous four-potential uses canonical boundary terms that cancel " *
                    "across the complete five-element block. Run `../support_codes/compare_wiggler_block.jl` " *
                    "for the physical block comparison. With the Bmad reference-time patch, " *
                    "the block affine/exit mismatch is about `7.12e-15`; the remaining " *
                    "`5.89e-6` `R12` difference is reproduced by Bmad Runge-Kutta field " *
                    "tracking and is absent only from Bmad's standard wiggler approximation.")
    end

    println("Wrote $OUTPUT")
    println("Wrote $SUMMARY")
    @printf("Maximum length error: %.6e m\n", max_length)
    @printf("Maximum local R error: %.6e\n", first(worst_local).local_R_max)
    @printf("Maximum vec0 error: %.6e\n", first(worst_vec0).vec0_max)
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end

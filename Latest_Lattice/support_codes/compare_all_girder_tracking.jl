#!/usr/bin/env julia

"""Compare every reconstructed SciBmad girder parameter with Bmad tracking."""

using LinearAlgebra
using Printf
using SciBmad

const LATTICE_DIR = normpath(joinpath(@__DIR__, ".."))
const BMAD_REFERENCE_DIR = joinpath(LATTICE_DIR, "bmad_reference")
const SCIBMAD_VALIDATION_DIR = joinpath(LATTICE_DIR, "scibmad_validation")

include(joinpath(LATTICE_DIR, "latest_cesr_scibmad_repaired.jl"))

const REFERENCE = joinpath(BMAD_REFERENCE_DIR, "girder", "tracking_response.csv")
const OUTPUT = joinpath(SCIBMAD_VALIDATION_DIR, "scibmad_girder_tracking_comparison.csv")
const SUMMARY = joinpath(SCIBMAD_VALIDATION_DIR, "GIRDER_VALIDATION.md")
const PARTICLE_START = zeros(6)
const COORDINATES = (:x, :px, :y, :py, :z, :pz)

function read_csv(path)
    lines = readlines(path)
    header = Symbol.(split(first(lines), ','))
    return [NamedTuple{Tuple(header)}(Tuple(split(line, ','))) for line in @view lines[2:end]]
end

function track_observations(indices::Set{Int})
    bunch = Bunch(v=reshape(copy(PARTICLE_START), 1, 6))
    SciBmad.BTBL.check_bl_bunch!(bunch, cesr, false)
    result = Dict{Int,Vector{Float64}}()
    for (index, element) in enumerate(cesr.line)
        track!(bunch, element)
        index in indices && (result[index] = vec(Float64.(bunch.coords.v)))
    end
    Set(keys(result)) == indices || error("Missing girder observation indices")
    return result
end

function csv_quote(value)
    text = replace(string(value), '"' => "\"\"")
    return '"' * text * '"'
end

function main()
    reference = read_csv(REFERENCE)
    grouped = Dict{Tuple{String,String},Vector{eltype(reference)}}()
    for row in reference
        push!(get!(grouped, (row.girder, row.parameter), eltype(reference)[]), row)
    end

    comparison = NamedTuple[]
    for (calculation, key) in enumerate(sort!(collect(keys(grouped))))
        girder, parameter_text = key
        parameter = Symbol(parameter_text)
        rows = grouped[key]
        indices = Set(parse(Int, row.observation_index) for row in rows)
        step = parse(Float64, first(rows).step)
        plus_setting = NamedTuple{(parameter,)}((step,))
        minus_setting = NamedTuple{(parameter,)}((-step,))
        set_latest_girder!(cesr, girder; plus_setting...)
        plus = track_observations(indices)
        set_latest_girder!(cesr, girder; minus_setting...)
        minus = track_observations(indices)
        set_latest_girder!(cesr, girder)

        for row in rows
            index = parse(Int, row.observation_index)
            candidate = (plus[index] - minus[index]) / (2step)
            bmad = [parse(Float64, getproperty(row, Symbol("d", name))) for name in COORDINATES]
            difference = candidate - bmad
            reference_norm = norm(bmad)
            push!(comparison, (
                girder,
                parameter=parameter_text,
                observation_index=index,
                bmad_norm=reference_norm,
                scibmad_norm=norm(candidate),
                absolute_max=maximum(abs, difference),
                relative_l2=norm(difference) / max(reference_norm, 1e-14),
                cosine_similarity=dot(candidate, bmad) / max(norm(candidate) * reference_norm, 1e-28),
                dx=candidate[1], dpx=candidate[2], dy=candidate[3],
                dpy=candidate[4], dz=candidate[5], dpz=candidate[6],
                bmad_dx=bmad[1], bmad_dpx=bmad[2], bmad_dy=bmad[3],
                bmad_dpy=bmad[4], bmad_dz=bmad[5], bmad_dpz=bmad[6],
            ))
        end
        println("SciBmad girder tracking $calculation/$(length(grouped)) $girder $parameter_text")
    end

    columns = propertynames(first(comparison))
    open(OUTPUT, "w") do io
        println(io, join(string.(columns), ','))
        for row in comparison
            println(io, join(csv_quote.(getproperty.(Ref(row), columns)), ','))
        end
    end

    informative = filter(row -> row.bmad_norm > 1e-10, comparison)
    worst_relative = sort(informative; by=row -> row.relative_l2, rev=true)
    worst_absolute = sort(comparison; by=row -> row.absolute_max, rev=true)
    open(SUMMARY, "w") do io
        println(io, "# Latest CESR girder tracking validation")
        println(io)
        println(io, "- Girders tested: `12`.")
        println(io, "- Parameters per girder: all six Bmad offsets/pitches/tilt.")
        println(io, "- Bmad member elements represented: `150` tracking elements/slices.")
        println(io, "- Tracking-response observations: `$(length(comparison))`.")
        @printf(io, "- Maximum absolute response-entry difference: `%.6e`.\n", first(worst_absolute).absolute_max)
        @printf(io, "- Maximum informative relative L2 difference: `%.6e`.\n", first(worst_relative).relative_l2)
        println(io)
        println(io, "## Largest relative discrepancies")
        println(io)
        println(io, "| Girder | Parameter | Observation | Bmad norm | Relative L2 | Max abs | Cosine |")
        println(io, "|---|---|---:|---:|---:|---:|---:|")
        for row in first(worst_relative, min(20, length(worst_relative)))
            @printf(
                io, "| `%s` | `%s` | %d | %.6e | %.6e | %.6e | %.9f |\n",
                row.girder, row.parameter, row.observation_index, row.bmad_norm,
                row.relative_l2, row.absolute_max, row.cosine_similarity,
            )
        end
    end
    println("Wrote $OUTPUT")
    println("Wrote $SUMMARY")
    @printf("Maximum absolute difference: %.6e\n", first(worst_absolute).absolute_max)
    @printf("Maximum informative relative L2: %.6e\n", first(worst_relative).relative_l2)
end

main()

#!/usr/bin/env julia

"""Compare all 124 live SciBmad controls with Bmad one-pass responses."""

using LinearAlgebra
using Printf
using SciBmad

include(joinpath(@__DIR__, "latest_cesr_scibmad_repaired.jl"))

const REFERENCE_DIR = joinpath(@__DIR__, "bmad_control_tracking_reference")
const REFERENCE_CSV = joinpath(REFERENCE_DIR, "control_tracking_response.csv")
const CONTROLS_CSV = joinpath(REFERENCE_DIR, "controls.csv")
const OUTPUT_CSV = joinpath(@__DIR__, "scibmad_control_tracking_comparison.csv")
const SUMMARY = joinpath(@__DIR__, "FULL_CONTROL_VALIDATION.md")
const PARTICLE_START = [1.0e-3, 2.0e-4, -8.0e-4, 1.5e-4, 5.0e-4, 2.0e-4]
const COORDINATES = (:x, :px, :y, :py, :z, :pz)

function read_csv(path)
    lines = readlines(path)
    header = Symbol.(split(first(lines), ','))
    return [
        NamedTuple{Tuple(header)}(Tuple(split(line, ',')))
        for line in @view lines[2:end]
    ]
end

function track_observations(indices::Set{Int})
    bunch = Bunch(v=reshape(copy(PARTICLE_START), 1, 6))
    SciBmad.BTBL.check_bl_bunch!(bunch, cesr, false)
    result = Dict{Int,Vector{Float64}}()
    for (index, element) in enumerate(cesr.line)
        track!(bunch, element)
        index in indices || continue
        result[index] = vec(Float64.(bunch.coords.v))
    end
    Set(keys(result)) == indices || error("Not all observation indices were tracked")
    return result
end

function csv_quote(value)
    text = replace(string(value), '"' => "\"\"")
    return '"' * text * '"'
end

function sample_median(values)
    ordered = sort(collect(values))
    n = length(ordered)
    isodd(n) && return ordered[(n + 1) ÷ 2]
    return (ordered[n ÷ 2] + ordered[n ÷ 2 + 1]) / 2
end

function main()
    reference = read_csv(REFERENCE_CSV)
    controls = read_csv(CONTROLS_CSV)
    by_control = Dict{String,Vector{eltype(reference)}}()
    for row in reference
        push!(get!(by_control, row.lord_name, eltype(reference)[]), row)
    end

    comparison_rows = NamedTuple[]
    for (control_index, control) in enumerate(controls)
        rows = by_control[control.lord_name]
        indices = Set(parse(Int, row.observation_index) for row in rows)
        step = parse(Float64, control.step)
        variable = Symbol(control.lord_name * "_" * lowercase(control.variable))
        isdefined(Main, variable) || error("Missing exported control variable $variable")

        Core.eval(Main, Expr(:(=), variable, step))
        plus = track_observations(indices)
        Core.eval(Main, Expr(:(=), variable, -step))
        minus = track_observations(indices)
        Core.eval(Main, Expr(:(=), variable, 0.0))

        for row in rows
            index = parse(Int, row.observation_index)
            candidate = (plus[index] - minus[index]) / (2step)
            bmad = [parse(Float64, getproperty(row, Symbol("d", name))) for name in COORDINATES]
            difference = candidate - bmad
            reference_norm = norm(bmad)
            push!(comparison_rows, (
                lord_name=control.lord_name,
                lord_key=control.lord_key,
                variable=control.variable,
                observation_index=index,
                relation_count=parse(Int, control.relation_count),
                bmad_norm=reference_norm,
                scibmad_norm=norm(candidate),
                absolute_max=maximum(abs, difference),
                relative_l2=norm(difference) / max(reference_norm, 1e-14),
                cosine_similarity=dot(candidate, bmad) /
                    max(norm(candidate) * reference_norm, 1e-28),
                dx=candidate[1], dpx=candidate[2], dy=candidate[3],
                dpy=candidate[4], dz=candidate[5], dpz=candidate[6],
                bmad_dx=bmad[1], bmad_dpx=bmad[2], bmad_dy=bmad[3],
                bmad_dpy=bmad[4], bmad_dz=bmad[5], bmad_dpz=bmad[6],
            ))
        end
        println(
            "SciBmad control tracking $control_index/$(length(controls)) " *
            "$(control.lord_name) ($(length(rows)) observations)",
        )
    end

    columns = propertynames(first(comparison_rows))
    open(OUTPUT_CSV, "w") do io
        println(io, join(string.(columns), ','))
        for row in comparison_rows
            println(io, join(csv_quote.(getproperty.(Ref(row), columns)), ','))
        end
    end

    informative = filter(row -> row.bmad_norm > 1e-10, comparison_rows)
    worst_relative = sort(informative; by=row -> row.relative_l2, rev=true)
    worst_absolute = sort(comparison_rows; by=row -> row.absolute_max, rev=true)
    control_max = Dict{String,Float64}()
    for row in informative
        control_max[row.lord_name] = max(
            get(control_max, row.lord_name, 0.0), row.relative_l2,
        )
    end

    open(SUMMARY, "w") do io
        println(io, "# Full latest-CESR control tracking validation")
        println(io)
        println(io, "- Bmad/Tao reference: `20260814-0`, branch 0 changed to open geometry for a one-pass map.")
        println(io, "- Controls tested: `$(length(controls))` (all 119 Overlay and 5 Group lords).")
        println(io, "- Bmad lord-to-slave relationships covered: `347`.")
        println(io, "- Response observations compared: `$(length(comparison_rows))`.")
        println(io, "- Fixed nonzero particle start: `$(PARTICLE_START)`." )
        @printf(io, "- Maximum absolute six-vector entry difference: `%.6e`.\n", first(worst_absolute).absolute_max)
        @printf(io, "- Maximum relative L2 difference for informative responses: `%.6e`.\n", first(worst_relative).relative_l2)
        println(io)
        println(io, "## Largest relative discrepancies")
        println(io)
        println(io, "| Control | Type | Observation | Bmad norm | Relative L2 | Max abs | Cosine |")
        println(io, "|---|---|---:|---:|---:|---:|---:|")
        for row in first(worst_relative, min(20, length(worst_relative)))
            @printf(
                io, "| `%s` | `%s` | %d | %.6e | %.6e | %.6e | %.9f |\n",
                row.lord_name, row.lord_key, row.observation_index,
                row.bmad_norm, row.relative_l2, row.absolute_max,
                row.cosine_similarity,
            )
        end
        println(io)
        println(io, "The CSV contains both six-component response vectors for every observation.")
    end

    println("Wrote $OUTPUT_CSV")
    println("Wrote $SUMMARY")
    @printf("Maximum absolute difference: %.6e\n", first(worst_absolute).absolute_max)
    @printf("Maximum informative relative L2: %.6e\n", first(worst_relative).relative_l2)
    @printf("Median control maximum relative L2: %.6e\n", sample_median(values(control_max)))
end

main()

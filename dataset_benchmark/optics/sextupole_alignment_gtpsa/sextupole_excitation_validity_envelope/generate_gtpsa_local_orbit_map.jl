#!/usr/bin/env julia

"""Generate one target's GTPSA local-orbit map versus two bump knobs and delta K2.

The normalized parameters are

    qx = requested x-bump / bump_scale
    qy = requested y-bump / bump_scale
    qk = delta_Kn2 / k2_scale

The periodic closed orbit is solved order by order.  Taylor coefficients are
saved at the entrance and exit of the target sextupole for all six canonical
coordinates.  This script intentionally handles one target per process: a
lower-level GTPSA failure at one magnet must not discard maps already produced
for the other targets.
"""

using Beamlines
using Dates
using GTPSA
using LinearAlgebra
using Printf
using SciBmad
using TOML

# An override permits a read-only frozen lattice copy to be used if the active
# workspace lattice is being edited concurrently.  The selected path is saved
# in every map's metadata.
const LATEST_LATTICE = get(
    ENV,
    "CESR_VALIDITY_LATTICE",
    normpath(joinpath(@__DIR__, "..", "..", "..", "..", "Latest_Lattice", "latest_cesr_scibmad_repaired.jl")),
)
include(LATEST_LATTICE)

constant_term(value) = Float64(GTPSA.scalar(Beamlines.deval(value)))
base_name(element) = first(split(uppercase(String(element.name)), '!'))

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

function active_sextupole_inventory(ring)
    entries = NamedTuple[]
    s_m = 0.0
    for (index, element) in enumerate(ring.line)
        length_m = constant_term(element.L)
        if string(element.kind) == "Sextupole"
            kn2_m3 = constant_term(element.Kn2)
            if !iszero(kn2_m3)
                push!(entries, (;
                    index,
                    name=base_name(element),
                    s_m=s_m + length_m/2,
                    kn2_m3,
                    length_m,
                ))
            end
        end
        s_m += length_m
    end
    length(entries) == 76 || error("Expected 76 active sextupoles")
    return entries
end

function corrector_axis(name, label, length_m)
    length_m > 0 || return nothing
    upper_name, upper_label = uppercase(name), uppercase(label)
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
        axis = corrector_axis(name, String(element.label), constant_term(element.L))
        isnothing(axis) && continue
        key = (name, axis)
        haskey(grouped, key) || push!(order, key)
        push!(get!(grouped, key, NamedTuple[]), (;
            index,
            original=axis == :Kn0 ? element.Kn0 : element.Ks0,
        ))
    end
    return [
        (;
            name=key[1],
            axis=key[2],
            indices=getproperty.(grouped[key], :index),
            originals=getproperty.(grouped[key], :original),
        )
        for key in order
    ]
end

function solve_closed_orbit(ring)
    solution = find_closed_orbit(ring; coasting_beam=false, batch=Val{false}(), warn=false)
    all(solution.sol.retcode .== SciBmad.BatchSolve.RETCODE_SUCCESS) ||
        error("RF-on closed-orbit solve failed")
    return solution
end

function calculate_direct_transport(ring, at, closed; descriptor)
    _, names, _, step_save = SciBmad._twiss_1(ring, at)
    _, eye, _, _, _, _ = SciBmad._twiss_2(
        step_save, (closed.v0, closed.coasting_beam), descriptor,
        Val{false}(), Val{false}(),
    )
    maps = [zero(eye) for _ in step_save]
    first_saved = 1
    if !isempty(step_save) && first(step_save) == 0
        maps[1] = one(eye)
        SciBmad.NNF.setscalar!(maps[1], eye.v0)
        first_saved = 2
    end
    let saved_steps=step_save, saved_maps=maps,
        current_step=Ref{Int}(0), saved_index=Ref{Int}(first_saved)
        callback = (i, coordinates, current_s, current_t_ref, last_ds_step, last_g,
                    transforms_out!, transforms_in!) -> begin
            current_step[] += 1
            if saved_index[] <= length(saved_steps) && current_step[] == saved_steps[saved_index[]]
                transforms_out!(i, coordinates, current_s, current_t_ref)
                SciBmad.NNF.setray!(saved_maps[saved_index[]].v; v=reshape(coordinates.v, :))
                transforms_in!(i, coordinates, current_s, current_t_ref)
                saved_index[] += 1
            end
            nothing
        end
        SciBmad._twiss_4(eye, callback, ring)
    end
    return (; names, maps)
end

function read_simple_csv(path)
    lines = readlines(path)
    header = split(first(lines), ',')
    return [
        Dict(header[index] => fields[index] for index in eachindex(header))
        for line in lines[2:end] for fields in (split(line, ','),)
    ]
end

function read_bump_knobs(path, target)
    selected = [row for row in read_simple_csv(path) if uppercase(row["target_sextupole"]) == uppercase(target)]
    isempty(selected) && error("No bump knob for $target")
    return selected
end

csv_value(value) = occursin(r"[,\"\n]", string(value)) ?
    "\"" * replace(string(value), "\"" => "\"\"") * "\"" : string(value)

function write_rows(path, rows)
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

function write_metadata(path, metadata)
    mkpath(dirname(path))
    open(path, "w") do io
        TOML.print(io, metadata; sorted=true)
    end
    return path
end

parameter_only(map) = map ∘ zero(map)

function fixed_point_residual(one_turn, fixed_point)
    image = parameter_only(one_turn ∘ fixed_point)
    origin = parameter_only(fixed_point)
    nv = SciBmad.NNF.nvars(one_turn)
    return image.v[1:nv] .- origin.v[1:nv]
end

residual_norm(residual) = sum(GTPSA.normTPS(value) for value in residual)

function parameter_dependent_fixed_point(one_turn, parameter_order)
    fixed_point = normal(one_turn, 0)
    linear = Matrix(SciBmad.NNF.jacobian(one_turn, SciBmad.NNF.HVARS))
    correction_operator = inv(I - linear)
    history = Float64[]
    for _ in 1:(parameter_order + 1)
        residual = fixed_point_residual(one_turn, fixed_point)
        push!(history, residual_norm(residual))
        correction = correction_operator * residual
        for coordinate in eachindex(correction)
            GTPSA.add!(fixed_point.v[coordinate], fixed_point.v[coordinate], correction[coordinate])
        end
    end
    push!(history, residual_norm(fixed_point_residual(one_turn, fixed_point)))
    return fixed_point, history
end

const DEFAULT_BUMP_KNOBS = joinpath(
    @__DIR__, "..", "quadrupole_affinity", "exact_11_triplet_validation",
    "results", "bump_knobs", "local_bump_knobs.csv",
)

function monomial_powers(maximum_order)
    powers = NTuple{3,Int}[]
    for total in 0:maximum_order
        for px in 0:total
            for py in 0:(total - px)
                pk = total - px - py
                push!(powers, (px, py, pk))
            end
        end
    end
    return powers
end

function taylor_coefficient(value, powers)
    result = value
    for (parameter_offset, power) in enumerate(powers)
        for _ in 1:power
            result = GTPSA.deriv(result, 6 + parameter_offset)
        end
    end
    divisor = prod(factorial(power) for power in powers)
    return Float64(constant_term(result)) / divisor
end

function normalized_parameter_map(map)
    return map ∘ zero(map)
end

function main(args=ARGS)
    defaults = Dict(
        "target" => "SEX_09AW",
        "order" => "4",
        "bump-scale-m" => "1.0e-3",
        "k2-scale-m3" => "0.1",
        "bump-knobs-csv" => DEFAULT_BUMP_KNOBS,
        "output-root" => joinpath(@__DIR__, "results", "gtpsa_maps"),
        "overwrite" => "false",
    )
    options = parse_exact11_options(defaults, args)
    target_name = uppercase(strip(options["target"]))
    maximum_order = parse(Int, options["order"])
    bump_scale = parse(Float64, options["bump-scale-m"])
    k2_scale = parse(Float64, options["k2-scale-m3"])
    maximum_order in 2:5 || error("--order must be between 2 and 5")
    all(isfinite.((bump_scale, k2_scale))) || error("Scales must be finite")
    min(bump_scale, k2_scale) > 0 || error("Scales must be positive")

    output_dir = joinpath(abspath(options["output-root"]), lowercase(target_name))
    metadata_path = joinpath(output_dir, "map_metadata.toml")
    isfile(metadata_path) && lowercase(options["overwrite"]) != "true" && begin
        println("Existing map retained: $metadata_path")
        return 0
    end

    ring = cesr
    sextupoles = active_sextupole_inventory(ring)
    controls = independent_corrector_inventory(ring)
    matches = [entry for entry in sextupoles if entry.name == target_name]
    length(matches) == 1 || error("Unknown or duplicate target: $target_name")
    target = only(matches)
    target.index < length(ring.line) || error("Target cannot be the last ring element")

    knob_rows = read_bump_knobs(abspath(options["bump-knobs-csv"]), target_name)
    knob_by_control = Dict(
        (row["corrector"], row["field"]) => (
            parse(Float64, row["field_per_x_bump_m"]),
            parse(Float64, row["field_per_y_bump_m"]),
        )
        for row in knob_rows
    )
    control_keys = Set((control.name, String(control.axis)) for control in controls)
    Set(keys(knob_by_control)) == control_keys || error("Bump-knob control mismatch")

    closed_timed = @timed solve_closed_orbit(ring)
    descriptor = Descriptor(6, maximum_order, 3, maximum_order)
    qx, qy, qk = params(descriptor)
    for control in controls
        cx, cy = knob_by_control[(control.name, String(control.axis))]
        for (element_index, original) in zip(control.indices, control.originals)
            value = constant_term(original) + bump_scale * (cx*qx + cy*qy)
            if control.axis == :Kn0
                ring.line[element_index].Kn0 = value
            else
                ring.line[element_index].Ks0 = value
            end
        end
    end
    ring.line[target.index].Kn2 = target.kn2_m3 + k2_scale*qk

    calculation_start = time()
    transport = calculate_direct_transport(
        ring,
        [ring.line[target.index], ring.line[target.index + 1], ring.line[end]],
        closed_timed.value;
        descriptor,
    )
    fixed_point, residual_history = parameter_dependent_fixed_point(
        transport.maps[end],
        maximum_order,
    )
    entry_map = normalized_parameter_map(transport.maps[1] ∘ fixed_point)
    exit_map = normalized_parameter_map(transport.maps[2] ∘ fixed_point)
    maps = (entry=entry_map, exit=exit_map)
    powers = monomial_powers(maximum_order)
    coordinate_names = ("x", "px", "y", "py", "z", "pz")
    coefficient_rows = NamedTuple[]
    for (location, local_map) in pairs(maps)
        for (coordinate_index, coordinate) in enumerate(coordinate_names)
            value = local_map.v[coordinate_index]
            for power in powers
                push!(coefficient_rows, (;
                    target=target_name,
                    location=String(location),
                    coordinate,
                    qx_power=power[1],
                    qy_power=power[2],
                    qk_power=power[3],
                    total_order=sum(power),
                    coefficient=taylor_coefficient(value, power),
                ))
            end
        end
    end
    all(isfinite(row.coefficient) for row in coefficient_rows) ||
        error("Non-finite local-orbit Taylor coefficient")

    mkpath(output_dir)
    write_rows(joinpath(output_dir, "local_orbit_taylor_coefficients.csv"), coefficient_rows)
    write_rows(joinpath(output_dir, "fixed_point_residual_history.csv"), [
        (; iteration=index - 1, residual_norm=value)
        for (index, value) in enumerate(residual_history)
    ])
    write_metadata(metadata_path, Dict(
        "format" => "cesr-sextupole-local-orbit-bump-k2-gtpsa-v1",
        "date" => string(Dates.today()),
        "engine" => "SciBmad/GTPSA direct transport with order-by-order parameter-dependent fixed point",
        "lattice" => LATEST_LATTICE,
        "target" => target_name,
        "target_s_m" => target.s_m,
        "target_length_m" => target.length_m,
        "target_nominal_kn2_m3" => target.kn2_m3,
        "descriptor" => "Descriptor(6, $maximum_order, 3, $maximum_order)",
        "maximum_total_parameter_order" => maximum_order,
        "parameters" => ["normalized_x_bump", "normalized_y_bump", "normalized_delta_Kn2"],
        "bump_scale_m" => bump_scale,
        "k2_scale_m3" => k2_scale,
        "coordinate_locations" => ["target entrance", "target exit"],
        "coefficient_semantics" => "Taylor coefficients in normalized qx,qy,qk; coefficient already includes inverse factorial",
        "nominal_closed_orbit_seconds" => closed_timed.time,
        "map_seconds" => time() - calculation_start,
        "initial_fixed_point_residual" => first(residual_history),
        "final_fixed_point_residual" => last(residual_history),
        "interpretation_boundary" => "model-validity map only; no hardware, aperture, or operator limit is inferred",
    ))
    @printf(
        "%s order-%d local map complete in %.3f s; residual %.3e\n",
        target_name,
        maximum_order,
        time() - calculation_start,
        last(residual_history),
    )
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end

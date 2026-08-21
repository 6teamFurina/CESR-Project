#!/usr/bin/env julia

"""
Compute per-sextupole GTPSA response coefficients for CESR beam-based
sextupole alignment.

Each selected active normal sextupole is assigned three GTPSA parameters:

    p1 = delta Kn2       [m^-3]
    p2 = delta x_offset  [m]
    p3 = delta y_offset  [m]

The descriptor uses third total/phase-space order and second parameter order:

    Descriptor(6, 3, 3, 2)

Third total order is required because linear-optics feed-down contains a
phase-space variable multiplied by both Kn2 and offset parameters. The saved
second derivatives are true derivatives (GTPSA Hessian entries), not Taylor
coefficients; diagonal quadratic Taylor coefficients are one half of the
corresponding saved second derivative.
"""

using Beamlines
using Dates
using GTPSA
using Printf
using SciBmad
using Statistics
using TOML

const HERE = @__DIR__
const PROJECT_DIR = normpath(joinpath(HERE, "..", "..", "..", ".."))
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

function parse_options(args)
    options = Dict{String,String}(
        "start" => "1",
        "stop" => "76",
        "output-dir" => joinpath(HERE, "results", "full"),
        "part-label" => "all",
        "warmup" => "true",
    )
    for argument in args
        startswith(argument, "--") ||
            error("Arguments must have --name=value form: $argument")
        fields = split(argument[3:end], '='; limit=2)
        length(fields) == 2 || error("Missing value in argument: $argument")
        haskey(options, fields[1]) || error("Unknown option: --$(fields[1])")
        options[fields[1]] = fields[2]
    end
    lowercase(options["warmup"]) in ("true", "false") ||
        error("--warmup must be true or false")
    return options
end

function active_sextupole_inventory(ring)
    entries = NamedTuple[]
    for (index, element) in enumerate(ring.line)
        string(element.kind) == "Sextupole" || continue
        kn2 = constant_term(Beamlines.deval(element.Kn2))
        iszero(kn2) && continue
        push!(entries, (;
            index,
            name=base_name(element),
            runtime_name=String(element.name),
            kn2_m3=kn2,
            length_m=constant_term(Beamlines.deval(element.L)),
            x_offset_m=constant_term(Beamlines.deval(element.x_offset)),
            y_offset_m=constant_term(Beamlines.deval(element.y_offset)),
        ))
    end
    sort!(entries; by=entry -> entry.name)
    names = getproperty.(entries, :name)
    length(unique(names)) == length(names) ||
        error("Active sextupole base names are not unique")
    length(entries) == 76 ||
        error("Expected 76 active normal sextupoles, found $(length(entries))")
    return entries
end

function detector_elements(ring)
    detectors = [
        element for element in ring.line
        if startswith(uppercase(String(element.name)), "DET_")
    ]
    length(detectors) == 99 ||
        error("Expected 99 detector markers, found $(length(detectors))")
    return detectors
end

function nominal_closed_orbit()
    model = load_cesr_model(zero_value=0.0, rf_on=true)
    solution = find_closed_orbit(
        model.ring;
        coasting_beam=false,
        batch=Val{false}(),
        warn=true,
    )
    all(solution.sol.retcode .== SciBmad.BatchSolve.RETCODE_SUCCESS) ||
        error("Nominal RF-on closed-orbit solve failed: $(solution.sol.retcode)")
    return (; v0=solution.v0, coasting_beam=solution.coasting_beam)
end

function prepare_parameterized_model(nominal_inventory)
    descriptor = Descriptor(6, 3, 3, 2)
    parameters = params(descriptor)
    typed_zero = zero(parameters[1])
    model = load_cesr_model(zero_value=typed_zero, rf_on=true)
    inventory = active_sextupole_inventory(model.ring)
    getproperty.(inventory, :name) == getproperty.(nominal_inventory, :name) ||
        error("Typed and nominal sextupole inventories differ")
    return (;
        descriptor,
        parameters,
        typed_zero,
        model,
        inventory,
        detectors=detector_elements(model.ring),
    )
end

function set_baseline!(setup, inventory_index)
    entry = setup.inventory[inventory_index]
    element = setup.model.ring.line[entry.index]
    element.Kn2 = entry.kn2_m3 + setup.typed_zero
    element.x_offset = entry.x_offset_m + setup.typed_zero
    element.y_offset = entry.y_offset_m + setup.typed_zero
    return element
end

function activate!(setup, inventory_index, previous_index=nothing)
    isnothing(previous_index) || set_baseline!(setup, previous_index)
    entry = setup.inventory[inventory_index]
    element = setup.model.ring.line[entry.index]
    element.Kn2 = entry.kn2_m3 + setup.parameters[1]
    element.x_offset = entry.x_offset_m + setup.parameters[2]
    element.y_offset = entry.y_offset_m + setup.parameters[3]
    return element
end

function parameterized_twiss(setup, closed)
    return twiss(
        setup.model.ring;
        GTPSA_descriptor=setup.descriptor,
        at=setup.detectors,
        v0=closed.v0,
        v0_and_coast=(closed.v0, closed.coasting_beam),
        spin=false,
        RDTs=false,
        normalizing_map=false,
    )
end

function coefficient_tuple(value)
    jacobian = GTPSA.jacobian([value]; include_params=true)
    hessian = GTPSA.hessian(value; include_params=true)
    size(jacobian, 2) == 9 ||
        error("Expected six variables plus three parameters in Jacobian")
    size(hessian) == (9, 9) || error("Unexpected Hessian size $(size(hessian))")
    all(isfinite, jacobian) || error("Non-finite Jacobian coefficient")
    all(isfinite, hessian) || error("Non-finite Hessian coefficient")
    return (;
        value=constant_term(value),
        d_k2=Float64(jacobian[1, 7]),
        d_x=Float64(jacobian[1, 8]),
        d_y=Float64(jacobian[1, 9]),
        d2_k2_k2=Float64(hessian[7, 7]),
        d2_k2_x=Float64(hessian[7, 8]),
        d2_k2_y=Float64(hessian[7, 9]),
        d2_x_x=Float64(hessian[8, 8]),
        d2_x_y=Float64(hessian[8, 9]),
        d2_y_y=Float64(hessian[9, 9]),
    )
end

function coefficient_rows(sextupole, optics)
    rows = NamedTuple[]
    length(optics.table.name) == 99 ||
        error("Expected 99 detector Twiss rows, found $(length(optics.table.name))")
    for detector in eachindex(optics.table.name)
        detector_name = lowercase(String(optics.table.name[detector]))
        for column in DETECTOR_COLUMNS
            value = getproperty(optics.table, column)[detector]
            # Accumulated eigenphases have a parameter-dependent additive
            # gauge in parameterized Twiss. Real difference-phase data use a
            # fixed BPM reference. Remove the first-detector phase for every
            # TPS order before extracting coefficients.
            if column in (:phi_1, :phi_2)
                value -= getproperty(optics.table, column)[1]
            end
            coefficients = coefficient_tuple(value)
            push!(rows, (;
                sextupole=sextupole.name,
                sextupole_index=sextupole.index,
                sextupole_s_m=sextupole.s_m,
                nominal_kn2_m3=sextupole.kn2_m3,
                nominal_k2l_m2=sextupole.k2l_m2,
                observation_scope="detector",
                observation_name=detector_name,
                observable=String(column),
                coefficients...,
            ))
        end
    end
    for tune_index in eachindex(optics.tunes)
        coefficients = coefficient_tuple(optics.tunes[tune_index])
        push!(rows, (;
            sextupole=sextupole.name,
            sextupole_index=sextupole.index,
            sextupole_s_m=sextupole.s_m,
            nominal_kn2_m3=sextupole.kn2_m3,
            nominal_k2l_m2=sextupole.k2l_m2,
            observation_scope="ring",
            observation_name="ring",
            observable="tune_$(tune_index)",
            coefficients...,
        ))
    end
    return rows
end

function with_s_positions(ring, inventory)
    s_m = 0.0
    center_by_index = Dict{Int,Float64}()
    for (index, element) in enumerate(ring.line)
        length_m = constant_term(Beamlines.deval(element.L))
        center_by_index[index] = s_m + length_m / 2
        s_m += length_m
    end
    return [
        (;
            entry...,
            s_m=center_by_index[entry.index],
            k2l_m2=entry.kn2_m3 * entry.length_m,
        )
        for entry in inventory
    ]
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

function distribution(values)
    return Dict(
        "minimum" => minimum(values),
        "median" => median(values),
        "maximum" => maximum(values),
        "mean" => mean(values),
    )
end

function main(args=ARGS)
    options = parse_options(args)
    output_dir = abspath(options["output-dir"])
    part_label = options["part-label"]
    first_index = parse(Int, options["start"])
    last_index = parse(Int, options["stop"])
    run_warmup = lowercase(options["warmup"]) == "true"

    baseline_model = load_cesr_model(zero_value=0.0, rf_on=true)
    nominal_inventory = with_s_positions(
        baseline_model.ring,
        active_sextupole_inventory(baseline_model.ring),
    )
    1 <= first_index <= last_index <= length(nominal_inventory) ||
        error("Requested inventory range $first_index:$last_index is invalid")

    closed_timed = @timed nominal_closed_orbit()
    closed = closed_timed.value
    setup_timed = @timed prepare_parameterized_model(nominal_inventory)
    setup = setup_timed.value

    warmup_seconds = 0.0
    if run_warmup
        activate!(setup, first_index)
        warmup_seconds = @elapsed parameterized_twiss(setup, closed)
        set_baseline!(setup, first_index)
        @printf("Part %s warmup: %.3f s\n", part_label, warmup_seconds)
        flush(stdout)
        GC.gc()
    end

    rows = NamedTuple[]
    timing_rows = NamedTuple[]
    previous_index = nothing
    for inventory_index in first_index:last_index
        entry = nominal_inventory[inventory_index]
        activation_seconds = @elapsed activate!(setup, inventory_index, previous_index)
        twiss_timed = @timed parameterized_twiss(setup, closed)
        extraction_timed = @timed coefficient_rows(entry, twiss_timed.value)
        append!(rows, extraction_timed.value)
        push!(timing_rows, (;
            part_label,
            inventory_index,
            sextupole=entry.name,
            activation_seconds,
            twiss_seconds=twiss_timed.time,
            extraction_seconds=extraction_timed.time,
            allocated_bytes=twiss_timed.bytes,
            coefficient_rows=length(extraction_timed.value),
        ))
        previous_index = inventory_index
        @printf(
            "Part %s sextupole %d/%d %-7s: Twiss %.3f s, extraction %.3f s\n",
            part_label,
            inventory_index,
            length(nominal_inventory),
            entry.name,
            twiss_timed.time,
            extraction_timed.time,
        )
        flush(stdout)
    end
    isnothing(previous_index) || set_baseline!(setup, previous_index)

    coefficients_path = write_rows(
        joinpath(output_dir, "alignment_coefficients_part_$(part_label).csv"),
        rows,
    )
    timings_path = write_rows(
        joinpath(output_dir, "alignment_timings_part_$(part_label).csv"),
        timing_rows,
    )
    twiss_seconds = getproperty.(timing_rows, :twiss_seconds)
    metadata = Dict(
        "format" => "cesr-sextupole-alignment-gtpsa-part-v1",
        "date" => string(Dates.today()),
        "engine" => "SciBmad/GTPSA",
        "rf_mode" => "on (six-dimensional periodic optics)",
        "descriptor" => "Descriptor(6, 3, 3, 2)",
        "phase_space_total_order" => 3,
        "parameter_count" => 3,
        "parameter_order" => 2,
        "parameters" => ["delta_Kn2_m3", "delta_x_offset_m", "delta_y_offset_m"],
        "parameter_scaling" => "unit physical increments: 1 m^-3, 1 m, 1 m",
        "hessian_semantics" => "saved d2 columns are true second derivatives; diagonal Taylor coefficient is d2/2",
        "part_label" => part_label,
        "inventory_start" => first_index,
        "inventory_stop" => last_index,
        "sextupole_count" => length(timing_rows),
        "detector_count" => 99,
        "detector_observables" => collect(String.(DETECTOR_COLUMNS)),
        "ring_tune_count" => 3,
        "coefficient_rows" => length(rows),
        "nominal_closed_orbit_seconds" => closed_timed.time,
        "model_setup_seconds" => setup_timed.time,
        "warmup_enabled" => run_warmup,
        "warmup_seconds" => warmup_seconds,
        "twiss_seconds" => sum(twiss_seconds),
        "twiss_seconds_distribution" => distribution(twiss_seconds),
        "coefficients_csv" => coefficients_path,
        "timings_csv" => timings_path,
        "julia_version" => string(VERSION),
        "julia_threads" => Threads.nthreads(),
    )
    metadata_path = joinpath(output_dir, "alignment_metadata_part_$(part_label).toml")
    open(metadata_path, "w") do io
        TOML.print(io, metadata; sorted=true)
    end
    println("Coefficients: $coefficients_path")
    println("Metadata: $metadata_path")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end

#!/usr/bin/env julia

"""Validate one sextupole's mixed GTPSA coefficients with four-corner differences."""

include(joinpath(@__DIR__, "run_sextupole_alignment_gtpsa.jl"))
using LinearAlgebra

function scalar_observables(entry, delta_k2, delta_x, delta_y)
    model = load_cesr_model(zero_value=0.0, rf_on=true)
    inventory = with_s_positions(model.ring, active_sextupole_inventory(model.ring))
    lookup = Dict(item.name => item for item in inventory)
    target = lookup[entry.name]
    element = model.ring.line[target.index]
    element.Kn2 = target.kn2_m3 + delta_k2
    element.x_offset = target.x_offset_m + delta_x
    element.y_offset = target.y_offset_m + delta_y
    closed = find_closed_orbit(
        model.ring;
        coasting_beam=false,
        batch=Val{false}(),
        warn=true,
    )
    all(closed.sol.retcode .== SciBmad.BatchSolve.RETCODE_SUCCESS) ||
        error("Validation closed-orbit solve failed")
    detectors = detector_elements(model.ring)
    optics = twiss(
        model.ring;
        GTPSA_descriptor=Descriptor(6, 3),
        at=detectors,
        v0=closed.v0,
        v0_and_coast=(closed.v0, closed.coasting_beam),
        spin=false,
        RDTs=false,
        normalizing_map=false,
    )
    values = Dict{Tuple{String,String,String},Float64}()
    for detector in eachindex(optics.table.name)
        name = lowercase(String(optics.table.name[detector]))
        for column in DETECTOR_COLUMNS
            value = getproperty(optics.table, column)[detector]
            if column in (:phi_1, :phi_2)
                value -= getproperty(optics.table, column)[1]
            end
            values[("detector", name, String(column))] = constant_term(value)
        end
    end
    for tune_index in eachindex(optics.tunes)
        values[("ring", "ring", "tune_$(tune_index)")] =
            constant_term(optics.tunes[tune_index])
    end
    return values
end

function mixed_four_corner(entry, k2_step, offset_step, plane)
    dx(sign) = plane == :x ? sign * offset_step : 0.0
    dy(sign) = plane == :y ? sign * offset_step : 0.0
    pp = scalar_observables(entry, k2_step, dx(1), dy(1))
    pm = scalar_observables(entry, k2_step, dx(-1), dy(-1))
    mp = scalar_observables(entry, -k2_step, dx(1), dy(1))
    mm = scalar_observables(entry, -k2_step, dx(-1), dy(-1))
    return Dict(
        key => (pp[key] - pm[key] - mp[key] + mm[key]) /
               (4 * k2_step * offset_step)
        for key in keys(pp)
    )
end

function main_validation(args=ARGS)
    options = Dict{String,String}(
        "sextupole" => "SEX_08W",
        "k2-step" => "0.01",
        "offset-step" => "0.0001",
        "output-dir" => joinpath(HERE, "results", "validation"),
    )
    for argument in args
        startswith(argument, "--") || error("Arguments must have --name=value form")
        fields = split(argument[3:end], '='; limit=2)
        haskey(options, fields[1]) || error("Unknown option: --$(fields[1])")
        options[fields[1]] = fields[2]
    end
    target_name = uppercase(options["sextupole"])
    k2_step = parse(Float64, options["k2-step"])
    offset_step = parse(Float64, options["offset-step"])
    output_dir = abspath(options["output-dir"])

    scalar_model = load_cesr_model(zero_value=0.0, rf_on=true)
    inventory = with_s_positions(
        scalar_model.ring,
        active_sextupole_inventory(scalar_model.ring),
    )
    entry = only(filter(item -> item.name == target_name, inventory))
    closed = nominal_closed_orbit()
    setup = prepare_parameterized_model(inventory)
    inventory_index = findfirst(item -> item.name == target_name, inventory)
    activate!(setup, inventory_index)
    optics_timed = @timed parameterized_twiss(setup, closed)
    gtpsa_rows = coefficient_rows(entry, optics_timed.value)
    gtpsa = Dict(
        (row.observation_scope, row.observation_name, row.observable) => row
        for row in gtpsa_rows
    )

    x_timed = @timed mixed_four_corner(entry, k2_step, offset_step, :x)
    y_timed = @timed mixed_four_corner(entry, k2_step, offset_step, :y)
    x_fd = x_timed.value
    y_fd = y_timed.value
    rows = NamedTuple[]
    for key in sort!(collect(keys(gtpsa)))
        row = gtpsa[key]
        push!(rows, (;
            sextupole=target_name,
            observation_scope=key[1],
            observation_name=key[2],
            observable=key[3],
            gtpsa_d2_k2_x=row.d2_k2_x,
            finite_difference_d2_k2_x=x_fd[key],
            residual_k2_x=x_fd[key] - row.d2_k2_x,
            gtpsa_d2_k2_y=row.d2_k2_y,
            finite_difference_d2_k2_y=y_fd[key],
            residual_k2_y=y_fd[key] - row.d2_k2_y,
        ))
    end
    validation_path = write_rows(
        joinpath(output_dir, "mixed_coefficient_validation.csv"), rows,
    )

    summaries = NamedTuple[]
    groups = sort!(unique(row.observable for row in rows))
    for observable in groups
        selected = filter(row -> row.observable == observable, rows)
        for plane in (:x, :y)
            gt = plane == :x ? getproperty.(selected, :gtpsa_d2_k2_x) :
                               getproperty.(selected, :gtpsa_d2_k2_y)
            fd = plane == :x ? getproperty.(selected, :finite_difference_d2_k2_x) :
                               getproperty.(selected, :finite_difference_d2_k2_y)
            residual = fd .- gt
            denominator = norm(gt)
            push!(summaries, (;
                observable,
                offset_plane=String(plane),
                count=length(selected),
                gtpsa_norm=norm(gt),
                finite_difference_norm=norm(fd),
                residual_norm=norm(residual),
                relative_l2=denominator == 0 ? NaN : norm(residual) / denominator,
                maximum_absolute_residual=maximum(abs, residual),
            ))
        end
    end
    summary_path = write_rows(
        joinpath(output_dir, "mixed_coefficient_validation_summary.csv"),
        summaries,
    )
    metadata = Dict(
        "format" => "cesr-sextupole-alignment-gtpsa-validation-v1",
        "date" => string(Dates.today()),
        "sextupole" => target_name,
        "k2_step_m3" => k2_step,
        "offset_step_m" => offset_step,
        "gtpsa_descriptor" => "Descriptor(6, 3, 3, 2)",
        "scalar_twiss_descriptor" => "Descriptor(6, 3)",
        "gtpsa_twiss_seconds" => optics_timed.time,
        "x_four_corner_seconds" => x_timed.time,
        "y_four_corner_seconds" => y_timed.time,
        "validation_csv" => validation_path,
        "summary_csv" => summary_path,
    )
    mkpath(output_dir)
    metadata_path = joinpath(output_dir, "validation_metadata.toml")
    open(metadata_path, "w") do io
        TOML.print(io, metadata; sorted=true)
    end
    println("Validation: $validation_path")
    println("Summary: $summary_path")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main_validation())
end

#!/usr/bin/env julia

"""Generate paired P1 and P2 inverse dictionaries for one saved exact scan.

P1 is an oracle-background-conditioned mixed response: the other sextupole
offsets are copied from the saved realization, while the target offset error is
set to zero. P2 additionally differentiates all saved observables with respect
to four integrated local multipoles on the target element.
"""

include(joinpath(@__DIR__, "run_exact_scan.jl"))

function read_simple_csv(path)
    lines = readlines(path)
    header = split(strip(first(lines)), ',')
    return [Dict(header .=> split(strip(line), ',')) for line in lines[2:end] if !isempty(strip(line))]
end

function set_background!(model, inventory, truth_rows, target_name, target_dx, target_dy; mode="saved")
    mode in ("saved", "nominal") || error("background mode must be saved or nominal")
    truth = Dict(uppercase(row["sextupole"]) => row for row in truth_rows)
    for entry in inventory
        row = truth[entry.name]
        dx = entry.name == target_name ? target_dx :
            (mode == "saved" ? parse(Float64, row["x_offset_error_m"]) : 0.0)
        dy = entry.name == target_name ? target_dy :
            (mode == "saved" ? parse(Float64, row["y_offset_error_m"]) : 0.0)
        element = model.ring.line[entry.index]
        element.x_offset = entry.x_offset_m + dx
        element.y_offset = entry.y_offset_m + dy
    end
end

function set_bump!(model, knobs, bump_x, bump_y)
    for knob in knobs
        model.controls[knob.corrector] =
            knob.kick_x_rad_per_m * bump_x + knob.kick_y_rad_per_m * bump_y
    end
end

function observation_vector(model, target)
    calculated = calculate_scalar_observables(model, target)
    keys = [(row.observation_scope, row.observation_name, row.observable) for row in calculated.rows]
    values = [row.value for row in calculated.rows]
    return keys, values
end

function k2_slope!(model, target, delta_k2)
    element = model.ring.line[target.index]
    element.Kn2 = target.kn2_m3 + delta_k2
    keys_plus, plus = observation_vector(model, target)
    element.Kn2 = target.kn2_m3 - delta_k2
    keys_minus, minus = observation_vector(model, target)
    element.Kn2 = target.kn2_m3
    keys_plus == keys_minus || error("Observation keys changed across K2 perturbation")
    return keys_plus, (plus .- minus) ./ (2delta_k2)
end

function offset_conditioned_response!(model, target, offset_step, k2_step)
    element = model.ring.line[target.index]
    x_base = target.x_offset_m
    y_base = target.y_offset_m
    element.x_offset = x_base
    element.y_offset = y_base
    keys, reference = k2_slope!(model, target, k2_step)

    element.x_offset = x_base + offset_step
    _, x_plus = k2_slope!(model, target, k2_step)
    element.x_offset = x_base - offset_step
    _, x_minus = k2_slope!(model, target, k2_step)
    element.x_offset = x_base

    element.y_offset = y_base + offset_step
    _, y_plus = k2_slope!(model, target, k2_step)
    element.y_offset = y_base - offset_step
    _, y_minus = k2_slope!(model, target, k2_step)
    element.y_offset = y_base
    return keys, reference, (x_plus .- x_minus) ./ (2offset_step),
        (y_plus .- y_minus) ./ (2offset_step)
end

function local_source_response!(model, target, source_steps)
    element = model.ring.line[target.index]
    columns = Vector{Vector{Float64}}()
    keys_reference = nothing
    for (property, step) in source_steps
        baseline = constant_term(Beamlines.deval(getproperty(element, property)))
        setproperty!(element, property, baseline + step)
        keys_plus, plus = observation_vector(model, target)
        setproperty!(element, property, baseline - step)
        keys_minus, minus = observation_vector(model, target)
        setproperty!(element, property, baseline)
        keys_plus == keys_minus || error("Observation keys changed for $property")
        isnothing(keys_reference) ? (keys_reference = keys_plus) :
            (keys_reference == keys_plus || error("Local-source keys differ"))
        push!(columns, (plus .- minus) ./ (2step))
    end
    return keys_reference, hcat(columns...)
end

function main(args=ARGS)
    defaults = Dict(
        "scan-dir" => joinpath(@__DIR__, "results", "smoke_background"),
        "output-dir" => joinpath(@__DIR__, "results", "paired_benchmark", "conditioned_models"),
        "knobs-csv" => joinpath(@__DIR__, "results", "bump_knobs", "local_bump_knobs.csv"),
        "offset-step-m" => "2.5e-5",
        "dipole-step" => "1e-7",
        "quadrupole-step" => "1e-5",
        "background-mode" => "saved",
        "nonlinear-calibration" => "false",
        "nonlinear-calibration-amplitude-m" => "5e-4",
    )
    options = parse_key_value_args(defaults, args)
    scan_dir = abspath(options["scan-dir"])
    output_dir = abspath(options["output-dir"])
    states = read_simple_csv(joinpath(scan_dir, "scan_states.csv"))
    truth_rows = read_simple_csv(joinpath(scan_dir, "sextupole_offset_truth.csv"))
    target_name = uppercase(first(states)["target_sextupole"])
    k2_step = abs(parse(Float64, first(row for row in states if parse(Float64, row["delta_k2_m3"]) != 0)["delta_k2_m3"]))
    offset_step = parse(Float64, options["offset-step-m"])
    background_mode = lowercase(options["background-mode"])
    nonlinear_calibration = lowercase(options["nonlinear-calibration"]) == "true"
    calibration_amplitude = parse(Float64, options["nonlinear-calibration-amplitude-m"])
    lowercase(options["nonlinear-calibration"]) in ("true", "false") ||
        error("--nonlinear-calibration must be true or false")

    model = load_cesr_model(zero_value=0.0, rf_on=true)
    inventory = active_sextupole_inventory(model.ring)
    target = find_inventory_entry(inventory, target_name)
    set_background!(model, inventory, truth_rows, target_name, 0.0, 0.0; mode=background_mode)
    knobs = read_target_knobs(abspath(options["knobs-csv"]), target_name)
    zero_states = sort(
        [row for row in states if abs(parse(Float64, row["delta_k2_m3"])) < 1e-18];
        by=row -> parse(Int, row["bump_index"]),
    )
    source_steps = [
        (:Kn0L, parse(Float64, options["dipole-step"])),
        (:Ks0L, parse(Float64, options["dipole-step"])),
        (:Kn1L, parse(Float64, options["quadrupole-step"])),
        (:Ks1L, parse(Float64, options["quadrupole-step"])),
    ]

    # Pay the Julia/SciBmad compilation cost once so the per-method timings
    # below represent warm exact-forward work rather than method ordering.
    set_bump!(model, knobs, 0.0, 0.0)
    warmup_seconds = @elapsed observation_vector(model, target)

    conditioned_rows = NamedTuple[]
    source_rows = NamedTuple[]
    nonlinear_rows = NamedTuple[]
    p1_seconds = 0.0
    p2_source_seconds = 0.0
    nonlinear_calibration_seconds = 0.0
    for (counter, row) in enumerate(zero_states)
        bump_index = parse(Int, row["bump_index"])
        bump_x = parse(Float64, row["bump_x_command_m"])
        bump_y = parse(Float64, row["bump_y_command_m"])
        set_bump!(model, knobs, bump_x, bump_y)
        local keys, reference, response_x, response_y
        p1_seconds += @elapsed begin
            keys, reference, response_x, response_y =
                offset_conditioned_response!(model, target, offset_step, k2_step)
        end
        local source_keys, source
        p2_source_seconds += @elapsed begin
            source_keys, source = local_source_response!(model, target, source_steps)
        end
        keys == source_keys || error("P1 and P2 observation keys differ")
        for index in eachindex(keys)
            scope, name, observable = keys[index]
            push!(conditioned_rows, (;
                target_sextupole=target.name, bump_index,
                observation_scope=scope, observation_name=name, observable,
                reference_k2_slope=reference[index],
                d2_k2_x=response_x[index], d2_k2_y=response_y[index],
            ))
            push!(source_rows, (;
                target_sextupole=target.name, bump_index,
                observation_scope=scope, observation_name=name, observable,
                d_kn0l=source[index, 1], d_ks0l=source[index, 2],
                d_kn1l=source[index, 3], d_ks1l=source[index, 4],
            ))
        end
        if nonlinear_calibration
            element = model.ring.line[target.index]
            calibration_points = bump_points("grid3", calibration_amplitude)
            nonlinear_calibration_seconds += @elapsed begin
                for (calibration_index, (offset_x, offset_y)) in enumerate(calibration_points)
                    element.x_offset = target.x_offset_m + offset_x
                    element.y_offset = target.y_offset_m + offset_y
                    calibration_keys, calibration_slope = k2_slope!(model, target, k2_step)
                    calibration_keys == keys || error("Nonlinear calibration keys differ")
                    for index in eachindex(keys)
                        scope, name, observable = keys[index]
                        push!(nonlinear_rows, (;
                            target_sextupole=target.name, bump_index, calibration_index,
                            calibration_x_offset_m=offset_x,
                            calibration_y_offset_m=offset_y,
                            observation_scope=scope, observation_name=name, observable,
                            k2_slope=calibration_slope[index],
                        ))
                    end
                end
                element.x_offset = target.x_offset_m
                element.y_offset = target.y_offset_m
            end
        end
        @printf("Conditioned model bump %d/%d complete\n", counter, length(zero_states))
        flush(stdout)
    end
    conditioned_path = write_rows(joinpath(output_dir, "conditioned_mixed_response.csv"), conditioned_rows)
    source_path = write_rows(joinpath(output_dir, "local_source_response.csv"), source_rows)
    nonlinear_path = nonlinear_calibration ?
        write_rows(joinpath(output_dir, "nonlinear_offset_calibration_slopes.csv"), nonlinear_rows) : ""
    metadata = Dict(
        "format" => "cesr-targeted-conditioned-inverse-models-v1",
        "date" => string(Dates.today()),
        "scan_dir" => scan_dir,
        "target_sextupole" => target.name,
        "background_mode" => background_mode,
        "background_semantics" => background_mode == "saved" ?
            "other 75 saved offsets known; target offset error fixed to zero" :
            "all other sextupole offset errors fixed to nominal zero; target offset error fixed to zero",
        "p1_derivative" => "exact central finite difference of exact K2 central-difference slopes",
        "p2_sources" => ["Kn0L", "Ks0L", "Kn1L", "Ks1L"],
        "offset_step_m" => offset_step,
        "k2_step_m3" => k2_step,
        "warmup_seconds" => warmup_seconds,
        "p1_conditioned_response_seconds" => p1_seconds,
        "p1_exact_forward_state_count" => 10 * length(zero_states),
        "p2_source_incremental_seconds" => p2_source_seconds,
        "p2_source_exact_forward_state_count" => 8 * length(zero_states),
        "p2_total_model_seconds_including_p1" => p1_seconds + p2_source_seconds,
        "nonlinear_calibration_enabled" => nonlinear_calibration,
        "nonlinear_calibration_amplitude_m" => calibration_amplitude,
        "nonlinear_calibration_seconds" => nonlinear_calibration_seconds,
        "nonlinear_calibration_state_count" => nonlinear_calibration ?
            2 * 9 * length(zero_states) : 0,
        "nonlinear_calibration_slopes_csv" => nonlinear_path,
        "conditioned_response_csv" => conditioned_path,
        "local_source_response_csv" => source_path,
    )
    mkpath(output_dir)
    open(joinpath(output_dir, "metadata.toml"), "w") do io
        TOML.print(io, metadata; sorted=true)
    end
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end

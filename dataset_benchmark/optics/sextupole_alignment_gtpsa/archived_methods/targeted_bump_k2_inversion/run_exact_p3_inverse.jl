#!/usr/bin/env julia

"""Small exact-SciBmad P3 inverse on a saved paired scan.

The other 75 sextupole offsets are fixed to the saved realization. The target
offset is updated by finite-difference Gauss--Newton against the complete
nonzero-K2 difference scan. This is an oracle-background diagnostic, not yet an
operational estimator for unknown nuisance offsets.
"""

include(joinpath(@__DIR__, "generate_conditioned_inverse_models.jl"))

function observable_family_name(name)
    startswith(name, "orbit_") && return "orbit"
    startswith(name, "phi_") && return "phase"
    startswith(name, "c") && return "coupling"
    startswith(name, "tune_") && return "tune"
    return "excluded"
end

function read_scales(path)
    rows = read_simple_csv(path)
    return Dict(row["observable"] => parse(Float64, row["global_mixed_response_rms"]) for row in rows)
end

function initial_estimate(path)
    rows = read_simple_csv(path)
    row = first(row for row in rows if row["method"] == "P2b_four_local_kicks")
    return [parse(Float64, row["estimated_x_offset_m"]), parse(Float64, row["estimated_y_offset_m"])]
end

function observed_difference_vector(states, observation_rows, scales)
    grouped = Dict{String,Vector{Dict{String,String}}}()
    for row in observation_rows
        push!(get!(grouped, row["state_id"], Dict{String,String}[]), row)
    end
    zero_id = Dict(
        parse(Int, row["bump_index"]) => row["state_id"]
        for row in states if abs(parse(Float64, row["delta_k2_m3"])) < 1e-18
    )
    selected = [
        index for (index, row) in enumerate(grouped[first(states)["state_id"]])
        if observable_family_name(row["observable"]) != "excluded"
    ]
    vector = Float64[]
    scale_vector = Float64[]
    for state in states
        abs(parse(Float64, state["delta_k2_m3"])) < 1e-18 && continue
        current = grouped[state["state_id"]]
        baseline = grouped[zero_id[parse(Int, state["bump_index"] )]]
        for index in selected
            push!(vector,
                parse(Float64, current[index]["observable_readback"]) -
                parse(Float64, baseline[index]["observable_readback"]))
            push!(scale_vector, scales[current[index]["observable"]])
        end
    end
    return vector ./ scale_vector, selected
end

function exact_difference_vector!(model, target, states, knobs, selected, scales, candidate)
    element = model.ring.line[target.index]
    element.x_offset = target.x_offset_m + candidate[1]
    element.y_offset = target.y_offset_m + candidate[2]
    values = Dict{String,Vector{Float64}}()
    zero_id = Dict{Int,String}()
    for state in states
        bump_index = parse(Int, state["bump_index"])
        set_bump!(
            model, knobs,
            parse(Float64, state["bump_x_command_m"]),
            parse(Float64, state["bump_y_command_m"]),
        )
        delta_k2 = parse(Float64, state["delta_k2_m3"])
        element.Kn2 = target.kn2_m3 + delta_k2
        keys, observable_values = observation_vector(model, target)
        values[state["state_id"]] = observable_values
        abs(delta_k2) < 1e-18 && (zero_id[bump_index] = state["state_id"])
    end
    element.Kn2 = target.kn2_m3
    vector = Float64[]
    for state in states
        delta_k2 = parse(Float64, state["delta_k2_m3"])
        abs(delta_k2) < 1e-18 && continue
        current = values[state["state_id"]]
        baseline = values[zero_id[parse(Int, state["bump_index"] )]]
        for index in selected
            observable = calculate_observable_name(index)
            push!(vector, (current[index] - baseline[index]) / scales[observable])
        end
    end
    return vector
end

# calculate_scalar_observables has a fixed order: 12 detector fields repeated
# for 99 detectors, followed by tune_1 and tune_2.
function calculate_observable_name(index)
    detector_count = 99 * length(DETECTOR_COLUMNS)
    if index <= detector_count
        return String(DETECTOR_COLUMNS[mod1(index, length(DETECTOR_COLUMNS))])
    end
    return "tune_$(index - detector_count)"
end

function main(args=ARGS)
    defaults = Dict(
        "scan-dir" => joinpath(@__DIR__, "results", "smoke_background"),
        "benchmark-dir" => joinpath(@__DIR__, "results", "paired_benchmark"),
        "knobs-csv" => joinpath(@__DIR__, "results", "bump_knobs", "local_bump_knobs.csv"),
        "iterations" => "3",
        "jacobian-step-m" => "1e-5",
        "max-update-m" => "2e-4",
    )
    options = parse_key_value_args(defaults, args)
    scan_dir = abspath(options["scan-dir"])
    benchmark_dir = abspath(options["benchmark-dir"])
    states = read_simple_csv(joinpath(scan_dir, "scan_states.csv"))
    observations = read_simple_csv(joinpath(scan_dir, "scan_observations.csv"))
    truth_rows = read_simple_csv(joinpath(scan_dir, "sextupole_offset_truth.csv"))
    target_name = uppercase(first(states)["target_sextupole"])
    scales = read_scales(joinpath(RESPONSE_MAP_DIR, "results", "full", "local_response_svd_scales.csv"))
    observed, selected = observed_difference_vector(states, observations, scales)

    model = load_cesr_model(zero_value=0.0, rf_on=true)
    inventory = active_sextupole_inventory(model.ring)
    target = find_inventory_entry(inventory, target_name)
    set_background!(model, inventory, truth_rows, target_name, 0.0, 0.0)
    knobs = read_target_knobs(abspath(options["knobs-csv"]), target_name)
    estimate = initial_estimate(joinpath(benchmark_dir, "p0_p2_offset_estimates.csv"))
    truth_row = first(row for row in truth_rows if uppercase(row["sextupole"]) == target_name)
    truth = [parse(Float64, truth_row["x_offset_error_m"]), parse(Float64, truth_row["y_offset_error_m"])]
    step = parse(Float64, options["jacobian-step-m"])
    max_update = parse(Float64, options["max-update-m"])
    iterations = parse(Int, options["iterations"])
    history = NamedTuple[]

    set_bump!(model, knobs, 0.0, 0.0)
    model.ring.line[target.index].Kn2 = target.kn2_m3
    warmup_seconds = @elapsed observation_vector(model, target)
    inverse_seconds = @elapsed begin
        for iteration in 0:iterations
            predicted = exact_difference_vector!(model, target, states, knobs, selected, scales, estimate)
            residual = observed - predicted
            rms = norm(residual) / sqrt(length(residual))
            push!(history, (;
                iteration,
                estimated_x_offset_m=estimate[1], estimated_y_offset_m=estimate[2],
                error_x_m=estimate[1] - truth[1], error_y_m=estimate[2] - truth[2],
                absolute_error_2d_m=norm(estimate - truth), weighted_residual_rms=rms,
            ))
            @printf("P3 iteration %d: x=%+.6f mm y=%+.6f mm error=%.6f um residual=%.6e\n",
                iteration, 1e3estimate[1], 1e3estimate[2], 1e6norm(estimate-truth), rms)
            iteration == iterations && break
            x_trial = copy(estimate); x_trial[1] += step
            y_trial = copy(estimate); y_trial[2] += step
            predicted_x = exact_difference_vector!(model, target, states, knobs, selected, scales, x_trial)
            predicted_y = exact_difference_vector!(model, target, states, knobs, selected, scales, y_trial)
            jacobian = hcat((predicted_x - predicted) ./ step, (predicted_y - predicted) ./ step)
            update = jacobian \ residual
            update = clamp.(update, -max_update, max_update)
            estimate .+= update
        end
    end

    history_path = write_rows(joinpath(benchmark_dir, "p3_exact_history.csv"), history)
    final = last(history)
    metadata = Dict(
        "format" => "cesr-exact-p3-inverse-v1",
        "date" => string(Dates.today()),
        "target" => target.name,
        "scan_dir" => scan_dir,
        "observable_view" => "orbit_phase_coupling_tune full finite scan differences",
        "background_semantics" => "other 75 saved offsets fixed to truth (oracle diagnostic)",
        "forward_engine" => "exact scalar SciBmad RF-on closed orbit and Twiss at every candidate",
        "iterations" => iterations,
        "jacobian_step_m" => step,
        "warmup_seconds" => warmup_seconds,
        "inverse_seconds" => inverse_seconds,
        "exact_forward_scan_evaluations" => 3 * iterations + 1,
        "exact_forward_state_count" => length(states) * (3 * iterations + 1),
        "estimated_x_offset_m" => final.estimated_x_offset_m,
        "estimated_y_offset_m" => final.estimated_y_offset_m,
        "absolute_error_2d_m" => final.absolute_error_2d_m,
        "history_csv" => history_path,
    )
    open(joinpath(benchmark_dir, "p3_exact_summary.toml"), "w") do io
        TOML.print(io, metadata; sorted=true)
    end
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end

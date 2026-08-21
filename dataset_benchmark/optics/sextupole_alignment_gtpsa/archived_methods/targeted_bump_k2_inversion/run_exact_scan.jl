#!/usr/bin/env julia

"""Generate one exact target-sextupole bump-by-K2 scan with known offset truth."""

include(joinpath(@__DIR__, "common.jl"))

function read_target_knobs(path, target)
    lines = readlines(path)
    header = split(strip(first(lines)), ',')
    lookup = Dict(name => index for (index, name) in enumerate(header))
    required = ("target_sextupole", "corrector", "kick_x_rad_per_m", "kick_y_rad_per_m")
    all(haskey(lookup, name) for name in required) || error("Unexpected knob CSV header")
    rows = NamedTuple[]
    for line in lines[2:end]
        isempty(strip(line)) && continue
        fields = split(strip(line), ',')
        uppercase(fields[lookup["target_sextupole"]]) == uppercase(target) || continue
        push!(rows, (;
            corrector=fields[lookup["corrector"]],
            kick_x_rad_per_m=parse(Float64, fields[lookup["kick_x_rad_per_m"]]),
            kick_y_rad_per_m=parse(Float64, fields[lookup["kick_y_rad_per_m"]]),
        ))
    end
    length(rows) == 119 || error("Expected 119 knob rows for $target, found $(length(rows))")
    return rows
end

function bump_points(protocol, amplitude)
    if protocol == "cross5"
        return [
            (0.0, 0.0),
            (-amplitude, 0.0),
            (amplitude, 0.0),
            (0.0, -amplitude),
            (0.0, amplitude),
        ]
    elseif protocol == "grid3"
        points = [(0.0, 0.0)]
        append!(points, [
            (x, y) for x in (-amplitude, 0.0, amplitude)
            for y in (-amplitude, 0.0, amplitude)
            if !(iszero(x) && iszero(y))
        ])
        return points
    end
    error("--bump-protocol must be cross5 or grid3")
end

function parse_levels(text)
    levels = parse.(Float64, strip.(split(text, ',')))
    0.0 in levels || error("K2 levels must contain zero")
    return levels
end

function main(args=ARGS)
    defaults = Dict(
        "target" => "SEX_08W",
        "true-x-offset-m" => "3.5e-4",
        "true-y-offset-m" => "-2.5e-4",
        "background-offset-rms-m" => "0.0",
        "seed" => "20260813",
        "bump-protocol" => "cross5",
        "bump-amplitude-m" => "5e-4",
        "k2-step-m3" => "0.01",
        "k2-levels" => "-1,0,1",
        "knobs-csv" => joinpath(@__DIR__, "results", "bump_knobs", "local_bump_knobs.csv"),
        "output-dir" => joinpath(@__DIR__, "results", "smoke_exact"),
    )
    options = parse_key_value_args(defaults, args)
    target_name = uppercase(options["target"])
    truth_x = parse(Float64, options["true-x-offset-m"])
    truth_y = parse(Float64, options["true-y-offset-m"])
    background_rms = parse(Float64, options["background-offset-rms-m"])
    seed = parse(Int, options["seed"])
    bump_amplitude = parse(Float64, options["bump-amplitude-m"])
    k2_step = parse(Float64, options["k2-step-m3"])
    k2_levels = parse_levels(options["k2-levels"])
    output_dir = abspath(options["output-dir"])
    knobs_path = abspath(options["knobs-csv"])
    knobs = read_target_knobs(knobs_path, target_name)
    points = bump_points(options["bump-protocol"], bump_amplitude)

    model = load_cesr_model(zero_value=0.0, rf_on=true)
    inventory = active_sextupole_inventory(model.ring)
    target = find_inventory_entry(inventory, target_name)
    rng = MersenneTwister(seed)
    background_rows = NamedTuple[]
    for entry in inventory
        dx = entry.name == target.name ? truth_x : background_rms * randn(rng)
        dy = entry.name == target.name ? truth_y : background_rms * randn(rng)
        element = model.ring.line[entry.index]
        element.x_offset = entry.x_offset_m + dx
        element.y_offset = entry.y_offset_m + dy
        push!(background_rows, (;
            sextupole=entry.name,
            x_offset_error_m=dx,
            y_offset_error_m=dy,
            is_target=entry.name == target.name,
        ))
    end

    observation_rows = NamedTuple[]
    state_rows = NamedTuple[]
    total_states = length(points) * length(k2_levels)
    state_counter = 0
    for (bump_index, (bump_x, bump_y)) in enumerate(points)
        for level in k2_levels
            state_counter += 1
            delta_k2 = level * k2_step
            for knob in knobs
                model.controls[knob.corrector] =
                    knob.kick_x_rad_per_m * bump_x + knob.kick_y_rad_per_m * bump_y
            end
            model.ring.line[target.index].Kn2 = target.kn2_m3 + delta_k2
            timed = @timed calculate_scalar_observables(model, target)
            calculated = timed.value
            state_id = @sprintf("%s_b%02d_k%+g", lowercase(target.name), bump_index, level)
            push!(state_rows, (;
                state_id,
                target_sextupole=target.name,
                bump_index,
                bump_x_command_m=bump_x,
                bump_y_command_m=bump_y,
                k2_level=level,
                delta_k2_m3=delta_k2,
                true_x_offset_m=truth_x,
                true_y_offset_m=truth_y,
                background_offset_rms_m=background_rms,
                target_orbit_x_m=calculated.target_orbit_x_m,
                target_orbit_y_m=calculated.target_orbit_y_m,
                calculation_seconds=timed.time,
            ))
            for row in calculated.rows
                push!(observation_rows, (;
                    state_id,
                    target_sextupole=target.name,
                    bump_index,
                    bump_x_command_m=bump_x,
                    bump_y_command_m=bump_y,
                    k2_level=level,
                    delta_k2_m3=delta_k2,
                    row.observation_scope,
                    row.observation_name,
                    row.observable,
                    observable_readback=row.value,
                ))
            end
            @printf(
                "State %d/%d bump=(%+.3f,%+.3f) mm dK2=%+.4g: %.3f s\n",
                state_counter,
                total_states,
                1e3 * bump_x,
                1e3 * bump_y,
                delta_k2,
                timed.time,
            )
            flush(stdout)
        end
    end

    observations_path = write_rows(joinpath(output_dir, "scan_observations.csv"), observation_rows)
    states_path = write_rows(joinpath(output_dir, "scan_states.csv"), state_rows)
    truth_path = write_rows(joinpath(output_dir, "sextupole_offset_truth.csv"), background_rows)
    metadata = Dict(
        "format" => "cesr-targeted-sextupole-exact-scan-v1",
        "date" => string(Dates.today()),
        "engine" => "SciBmad exact scalar RF-on closed orbit and Twiss",
        "target_sextupole" => target.name,
        "true_x_offset_m" => truth_x,
        "true_y_offset_m" => truth_y,
        "background_offset_rms_m" => background_rms,
        "seed" => seed,
        "bump_protocol" => options["bump-protocol"],
        "bump_amplitude_m" => bump_amplitude,
        "k2_step_m3" => k2_step,
        "k2_levels" => k2_levels,
        "state_count" => length(state_rows),
        "observation_rows" => length(observation_rows),
        "knobs_csv" => knobs_path,
        "observations_csv" => observations_path,
        "states_csv" => states_path,
        "truth_csv" => truth_path,
        "phase_reference" => "DET_00W",
    )
    open(joinpath(output_dir, "scan_metadata.toml"), "w") do io
        TOML.print(io, metadata; sorted=true)
    end
    println("Observations: $observations_path")
    println("States: $states_path")
    println("Truth: $truth_path")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end


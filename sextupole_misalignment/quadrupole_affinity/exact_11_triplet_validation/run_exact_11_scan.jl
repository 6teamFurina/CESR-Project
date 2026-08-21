#!/usr/bin/env julia

"""Run an exact repaired-lattice 11-K1 by 9-bump by 5-K2 orbit scan."""

include(joinpath(@__DIR__, "common.jl"))

function parse_levels(text)
    values = parse.(Float64, strip.(split(text, ',')))
    values == [-2.0, -1.0, 0.0, 1.0, 2.0] ||
        error("This validation requires K2 levels -2,-1,0,1,2")
    return values
end

function grid3_points(amplitude)
    return [(x, y) for x in (-amplitude, 0.0, amplitude) for y in (-amplitude, 0.0, amplitude)]
end

function solve_scalar_warm(ring, initial_v0)
    solution = find_closed_orbit(
        ring;
        v0=copy(initial_v0),
        coasting_beam=false,
        batch=Val{false}(),
        warn=false,
    )
    all(solution.sol.retcode .== SciBmad.BatchSolve.RETCODE_SUCCESS) ||
        error("RF-on closed-orbit solve failed: $(solution.sol.retcode)")
    return solution
end

function run_scan(args=ARGS)
    defaults = Dict(
        "target" => "SEX_09AW",
        "true-x-offset-m" => "3.5e-4",
        "true-y-offset-m" => "-2.5e-4",
        "offset-fd-step-m" => "5e-5",
        "bump-amplitude-m" => "5e-4",
        "k2-step-m3" => "0.01",
        "k2-levels" => "-2,-1,0,1,2",
        "quadrupole-fraction" => "0.001",
        "scenario-mode" => "truth-only",
        "candidate-csv" => joinpath(
            AFFINITY_HERE,
            "results",
            "scibmad_latest",
            "selection",
            "quadrupole_sets_by_sextupole.csv",
        ),
        "bump-knobs-csv" => joinpath(@__DIR__, "results", "bump_knobs", "local_bump_knobs.csv"),
        "output-dir" => joinpath(@__DIR__, "results", "scans"),
        "overwrite" => "false",
    )
    options = parse_exact11_options(defaults, args)
    target_name = uppercase(options["target"])
    truth_x = parse(Float64, options["true-x-offset-m"])
    truth_y = parse(Float64, options["true-y-offset-m"])
    fd_step = parse(Float64, options["offset-fd-step-m"])
    bump_amplitude = parse(Float64, options["bump-amplitude-m"])
    k2_step = parse(Float64, options["k2-step-m3"])
    k2_levels = parse_levels(options["k2-levels"])
    quadrupole_fraction = parse(Float64, options["quadrupole-fraction"])
    overwrite = lowercase(options["overwrite"]) == "true"
    output_dir = joinpath(abspath(options["output-dir"]), safe_name(target_name))
    metadata_path = joinpath(output_dir, "scan_metadata.toml")
    isfile(metadata_path) && !overwrite && error("Output exists; use --overwrite=true: $metadata_path")

    ring = cesr
    sextupoles = active_sextupole_inventory(ring)
    quadrupoles = active_quadrupole_inventory(ring)
    detectors = measurable_bpms(ring)
    target = findfirst(entry -> entry.name == target_name, sextupoles)
    isnothing(target) && error("Unknown target sextupole: $target_name")
    target_entry = sextupoles[target]
    quadrupole_by_name = Dict(entry.name => entry for entry in quadrupoles)
    candidate_sets = read_five_candidate_sets(abspath(options["candidate-csv"]))
    candidates = candidate_sets[target_name]
    length(candidates) == 5 || error("Expected five candidates for $target_name")
    candidate_entries = [quadrupole_by_name[name] for name in candidates]

    controls = independent_corrector_inventory(ring)
    control_by_key = Dict((control.name, String(control.axis)) => control for control in controls)
    knob_rows = read_bump_knobs(abspath(options["bump-knobs-csv"]), target_name)
    length(knob_rows) == length(controls) || error("Bump knob control count mismatch")
    knob_by_key = Dict(
        (row["corrector"], row["field"]) => (
            parse(Float64, row["field_per_x_bump_m"]),
            parse(Float64, row["field_per_y_bump_m"]),
        )
        for row in knob_rows
    )
    Set(keys(knob_by_key)) == Set(keys(control_by_key)) || error("Bump knob controls do not match lattice inventory")

    condition_rows = NamedTuple[(;
        condition_index=1,
        condition_name="nominal",
        quadrupole="nominal",
        sign=0,
        delta_k1_m2=0.0,
    )]
    for entry in candidate_entries
        step = quadrupole_fraction * abs(entry.kn1_m2)
        for sign in (1, -1)
            push!(condition_rows, (;
                condition_index=length(condition_rows) + 1,
                condition_name="$(lowercase(entry.name))_$(sign > 0 ? "plus" : "minus")",
                quadrupole=entry.name,
                sign,
                delta_k1_m2=sign * step,
            ))
        end
    end
    length(condition_rows) == 11 || error("Expected 11 K1 conditions")
    bump_points = grid3_points(bump_amplitude)
    scenario_mode = lowercase(options["scenario-mode"])
    scenario_mode in ("truth-only", "full") || error("--scenario-mode must be truth-only or full")
    scenario_rows = scenario_mode == "truth-only" ?
        [(name="truth", x=truth_x, y=truth_y)] :
        [
            (name="truth", x=truth_x, y=truth_y),
            (name="zero", x=0.0, y=0.0),
            (name="x_plus", x=fd_step, y=0.0),
            (name="x_minus", x=-fd_step, y=0.0),
            (name="y_plus", x=0.0, y=fd_step),
            (name="y_minus", x=0.0, y=-fd_step),
        ]
    state_count = length(condition_rows) * length(bump_points) * length(k2_levels)
    bpm_names = String.(base_name.(detectors))
    observations = zeros(
        length(scenario_rows),
        length(condition_rows),
        length(bump_points),
        length(k2_levels),
        length(detectors),
        2,
    )
    target_orbits = zeros(
        length(scenario_rows),
        length(condition_rows),
        length(bump_points),
        length(k2_levels),
        2,
    )
    timing_rows = NamedTuple[]

    initial_timed = @timed solve_closed_orbit(ring)
    nominal_v0 = copy(initial_timed.value.v0)
    target_element = ring.line[target_entry.index]
    for (scenario_index, scenario) in enumerate(scenario_rows)
        target_element.x_offset = target_entry.x_offset_m + scenario.x
        target_element.y_offset = target_entry.y_offset_m + scenario.y
        current_v0 = copy(nominal_v0)
        scenario_start = time()
        state_counter = 0
        for (condition_index, condition) in enumerate(condition_rows)
            for entry in candidate_entries
                set_quadrupole!(ring, entry, entry.kn1_m2)
            end
            if condition.sign != 0
                entry = quadrupole_by_name[condition.quadrupole]
                set_quadrupole!(ring, entry, entry.kn1_m2 + condition.delta_k1_m2)
            end
            for (bump_index, (bump_x, bump_y)) in enumerate(bump_points)
                for control in controls
                    key = (control.name, String(control.axis))
                    coefficient_x, coefficient_y = knob_by_key[key]
                    baseline = constant_term(first(control.originals))
                    set_corrector_scalar!(
                        ring,
                        control,
                        baseline + coefficient_x * bump_x + coefficient_y * bump_y,
                    )
                end
                for (k2_index, level) in enumerate(k2_levels)
                    state_counter += 1
                    target_element.Kn2 = target_entry.kn2_m3 + level * k2_step
                    solve_timed = @timed solve_scalar_warm(ring, current_v0)
                    closed = solve_timed.value
                    current_v0 .= closed.v0
                    track_timed = @timed track_orbits_at_names(
                        ring,
                        closed,
                        vcat(bpm_names, [target_name]),
                    )
                    tracked = track_timed.value
                    for (bpm_index, bpm_name) in enumerate(bpm_names)
                        observations[scenario_index, condition_index, bump_index, k2_index, bpm_index, 1] =
                            tracked.horizontal[bpm_name][1]
                        observations[scenario_index, condition_index, bump_index, k2_index, bpm_index, 2] =
                            tracked.vertical[bpm_name][1]
                    end
                    target_orbits[scenario_index, condition_index, bump_index, k2_index, 1] =
                        tracked.horizontal[target_name][1]
                    target_orbits[scenario_index, condition_index, bump_index, k2_index, 2] =
                        tracked.vertical[target_name][1]
                    push!(timing_rows, (;
                        scenario=scenario.name,
                        condition_index,
                        bump_index,
                        k2_index,
                        closed_orbit_seconds=solve_timed.time,
                        tracking_seconds=track_timed.time,
                    ))
                end
            end
            @printf(
                "%s %s: %d/%d states, %.1f s\n",
                target_name,
                scenario.name,
                state_counter,
                state_count,
                time() - scenario_start,
            )
            flush(stdout)
        end
    end

    mkpath(output_dir)
    write_npy(joinpath(output_dir, "bpm_orbits.npy"), observations)
    write_npy(joinpath(output_dir, "target_orbits.npy"), target_orbits)
    write_rows(joinpath(output_dir, "k1_conditions.csv"), condition_rows)
    write_rows(joinpath(output_dir, "timings.csv"), timing_rows)
    write_lines(joinpath(output_dir, "scenario_labels.txt"), getproperty.(scenario_rows, :name))
    write_lines(joinpath(output_dir, "bpm_names.txt"), bpm_names)
    write_rows(
        joinpath(output_dir, "bump_points.csv"),
        [(; bump_index=index, bump_x_command_m=point[1], bump_y_command_m=point[2]) for (index, point) in enumerate(bump_points)],
    )
    metadata = Dict(
        "format" => "cesr-repaired-lattice-exact-11-k1-scan-v1",
        "date" => string(Dates.today()),
        "engine" => "SciBmad exact scalar RF-on closed orbit and tracking",
        "lattice" => LATEST_LATTICE,
        "target_sextupole" => target_name,
        "candidate_quadrupoles" => candidates,
        "scenario_count" => length(scenario_rows),
        "k1_condition_count" => length(condition_rows),
        "bump_count" => length(bump_points),
        "k2_count" => length(k2_levels),
        "state_count_per_scenario" => state_count,
        "true_x_offset_m" => truth_x,
        "true_y_offset_m" => truth_y,
        "baseline_x_offset_m" => target_entry.x_offset_m,
        "baseline_y_offset_m" => target_entry.y_offset_m,
        "offset_fd_step_m" => fd_step,
        "scenario_mode" => scenario_mode,
        "bump_amplitude_m" => bump_amplitude,
        "k2_step_m3" => k2_step,
        "k2_levels" => k2_levels,
        "quadrupole_fraction" => quadrupole_fraction,
        "bpm_count" => length(detectors),
        "corrector_count" => length(controls),
        "initial_closed_orbit_seconds" => initial_timed.time,
        "observable_scope" => "BPM x/y closed orbit; direct phase/coupling ablation deferred",
        "nuisance_scope" => "nominal lattice and target offset only; measurement-noise covariance is added in analysis",
    )
    write_metadata(metadata_path, metadata)
    target_element.Kn2 = target_entry.original_kn2
    target_element.x_offset = target_entry.original_x_offset
    target_element.y_offset = target_entry.original_y_offset
    for entry in candidate_entries
        restore_quadrupole!(ring, entry)
    end
    for control in controls
        for (index, original) in zip(control.indices, control.originals)
            if control.axis == :Kn0
                ring.line[index].Kn0 = original
            else
                ring.line[index].Ks0 = original
            end
        end
    end
    println("Output: $output_dir")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(run_scan())
end

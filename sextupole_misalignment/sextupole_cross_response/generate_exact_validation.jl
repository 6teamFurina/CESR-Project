#!/usr/bin/env julia

"""Generate paired exact SciBmad scans for selected cross-response targets.

Each target is evaluated in an aligned and a deliberately misaligned scenario
on the same nine-state K2/bump stencil.  Subtracting the aligned gradient from
the misaligned gradient isolates the center-dependent response and removes the
finite-amplitude aligned background before comparison with the nominal GTPSA
alignment design.
"""

include(joinpath(
    @__DIR__, "..", "quadrupole_affinity", "exact_11_triplet_validation", "common.jl",
))

const DEFAULT_BUMP_KNOBS = joinpath(
    @__DIR__, "..", "quadrupole_affinity", "exact_11_triplet_validation",
    "results", "bump_knobs", "local_bump_knobs.csv",
)

function scan_states(bump_amplitude, k2_step)
    rows = NamedTuple[(;
        state_index=1,
        state="zero",
        bump_axis="zero",
        bump_sign=0,
        k2_sign=0,
        bump_x_m=0.0,
        bump_y_m=0.0,
        delta_k2_m3=0.0,
    )]
    for (axis, bx, by) in (("x", 1.0, 0.0), ("y", 0.0, 1.0))
        for bump_sign in (-1, 1), k2_sign in (-1, 1)
            push!(rows, (;
                state_index=length(rows) + 1,
                state="$(axis)_b$(bump_sign)_k$(k2_sign)",
                bump_axis=axis,
                bump_sign,
                k2_sign,
                bump_x_m=bump_sign * bx * bump_amplitude,
                bump_y_m=bump_sign * by * bump_amplitude,
                delta_k2_m3=k2_sign * k2_step,
            ))
        end
    end
    return rows
end

function main(args=ARGS)
    defaults = Dict(
        "targets" => "SEX_09AW,SEX_14W,SEX_18W,SEX_39W,SEX_44E",
        "bump-amplitude-m" => "5.0e-4",
        "k2-step-m3" => "2.0e-2",
        "offset-x-m" => "3.5e-4",
        "offset-y-m" => "-2.5e-4",
        "bump-knobs-csv" => DEFAULT_BUMP_KNOBS,
        "output-dir" => joinpath(@__DIR__, "results", "exact_validation"),
        "overwrite" => "false",
    )
    options = parse_exact11_options(defaults, args)
    output_dir = abspath(options["output-dir"])
    metadata_path = joinpath(output_dir, "scan_metadata.toml")
    if isfile(metadata_path) && lowercase(options["overwrite"]) != "true"
        println("Existing exact validation retained: $metadata_path")
        return 0
    end

    bump_amplitude = parse(Float64, options["bump-amplitude-m"])
    k2_step = parse(Float64, options["k2-step-m3"])
    offset_x = parse(Float64, options["offset-x-m"])
    offset_y = parse(Float64, options["offset-y-m"])
    all(isfinite.((bump_amplitude, k2_step, offset_x, offset_y))) ||
        error("Validation settings must be finite")
    min(bump_amplitude, k2_step) > 0 || error("Bump and K2 steps must be positive")

    ring = cesr
    sextupoles = active_sextupole_inventory(ring)
    controls = independent_corrector_inventory(ring)
    all_names = String.(getproperty.(sextupoles, :name))
    selected_names = [uppercase(strip(name)) for name in split(options["targets"], ',')]
    unknown = setdiff(selected_names, all_names)
    isempty(unknown) || error("Unknown targets: $(join(unknown, ", "))")
    selected = [target for target in sextupoles if target.name in selected_names]
    length(selected) == length(selected_names) || error("Duplicate requested targets")

    knob_rows = read_simple_csv(abspath(options["bump-knobs-csv"]))
    knob_by_target = Dict{String,Dict{Tuple{String,String},Tuple{Float64,Float64}}}()
    for target in selected
        rows = [row for row in knob_rows if uppercase(row["target_sextupole"]) == target.name]
        length(rows) == length(controls) || error("Incomplete bump knob for $(target.name)")
        knob_by_target[target.name] = Dict(
            (row["corrector"], row["field"]) => (
                parse(Float64, row["field_per_x_bump_m"]),
                parse(Float64, row["field_per_y_bump_m"]),
            )
            for row in rows
        )
    end

    states = scan_states(bump_amplitude, k2_step)
    scenario_offsets = ((0.0, 0.0), (offset_x, offset_y))
    scenario_names = ("aligned", "misaligned")
    state_count = length(states)
    batch_count = length(scenario_names) * state_count
    orbits = zeros(length(selected), length(scenario_names), state_count, length(sextupoles), 2)
    centers = zeros(length(selected), length(scenario_names), 2)
    timing_rows = NamedTuple[]

    # The exported lattice retains some scalar parameter blocks as DefExpr
    # values.  Materialize only the correctors and selected sextupoles that
    # will be promoted to BatchParam; otherwise mixed DefExpr{Any}/BatchParam
    # blocks are ambiguous to Beamlines.  This deliberately avoids touching
    # unrelated lattice elements.
    for element_index in Set(index for control in controls for index in control.indices)
        element = ring.line[element_index]
        element.BMultipoleParams = Beamlines.deval(element.BMultipoleParams)
    end
    for target in selected
        element = ring.line[target.index]
        element.BMultipoleParams = Beamlines.deval(element.BMultipoleParams)
        element.AlignmentParams = Beamlines.deval(element.AlignmentParams)
    end

    for (target_counter, target) in enumerate(selected)
        batch_bx = zeros(batch_count)
        batch_by = zeros(batch_count)
        batch_k2 = zeros(batch_count)
        batch_x_offset = zeros(batch_count)
        batch_y_offset = zeros(batch_count)
        for scenario in eachindex(scenario_names), state in eachindex(states)
            index = (scenario - 1) * state_count + state
            batch_bx[index] = states[state].bump_x_m
            batch_by[index] = states[state].bump_y_m
            batch_k2[index] = target.kn2_m3 + states[state].delta_k2_m3
            batch_x_offset[index] = target.x_offset_m + scenario_offsets[scenario][1]
            batch_y_offset[index] = target.y_offset_m + scenario_offsets[scenario][2]
        end

        knobs = knob_by_target[target.name]
        for control in controls
            cx, cy = knobs[(control.name, String(control.axis))]
            baseline = constant_term(first(control.originals))
            set_corrector_values!(ring, control, baseline .+ cx .* batch_bx .+ cy .* batch_by)
        end
        element = ring.line[target.index]
        if iszero(target.kn1_m2) && iszero(target.ks1_m2)
            element.tracking_method = BeamTracking.DriftKick()
        end
        element.Kn2 = BatchParam(batch_k2)
        element.x_offset = BatchParam(batch_x_offset)
        element.y_offset = BatchParam(batch_y_offset)

        solve_timed = @timed solve_batch_closed_orbit(ring, batch_count)
        track_timed = @timed track_orbits_at_names(ring, solve_timed.value, all_names)
        tracked = track_timed.value
        for scenario in eachindex(scenario_names), state in eachindex(states)
            batch_index = (scenario - 1) * state_count + state
            for (observation, name) in enumerate(all_names)
                orbits[target_counter, scenario, state, observation, 1] =
                    tracked.horizontal[name][batch_index]
                orbits[target_counter, scenario, state, observation, 2] =
                    tracked.vertical[name][batch_index]
            end
        end
        for scenario in eachindex(scenario_names)
            centers[target_counter, scenario, 1] = target.x_offset_m + scenario_offsets[scenario][1]
            centers[target_counter, scenario, 2] = target.y_offset_m + scenario_offsets[scenario][2]
        end

        for control in controls
            set_corrector_scalar!(ring, control, constant_term(first(control.originals)))
        end
        restore_sextupole!(ring, target)
        push!(timing_rows, (;
            target_index=target_counter,
            target=target.name,
            batch_states=batch_count,
            closed_orbit_seconds=solve_timed.time,
            tracking_seconds=track_timed.time,
        ))
        @printf(
            "Exact validation %d/%d %-10s closed %.3f s track %.3f s\n",
            target_counter, length(selected), target.name, solve_timed.time, track_timed.time,
        )
        flush(stdout)
    end

    target_rows = [
        (;
            selected_index=index,
            inventory_index=findfirst(==(target.name), all_names),
            target=target.name,
            s_m=target.s_m,
            length_m=target.length_m,
        )
        for (index, target) in enumerate(selected)
    ]
    scenario_rows = [
        (;
            scenario_index=index,
            scenario=scenario_names[index],
            added_x_offset_m=scenario_offsets[index][1],
            added_y_offset_m=scenario_offsets[index][2],
        )
        for index in eachindex(scenario_names)
    ]
    mkpath(output_dir)
    write_npy(joinpath(output_dir, "exact_sextupole_orbits.npy"), orbits)
    write_npy(joinpath(output_dir, "scenario_centers.npy"), centers)
    write_rows(joinpath(output_dir, "states.csv"), states)
    write_rows(joinpath(output_dir, "scenarios.csv"), scenario_rows)
    write_rows(joinpath(output_dir, "selected_targets.csv"), target_rows)
    write_rows(joinpath(output_dir, "timings.csv"), timing_rows)
    write_metadata(metadata_path, Dict(
        "format" => "cesr-sextupole-cross-response-exact-validation-v1",
        "date" => string(Dates.today()),
        "engine" => "exact scalar-field BatchParam SciBmad RF-on closed orbit and tracking",
        "lattice" => LATEST_LATTICE,
        "all_sextupole_count" => length(sextupoles),
        "selected_target_count" => length(selected),
        "selected_targets" => String.(getproperty.(selected, :name)),
        "scenario_order" => collect(scenario_names),
        "state_count_per_scenario" => state_count,
        "bump_amplitude_m" => bump_amplitude,
        "k2_step_m3" => k2_step,
        "misaligned_added_x_offset_m" => offset_x,
        "misaligned_added_y_offset_m" => offset_y,
        "observation_location" => "first element-entry occurrence of every active normal sextupole",
        "paired_subtraction" => "misaligned K2-odd/bump-odd gradient minus aligned gradient",
        "bump_knobs_csv" => abspath(options["bump-knobs-csv"]),
        "julia_version" => string(VERSION),
    ))
    println("Exact cross-response validation written to $output_dir")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end

#!/usr/bin/env julia

"""Generate exact latest-lattice scans with sextupole misalignment as the only error.

Every target/realization contains fixed x/y offsets on all 76 active normal
sextupoles.  Quadrupoles, correctors, BPMs, RF, and all other model parameters
remain at their validated latest-lattice values.  The dense 5 x 5 local-bump
grid and five symmetric target-K2 levels support derivative and high-order
Taylor-surface inversions without changing the underlying physical states.
"""

include(joinpath(
    @__DIR__, "..", "quadrupole_affinity", "exact_11_triplet_validation", "common.jl",
))

const DEFAULT_BUMP_KNOBS = joinpath(
    @__DIR__, "..", "quadrupole_affinity", "exact_11_triplet_validation",
    "results", "bump_knobs", "local_bump_knobs.csv",
)

function scalar_warm(ring, v0)
    solution = find_closed_orbit(
        ring; v0=copy(v0), coasting_beam=false, batch=Val{false}(), warn=false,
    )
    all(solution.sol.retcode .== SciBmad.BatchSolve.RETCODE_SUCCESS) ||
        error("RF-on closed-orbit solve failed: $(solution.sol.retcode)")
    return solution
end

function restore_correctors!(ring, controls)
    for control in controls
        for (index, original) in zip(control.indices, control.originals)
            if control.axis == :Kn0
                ring.line[index].Kn0 = original
            else
                ring.line[index].Ks0 = original
            end
        end
    end
    return nothing
end

function selected_target_entries(requested, sextupoles)
    text = lowercase(strip(requested))
    text == "all" && return collect(enumerate(sextupoles))
    wanted = Set(uppercase(strip(name)) for name in split(requested, ',') if !isempty(strip(name)))
    known = Set(getproperty.(sextupoles, :name))
    missing = setdiff(wanted, known)
    isempty(missing) || error("Unknown target sextupoles: $(join(sort!(collect(missing)), ", "))")
    return [(index, entry) for (index, entry) in enumerate(sextupoles) if entry.name in wanted]
end

function main(args=ARGS)
    defaults = Dict(
        "targets" => "all",
        "realizations" => "1",
        "seed" => "20260819",
        "target-halfwidth-m" => "3.5e-4",
        "other-sext-rms-m" => "3.0e-4",
        "bump-amplitude-m" => "5.0e-4",
        "k2-step-m3" => "0.01",
        "bump-knobs-csv" => DEFAULT_BUMP_KNOBS,
        "output-dir" => joinpath(@__DIR__, "results", "exact_scans"),
        "overwrite" => "false",
    )
    options = parse_exact11_options(defaults, args)
    realization_count = parse(Int, options["realizations"])
    seed = parse(Int, options["seed"])
    target_halfwidth = parse(Float64, options["target-halfwidth-m"])
    other_sext_rms = parse(Float64, options["other-sext-rms-m"])
    bump_amplitude = parse(Float64, options["bump-amplitude-m"])
    k2_step = parse(Float64, options["k2-step-m3"])
    output_dir = abspath(options["output-dir"])
    metadata_path = joinpath(output_dir, "scan_metadata.toml")
    isfile(metadata_path) && lowercase(options["overwrite"]) != "true" &&
        error("Output exists; use --overwrite=true: $metadata_path")
    realization_count > 0 || error("--realizations must be positive")

    ring = cesr
    sextupoles = active_sextupole_inventory(ring)
    detectors = measurable_bpms(ring)
    controls = independent_corrector_inventory(ring)
    length(sextupoles) == 76 || error("Expected 76 active sextupoles")
    selected = selected_target_entries(options["targets"], sextupoles)
    isempty(selected) && error("No target sextupoles selected")
    target_names = [entry.name for (_, entry) in selected]
    target_inventory_indices = [index for (index, _) in selected]
    bpm_names = String.(base_name.(detectors))

    knob_rows = read_simple_csv(abspath(options["bump-knobs-csv"]))
    knobs_by_target = Dict{String,Dict{Tuple{String,String},Tuple{Float64,Float64}}}()
    for row in knob_rows
        target = uppercase(row["target_sextupole"])
        target_knobs = get!(
            knobs_by_target,
            target,
            Dict{Tuple{String,String},Tuple{Float64,Float64}}(),
        )
        target_knobs[(row["corrector"], row["field"])] = (
            parse(Float64, row["field_per_x_bump_m"]),
            parse(Float64, row["field_per_y_bump_m"]),
        )
    end
    control_keys = Set((control.name, String(control.axis)) for control in controls)
    all(haskey(knobs_by_target, name) for name in target_names) ||
        error("At least one selected target has no bump knob")
    all(Set(keys(knobs_by_target[name])) == control_keys for name in target_names) ||
        error("At least one selected target bump knob has a control mismatch")

    grid = bump_amplitude .* (-1.0, -0.5, 0.0, 0.5, 1.0)
    bumps = [(x, y) for x in grid for y in grid]
    k2_levels = [-2.0, -1.0, 0.0, 1.0, 2.0]
    nt, nr = length(selected), realization_count
    nb, nk, nd, ns = length(bumps), length(k2_levels), length(detectors), length(sextupoles)
    bpm_orbits = zeros(nt, nr, nb, nk, nd, 2)
    target_orbits = zeros(nt, nr, nb, nk, 2)
    target_truth = zeros(nt, nr, 2)
    latent_offsets = zeros(nt, nr, ns, 2)
    realization_seconds = zeros(nt, nr)

    nominal_closed = solve_closed_orbit(ring)
    nominal_v0 = copy(nominal_closed.v0)
    all_sextupole_names = getproperty.(sextupoles, :name)
    nominal_tracked = track_orbits_at_names(
        ring,
        nominal_closed,
        vcat(bpm_names, all_sextupole_names),
    )
    nominal_bpm_orbits = zeros(nd, 2)
    for (index, name) in enumerate(bpm_names)
        nominal_bpm_orbits[index, :] .= (
            nominal_tracked.horizontal[name][1], nominal_tracked.vertical[name][1],
        )
    end
    nominal_target_orbits = zeros(nt, 2)
    nominal_target_centers = zeros(nt, 2)
    for (index, ((_, entry), name)) in enumerate(zip(selected, target_names))
        nominal_target_orbits[index, :] .= (
            nominal_tracked.horizontal[name][1], nominal_tracked.vertical[name][1],
        )
        nominal_target_centers[index, :] .= (entry.x_offset_m, entry.y_offset_m)
    end

    calculation_start = time()
    try
        for (target_counter, (target_inventory_index, target_entry)) in enumerate(selected)
            target_start = time()
            target_element = ring.line[target_entry.index]
            target_knobs = knobs_by_target[target_entry.name]
            rng = MersenneTwister(seed + target_inventory_index)
            for realization in 1:nr
                realization_start = time()
                for entry in sextupoles
                    restore_sextupole!(ring, entry)
                end
                truth_x = (2rand(rng) - 1) * target_halfwidth
                truth_y = (2rand(rng) - 1) * target_halfwidth
                target_truth[target_counter, realization, :] .= (truth_x, truth_y)
                for (sextupole_index, entry) in enumerate(sextupoles)
                    dx, dy = if sextupole_index == target_inventory_index
                        truth_x, truth_y
                    else
                        other_sext_rms * randn(rng), other_sext_rms * randn(rng)
                    end
                    latent_offsets[target_counter, realization, sextupole_index, :] .= (dx, dy)
                    element = ring.line[entry.index]
                    element.x_offset = entry.x_offset_m + dx
                    element.y_offset = entry.y_offset_m + dy
                end

                current_v0 = copy(nominal_v0)
                for (bump_index, (bump_x, bump_y)) in enumerate(bumps)
                    for control in controls
                        cx, cy = target_knobs[(control.name, String(control.axis))]
                        set_corrector_scalar!(
                            ring,
                            control,
                            Float64(constant_term(first(control.originals))) +
                            cx*bump_x + cy*bump_y,
                        )
                    end
                    for (k2_index, level) in enumerate(k2_levels)
                        target_element.Kn2 = target_entry.kn2_m3 + level*k2_step
                        closed = scalar_warm(ring, current_v0)
                        current_v0 .= closed.v0
                        tracked = track_orbits_at_names(
                            ring,
                            closed,
                            vcat(bpm_names, [target_entry.name]),
                        )
                        for (bpm_index, bpm_name) in enumerate(bpm_names)
                            bpm_orbits[
                                target_counter, realization, bump_index, k2_index,
                                bpm_index, :,
                            ] .= (
                                tracked.horizontal[bpm_name][1],
                                tracked.vertical[bpm_name][1],
                            )
                        end
                        target_orbits[
                            target_counter, realization, bump_index, k2_index, :,
                        ] .= (
                            tracked.horizontal[target_entry.name][1],
                            tracked.vertical[target_entry.name][1],
                        )
                    end
                end
                realization_seconds[target_counter, realization] = time() - realization_start
            end
            @printf(
                "%s target %d/%d complete in %.1f s (elapsed %.1f s)\n",
                target_entry.name,
                target_counter,
                nt,
                time() - target_start,
                time() - calculation_start,
            )
            flush(stdout)
        end
    finally
        restore_correctors!(ring, controls)
        for entry in sextupoles
            restore_sextupole!(ring, entry)
        end
    end

    mkpath(output_dir)
    write_npy(joinpath(output_dir, "bpm_orbits.npy"), bpm_orbits)
    write_npy(joinpath(output_dir, "target_orbits.npy"), target_orbits)
    write_npy(joinpath(output_dir, "target_truth.npy"), target_truth)
    write_npy(joinpath(output_dir, "latent_sextupole_offsets.npy"), latent_offsets)
    write_npy(joinpath(output_dir, "nominal_bpm_orbits.npy"), nominal_bpm_orbits)
    write_npy(joinpath(output_dir, "nominal_target_orbits.npy"), nominal_target_orbits)
    write_npy(joinpath(output_dir, "nominal_target_centers.npy"), nominal_target_centers)
    write_npy(joinpath(output_dir, "realization_seconds.npy"), realization_seconds)
    write_lines(joinpath(output_dir, "target_names.txt"), target_names)
    write_lines(joinpath(output_dir, "bpm_names.txt"), bpm_names)
    write_lines(joinpath(output_dir, "sextupole_names.txt"), all_sextupole_names)
    write_rows(joinpath(output_dir, "target_inventory.csv"), [
        (; target_index=index, sextupole_inventory_index=inventory_index, target=name)
        for (index, (inventory_index, name)) in
            enumerate(zip(target_inventory_indices, target_names))
    ])
    write_rows(joinpath(output_dir, "bump_points.csv"), [
        (; bump_index=index, bump_x_command_m=x, bump_y_command_m=y)
        for (index, (x, y)) in enumerate(bumps)
    ])
    write_rows(joinpath(output_dir, "realization_timings.csv"), [
        (;
            target_index=target,
            target=target_names[target],
            realization,
            seconds=realization_seconds[target, realization],
        )
        for target in 1:nt for realization in 1:nr
    ])
    wall_seconds = time() - calculation_start
    write_metadata(metadata_path, Dict(
        "format" => "cesr-sextupole-misalignment-only-bpm-taylor-map-v1",
        "date" => string(Dates.today()),
        "engine" => "SciBmad exact scalar RF-on closed orbit and tracking",
        "lattice" => LATEST_LATTICE,
        "target_count" => nt,
        "target_inventory" => target_names,
        "realization_count_per_target" => nr,
        "random_seed_base" => seed,
        "target_seed_rule" => "base seed plus one-based full sextupole inventory index",
        "only_machine_error" => "fixed x/y misalignment of all 76 active normal sextupoles",
        "omitted_errors" => [
            "BPM noise/offset/gain/roll/missing channels",
            "time drift",
            "corrector calibration error",
            "target-K2 calibration error",
            "quadrupole strength/roll/misalignment error",
            "additional RF or lattice errors",
        ],
        "target_offset_distribution" => "independent uniform[-halfwidth,+halfwidth] in x and y",
        "target_offset_halfwidth_m" => target_halfwidth,
        "other_sextupole_offset_distribution" => "independent Gaussian in x and y, fixed within each scan tensor",
        "other_sextupole_offset_rms_m" => other_sext_rms,
        "sextupole_count" => ns,
        "bump_protocol" => "5 x 5 Cartesian grid at normalized coordinates -1,-0.5,0,0.5,1",
        "bump_count" => nb,
        "bump_amplitude_m" => bump_amplitude,
        "k2_levels" => k2_levels,
        "k2_count" => nk,
        "k2_step_m3" => k2_step,
        "state_count_per_realization" => nb*nk,
        "total_state_count" => nt*nr*nb*nk,
        "bpm_count" => nd,
        "calculation_wall_seconds" => wall_seconds,
        "truth_semantics" => "incremental physical target offset relative to the validated latest-lattice element offset",
    ))
    println("Wrote only-sextupole-misalignment scans to $output_dir in $(round(wall_seconds; digits=1)) s")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end

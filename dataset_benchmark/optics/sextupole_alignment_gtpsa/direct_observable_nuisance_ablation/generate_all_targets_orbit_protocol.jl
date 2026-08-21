#!/usr/bin/env julia

"""Generate the 76-target, five-axial-bump, three-K2, nominal-K1 protocol.

Each target is evaluated in eight independent latent machines by default. The
other 75 sextupoles carry unknown Gaussian x/y offsets and every active
quadrupole carries a fixed uniform ±1% physical strength error within one scan
tensor. Only the selected target K2 and the local orbit bump are intervened on.
"""

include(joinpath(
    @__DIR__, "..", "quadrupole_affinity", "exact_11_triplet_validation", "common.jl",
))

const ALL_TARGET_BUMP_KNOBS = joinpath(
    @__DIR__, "..", "quadrupole_affinity", "exact_11_triplet_validation",
    "results", "bump_knobs", "local_bump_knobs.csv",
)

function scalar_warm_all_targets(ring, v0)
    solution = find_closed_orbit(
        ring; v0=copy(v0), coasting_beam=false, batch=Val{false}(), warn=false,
    )
    all(solution.sol.retcode .== SciBmad.BatchSolve.RETCODE_SUCCESS) ||
        error("RF-on closed-orbit solve failed: $(solution.sol.retcode)")
    return solution
end

function restore_all_target_controls!(ring, controls)
    for control in controls
        for (index, original) in zip(control.indices, control.originals)
            if control.axis == :Kn0
                ring.line[index].Kn0 = original
            else
                ring.line[index].Ks0 = original
            end
        end
    end
end

function main(args=ARGS)
    defaults = Dict(
        "realizations" => "8",
        "seed" => "20260818",
        "target-halfwidth-m" => "3.5e-4",
        "other-sext-rms-m" => "3.0e-4",
        "quadrupole-error-fraction" => "0.01",
        "bump-amplitude-m" => "5.0e-4",
        "k2-step-m3" => "0.01",
        "bump-knobs-csv" => ALL_TARGET_BUMP_KNOBS,
        "output-dir" => joinpath(@__DIR__, "results", "all_76_orbit_protocol"),
        "overwrite" => "false",
    )
    options = parse_exact11_options(defaults, args)
    nr = parse(Int, options["realizations"])
    seed = parse(Int, options["seed"])
    target_halfwidth = parse(Float64, options["target-halfwidth-m"])
    other_sext_rms = parse(Float64, options["other-sext-rms-m"])
    quadrupole_error_fraction = parse(Float64, options["quadrupole-error-fraction"])
    bump_amplitude = parse(Float64, options["bump-amplitude-m"])
    k2_step = parse(Float64, options["k2-step-m3"])
    output_dir = abspath(options["output-dir"])
    metadata_path = joinpath(output_dir, "scan_metadata.toml")
    isfile(metadata_path) && lowercase(options["overwrite"]) != "true" &&
        error("Output exists; use --overwrite=true: $metadata_path")

    ring = cesr
    sextupoles = active_sextupole_inventory(ring)
    quadrupoles = active_quadrupole_inventory(ring)
    detectors = measurable_bpms(ring)
    controls = independent_corrector_inventory(ring)
    length(sextupoles) == 76 || error("Expected 76 active sextupoles")
    bpm_names = String.(base_name.(detectors))
    target_names = getproperty.(sextupoles, :name)

    knob_rows = read_simple_csv(abspath(options["bump-knobs-csv"]))
    knobs_by_target = Dict{String,Dict{Tuple{String,String},Tuple{Float64,Float64}}}()
    for row in knob_rows
        target = uppercase(row["target_sextupole"])
        target_knobs = get!(knobs_by_target, target, Dict{Tuple{String,String},Tuple{Float64,Float64}}())
        target_knobs[(row["corrector"], row["field"])] = (
            parse(Float64, row["field_per_x_bump_m"]),
            parse(Float64, row["field_per_y_bump_m"]),
        )
    end
    Set(keys(knobs_by_target)) == Set(target_names) ||
        error("Bump knob target inventory does not cover all 76 sextupoles")
    control_keys = Set((control.name, String(control.axis)) for control in controls)
    all(Set(keys(knobs_by_target[name])) == control_keys for name in target_names) ||
        error("At least one target bump knob has a control mismatch")

    bumps = [
        (-bump_amplitude, 0.0),
        (0.0, -bump_amplitude),
        (0.0, 0.0),
        (0.0, bump_amplitude),
        (bump_amplitude, 0.0),
    ]
    k2_levels = [-2.0, 0.0, 2.0]
    nt, nb, nk, nd, nq = length(sextupoles), length(bumps), length(k2_levels), length(detectors), length(quadrupoles)
    bpm_orbits = zeros(nt, nr, nb, nk, nd, 2)
    target_orbits = zeros(nt, nr, nb, nk, 2)
    target_truth = zeros(nt, nr, 2)
    latent_sextupole_offsets = zeros(nt, nr, nt, 2)
    latent_quadrupole_errors = zeros(nt, nr, nq)
    realization_seconds = zeros(nt, nr)
    nominal_closed = solve_closed_orbit(ring)
    nominal_v0 = copy(nominal_closed.v0)
    calculation_start = time()

    try
        for (target_index, target_entry) in enumerate(sextupoles)
            target_start = time()
            target_element = ring.line[target_entry.index]
            target_knobs = knobs_by_target[target_entry.name]
            rng = MersenneTwister(seed + target_index)
            for realization in 1:nr
                realization_start = time()
                for entry in sextupoles
                    ring.line[entry.index].Kn2 = entry.kn2_m3
                end
                truth_x = (2rand(rng) - 1) * target_halfwidth
                truth_y = (2rand(rng) - 1) * target_halfwidth
                target_truth[target_index, realization, :] .= (truth_x, truth_y)
                for (sextupole_index, entry) in enumerate(sextupoles)
                    dx, dy = if sextupole_index == target_index
                        truth_x, truth_y
                    else
                        other_sext_rms * randn(rng), other_sext_rms * randn(rng)
                    end
                    latent_sextupole_offsets[target_index, realization, sextupole_index, :] .= (dx, dy)
                    element = ring.line[entry.index]
                    element.x_offset = entry.x_offset_m + dx
                    element.y_offset = entry.y_offset_m + dy
                end
                for (quadrupole_index, entry) in enumerate(quadrupoles)
                    relative_error = (2rand(rng) - 1) * quadrupole_error_fraction
                    latent_quadrupole_errors[target_index, realization, quadrupole_index] = relative_error
                    set_quadrupole!(ring, entry, entry.kn1_m2 * (1 + relative_error))
                end

                current_v0 = copy(nominal_v0)
                for (bump_index, (bump_x, bump_y)) in enumerate(bumps)
                    for control in controls
                        cx, cy = target_knobs[(control.name, String(control.axis))]
                        set_corrector_scalar!(
                            ring,
                            control,
                            Float64(constant_term(first(control.originals))) + cx*bump_x + cy*bump_y,
                        )
                    end
                    for (k2_index, level) in enumerate(k2_levels)
                        target_element.Kn2 = target_entry.kn2_m3 + level*k2_step
                        closed = scalar_warm_all_targets(ring, current_v0)
                        current_v0 .= closed.v0
                        tracked = track_orbits_at_names(
                            ring, closed, vcat(bpm_names, [target_entry.name]),
                        )
                        for (bpm_index, bpm_name) in enumerate(bpm_names)
                            bpm_orbits[target_index, realization, bump_index, k2_index, bpm_index, :] .= (
                                tracked.horizontal[bpm_name][1],
                                tracked.vertical[bpm_name][1],
                            )
                        end
                        target_orbits[target_index, realization, bump_index, k2_index, :] .= (
                            tracked.horizontal[target_entry.name][1],
                            tracked.vertical[target_entry.name][1],
                        )
                    end
                end
                realization_seconds[target_index, realization] = time() - realization_start
            end
            @printf(
                "%s target %d/%d complete in %.1f s (elapsed %.1f s)\n",
                target_entry.name, target_index, nt, time()-target_start, time()-calculation_start,
            )
            flush(stdout)
        end
    finally
        restore_all_target_controls!(ring, controls)
        for entry in quadrupoles
            restore_quadrupole!(ring, entry)
        end
        for entry in sextupoles
            restore_sextupole!(ring, entry)
        end
    end

    mkpath(output_dir)
    write_npy(joinpath(output_dir, "bpm_orbits.npy"), bpm_orbits)
    write_npy(joinpath(output_dir, "target_orbits.npy"), target_orbits)
    write_npy(joinpath(output_dir, "target_truth.npy"), target_truth)
    write_npy(joinpath(output_dir, "latent_sextupole_offsets.npy"), latent_sextupole_offsets)
    write_npy(joinpath(output_dir, "latent_quadrupole_relative_errors.npy"), latent_quadrupole_errors)
    write_npy(joinpath(output_dir, "realization_seconds.npy"), realization_seconds)
    write_lines(joinpath(output_dir, "target_names.txt"), target_names)
    write_lines(joinpath(output_dir, "bpm_names.txt"), bpm_names)
    write_lines(joinpath(output_dir, "quadrupole_names.txt"), getproperty.(quadrupoles, :name))
    write_rows(joinpath(output_dir, "bump_points.csv"), [
        (; bump_index=index, bump_x_command_m=x, bump_y_command_m=y)
        for (index, (x, y)) in enumerate(bumps)
    ])
    wall_seconds = time() - calculation_start
    write_metadata(metadata_path, Dict(
        "format" => "cesr-all-76-orbit-protocol-v1",
        "date" => string(Dates.today()),
        "engine" => "SciBmad exact scalar RF-on closed orbit and tracking",
        "lattice" => LATEST_LATTICE,
        "target_count" => nt,
        "realization_count_per_target" => nr,
        "state_count_per_realization" => nb*nk,
        "total_state_count" => nt*nr*nb*nk,
        "random_seed_base" => seed,
        "target_seed_rule" => "base seed plus one-based target inventory index",
        "target_offset_distribution" => "independent uniform[-halfwidth,+halfwidth] in x and y",
        "target_offset_halfwidth_m" => target_halfwidth,
        "other_sextupole_count" => nt-1,
        "other_sextupole_offset_distribution" => "independent Gaussian in x and y, fixed within each scan tensor",
        "other_sextupole_offset_rms_m" => other_sext_rms,
        "quadrupole_count" => nq,
        "quadrupole_relative_error_distribution" => "independent uniform[-fraction,+fraction], fixed within each scan tensor",
        "quadrupole_error_fraction" => quadrupole_error_fraction,
        "k1_protocol" => "nominal commands only; random physical quadrupole errors remain active",
        "bump_protocol" => "five-point axial cross: (-x,0),(0,-y),(0,0),(0,+y),(+x,0)",
        "bump_count" => nb,
        "bump_amplitude_m" => bump_amplitude,
        "k2_count" => nk,
        "k2_levels" => k2_levels,
        "k2_step_m3" => k2_step,
        "bpm_count" => nd,
        "local_orbit_measurement_noise" => "none",
        "calculation_wall_seconds" => wall_seconds,
    ))
    println("Wrote all-target protocol to $output_dir in $(round(wall_seconds; digits=1)) s")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end

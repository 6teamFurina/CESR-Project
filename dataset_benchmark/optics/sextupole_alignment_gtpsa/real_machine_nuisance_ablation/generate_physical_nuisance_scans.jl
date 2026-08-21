#!/usr/bin/env julia

"""Generate paired exact-SciBmad scans for physical nuisance studies.

Every case shares the same target truth and all-sextupole alignment realization.
The baseline deliberately contains none of the nuisances studied here.  Each
legacy non-baseline case activates exactly one additional physical nuisance.
The two combined cases activate every maintained nuisance except quadrupole
misalignment, with and without the paired drift secant needed by the repeated
eight-state time-series inverse.
"""

include(joinpath(
    @__DIR__, "..", "quadrupole_affinity", "exact_11_triplet_validation", "common.jl",
))

const BUMP_KNOBS = joinpath(
    @__DIR__, "..", "quadrupole_affinity", "exact_11_triplet_validation",
    "results", "bump_knobs", "local_bump_knobs.csv",
)

const PHYSICAL_CASES = [
    "baseline",
    "corrector_gain",
    "k2_calibration",
    "quadrupole_strength",
    "quadrupole_roll",
    "quadrupole_misalignment",
    "time_drift",
    "combined_without_quadrupole_misalignment",
    "combined_without_quadrupole_misalignment_time_drift",
]

const CASE_SEED_OFFSETS = Dict(
    "baseline" => 0,
    "corrector_gain" => 10_000_000,
    "k2_calibration" => 20_000_000,
    "quadrupole_strength" => 30_000_000,
    "quadrupole_roll" => 40_000_000,
    "quadrupole_misalignment" => 50_000_000,
    "time_drift" => 60_000_000,
)

function scalar_closed_orbit(ring, v0)
    solution = find_closed_orbit(
        ring; v0=copy(v0), coasting_beam=false, batch=Val{false}(), warn=false,
    )
    all(solution.sol.retcode .== SciBmad.BatchSolve.RETCODE_SUCCESS) ||
        error("RF-on closed-orbit solve failed: $(solution.sol.retcode)")
    return solution
end

function restore_controls!(ring, controls)
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

function quadrupole_geometry_inventory(ring, quadrupoles)
    return Dict(
        index => (;
            x_offset=ring.line[index].x_offset,
            y_offset=ring.line[index].y_offset,
            tilt=ring.line[index].tilt,
        )
        for entry in quadrupoles for index in entry.indices
    )
end

function restore_quadrupole_geometry!(ring, geometry)
    for (index, original) in geometry
        ring.line[index].x_offset = original.x_offset
        ring.line[index].y_offset = original.y_offset
        ring.line[index].tilt = original.tilt
    end
    return nothing
end

function restore_machine!(ring, sextupoles, quadrupoles, controls, geometry)
    restore_controls!(ring, controls)
    for entry in quadrupoles
        restore_quadrupole!(ring, entry)
    end
    restore_quadrupole_geometry!(ring, geometry)
    for entry in sextupoles
        restore_sextupole!(ring, entry)
    end
    return nothing
end

function read_knobs(path, target_names, controls)
    rows = read_simple_csv(abspath(path))
    result = Dict{String,Dict{Tuple{String,String},Tuple{Float64,Float64}}}()
    for row in rows
        target = uppercase(row["target_sextupole"])
        target_rows = get!(
            result, target, Dict{Tuple{String,String},Tuple{Float64,Float64}}(),
        )
        target_rows[(row["corrector"], row["field"])] = (
            parse(Float64, row["field_per_x_bump_m"]),
            parse(Float64, row["field_per_y_bump_m"]),
        )
    end
    Set(keys(result)) == Set(target_names) || error("Bump knobs do not cover target inventory")
    control_keys = Set((control.name, String(control.axis)) for control in controls)
    all(Set(keys(result[name])) == control_keys for name in target_names) ||
        error("At least one target bump knob has a control mismatch")
    return result
end

function set_base_latent_machine!(
    ring, sextupoles, target_index, realization, seed, target_halfwidth,
    other_sext_rms, target_truth, latent_sextupole_offsets,
)
    rng = MersenneTwister(seed + 10_000 * target_index + realization)
    truth_x = (2rand(rng) - 1) * target_halfwidth
    truth_y = (2rand(rng) - 1) * target_halfwidth
    target_truth[target_index, realization, :] .= (truth_x, truth_y)
    for (sextupole_index, entry) in enumerate(sextupoles)
        dx, dy = if sextupole_index == target_index
            truth_x, truth_y
        else
            other_sext_rms * randn(rng), other_sext_rms * randn(rng)
        end
        latent_sextupole_offsets[
            target_index, realization, sextupole_index, :,
        ] .= (dx, dy)
        element = ring.line[entry.index]
        element.Kn2 = entry.kn2_m3
        element.x_offset = entry.x_offset_m + dx
        element.y_offset = entry.y_offset_m + dy
    end
    return nothing
end

function apply_physical_nuisance!(
    case_name, ring, quadrupoles, controls, geometry, target_index, realization,
    seed, corrector_gain_rms, k2_gain_rms, quadrupole_strength_fraction,
    quadrupole_roll_rms, quadrupole_alignment_rms,
    latent_corrector_gain_errors, latent_k2_gain_errors,
    latent_quadrupole_relative_errors, latent_quadrupole_rolls,
    latent_quadrupole_offsets, latent_drift_directions,
)
    if case_name == "combined_without_quadrupole_misalignment" ||
       case_name == "combined_without_quadrupole_misalignment_time_drift"
        components = [
            "corrector_gain",
            "k2_calibration",
            "quadrupole_strength",
            "quadrupole_roll",
        ]
        case_name == "combined_without_quadrupole_misalignment_time_drift" &&
            push!(components, "time_drift")
        for component in components
            apply_physical_nuisance!(
                component, ring, quadrupoles, controls, geometry, target_index,
                realization, seed, corrector_gain_rms, k2_gain_rms,
                quadrupole_strength_fraction, quadrupole_roll_rms,
                quadrupole_alignment_rms, latent_corrector_gain_errors,
                latent_k2_gain_errors, latent_quadrupole_relative_errors,
                latent_quadrupole_rolls, latent_quadrupole_offsets,
                latent_drift_directions,
            )
        end
        return nothing
    end
    nuisance_seed = seed + CASE_SEED_OFFSETS[case_name] + 10_000 * target_index + realization
    rng = MersenneTwister(nuisance_seed)
    if case_name == "corrector_gain"
        for (control_index, _) in enumerate(controls)
            latent_corrector_gain_errors[target_index, realization, control_index] =
                corrector_gain_rms * randn(rng)
        end
    elseif case_name == "k2_calibration"
        latent_k2_gain_errors[target_index, realization] = k2_gain_rms * randn(rng)
    elseif case_name == "quadrupole_strength"
        for (quadrupole_index, entry) in enumerate(quadrupoles)
            relative_error = (2rand(rng) - 1) * quadrupole_strength_fraction
            latent_quadrupole_relative_errors[
                target_index, realization, quadrupole_index,
            ] = relative_error
            set_quadrupole!(ring, entry, entry.kn1_m2 * (1 + relative_error))
        end
    elseif case_name == "quadrupole_roll"
        for (quadrupole_index, entry) in enumerate(quadrupoles)
            roll = quadrupole_roll_rms * randn(rng)
            latent_quadrupole_rolls[target_index, realization, quadrupole_index] = roll
            for index in entry.indices
                ring.line[index].tilt = geometry[index].tilt + roll
            end
        end
    elseif case_name == "quadrupole_misalignment"
        for (quadrupole_index, entry) in enumerate(quadrupoles)
            dx = quadrupole_alignment_rms * randn(rng)
            dy = quadrupole_alignment_rms * randn(rng)
            latent_quadrupole_offsets[
                target_index, realization, quadrupole_index, :,
            ] .= (dx, dy)
            for index in entry.indices
                ring.line[index].x_offset = geometry[index].x_offset + dx
                ring.line[index].y_offset = geometry[index].y_offset + dy
            end
        end
    elseif case_name == "time_drift"
        angle = 2pi * rand(rng)
        latent_drift_directions[target_index, realization, :] .= (cos(angle), sin(angle))
    elseif case_name != "baseline"
        error("Unsupported physical nuisance case: $case_name")
    end
    return nothing
end

function generate_case(options, case_name, ring, sextupoles, quadrupoles, detectors,
                       controls, geometry, knobs_by_target, nominal_v0)
    nr = parse(Int, options["realizations"])
    seed = parse(Int, options["seed"])
    target_limit = parse(Int, options["target-limit"])
    nt_all = length(sextupoles)
    nt = target_limit == 0 ? nt_all : min(target_limit, nt_all)
    selected_targets = sextupoles[1:nt]
    target_halfwidth = parse(Float64, options["target-halfwidth-m"])
    other_sext_rms = parse(Float64, options["other-sext-rms-m"])
    bump_amplitude = parse(Float64, options["bump-amplitude-m"])
    k2_step = parse(Float64, options["k2-step-m3"])
    corrector_gain_rms = parse(Float64, options["corrector-gain-rms"])
    k2_gain_rms = parse(Float64, options["k2-gain-rms"])
    quadrupole_strength_fraction = parse(Float64, options["quadrupole-strength-fraction"])
    quadrupole_roll_rms = parse(Float64, options["quadrupole-roll-rms-rad"])
    quadrupole_alignment_rms = parse(Float64, options["quadrupole-alignment-rms-m"])
    drift_halfwidth = parse(Float64, options["drift-halfwidth-m"])
    output_dir = joinpath(abspath(options["output-root"]), case_name)
    metadata_path = joinpath(output_dir, "scan_metadata.toml")
    isfile(metadata_path) && lowercase(options["overwrite"]) != "true" &&
        error("Output exists; use --overwrite=true: $metadata_path")

    target_names = String.(getproperty.(selected_targets, :name))
    bpm_names = String.(base_name.(detectors))
    bumps = [
        (-bump_amplitude, 0.0),
        (0.0, -bump_amplitude),
        (0.0, 0.0),
        (0.0, bump_amplitude),
        (bump_amplitude, 0.0),
    ]
    k2_levels = [-2.0, 0.0, 2.0]
    nb, nk, nd = length(bumps), length(k2_levels), length(detectors)
    nq, nc = length(quadrupoles), length(controls)
    state_count = nb * nk
    has_time_drift = case_name == "time_drift" ||
                     case_name == "combined_without_quadrupole_misalignment_time_drift"

    bpm_orbits = zeros(nt, nr, nb, nk, nd, 2)
    target_orbits = zeros(nt, nr, nb, nk, 2)
    target_truth = zeros(nt, nr, 2)
    latent_sextupole_offsets = zeros(nt, nr, nt_all, 2)
    latent_corrector_gain_errors = zeros(nt, nr, nc)
    latent_k2_gain_errors = zeros(nt, nr)
    latent_quadrupole_relative_errors = zeros(nt, nr, nq)
    latent_quadrupole_rolls = zeros(nt, nr, nq)
    latent_quadrupole_offsets = zeros(nt, nr, nq, 2)
    latent_drift_directions = zeros(nt, nr, 2)
    realization_seconds = zeros(nt, nr)
    calculation_start = time()

    try
        for (target_index, target_entry) in enumerate(selected_targets)
            target_start = time()
            target_knobs = knobs_by_target[target_entry.name]
            for realization in 1:nr
                realization_start = time()
                restore_machine!(ring, sextupoles, quadrupoles, controls, geometry)
                set_base_latent_machine!(
                    ring, sextupoles, target_index, realization, seed,
                    target_halfwidth, other_sext_rms, target_truth,
                    latent_sextupole_offsets,
                )
                apply_physical_nuisance!(
                    case_name, ring, quadrupoles, controls, geometry,
                    target_index, realization, seed, corrector_gain_rms,
                    k2_gain_rms, quadrupole_strength_fraction,
                    quadrupole_roll_rms, quadrupole_alignment_rms,
                    latent_corrector_gain_errors, latent_k2_gain_errors,
                    latent_quadrupole_relative_errors, latent_quadrupole_rolls,
                    latent_quadrupole_offsets, latent_drift_directions,
                )
                current_v0 = copy(nominal_v0)
                state_index = 0
                for (bump_index, (bump_x, bump_y)) in enumerate(bumps)
                    for (k2_index, level) in enumerate(k2_levels)
                        state_index += 1
                        drift_fraction = (state_index - (state_count + 1) / 2) /
                                         ((state_count - 1) / 2)
                        drift_x = has_time_drift ?
                            drift_halfwidth * drift_fraction *
                            latent_drift_directions[target_index, realization, 1] : 0.0
                        drift_y = has_time_drift ?
                            drift_halfwidth * drift_fraction *
                            latent_drift_directions[target_index, realization, 2] : 0.0
                        for (control_index, control) in enumerate(controls)
                            cx, cy = target_knobs[(control.name, String(control.axis))]
                            commanded_delta = cx * bump_x + cy * bump_y
                            gain = 1 + latent_corrector_gain_errors[
                                target_index, realization, control_index,
                            ]
                            drift_delta = cx * drift_x + cy * drift_y
                            set_corrector_scalar!(
                                ring, control,
                                Float64(constant_term(first(control.originals))) +
                                gain * commanded_delta + drift_delta,
                            )
                        end
                        k2_gain = 1 + latent_k2_gain_errors[target_index, realization]
                        ring.line[target_entry.index].Kn2 =
                            target_entry.kn2_m3 + k2_gain * level * k2_step
                        closed = scalar_closed_orbit(ring, current_v0)
                        current_v0 .= closed.v0
                        tracked = track_orbits_at_names(
                            ring, closed, vcat(bpm_names, [target_entry.name]),
                        )
                        for (bpm_index, bpm_name) in enumerate(bpm_names)
                            bpm_orbits[
                                target_index, realization, bump_index, k2_index,
                                bpm_index, :,
                            ] .= (
                                tracked.horizontal[bpm_name][1],
                                tracked.vertical[bpm_name][1],
                            )
                        end
                        target_orbits[
                            target_index, realization, bump_index, k2_index, :,
                        ] .= (
                            tracked.horizontal[target_entry.name][1],
                            tracked.vertical[target_entry.name][1],
                        )
                    end
                end
                realization_seconds[target_index, realization] = time() - realization_start
            end
            @printf(
                "%s %s target %d/%d complete in %.1f s (elapsed %.1f s)\n",
                case_name, target_entry.name, target_index, nt,
                time() - target_start, time() - calculation_start,
            )
            flush(stdout)
        end
    finally
        restore_machine!(ring, sextupoles, quadrupoles, controls, geometry)
    end

    mkpath(output_dir)
    write_npy(joinpath(output_dir, "bpm_orbits.npy"), bpm_orbits)
    write_npy(joinpath(output_dir, "target_orbits.npy"), target_orbits)
    write_npy(joinpath(output_dir, "target_truth.npy"), target_truth)
    write_npy(joinpath(output_dir, "latent_sextupole_offsets.npy"), latent_sextupole_offsets)
    write_npy(joinpath(output_dir, "latent_corrector_gain_errors.npy"), latent_corrector_gain_errors)
    write_npy(joinpath(output_dir, "latent_k2_gain_errors.npy"), latent_k2_gain_errors)
    write_npy(joinpath(output_dir, "latent_quadrupole_relative_errors.npy"), latent_quadrupole_relative_errors)
    write_npy(joinpath(output_dir, "latent_quadrupole_rolls.npy"), latent_quadrupole_rolls)
    write_npy(joinpath(output_dir, "latent_quadrupole_offsets.npy"), latent_quadrupole_offsets)
    write_npy(joinpath(output_dir, "latent_drift_directions.npy"), latent_drift_directions)
    write_npy(joinpath(output_dir, "realization_seconds.npy"), realization_seconds)
    write_lines(joinpath(output_dir, "target_names.txt"), target_names)
    write_lines(joinpath(output_dir, "bpm_names.txt"), bpm_names)
    write_lines(joinpath(output_dir, "quadrupole_names.txt"), String.(getproperty.(quadrupoles, :name)))
    write_lines(joinpath(output_dir, "corrector_labels.txt"), [
        "$(control.name):$(control.axis)" for control in controls
    ])
    write_rows(joinpath(output_dir, "bump_points.csv"), [
        (; bump_index=index, bump_x_command_m=x, bump_y_command_m=y)
        for (index, (x, y)) in enumerate(bumps)
    ])
    wall_seconds = time() - calculation_start
    write_metadata(metadata_path, Dict(
        "format" => "cesr-real-machine-nuisance-physical-scan-v1",
        "date" => string(Dates.today()),
        "engine" => "SciBmad exact scalar RF-on closed orbit and tracking",
        "lattice" => LATEST_LATTICE,
        "nuisance_case" => case_name,
        "target_count" => nt,
        "full_sextupole_count" => nt_all,
        "realization_count_per_target" => nr,
        "total_state_count" => nt * nr * state_count,
        "random_seed_base" => seed,
        "target_offset_halfwidth_m" => target_halfwidth,
        "other_sextupole_offset_rms_m" => other_sext_rms,
        "bump_count" => nb,
        "bump_amplitude_m" => bump_amplitude,
        "k2_count" => nk,
        "k2_levels" => k2_levels,
        "k2_step_m3" => k2_step,
        "bpm_count" => nd,
        "corrector_count" => nc,
        "quadrupole_count" => nq,
        "corrector_gain_rms" => corrector_gain_rms,
        "k2_gain_rms" => k2_gain_rms,
        "quadrupole_strength_fraction" => quadrupole_strength_fraction,
        "quadrupole_roll_rms_rad" => quadrupole_roll_rms,
        "quadrupole_alignment_rms_m" => quadrupole_alignment_rms,
        "drift_halfwidth_m" => drift_halfwidth,
        "drift_protocol" => "linear in acquisition order; zero at zero-bump nominal-K2 state",
        "measurement_nuisance" => "none in physical tensors",
        "calculation_wall_seconds" => wall_seconds,
    ))
    println("Wrote $case_name physical scan to $output_dir in $(round(wall_seconds; digits=1)) s")
    return nothing
end

function main(args=ARGS)
    defaults = Dict(
        "cases" => join(PHYSICAL_CASES, ','),
        "realizations" => "4",
        "target-limit" => "0",
        "seed" => "20260819",
        "target-halfwidth-m" => "3.5e-4",
        "other-sext-rms-m" => "3.0e-4",
        "bump-amplitude-m" => "5.0e-4",
        "k2-step-m3" => "0.01",
        "corrector-gain-rms" => "0.01",
        "k2-gain-rms" => "0.01",
        "quadrupole-strength-fraction" => "0.01",
        "quadrupole-roll-rms-rad" => "1.0e-3",
        "quadrupole-alignment-rms-m" => "1.0e-4",
        "drift-halfwidth-m" => "5.0e-6",
        "bump-knobs-csv" => BUMP_KNOBS,
        "output-root" => joinpath(@__DIR__, "results", "physical_scans"),
        "overwrite" => "false",
    )
    options = parse_exact11_options(defaults, args)
    cases = strip.(split(options["cases"], ','))
    isempty(cases) && error("At least one nuisance case is required")
    all(case in PHYSICAL_CASES for case in cases) ||
        error("Unknown case; allowed cases are $(join(PHYSICAL_CASES, ", "))")

    ring = cesr
    sextupoles = active_sextupole_inventory(ring)
    quadrupoles = active_quadrupole_inventory(ring)
    detectors = measurable_bpms(ring)
    controls = independent_corrector_inventory(ring)
    geometry = quadrupole_geometry_inventory(ring, quadrupoles)
    target_names = String.(getproperty.(sextupoles, :name))
    knobs_by_target = read_knobs(options["bump-knobs-csv"], target_names, controls)
    nominal_v0 = copy(solve_closed_orbit(ring).v0)
    for case_name in cases
        generate_case(
            options, case_name, ring, sextupoles, quadrupoles, detectors,
            controls, geometry, knobs_by_target, nominal_v0,
        )
    end
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end

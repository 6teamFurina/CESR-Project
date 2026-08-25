#!/usr/bin/env julia

"""Generate machine-indexed, one-target-at-a-time sextupole scan tensors.

One latent machine realization contains all 76 sextupole offsets and all
static nuisance draws.  Those values remain fixed while every target
sextupole is scanned in sequence.  The two maintained cases share the same
draws; the second case additionally applies a fixed quadrupole alignment drift.

The primary 50 micrometer setting is the per-plane RMS of an independent
Gaussian x/y displacement for each physical quadrupole, coherent across all
tracking slices belonging to that quadrupole.  It represents residual annual
drift relative to the latest measured nominal geometry, not survey uncertainty
within one scan and not a time-varying offset during acquisition.
"""

include(joinpath(
    @__DIR__, "..", "quadrupole_affinity", "exact_11_triplet_validation", "common.jl",
))

const JOINT_BUMP_KNOBS = joinpath(
    @__DIR__, "..", "quadrupole_affinity", "exact_11_triplet_validation",
    "results", "bump_knobs", "local_bump_knobs.csv",
)

const JOINT_CASES = [
    "without_quadrupole_misalignment",
    "with_quadrupole_misalignment",
]

function scalar_closed_orbit_joint(ring, v0)
    solution = find_closed_orbit(
        ring;
        v0=copy(v0),
        coasting_beam=false,
        batch=Val{false}(),
        warn=false,
    )
    all(solution.sol.retcode .== SciBmad.BatchSolve.RETCODE_SUCCESS) ||
        error("RF-on closed-orbit solve failed: $(solution.sol.retcode)")
    return solution
end

function quadrupole_geometry_joint(ring, quadrupoles)
    return Dict(
        index => (;
            x_offset=ring.line[index].x_offset,
            y_offset=ring.line[index].y_offset,
            tilt=ring.line[index].tilt,
        )
        for entry in quadrupoles for index in entry.indices
    )
end

function restore_joint_machine!(ring, sextupoles, quadrupoles, controls, geometry)
    for control in controls
        for (index, original) in zip(control.indices, control.originals)
            if control.axis == :Kn0
                ring.line[index].Kn0 = original
            else
                ring.line[index].Ks0 = original
            end
        end
    end
    for entry in quadrupoles
        restore_quadrupole!(ring, entry)
        for index in entry.indices
            ring.line[index].x_offset = geometry[index].x_offset
            ring.line[index].y_offset = geometry[index].y_offset
            ring.line[index].tilt = geometry[index].tilt
        end
    end
    for entry in sextupoles
        restore_sextupole!(ring, entry)
    end
    return nothing
end

function materialize_multipole_group!(element)
    multipoles = element.BMultipoleParams
    isnothing(multipoles) && return nothing
    scalar_value(value) = value isa BatchParam ?
        Float64(value.batch isa AbstractArray ? first(value.batch) : value.batch) :
        Float64(constant_term(value))
    element.BMultipoleParams = BMultipoleParams(
        scalar_value.(multipoles.n),
        scalar_value.(multipoles.s),
        scalar_value.(multipoles.tilt),
        multipoles.order,
        multipoles.normalized,
        multipoles.integrated,
    )
    return nothing
end

function materialize_batch_elements!(ring, sextupoles, controls)
    indices = Set{Int}()
    for control in controls
        union!(indices, control.indices)
    end
    union!(indices, getproperty.(sextupoles, :index))
    for index in indices
        materialize_multipole_group!(ring.line[index])
    end
    return nothing
end

function read_joint_knobs(path, target_names, controls)
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
    Set(keys(result)) == Set(target_names) ||
        error("Bump knobs do not cover the selected sextupole inventory")
    control_keys = Set((control.name, String(control.axis)) for control in controls)
    all(Set(keys(result[name])) == control_keys for name in target_names) ||
        error("At least one target bump knob has a control mismatch")
    return result
end

function sample_joint_latents(
    machine_count, sextupole_count, quadrupole_count, corrector_count, bpm_count,
    seed, sextupole_rms, corrector_gain_rms, k2_gain_rms, bpm_gain_rms,
    quadrupole_strength_fraction, quadrupole_roll_rms,
)
    sextupole_offsets = zeros(machine_count, sextupole_count, 2)
    corrector_gain_errors = zeros(machine_count, corrector_count)
    k2_gain_errors = zeros(machine_count, sextupole_count)
    quadrupole_relative_errors = zeros(machine_count, quadrupole_count)
    quadrupole_rolls = zeros(machine_count, quadrupole_count)
    quadrupole_alignment_standard_normals = zeros(machine_count, quadrupole_count, 2)
    bpm_gain_errors = zeros(machine_count, bpm_count, 2)
    drift_directions = zeros(machine_count, sextupole_count, 2)

    for machine in 1:machine_count
        sext_rng = MersenneTwister(seed + 1_000_000 + machine)
        corrector_rng = MersenneTwister(seed + 2_000_000 + machine)
        k2_rng = MersenneTwister(seed + 3_000_000 + machine)
        strength_rng = MersenneTwister(seed + 4_000_000 + machine)
        roll_rng = MersenneTwister(seed + 5_000_000 + machine)
        alignment_rng = MersenneTwister(seed + 6_000_000 + machine)
        bpm_rng = MersenneTwister(seed + 7_000_000 + machine)
        drift_rng = MersenneTwister(seed + 8_000_000 + machine)

        sextupole_offsets[machine, :, :] .=
            sextupole_rms .* randn(sext_rng, sextupole_count, 2)
        corrector_gain_errors[machine, :] .=
            corrector_gain_rms .* randn(corrector_rng, corrector_count)
        k2_gain_errors[machine, :] .=
            k2_gain_rms .* randn(k2_rng, sextupole_count)
        quadrupole_relative_errors[machine, :] .=
            quadrupole_strength_fraction .* (2rand(strength_rng, quadrupole_count) .- 1)
        quadrupole_rolls[machine, :] .=
            quadrupole_roll_rms .* randn(roll_rng, quadrupole_count)
        quadrupole_alignment_standard_normals[machine, :, :] .=
            randn(alignment_rng, quadrupole_count, 2)
        bpm_gain_errors[machine, :, :] .= bpm_gain_rms .* randn(bpm_rng, bpm_count, 2)
        for target in 1:sextupole_count
            angle = 2pi * rand(drift_rng)
            drift_directions[machine, target, :] .= (cos(angle), sin(angle))
        end
    end
    return (;
        sextupole_offsets,
        corrector_gain_errors,
        k2_gain_errors,
        quadrupole_relative_errors,
        quadrupole_rolls,
        quadrupole_alignment_standard_normals,
        bpm_gain_errors,
        drift_directions,
    )
end

function apply_joint_machine!(
    ring, case_name, machine, sextupoles, quadrupoles, geometry, latents,
    quadrupole_alignment_rms,
)
    for (index, entry) in enumerate(sextupoles)
        dx, dy = latents.sextupole_offsets[machine, index, :]
        element = ring.line[entry.index]
        element.Kn2 = entry.kn2_m3
        element.x_offset = entry.x_offset_m + dx
        element.y_offset = entry.y_offset_m + dy
    end
    include_alignment = case_name == "with_quadrupole_misalignment"
    for (index, entry) in enumerate(quadrupoles)
        relative_error = latents.quadrupole_relative_errors[machine, index]
        roll = latents.quadrupole_rolls[machine, index]
        dx = include_alignment ?
            quadrupole_alignment_rms * latents.quadrupole_alignment_standard_normals[machine, index, 1] : 0.0
        dy = include_alignment ?
            quadrupole_alignment_rms * latents.quadrupole_alignment_standard_normals[machine, index, 2] : 0.0
        set_quadrupole!(ring, entry, entry.kn1_m2 * (1 + relative_error))
        for lattice_index in entry.indices
            ring.line[lattice_index].tilt = geometry[lattice_index].tilt + roll
            ring.line[lattice_index].x_offset = geometry[lattice_index].x_offset + dx
            ring.line[lattice_index].y_offset = geometry[lattice_index].y_offset + dy
        end
    end
    return nothing
end

function track_joint_machine_reference(ring, closed, bpm_names, target_names)
    tracked = track_orbits_at_names(ring, closed, vcat(bpm_names, target_names))
    bpm = zeros(length(bpm_names), 2)
    target = zeros(length(target_names), 2)
    for (index, name) in enumerate(bpm_names)
        bpm[index, :] .= (tracked.horizontal[name][1], tracked.vertical[name][1])
    end
    for (index, name) in enumerate(target_names)
        target[index, :] .= (tracked.horizontal[name][1], tracked.vertical[name][1])
    end
    return bpm, target
end

function scan_joint_target!(
    ring, target_index, target_entry, target_knobs, controls, bpm_names,
    bumps, delta_k2, base_closed, corrector_gain_errors, k2_gain_error,
    drift_direction, drift_halfwidth,
)
    nb, nk = length(bumps), length(delta_k2)
    physical_state_count = nb * nk
    total_state_count = 2physical_state_count
    bump_x = zeros(total_state_count)
    bump_y = zeros(total_state_count)
    k2_values = zeros(total_state_count)
    drift_fraction = zeros(total_state_count)
    state = 0
    for copy_index in 1:2
        for (bump_index, (x, y)) in enumerate(bumps)
            for (k2_index, k2) in enumerate(delta_k2)
                state += 1
                bump_x[state] = x
                bump_y[state] = y
                k2_values[state] = target_entry.kn2_m3 + (1 + k2_gain_error) * k2
                if copy_index == 2
                    base_index = (bump_index - 1) * nk + k2_index
                    drift_fraction[state] =
                        (base_index - (physical_state_count + 1) / 2) /
                        ((physical_state_count - 1) / 2)
                end
            end
        end
    end
    for (control_index, control) in enumerate(controls)
        cx, cy = target_knobs[(control.name, String(control.axis))]
        baseline = Float64(constant_term(first(control.originals)))
        command = cx .* bump_x .+ cy .* bump_y
        drift_command = cx .* drift_fraction .* drift_halfwidth .* drift_direction[1] .+
                        cy .* drift_fraction .* drift_halfwidth .* drift_direction[2]
        values = baseline .+ (1 + corrector_gain_errors[control_index]) .* command .+
                 drift_command
        set_corrector_values!(ring, control, values)
    end
    ring.line[target_entry.index].Kn2 = BatchParam(k2_values)
    initial = repeat(reshape(Array(base_closed.v0), 1, 6), total_state_count, 1)
    closed = solve_batch_closed_orbit(ring, total_state_count; initial_v0=initial)
    tracked = track_orbits_at_names(ring, closed, vcat(bpm_names, [target_entry.name]))

    reshape_states(values) = permutedims(reshape(values, nk, nb), (2, 1))

    static_bpm = zeros(nb, nk, length(bpm_names), 2)
    drift_bpm = similar(static_bpm)
    static_target = zeros(nb, nk, 2)
    drift_target = similar(static_target)
    for (bpm_index, name) in enumerate(bpm_names)
        x = tracked.horizontal[name]
        y = tracked.vertical[name]
        static_bpm[:, :, bpm_index, 1] .= reshape_states(x[1:physical_state_count])
        static_bpm[:, :, bpm_index, 2] .= reshape_states(y[1:physical_state_count])
        drift_bpm[:, :, bpm_index, 1] .= reshape_states(x[physical_state_count+1:end])
        drift_bpm[:, :, bpm_index, 2] .= reshape_states(y[physical_state_count+1:end])
    end
    target_x = tracked.horizontal[target_entry.name]
    target_y = tracked.vertical[target_entry.name]
    static_target[:, :, 1] .= reshape_states(target_x[1:physical_state_count])
    static_target[:, :, 2] .= reshape_states(target_y[1:physical_state_count])
    drift_target[:, :, 1] .= reshape_states(target_x[physical_state_count+1:end])
    drift_target[:, :, 2] .= reshape_states(target_y[physical_state_count+1:end])
    return static_bpm, drift_bpm, static_target, drift_target
end

function generate_joint_case(
    options, case_name, ring, sextupoles, quadrupoles, detectors, controls,
    geometry, knobs_by_target, nominal_v0, latents,
)
    machine_count = parse(Int, options["machines"])
    target_limit = parse(Int, options["target-limit"])
    selected = target_limit == 0 ? sextupoles : sextupoles[1:min(target_limit, length(sextupoles))]
    nt, nd = length(selected), length(detectors)
    full_nt = length(sextupoles)
    target_names = String.(getproperty.(selected, :name))
    bpm_names = String.(base_name.(detectors))
    bump_amplitude = parse(Float64, options["bump-amplitude-m"])
    k2_amplitude = parse(Float64, options["k2-amplitude-m3"])
    drift_halfwidth = parse(Float64, options["drift-halfwidth-m"])
    quadrupole_alignment_rms = parse(Float64, options["quadrupole-alignment-rms-m"])
    bumps = [
        (-bump_amplitude, 0.0),
        (0.0, -bump_amplitude),
        (0.0, 0.0),
        (0.0, bump_amplitude),
        (bump_amplitude, 0.0),
    ]
    delta_k2 = [-k2_amplitude, 0.0, k2_amplitude]
    nb, nk = length(bumps), length(delta_k2)

    bpm_orbits = zeros(machine_count, nt, nb, nk, nd, 2)
    drift_bpm_orbits = similar(bpm_orbits)
    target_orbits = zeros(machine_count, nt, nb, nk, 2)
    drift_target_orbits = similar(target_orbits)
    reference_bpm_orbits = zeros(machine_count, nd, 2)
    reference_target_orbits = zeros(machine_count, nt, 2)
    scan_seconds = zeros(machine_count, nt)
    output_dir = joinpath(abspath(options["output-root"]), case_name)
    metadata_path = joinpath(output_dir, "scan_metadata.toml")
    isfile(metadata_path) && lowercase(options["overwrite"]) != "true" &&
        error("Output exists; use --overwrite=true: $metadata_path")
    started = time()

    try
        for machine in 1:machine_count
            machine_start = time()
            restore_joint_machine!(ring, sextupoles, quadrupoles, controls, geometry)
            apply_joint_machine!(
                ring, case_name, machine, sextupoles, quadrupoles, geometry,
                latents, quadrupole_alignment_rms,
            )
            # Some latest-lattice correctors and combined-function sextupoles
            # carry zero-valued DefExpr multipole components.  Promoting those
            # containers directly to BatchParam attempts `Any(BatchParam)`.
            # Materialize only the elements parameterized in this scan, after
            # all physical nuisance values have been applied.
            materialize_batch_elements!(ring, sextupoles, controls)
            base_closed = scalar_closed_orbit_joint(ring, nominal_v0)
            bpm_reference, target_reference = track_joint_machine_reference(
                ring, base_closed, bpm_names, target_names,
            )
            reference_bpm_orbits[machine, :, :] .= bpm_reference
            reference_target_orbits[machine, :, :] .= target_reference

            for (target_index, target_entry) in enumerate(selected)
                target_start = time()
                target_knobs = knobs_by_target[target_entry.name]
                try
                    static_bpm, drift_bpm, static_target, drift_target = scan_joint_target!(
                        ring, target_index, target_entry, target_knobs, controls,
                        bpm_names, bumps, delta_k2, base_closed,
                        latents.corrector_gain_errors[machine, :],
                        latents.k2_gain_errors[machine, target_index],
                        latents.drift_directions[machine, target_index, :],
                        drift_halfwidth,
                    )
                    bpm_orbits[machine, target_index, :, :, :, :] .= static_bpm
                    drift_bpm_orbits[machine, target_index, :, :, :, :] .= drift_bpm
                    target_orbits[machine, target_index, :, :, :] .= static_target
                    drift_target_orbits[machine, target_index, :, :, :] .= drift_target
                finally
                    for control in controls
                        baseline = Float64(constant_term(first(control.originals)))
                        set_corrector_scalar!(ring, control, baseline)
                    end
                    ring.line[target_entry.index].Kn2 = target_entry.kn2_m3
                end
                scan_seconds[machine, target_index] = time() - target_start
                if target_index % 10 == 0 || target_index == nt
                    @printf(
                        "%s machine %d/%d target %d/%d complete (%.1f s machine elapsed)\n",
                        case_name, machine, machine_count, target_index, nt,
                        time() - machine_start,
                    )
                    flush(stdout)
                end
            end
            @printf(
                "%s machine %d/%d complete in %.1f s (elapsed %.1f s)\n",
                case_name, machine, machine_count, time() - machine_start, time() - started,
            )
            flush(stdout)
        end
    finally
        restore_joint_machine!(ring, sextupoles, quadrupoles, controls, geometry)
    end

    mkpath(output_dir)
    write_npy(joinpath(output_dir, "bpm_orbits.npy"), bpm_orbits)
    write_npy(joinpath(output_dir, "drift_bpm_orbits.npy"), drift_bpm_orbits)
    write_npy(joinpath(output_dir, "target_orbits.npy"), target_orbits)
    write_npy(joinpath(output_dir, "drift_target_orbits.npy"), drift_target_orbits)
    write_npy(joinpath(output_dir, "reference_bpm_orbits.npy"), reference_bpm_orbits)
    write_npy(joinpath(output_dir, "reference_target_orbits.npy"), reference_target_orbits)
    write_npy(joinpath(output_dir, "scan_seconds.npy"), scan_seconds)
    write_lines(joinpath(output_dir, "target_names.txt"), target_names)
    write_lines(joinpath(output_dir, "bpm_names.txt"), bpm_names)
    write_rows(joinpath(output_dir, "bump_points.csv"), [
        (; bump_index=index, bump_x_command_m=x, bump_y_command_m=y)
        for (index, (x, y)) in enumerate(bumps)
    ])
    wall_seconds = time() - started
    write_metadata(metadata_path, Dict(
        "format" => "cesr-sequential-joint-machine-scan-v1",
        "date" => string(Dates.today()),
        "engine" => "SciBmad exact RF-on BatchParam closed orbit and tracking",
        "lattice" => LATEST_LATTICE,
        "case" => case_name,
        "machine_count" => machine_count,
        "target_count" => nt,
        "full_sextupole_count" => full_nt,
        "state_count_per_target" => nb * nk,
        "batch_lane_count_per_target" => 2nb * nk,
        "total_exact_states" => 2machine_count * nt * nb * nk,
        "machine_atomic_unit" => "all latent errors fixed across all target scans",
        "paired_case_latents" => "all draws shared; quadrupole alignment is zero versus enabled",
        "sextupole_offset_distribution" => "independent Gaussian x/y, fixed across the machine realization",
        "sextupole_offset_rms_m" => parse(Float64, options["sextupole-offset-rms-m"]),
        "corrector_gain_rms" => parse(Float64, options["corrector-gain-rms"]),
        "k2_gain_rms" => parse(Float64, options["k2-gain-rms"]),
        "quadrupole_strength_distribution" => "independent uniform +/-fraction, fixed across all target scans",
        "quadrupole_strength_fraction" => parse(Float64, options["quadrupole-strength-fraction"]),
        "quadrupole_roll_rms_rad" => parse(Float64, options["quadrupole-roll-rms-rad"]),
        "quadrupole_alignment_distribution" => "independent Gaussian x/y per physical quadrupole, coherent across slices",
        "quadrupole_alignment_rms_m_per_plane" => quadrupole_alignment_rms,
        "quadrupole_alignment_semantics" => "annual residual drift relative to measured nominal; fixed during all scans",
        "bpm_gain_rms" => parse(Float64, options["bpm-gain-rms"]),
        "bump_amplitude_m" => bump_amplitude,
        "bump_count" => nb,
        "k2_delta_m3" => delta_k2,
        "k2_count" => nk,
        "drift_halfwidth_m" => drift_halfwidth,
        "drift_secant" => "same target bump knobs, linear state-order coefficient from -1 to +1",
        "bpm_count" => nd,
        "quadrupole_count" => length(quadrupoles),
        "corrector_count" => length(controls),
        "calculation_wall_seconds" => wall_seconds,
    ))
    println("Wrote $case_name joint-machine scans to $output_dir in $(round(wall_seconds; digits=1)) s")
    return nothing
end

function main(args=ARGS)
    defaults = Dict(
        "cases" => join(JOINT_CASES, ','),
        "machines" => "16",
        "target-limit" => "0",
        "seed" => "20260823",
        "sextupole-offset-rms-m" => "3.0e-4",
        "corrector-gain-rms" => "0.01",
        "k2-gain-rms" => "0.01",
        "bpm-gain-rms" => "0.01",
        "quadrupole-strength-fraction" => "0.01",
        "quadrupole-roll-rms-rad" => "1.0e-3",
        "quadrupole-alignment-rms-m" => "5.0e-5",
        "bump-amplitude-m" => "1.5e-3",
        "k2-amplitude-m3" => "1.0e-1",
        "drift-halfwidth-m" => "5.0e-6",
        "bump-knobs-csv" => JOINT_BUMP_KNOBS,
        "output-root" => joinpath(@__DIR__, "results", "exact_joint_machines"),
        "overwrite" => "false",
    )
    options = parse_exact11_options(defaults, args)
    cases = strip.(split(options["cases"], ','))
    isempty(cases) && error("At least one case is required")
    all(case in JOINT_CASES for case in cases) ||
        error("Unknown case; allowed cases are $(join(JOINT_CASES, ", "))")

    ring = cesr
    sextupoles = active_sextupole_inventory(ring)
    quadrupoles = active_quadrupole_inventory(ring)
    detectors = measurable_bpms(ring)
    controls = independent_corrector_inventory(ring)
    length(sextupoles) == 76 || error("Expected 76 active normal sextupoles")
    geometry = quadrupole_geometry_joint(ring, quadrupoles)
    target_names = String.(getproperty.(sextupoles, :name))
    knobs_by_target = read_joint_knobs(options["bump-knobs-csv"], target_names, controls)
    nominal_v0 = copy(solve_closed_orbit(ring).v0)
    latents = sample_joint_latents(
        parse(Int, options["machines"]), length(sextupoles), length(quadrupoles),
        length(controls), length(detectors), parse(Int, options["seed"]),
        parse(Float64, options["sextupole-offset-rms-m"]),
        parse(Float64, options["corrector-gain-rms"]),
        parse(Float64, options["k2-gain-rms"]),
        parse(Float64, options["bpm-gain-rms"]),
        parse(Float64, options["quadrupole-strength-fraction"]),
        parse(Float64, options["quadrupole-roll-rms-rad"]),
    )

    latent_root = joinpath(abspath(options["output-root"]), "paired_latents")
    mkpath(latent_root)
    write_npy(joinpath(latent_root, "sextupole_offsets.npy"), latents.sextupole_offsets)
    write_npy(joinpath(latent_root, "corrector_gain_errors.npy"), latents.corrector_gain_errors)
    write_npy(joinpath(latent_root, "k2_gain_errors.npy"), latents.k2_gain_errors)
    write_npy(joinpath(latent_root, "quadrupole_relative_errors.npy"), latents.quadrupole_relative_errors)
    write_npy(joinpath(latent_root, "quadrupole_rolls.npy"), latents.quadrupole_rolls)
    write_npy(
        joinpath(latent_root, "quadrupole_alignment_standard_normals.npy"),
        latents.quadrupole_alignment_standard_normals,
    )
    write_npy(joinpath(latent_root, "bpm_gain_errors.npy"), latents.bpm_gain_errors)
    write_npy(joinpath(latent_root, "drift_directions.npy"), latents.drift_directions)
    write_lines(joinpath(latent_root, "target_names.txt"), target_names)
    write_lines(joinpath(latent_root, "quadrupole_names.txt"), String.(getproperty.(quadrupoles, :name)))
    write_lines(joinpath(latent_root, "bpm_names.txt"), String.(base_name.(detectors)))
    write_metadata(joinpath(latent_root, "latent_metadata.toml"), Dict(
        "format" => "cesr-sequential-joint-paired-latents-v1",
        "date" => string(Dates.today()),
        "random_seed_base" => parse(Int, options["seed"]),
        "machine_count" => parse(Int, options["machines"]),
        "target_count" => length(sextupoles),
        "quadrupole_count" => length(quadrupoles),
        "corrector_count" => length(controls),
        "bpm_count" => length(detectors),
        "sextupole_offset_rms_m_per_plane" => parse(Float64, options["sextupole-offset-rms-m"]),
        "corrector_gain_rms" => parse(Float64, options["corrector-gain-rms"]),
        "k2_gain_rms" => parse(Float64, options["k2-gain-rms"]),
        "bpm_gain_rms" => parse(Float64, options["bpm-gain-rms"]),
        "quadrupole_strength_halfwidth_fraction" => parse(Float64, options["quadrupole-strength-fraction"]),
        "quadrupole_roll_rms_rad" => parse(Float64, options["quadrupole-roll-rms-rad"]),
        "quadrupole_alignment_standard_normals" => "unit Gaussian x/y draws; multiply by the case metadata RMS",
        "pairing" => "all latent arrays are shared by both physical cases",
    ))

    for case_name in cases
        generate_joint_case(
            options, case_name, ring, sextupoles, quadrupoles, detectors, controls,
            geometry, knobs_by_target, nominal_v0, latents,
        )
    end
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end

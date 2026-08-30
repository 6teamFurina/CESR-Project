#!/usr/bin/env julia

"""Generate sextupole scans after a fixed BPM-reference orbit correction.

For every latent machine, this script first measures a corrector-to-BPM
response matrix in the paired zero-quadrupole-offset state.  It then applies
the 50-micrometer/plane quadrupole offsets, solves one baseline correction from
the stored reference BPM readback, and holds that baseline command fixed while
all selected sextupoles are scanned one at a time.  Local bump and K2 commands
are superposed on the fixed corrected machine; no latent quadrupole offset or
target-local orbit is supplied to the correction or inverse.
"""

include(joinpath(@__DIR__, "run_quadrupole_orbit_correction.jl"))
include(joinpath(@__DIR__, "gtpsa_noisy_response.jl"))

const CORRECTED_JOINT_CASE = "with_quadrupole_misalignment_corrected"
const CORRECTED_SCAN_ROOT = joinpath(
    @__DIR__, "..", "sequential_joint_inverse", "results", "exact_joint_machines",
)

function parse_boolean_option(name, value)
    lowercase(value) in ("true", "false") ||
        error("--$name must be true or false")
    return lowercase(value) == "true"
end

function sample_baseline_corrector_gain_errors(
    machine_count, corrector_count, seed, rms,
)
    errors = zeros(machine_count, corrector_count)
    for machine in 1:machine_count
        rng = MersenneTwister(seed + 2_000_000 + machine)
        errors[machine, :] .= rms .* randn(rng, corrector_count)
    end
    return errors
end

function batch_element_indices(sextupoles, controls)
    indices = Set{Int}(getproperty.(sextupoles, :index))
    for control in controls
        union!(indices, control.indices)
    end
    return sort!(collect(indices))
end

function restore_batch_multipoles!(ring, originals)
    for (index, multipoles) in originals
        ring.line[index].BMultipoleParams = multipoles
    end
    return nothing
end

function physical_bump_control_baselines(ring, controls)
    baselines = zeros(length(controls))
    for (control_index, control) in enumerate(controls)
        values = [
            Float64(constant_term(
                control.axis == :Kn0 ? ring.line[index].Kn0 : ring.line[index].Ks0,
            ))
            for index in control.indices
        ]
        maximum(abs.(values .- first(values))) <= 1e-12 ||
            error("Corrected split-element field mismatch for $(control.name)")
        baselines[control_index] = first(values)
    end
    return baselines
end

"""Create one independent latest-lattice scan context for a Julia worker.

Each context owns its mutable ring, control registry, element inventories, and
saved multipole containers.  No SciBmad lattice object is shared between
simultaneous target scans.
"""
function build_target_scan_context(options)
    model = load_ring_model(; ring=:latest, zero_value=0.0, rf_on=true)
    ring = model.ring
    sextupoles = active_sextupole_inventory(ring)
    quadrupoles = active_quadrupole_inventory(ring)
    detectors = measurable_bpms(ring)
    local_controls = independent_corrector_inventory(ring)
    geometry = quadrupole_geometry_joint(ring, quadrupoles)
    target_names = String.(getproperty.(sextupoles, :name))
    knobs_by_target = read_joint_knobs(
        options["bump-knobs-csv"], target_names, local_controls,
    )
    indices = batch_element_indices(sextupoles, local_controls)
    batch_originals = Dict(
        index => ring.line[index].BMultipoleParams for index in indices
    )
    return (;
        model,
        ring,
        sextupoles,
        quadrupoles,
        detectors,
        local_controls,
        geometry,
        knobs_by_target,
        target_names,
        bpm_names=String.(base_name.(detectors)),
        control_names=String.(model.metadata.steering_control_names),
        batch_originals,
    )
end


"""Restore one worker to the fixed corrected machine for a target scan."""
function prepare_target_scan_context!(
    context,
    machine,
    latents,
    quadrupole_alignment_rms,
    physical_baseline_gains,
    baseline_commands,
    corrected_initial_v0,
)
    restore_batch_multipoles!(context.ring, context.batch_originals)
    prepare_machine!(
        context.model,
        "with_quadrupole_misalignment",
        machine,
        context.sextupoles,
        context.quadrupoles,
        context.geometry,
        latents,
        quadrupole_alignment_rms,
    )
    set_commanded_correctors!(
        context.model,
        context.control_names,
        physical_baseline_gains,
        baseline_commands,
    )
    closed = scalar_closed_orbit_joint(context.ring, corrected_initial_v0)
    local_baseline = physical_bump_control_baselines(
        context.ring, context.local_controls,
    )
    materialize_batch_elements!(
        context.ring, context.sextupoles, context.local_controls,
    )
    closed = scalar_closed_orbit_joint(context.ring, closed.v0)
    return closed, local_baseline
end

function write_or_validate_joint_latents!(
    output_root,
    options,
    latents,
    sextupoles,
    quadrupoles,
    detectors,
    local_controls,
)
    latent_root = joinpath(abspath(output_root), "paired_latents")
    metadata_path = joinpath(latent_root, "latent_metadata.toml")
    if isfile(metadata_path)
        metadata = TOML.parsefile(metadata_path)
        checks = (
            ("random_seed_base", parse(Int, options["seed"])),
            ("machine_count", parse(Int, options["machines"])),
            ("target_count", length(sextupoles)),
            ("quadrupole_count", length(quadrupoles)),
            ("corrector_count", length(local_controls)),
            ("bpm_count", length(detectors)),
        )
        for (key, expected) in checks
            get(metadata, key, nothing) == expected ||
                error("Existing paired latent metadata has incompatible $key")
        end
        return latent_root
    end

    mkpath(latent_root)
    write_npy(joinpath(latent_root, "sextupole_offsets.npy"), latents.sextupole_offsets)
    write_npy(
        joinpath(latent_root, "corrector_gain_errors.npy"),
        latents.corrector_gain_errors,
    )
    write_npy(joinpath(latent_root, "k2_gain_errors.npy"), latents.k2_gain_errors)
    write_npy(
        joinpath(latent_root, "quadrupole_relative_errors.npy"),
        latents.quadrupole_relative_errors,
    )
    write_npy(joinpath(latent_root, "quadrupole_rolls.npy"), latents.quadrupole_rolls)
    write_npy(
        joinpath(latent_root, "quadrupole_alignment_standard_normals.npy"),
        latents.quadrupole_alignment_standard_normals,
    )
    write_npy(joinpath(latent_root, "bpm_gain_errors.npy"), latents.bpm_gain_errors)
    write_npy(joinpath(latent_root, "drift_directions.npy"), latents.drift_directions)
    write_lines(
        joinpath(latent_root, "target_names.txt"),
        String.(getproperty.(sextupoles, :name)),
    )
    write_lines(
        joinpath(latent_root, "quadrupole_names.txt"),
        String.(getproperty.(quadrupoles, :name)),
    )
    write_lines(
        joinpath(latent_root, "bpm_names.txt"),
        String.(base_name.(detectors)),
    )
    write_metadata(metadata_path, Dict(
        "format" => "cesr-sequential-joint-paired-latents-v1",
        "date" => string(Dates.today()),
        "random_seed_base" => parse(Int, options["seed"]),
        "machine_count" => parse(Int, options["machines"]),
        "target_count" => length(sextupoles),
        "quadrupole_count" => length(quadrupoles),
        "corrector_count" => length(local_controls),
        "bpm_count" => length(detectors),
        "sextupole_offset_rms_m_per_plane" =>
            parse(Float64, options["sextupole-offset-rms-m"]),
        "corrector_gain_rms" => parse(Float64, options["corrector-gain-rms"]),
        "k2_gain_rms" => parse(Float64, options["k2-gain-rms"]),
        "bpm_gain_rms" => parse(Float64, options["bpm-gain-rms"]),
        "quadrupole_strength_halfwidth_fraction" =>
            parse(Float64, options["quadrupole-strength-fraction"]),
        "quadrupole_roll_rms_rad" =>
            parse(Float64, options["quadrupole-roll-rms-rad"]),
        "quadrupole_alignment_standard_normals" =>
            "unit Gaussian x/y draws; multiply by the case metadata RMS",
        "pairing" => "all latent arrays are shared by both physical cases",
    ))
    return latent_root
end

function generate_corrected_joint_case(
    options,
    model,
    sextupoles,
    quadrupoles,
    detectors,
    local_controls,
    geometry,
    knobs_by_target,
    nominal_v0,
    latents,
    baseline_gain_errors,
)
    ring = model.ring
    machine_count = parse(Int, options["machines"])
    target_limit = parse(Int, options["target-limit"])
    selected = target_limit == 0 ?
        sextupoles : sextupoles[1:min(target_limit, length(sextupoles))]
    target_names = String.(getproperty.(selected, :name))
    full_target_names = String.(getproperty.(sextupoles, :name))
    bpm_names = String.(base_name.(detectors))
    control_names = String.(model.metadata.steering_control_names)
    nt, nd, nc = length(selected), length(detectors), length(control_names)
    nl = length(local_controls)
    nc == 103 || error("Expected 103 latest-lattice steering controls, found $nc")
    nl == 62 || error("Expected 62 local-bump physical controls, found $nl")

    bump_amplitude = parse(Float64, options["bump-amplitude-m"])
    k2_amplitude = parse(Float64, options["k2-amplitude-m3"])
    drift_halfwidth = parse(Float64, options["drift-halfwidth-m"])
    quadrupole_alignment_rms =
        parse(Float64, options["quadrupole-alignment-rms-m"])
    response_step = parse(Float64, options["response-step"])
    ridge_ratio = parse(Float64, options["ridge-ratio"])
    relative_cutoff = parse(Float64, options["relative-svd-cutoff"])
    max_iterations = parse(Int, options["iterations"])
    tolerance_m = parse(Float64, options["tolerance-m"])
    max_update = parse(Float64, options["max-update"])
    line_search_steps = parse(Int, options["line-search-steps"])
    baseline_response_method = options["baseline-response-method"]
    baseline_response_method in ("finite_difference", "gtpsa") ||
        error("--baseline-response-method must be finite_difference or gtpsa")
    gtpsa_response_model = options["gtpsa-response-model"]
    gtpsa_response_model in ("realized", "nominal") ||
        error("--gtpsa-response-model must be realized or nominal")
    target_parallelism = options["target-parallelism"]
    target_parallelism in ("serial", "threads") ||
        error("--target-parallelism must be serial or threads")
    requested_scan_threads = parse(Int, options["scan-thread-count"])
    requested_scan_threads >= 0 || error("--scan-thread-count must be nonnegative")
    thread_equivalence_check = parse_boolean_option(
        "thread-equivalence-check", options["thread-equivalence-check"],
    )
    correction_bpm_noise_rms =
        parse(Float64, options["correction-bpm-noise-rms-m"])
    correction_measurement_repeats =
        parse(Int, options["correction-measurement-repeats"])
    correction_noise_seed = parse(Int, options["correction-noise-seed"])
    validate_gtpsa_with_finite_difference = parse_boolean_option(
        "validate-gtpsa-with-finite-difference",
        options["validate-gtpsa-with-finite-difference"],
    )
    correction_noise_std = measurement_mean_noise_std(
        correction_bpm_noise_rms, correction_measurement_repeats,
    )
    bumps = [
        (-bump_amplitude, 0.0),
        (0.0, -bump_amplitude),
        (0.0, 0.0),
        (0.0, bump_amplitude),
        (bump_amplitude, 0.0),
    ]
    delta_k2 = [-k2_amplitude, 0.0, k2_amplitude]
    nb, nk = length(bumps), length(delta_k2)

    bpm_orbits = fill(NaN, machine_count, nt, nb, nk, nd, 2)
    drift_bpm_orbits = similar(bpm_orbits)
    target_orbits = fill(NaN, machine_count, nt, nb, nk, 2)
    drift_target_orbits = similar(target_orbits)
    reference_bpm_orbits = fill(NaN, machine_count, nd, 2)
    reference_target_orbits = fill(NaN, machine_count, nt, 2)
    zero_offset_reference_bpm_orbits = similar(reference_bpm_orbits)
    zero_offset_reference_target_orbits = similar(reference_target_orbits)
    uncorrected_bpm_orbits = similar(reference_bpm_orbits)
    uncorrected_target_orbits = similar(reference_target_orbits)
    baseline_commands = fill(NaN, machine_count, nc)
    baseline_local_fields = fill(NaN, machine_count, nl)
    correction_histories = fill(NaN, machine_count, max_iterations + 1)
    physical_correction_histories = fill(NaN, machine_count, max_iterations + 1)
    response_singular_values = fill(NaN, machine_count, nc)
    reference_bpm_noise = zeros(machine_count, nd, 2)
    final_correction_bpm_noise = zeros(machine_count, nd, 2)
    validation_correction_bpm_noise = zeros(machine_count, nd, 2)
    validation_correction_residual_rms = fill(NaN, machine_count)
    response_seconds = fill(NaN, machine_count)
    response_closure_norm = fill(NaN, machine_count)
    response_fd_relative_l2 = fill(NaN, machine_count)
    response_fd_max_abs = fill(NaN, machine_count)
    scan_seconds = fill(NaN, machine_count, nt)
    correction_rows = NamedTuple[]
    thread_equivalence_bpm_max_abs_m = 0.0
    thread_equivalence_drift_bpm_max_abs_m = 0.0
    thread_equivalence_target_max_abs_m = 0.0
    thread_equivalence_drift_target_max_abs_m = 0.0
    thread_equivalence_checked = false

    scan_worker_count = if target_parallelism == "threads"
        requested = requested_scan_threads == 0 ? Threads.nthreads() :
            requested_scan_threads
        min(requested, Threads.nthreads(), nt)
    else
        1
    end
    target_parallelism == "threads" && scan_worker_count < 2 &&
        error("Threaded target scans require at least two Julia threads")
    target_parallelism == "threads" && BLAS.set_num_threads(1)
    scan_contexts = target_parallelism == "threads" ?
        [build_target_scan_context(options) for _ in 1:scan_worker_count] : Any[]
    for context in scan_contexts
        context.target_names == full_target_names ||
            error("Thread scan context sextupole inventory mismatch")
        context.bpm_names == bpm_names ||
            error("Thread scan context BPM inventory mismatch")
        context.control_names == control_names ||
            error("Thread scan context steering-control inventory mismatch")
    end

    nominal_gtpsa = if baseline_response_method == "gtpsa" &&
                       gtpsa_response_model == "nominal"
        validate_gtpsa_with_finite_difference && error(
            "Nominal GTPSA production mode does not reapply a latent-machine " *
            "central-difference ORM validation",
        )
        nominal_gtpsa_bpm_orm(control_names, bpm_names, nominal_v0)
    else
        nothing
    end

    indices = batch_element_indices(sextupoles, local_controls)
    batch_originals = Dict(
        index => ring.line[index].BMultipoleParams for index in indices
    )
    output_dir = joinpath(
        abspath(options["output-root"]), options["corrected-case-name"],
    )
    metadata_path = joinpath(output_dir, "scan_metadata.toml")
    isfile(metadata_path) && !parse_boolean_option("overwrite", options["overwrite"]) &&
        error("Output exists; use --overwrite=true: $metadata_path")
    started = time()

    try
        for machine in 1:machine_count
            machine_start = time()
            restore_batch_multipoles!(ring, batch_originals)
            restore_joint_machine!(
                ring, sextupoles, quadrupoles, NamedTuple[], geometry,
            )
            set_all_controls!(model.controls, 0.0)

            bpm_gain_errors = Matrix(latents.bpm_gain_errors[machine, :, :])
            physical_baseline_gains = 1 .+ baseline_gain_errors[machine, :]

            prepare_machine!(
                model,
                "without_quadrupole_misalignment",
                machine,
                sextupoles,
                quadrupoles,
                geometry,
                latents,
                quadrupole_alignment_rms,
            )
            zero_closed = scalar_closed_orbit_joint(ring, nominal_v0)
            zero_bpm, zero_target_full = read_machine_orbits(
                ring, zero_closed, bpm_names, full_target_names,
            )
            zero_measured = measured_bpm(zero_bpm, bpm_gain_errors)
            reference_rng = MersenneTwister(
                correction_noise_seed + 1_000_000 + machine,
            )
            noisy_zero_measured, zero_noise = noisy_measured_bpm(
                zero_bpm,
                bpm_gain_errors,
                reference_rng,
                correction_noise_std,
            )
            reference_bpm_noise[machine, :, :] .= zero_noise
            zero_offset_reference_bpm_orbits[machine, :, :] .= zero_bpm
            zero_offset_reference_target_orbits[machine, :, :] .=
                zero_target_full[1:nt, :]
            reference_response = if baseline_response_method == "gtpsa"
                gtpsa = if gtpsa_response_model == "nominal"
                    nominal_gtpsa
                else
                    measured_gtpsa_orm(
                        control_names,
                        physical_baseline_gains,
                        bpm_names,
                        bpm_gain_errors,
                        zero_closed.v0,
                        machine,
                        "without_quadrupole_misalignment",
                        latents,
                        quadrupole_alignment_rms,
                    )
                end
                response_seconds[machine] =
                    gtpsa_response_model == "nominal" ?
                        (machine == 1 ? gtpsa.seconds : 0.0) : gtpsa.seconds
                response_closure_norm[machine] = gtpsa.closure_norm_max
                if validate_gtpsa_with_finite_difference
                    finite_difference = finite_difference_measured_orm!(
                        model,
                        control_names,
                        physical_baseline_gains,
                        zeros(nc),
                        bpm_names,
                        bpm_gain_errors,
                        zero_closed,
                        response_step,
                    )
                    difference = gtpsa.response - finite_difference
                    response_fd_relative_l2[machine] =
                        norm(difference) / norm(finite_difference)
                    response_fd_max_abs[machine] = maximum(abs, difference)
                end
                gtpsa.response
            else
                response = nothing
                response_seconds[machine] = @elapsed response =
                    finite_difference_measured_orm!(
                        model,
                        control_names,
                        physical_baseline_gains,
                        zeros(nc),
                        bpm_names,
                        bpm_gain_errors,
                        zero_closed,
                        response_step,
                    )
                response_closure_norm[machine] = 0.0
                response
            end

            prepare_machine!(
                model,
                "with_quadrupole_misalignment",
                machine,
                sextupoles,
                quadrupoles,
                geometry,
                latents,
                quadrupole_alignment_rms,
            )
            uncorrected_closed = scalar_closed_orbit_joint(ring, nominal_v0)
            uncorrected_bpm, uncorrected_target_full = read_machine_orbits(
                ring, uncorrected_closed, bpm_names, full_target_names,
            )
            uncorrected_bpm_orbits[machine, :, :] .= uncorrected_bpm
            uncorrected_target_orbits[machine, :, :] .=
                uncorrected_target_full[1:nt, :]
            if correction_bpm_noise_rms > 0
                measurement_rng = MersenneTwister(
                    correction_noise_seed + 2_000_000 + machine,
                )
                validation_rng = MersenneTwister(
                    correction_noise_seed + 3_000_000 + machine,
                )
                result = solve_noisy_corrected_machine!(
                    model,
                    control_names,
                    physical_baseline_gains,
                    bpm_names,
                    full_target_names,
                    bpm_gain_errors,
                    noisy_zero_measured,
                    zero_measured,
                    reference_response,
                    uncorrected_closed,
                    measurement_rng,
                    validation_rng,
                    correction_noise_std;
                    ridge_ratio,
                    relative_cutoff,
                    max_iterations,
                    tolerance_m=max(tolerance_m, sqrt(2) * correction_noise_std),
                    max_update,
                    line_search_steps,
                )
                physical_correction_histories[machine, :] .= result.physical_history
                final_correction_bpm_noise[machine, :, :] .=
                    result.final_measurement_noise
                validation_correction_bpm_noise[machine, :, :] .=
                    result.validation_noise
                validation_correction_residual_rms[machine] =
                    coordinate_rms(result.validation_residual)
            else
                result = solve_corrected_machine!(
                    model,
                    control_names,
                    physical_baseline_gains,
                    bpm_names,
                    full_target_names,
                    bpm_gain_errors,
                    zero_measured,
                    reference_response,
                    uncorrected_closed;
                    ridge_ratio,
                    relative_cutoff,
                    max_iterations,
                    tolerance_m,
                    max_update,
                    line_search_steps,
                )
                physical_correction_histories[machine, :] .= result.history
                validation_correction_residual_rms[machine] =
                    coordinate_rms(result.residual)
            end
            baseline_commands[machine, :] .= result.commands
            correction_histories[machine, :] .= result.history
            response_singular_values[machine, :] .= result.singular
            corrected_local_baseline =
                physical_bump_control_baselines(ring, local_controls)
            baseline_local_fields[machine, :] .= corrected_local_baseline

            materialize_batch_elements!(ring, sextupoles, local_controls)
            corrected_closed = scalar_closed_orbit_joint(ring, result.closed.v0)
            corrected_bpm, corrected_target_full = read_machine_orbits(
                ring, corrected_closed, bpm_names, full_target_names,
            )
            reference_bpm_orbits[machine, :, :] .= corrected_bpm
            reference_target_orbits[machine, :, :] .= corrected_target_full[1:nt, :]
            materialization_error = maximum(abs, corrected_bpm - result.physical_bpm)
            materialization_error <= 1e-10 ||
                error("Materializing corrected fields changed BPM orbit by $materialization_error m")

            before_measured = measured_bpm(uncorrected_bpm, bpm_gain_errors) - zero_measured
            after_measured = measured_bpm(corrected_bpm, bpm_gain_errors) - zero_measured
            before_target = uncorrected_target_full - zero_target_full
            after_target = corrected_target_full - zero_target_full
            push!(correction_rows, (;
                machine,
                before_bpm_rms_um=coordinate_rms(before_measured) * 1e6,
                after_bpm_rms_um=coordinate_rms(after_measured) * 1e6,
                before_target_2d_rms_um=two_plane_rms(before_target) * 1e6,
                after_target_2d_rms_um=two_plane_rms(after_target) * 1e6,
                command_rms=coordinate_rms(result.commands),
                max_abs_command=maximum(abs, result.commands),
                iterations=result.iterations,
                retained_response_rank=result.retained_rank,
                retained_response_condition=result.condition,
                response_method=baseline_response_method,
                response_seconds=response_seconds[machine],
                gtpsa_closure_norm_max=response_closure_norm[machine],
                gtpsa_vs_finite_difference_relative_l2=
                    response_fd_relative_l2[machine],
                gtpsa_vs_finite_difference_max_abs=
                    response_fd_max_abs[machine],
                correction_measurement_residual_rms_um=
                    last(filter(isfinite, result.history)) * 1e6,
                correction_independent_validation_residual_rms_um=
                    validation_correction_residual_rms[machine] * 1e6,
                materialization_bpm_max_abs_error_m=materialization_error,
            ))

            function scan_one_target!(
                scan_ring,
                scan_sextupoles,
                scan_controls,
                scan_knobs_by_target,
                scan_closed,
                scan_local_baseline,
                target_index,
            )
                target_entry = scan_sextupoles[target_index]
                target_start = time()
                target_knobs = scan_knobs_by_target[target_entry.name]
                try
                    static_bpm, drift_bpm, static_target, drift_target =
                        scan_joint_target!(
                            scan_ring,
                            target_index,
                            target_entry,
                            target_knobs,
                            scan_controls,
                            bpm_names,
                            bumps,
                            delta_k2,
                            scan_closed,
                            latents.corrector_gain_errors[machine, :],
                            latents.k2_gain_errors[machine, target_index],
                            latents.drift_directions[machine, target_index, :],
                            drift_halfwidth;
                            baseline_corrector_fields=scan_local_baseline,
                        )
                    bpm_orbits[machine, target_index, :, :, :, :] .= static_bpm
                    drift_bpm_orbits[machine, target_index, :, :, :, :] .= drift_bpm
                    target_orbits[machine, target_index, :, :, :] .= static_target
                    drift_target_orbits[machine, target_index, :, :, :] .= drift_target
                finally
                    for (control_index, control) in enumerate(scan_controls)
                        set_corrector_scalar!(
                            scan_ring, control, scan_local_baseline[control_index],
                        )
                    end
                    scan_ring.line[target_entry.index].Kn2 = target_entry.kn2_m3
                end
                scan_seconds[machine, target_index] = time() - target_start
                return nothing
            end

            if target_parallelism == "threads"
                worker_closed = Vector{Any}(undef, scan_worker_count)
                worker_local_baseline = Vector{Vector{Float64}}(
                    undef, scan_worker_count,
                )
                Threads.@threads :static for worker in 1:scan_worker_count
                    worker_closed[worker], worker_local_baseline[worker] =
                        prepare_target_scan_context!(
                            scan_contexts[worker],
                            machine,
                            latents,
                            quadrupole_alignment_rms,
                            physical_baseline_gains,
                            result.commands,
                            corrected_closed.v0,
                        )
                end
                worker_bpm, _ = read_machine_orbits(
                    scan_contexts[1].ring,
                    worker_closed[1],
                    bpm_names,
                    full_target_names,
                )
                worker_reference_error = maximum(abs, worker_bpm - corrected_bpm)
                worker_reference_error <= 1e-10 || error(
                    "Thread scan context changed corrected BPM orbit by " *
                    "$worker_reference_error m",
                )
                progress_lock = ReentrantLock()
                completed = Threads.Atomic{Int}(0)
                Threads.@threads :static for worker in 1:scan_worker_count
                    context = scan_contexts[worker]
                    for target_index in worker:scan_worker_count:nt
                        scan_one_target!(
                            context.ring,
                            context.sextupoles,
                            context.local_controls,
                            context.knobs_by_target,
                            worker_closed[worker],
                            worker_local_baseline[worker],
                            target_index,
                        )
                        finished = Threads.atomic_add!(completed, 1) + 1
                        if finished % 10 == 0 || finished == nt
                            lock(progress_lock) do
                                @printf(
                                    "%s machine %d/%d targets %d/%d complete with %d workers (%.1f s machine elapsed)\n",
                                    options["corrected-case-name"],
                                    machine,
                                    machine_count,
                                    finished,
                                    nt,
                                    scan_worker_count,
                                    time() - machine_start,
                                )
                                flush(stdout)
                            end
                        end
                    end
                end
                if thread_equivalence_check && machine == 1
                    parallel_bpm = copy(bpm_orbits[machine, 1, :, :, :, :])
                    parallel_drift_bpm = copy(
                        drift_bpm_orbits[machine, 1, :, :, :, :],
                    )
                    parallel_target = copy(target_orbits[machine, 1, :, :, :])
                    parallel_drift_target = copy(
                        drift_target_orbits[machine, 1, :, :, :],
                    )
                    scan_one_target!(
                        ring,
                        sextupoles,
                        local_controls,
                        knobs_by_target,
                        corrected_closed,
                        corrected_local_baseline,
                        1,
                    )
                    thread_equivalence_bpm_max_abs_m = maximum(
                        abs, bpm_orbits[machine, 1, :, :, :, :] - parallel_bpm,
                    )
                    thread_equivalence_drift_bpm_max_abs_m = maximum(
                        abs,
                        drift_bpm_orbits[machine, 1, :, :, :, :] -
                        parallel_drift_bpm,
                    )
                    thread_equivalence_target_max_abs_m = maximum(
                        abs, target_orbits[machine, 1, :, :, :] - parallel_target,
                    )
                    thread_equivalence_drift_target_max_abs_m = maximum(
                        abs,
                        drift_target_orbits[machine, 1, :, :, :] -
                        parallel_drift_target,
                    )
                    maximum((
                        thread_equivalence_bpm_max_abs_m,
                        thread_equivalence_drift_bpm_max_abs_m,
                        thread_equivalence_target_max_abs_m,
                        thread_equivalence_drift_target_max_abs_m,
                    )) <= 1e-13 || error(
                        "Threaded and serial target scans are not equivalent",
                    )
                    thread_equivalence_checked = true
                end
            else
                for target_index in eachindex(selected)
                    scan_one_target!(
                        ring,
                        sextupoles,
                        local_controls,
                        knobs_by_target,
                        corrected_closed,
                        corrected_local_baseline,
                        target_index,
                    )
                    if target_index % 10 == 0 || target_index == nt
                        @printf(
                            "%s machine %d/%d target %d/%d complete (%.1f s machine elapsed)\n",
                            options["corrected-case-name"],
                            machine,
                            machine_count,
                            target_index,
                            nt,
                            time() - machine_start,
                        )
                        flush(stdout)
                    end
                end
            end
            @printf(
                "%s machine %d/%d complete in %.1f s (elapsed %.1f s)\n",
                options["corrected-case-name"],
                machine,
                machine_count,
                time() - machine_start,
                time() - started,
            )
            flush(stdout)
        end
    finally
        restore_batch_multipoles!(ring, batch_originals)
        restore_joint_machine!(ring, sextupoles, quadrupoles, NamedTuple[], geometry)
        set_all_controls!(model.controls, 0.0)
    end

    mkpath(output_dir)
    write_npy(joinpath(output_dir, "bpm_orbits.npy"), bpm_orbits)
    write_npy(joinpath(output_dir, "drift_bpm_orbits.npy"), drift_bpm_orbits)
    bpm_gain_factors = 1 .+ latents.bpm_gain_errors
    observable_gain_shape = reshape(
        bpm_gain_factors, machine_count, 1, 1, 1, nd, 2,
    )
    write_npy(
        joinpath(output_dir, "observable_bpm_readbacks.npy"),
        bpm_orbits .* observable_gain_shape,
    )
    write_npy(
        joinpath(output_dir, "observable_drift_bpm_readbacks.npy"),
        drift_bpm_orbits .* observable_gain_shape,
    )
    write_npy(joinpath(output_dir, "target_orbits.npy"), target_orbits)
    write_npy(joinpath(output_dir, "drift_target_orbits.npy"), drift_target_orbits)
    write_npy(joinpath(output_dir, "reference_bpm_orbits.npy"), reference_bpm_orbits)
    write_npy(
        joinpath(output_dir, "reference_target_orbits.npy"),
        reference_target_orbits,
    )
    write_npy(
        joinpath(output_dir, "zero_offset_reference_bpm_orbits.npy"),
        zero_offset_reference_bpm_orbits,
    )
    write_npy(
        joinpath(output_dir, "zero_offset_reference_target_orbits.npy"),
        zero_offset_reference_target_orbits,
    )
    write_npy(
        joinpath(output_dir, "uncorrected_bpm_orbits.npy"), uncorrected_bpm_orbits,
    )
    write_npy(
        joinpath(output_dir, "uncorrected_target_orbits.npy"),
        uncorrected_target_orbits,
    )
    write_npy(
        joinpath(output_dir, "baseline_corrector_commands.npy"), baseline_commands,
    )
    write_npy(
        joinpath(output_dir, "baseline_local_corrector_fields.npy"),
        baseline_local_fields,
    )
    write_npy(
        joinpath(output_dir, "baseline_corrector_gain_errors.npy"),
        baseline_gain_errors,
    )
    write_npy(
        joinpath(output_dir, "baseline_correction_history_m.npy"),
        correction_histories,
    )
    write_npy(
        joinpath(output_dir, "baseline_physical_correction_history_m.npy"),
        physical_correction_histories,
    )
    write_npy(
        joinpath(output_dir, "baseline_response_singular_values.npy"),
        response_singular_values,
    )
    write_npy(
        joinpath(output_dir, "baseline_reference_bpm_noise_m.npy"),
        reference_bpm_noise,
    )
    write_npy(
        joinpath(output_dir, "baseline_final_correction_bpm_noise_m.npy"),
        final_correction_bpm_noise,
    )
    write_npy(
        joinpath(output_dir, "baseline_validation_bpm_noise_m.npy"),
        validation_correction_bpm_noise,
    )
    write_npy(
        joinpath(output_dir, "baseline_validation_residual_rms_m.npy"),
        validation_correction_residual_rms,
    )
    write_npy(
        joinpath(output_dir, "baseline_response_seconds.npy"), response_seconds,
    )
    write_npy(
        joinpath(output_dir, "baseline_response_closure_norm.npy"),
        response_closure_norm,
    )
    write_npy(
        joinpath(output_dir, "baseline_gtpsa_vs_fd_relative_l2.npy"),
        response_fd_relative_l2,
    )
    write_npy(
        joinpath(output_dir, "baseline_gtpsa_vs_fd_max_abs.npy"),
        response_fd_max_abs,
    )
    write_npy(joinpath(output_dir, "scan_seconds.npy"), scan_seconds)
    if !isnothing(nominal_gtpsa)
        write_npy(
            joinpath(output_dir, "baseline_nominal_gtpsa_orm.npy"),
            nominal_gtpsa.response,
        )
    end
    write_lines(joinpath(output_dir, "target_names.txt"), target_names)
    write_lines(joinpath(output_dir, "bpm_names.txt"), bpm_names)
    write_lines(joinpath(output_dir, "baseline_corrector_names.txt"), control_names)
    write_lines(
        joinpath(output_dir, "local_bump_corrector_names.txt"),
        ["$(control.name):$(String(control.axis))" for control in local_controls],
    )
    write_rows(joinpath(output_dir, "baseline_correction.csv"), correction_rows)
    write_rows(joinpath(output_dir, "bump_points.csv"), [
        (; bump_index=index, bump_x_command_m=x, bump_y_command_m=y)
        for (index, (x, y)) in enumerate(bumps)
    ])

    before_bpm = coordinate_rms(
        uncorrected_bpm_orbits .* bpm_gain_factors -
        zero_offset_reference_bpm_orbits .* bpm_gain_factors,
    ) * 1e6
    after_bpm = coordinate_rms(
        reference_bpm_orbits .* bpm_gain_factors -
        zero_offset_reference_bpm_orbits .* bpm_gain_factors,
    ) * 1e6
    before_target = two_plane_rms(reshape(
        uncorrected_target_orbits - zero_offset_reference_target_orbits, :, 2,
    )) * 1e6
    after_target = two_plane_rms(reshape(
        reference_target_orbits - zero_offset_reference_target_orbits, :, 2,
    )) * 1e6
    wall_seconds = time() - started
    scan_format = baseline_response_method == "gtpsa" || correction_bpm_noise_rms > 0 ?
        "cesr-sequential-joint-machine-scan-v3" :
        "cesr-sequential-joint-machine-scan-v2"
    response_engine = baseline_response_method == "gtpsa" ?
        "SciBmad/GTPSA implicit periodic closed-orbit Jacobian" :
        "SciBmad exact central finite-difference closed-orbit response"
    write_metadata(metadata_path, Dict(
        "format" => scan_format,
        "date" => string(Dates.today()),
        "engine" => "SciBmad exact RF-on correction, $response_engine, and BatchParam closed-orbit scans",
        "lattice" => LATEST_LATTICE,
        "case" => options["corrected-case-name"],
        "machine_count" => machine_count,
        "target_count" => nt,
        "full_sextupole_count" => length(sextupoles),
        "state_count_per_target" => nb * nk,
        "batch_lane_count_per_target" => 2nb * nk,
        "total_exact_states" => 2machine_count * nt * nb * nk,
        "machine_atomic_unit" => "all latent errors and baseline correction fixed across all target scans",
        "paired_case_latents" => "same draws as zero-offset and uncorrected cases; only quadrupole alignment and recorded baseline correction differ",
        "baseline_orbit_correction_applied" => true,
        "baseline_reference_semantics" => "paired BPM closed orbit with quadrupole alignment offsets disabled; not BPM zero",
        "baseline_solver_inputs" => "reference/current BPM readbacks and zero-offset response only; latent offsets and target orbit excluded",
        "baseline_response_method" => baseline_response_method == "gtpsa" ?
            "reference_gtpsa_orm" : "reference_finite_difference_orm",
        "baseline_gtpsa_response_model" => gtpsa_response_model,
        "baseline_response_engine" => response_engine,
        "baseline_gtpsa_validation_enabled" =>
            validate_gtpsa_with_finite_difference,
        "baseline_gtpsa_vs_finite_difference_relative_l2_max" =>
            any(isfinite, response_fd_relative_l2) ?
                maximum(filter(isfinite, response_fd_relative_l2)) : NaN,
        "baseline_gtpsa_vs_finite_difference_max_abs" =>
            any(isfinite, response_fd_max_abs) ?
                maximum(filter(isfinite, response_fd_max_abs)) : NaN,
        "baseline_response_closure_norm_max" =>
            maximum(response_closure_norm),
        "baseline_response_seconds_total" => sum(response_seconds),
        "baseline_corrector_count" => nc,
        "baseline_corrector_registry" => "latest-lattice H/V steering Overlay controls",
        "baseline_corrector_units" => "HKICK/VKICK command radians; no CESR hardware limit asserted",
        "baseline_model_calibration_semantics" =>
            gtpsa_response_model == "nominal" ?
                "nominal theoretical GTPSA ORM; all BPM/corrector gains and static magnet/alignment errors remain unknown to the correction response model" :
                "response rows/columns use the realized BPM/corrector gains; this run has no unknown gain mismatch between the GTPSA/finite-difference ORM and the simulated machine",
        "baseline_vs_local_corrector_gain_semantics" =>
            "the 103 baseline-control gains and 62 local-bump-control gains are separate deterministic draws from the same prior; shared physical-device calibration is not yet unified",
        "baseline_response_step" => response_step,
        "baseline_ridge_ratio" => ridge_ratio,
        "baseline_relative_svd_cutoff" => relative_cutoff,
        "baseline_maximum_iterations" => max_iterations,
        "baseline_bpm_tolerance_m" => tolerance_m,
        "baseline_effective_bpm_tolerance_m" =>
            max(tolerance_m, sqrt(2) * correction_noise_std),
        "baseline_maximum_update" => max_update,
        "baseline_line_search_steps" => line_search_steps,
        "baseline_bpm_noise_rms_m_per_read" => correction_bpm_noise_rms,
        "baseline_measurement_repeats" => correction_measurement_repeats,
        "baseline_bpm_mean_noise_std_m" => correction_noise_std,
        "baseline_reference_current_noise_semantics" =>
            "same independent Gaussian per-read model and shared static BPM gains; distinct reference, iteration, and validation draws",
        "baseline_noise_seed" => correction_noise_seed,
        "baseline_independent_validation_residual_rms_um" =>
            coordinate_rms(validation_correction_residual_rms) * 1e6,
        "baseline_bpm_rms_before_um" => before_bpm,
        "baseline_bpm_rms_after_um" => after_bpm,
        "baseline_target_2d_rms_before_um" => before_target,
        "baseline_target_2d_rms_after_um" => after_target,
        "baseline_hold_semantics" => "fixed during every target scan; local bump commands are additive",
        "observable_bpm_readback_semantics" =>
            "forward-generated physical BPM coordinates multiplied by the realized fixed BPM gain factors; the inverse consumes these saved readbacks and receives no gain realization",
        "target_scan_parallelism" => target_parallelism,
        "target_scan_worker_count" => scan_worker_count,
        "julia_thread_count" => Threads.nthreads(),
        "blas_thread_count_during_scan" => BLAS.get_num_threads(),
        "thread_safety_semantics" => "each worker owns an independently loaded mutable latest-lattice ring/model; latent draws and fixed baseline commands are shared only as immutable numeric inputs",
        "thread_equivalence_checked" => thread_equivalence_checked,
        "thread_equivalence_bpm_max_abs_m" =>
            thread_equivalence_bpm_max_abs_m,
        "thread_equivalence_drift_bpm_max_abs_m" =>
            thread_equivalence_drift_bpm_max_abs_m,
        "thread_equivalence_target_max_abs_m" =>
            thread_equivalence_target_max_abs_m,
        "thread_equivalence_drift_target_max_abs_m" =>
            thread_equivalence_drift_target_max_abs_m,
        "local_bump_corrector_count" => nl,
        "sextupole_offset_distribution" => "independent Gaussian x/y, fixed across the machine realization",
        "sextupole_offset_rms_m" =>
            parse(Float64, options["sextupole-offset-rms-m"]),
        "corrector_gain_rms" => parse(Float64, options["corrector-gain-rms"]),
        "k2_gain_rms" => parse(Float64, options["k2-gain-rms"]),
        "quadrupole_strength_distribution" => "independent uniform +/-fraction, fixed across all target scans",
        "quadrupole_strength_fraction" =>
            parse(Float64, options["quadrupole-strength-fraction"]),
        "quadrupole_roll_rms_rad" =>
            parse(Float64, options["quadrupole-roll-rms-rad"]),
        "quadrupole_alignment_distribution" => "independent Gaussian x/y per physical quadrupole, coherent across slices",
        "quadrupole_alignment_rms_m_per_plane" => quadrupole_alignment_rms,
        "quadrupole_alignment_semantics" => "annual residual drift relative to measured nominal; fixed during correction and all scans",
        "bpm_gain_rms" => parse(Float64, options["bpm-gain-rms"]),
        "bump_amplitude_m" => bump_amplitude,
        "bump_count" => nb,
        "k2_delta_m3" => delta_k2,
        "k2_count" => nk,
        "drift_halfwidth_m" => drift_halfwidth,
        "drift_secant" => "same target bump knobs, linear state-order coefficient from -1 to +1",
        "bpm_count" => nd,
        "quadrupole_count" => length(quadrupoles),
        "calculation_wall_seconds" => wall_seconds,
    ))
    @printf(
        "Corrected scans complete: BPM %.3f -> %.3f um, target %.3f -> %.3f um 2D\n",
        before_bpm,
        after_bpm,
        before_target,
        after_target,
    )
    println("Wrote corrected joint-machine scans to $output_dir")
    return output_dir
end

function main_corrected_scans(args=ARGS)
    defaults = Dict(
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
        "baseline-response-method" => "finite_difference",
        "gtpsa-response-model" => "realized",
        "target-parallelism" => "serial",
        "scan-thread-count" => "0",
        "thread-equivalence-check" => "true",
        "response-step" => "1.0e-6",
        "validate-gtpsa-with-finite-difference" => "false",
        "correction-bpm-noise-rms-m" => "0.0",
        "correction-measurement-repeats" => "1",
        "correction-noise-seed" => "20261123",
        "ridge-ratio" => "1.0e-2",
        "relative-svd-cutoff" => "1.0e-8",
        "iterations" => "3",
        "tolerance-m" => "1.0e-7",
        "max-update" => "0.0",
        "line-search-steps" => "8",
        "corrected-case-name" => CORRECTED_JOINT_CASE,
        "include-reference-case" => "false",
        "output-root" => CORRECTED_SCAN_ROOT,
        "overwrite" => "false",
    )
    options = parse_exact11_options(defaults, args)
    machine_count = parse(Int, options["machines"])
    target_limit = parse(Int, options["target-limit"])
    machine_count > 0 || error("--machines must be positive")
    target_limit >= 0 || error("--target-limit must be nonnegative")
    parse_boolean_option("include-reference-case", options["include-reference-case"])
    parse_boolean_option("overwrite", options["overwrite"])
    options["baseline-response-method"] in ("finite_difference", "gtpsa") ||
        error("--baseline-response-method must be finite_difference or gtpsa")
    options["gtpsa-response-model"] in ("realized", "nominal") ||
        error("--gtpsa-response-model must be realized or nominal")
    options["target-parallelism"] in ("serial", "threads") ||
        error("--target-parallelism must be serial or threads")
    parse(Int, options["scan-thread-count"]) >= 0 ||
        error("--scan-thread-count must be nonnegative")
    parse_boolean_option(
        "thread-equivalence-check", options["thread-equivalence-check"],
    )
    parse_boolean_option(
        "validate-gtpsa-with-finite-difference",
        options["validate-gtpsa-with-finite-difference"],
    )
    parse(Float64, options["correction-bpm-noise-rms-m"]) >= 0 ||
        error("--correction-bpm-noise-rms-m must be nonnegative")
    parse(Int, options["correction-measurement-repeats"]) > 0 ||
        error("--correction-measurement-repeats must be positive")
    isempty(strip(options["corrected-case-name"])) &&
        error("--corrected-case-name must be nonempty")

    model = load_ring_model(; ring=:latest, zero_value=0.0, rf_on=true)
    ring = model.ring
    sextupoles = active_sextupole_inventory(ring)
    quadrupoles = active_quadrupole_inventory(ring)
    detectors = measurable_bpms(ring)
    local_controls = independent_corrector_inventory(ring)
    geometry = quadrupole_geometry_joint(ring, quadrupoles)
    target_names = String.(getproperty.(sextupoles, :name))
    knobs_by_target = read_joint_knobs(
        options["bump-knobs-csv"], target_names, local_controls,
    )
    nominal_v0 = copy(solve_closed_orbit(ring).v0)
    latents = sample_joint_latents(
        machine_count,
        length(sextupoles),
        length(quadrupoles),
        length(local_controls),
        length(detectors),
        parse(Int, options["seed"]),
        parse(Float64, options["sextupole-offset-rms-m"]),
        parse(Float64, options["corrector-gain-rms"]),
        parse(Float64, options["k2-gain-rms"]),
        parse(Float64, options["bpm-gain-rms"]),
        parse(Float64, options["quadrupole-strength-fraction"]),
        parse(Float64, options["quadrupole-roll-rms-rad"]),
    )
    baseline_gain_errors = sample_baseline_corrector_gain_errors(
        machine_count,
        length(model.metadata.steering_control_names),
        parse(Int, options["seed"]),
        parse(Float64, options["corrector-gain-rms"]),
    )
    write_or_validate_joint_latents!(
        options["output-root"],
        options,
        latents,
        sextupoles,
        quadrupoles,
        detectors,
        local_controls,
    )

    generate_corrected_joint_case(
        options,
        model,
        sextupoles,
        quadrupoles,
        detectors,
        local_controls,
        geometry,
        knobs_by_target,
        nominal_v0,
        latents,
        baseline_gain_errors,
    )
    if parse_boolean_option("include-reference-case", options["include-reference-case"])
        reference_model = load_ring_model(; ring=:latest, zero_value=0.0, rf_on=true)
        reference_ring = reference_model.ring
        reference_sextupoles = active_sextupole_inventory(reference_ring)
        reference_quadrupoles = active_quadrupole_inventory(reference_ring)
        reference_detectors = measurable_bpms(reference_ring)
        reference_controls = independent_corrector_inventory(reference_ring)
        reference_geometry = quadrupole_geometry_joint(
            reference_ring, reference_quadrupoles,
        )
        reference_knobs = read_joint_knobs(
            options["bump-knobs-csv"], target_names, reference_controls,
        )
        reference_v0 = copy(solve_closed_orbit(reference_ring).v0)
        generate_joint_case(
            options,
            "without_quadrupole_misalignment",
            reference_ring,
            reference_sextupoles,
            reference_quadrupoles,
            reference_detectors,
            reference_controls,
            reference_geometry,
            reference_knobs,
            reference_v0,
            latents,
        )
    end
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main_corrected_scans())
end

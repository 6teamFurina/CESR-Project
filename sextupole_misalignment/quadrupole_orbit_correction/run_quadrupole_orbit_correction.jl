#!/usr/bin/env julia

"""Restore the BPM orbit after a paired 50-um/plane quadrupole alignment drift.

For each latent machine, all maintained static nuisance draws are identical in
the reference and drifted cases.  The only switched quantity is the x/y
alignment offset of every active quadrupole.  Correction consumes only the
measured reference BPM vector, the measured current BPM vector, and a
corrector-to-BPM response matrix.  Latent quadrupole offsets and target-local
orbits are retained for scoring only.
"""

include(joinpath(
    @__DIR__, "..", "sequential_joint_inverse", "generate_joint_machine_scans.jl",
))
include(joinpath(
    @__DIR__, "..", "..", "orbit", "Orbit_Calculation", "ring_model_adapter.jl",
))

const ORBIT_CORRECTION_METHODS = ("reference_orm", "current_orm")

function parse_method_names(text)
    methods = strip.(split(text, ','))
    isempty(methods) && error("At least one response method is required")
    all(method in ORBIT_CORRECTION_METHODS for method in methods) ||
        error("Unknown response method; allowed methods are $(join(ORBIT_CORRECTION_METHODS, ", "))")
    length(unique(methods)) == length(methods) || error("Response methods must be unique")
    return methods
end

function set_commanded_correctors!(model, control_names, physical_gains, commands)
    length(control_names) == length(physical_gains) == length(commands) ||
        error("Corrector command vectors have inconsistent lengths")
    for index in eachindex(control_names)
        model.controls[control_names[index]] = physical_gains[index] * commands[index]
    end
    return nothing
end

function flatten_bpm(values)
    size(values, 2) == 2 || error("BPM orbit must have x/y columns")
    return vcat(values[:, 1], values[:, 2])
end

function measured_bpm(physical, gain_errors)
    size(physical) == size(gain_errors) || error("BPM gain array has the wrong shape")
    return physical .* (1 .+ gain_errors)
end

coordinate_rms(values) = sqrt(mean(abs2, values))
two_plane_rms(values) = sqrt(mean(sum(abs2, values; dims=2)))

function read_machine_orbits(ring, closed, bpm_names, target_names)
    return track_joint_machine_reference(ring, closed, bpm_names, target_names)
end

function finite_difference_measured_orm!(
    model,
    control_names,
    physical_corrector_gains,
    commands,
    bpm_names,
    bpm_gain_errors,
    closed,
    step,
)
    ring = model.ring
    nc = length(control_names)
    nd = length(bpm_names)
    state_count = 2nc
    response = zeros(2nd, nc)
    try
        for control_index in eachindex(control_names)
            values = fill(
                physical_corrector_gains[control_index] * commands[control_index],
                state_count,
            )
            values[2control_index - 1] += physical_corrector_gains[control_index] * step
            values[2control_index] -= physical_corrector_gains[control_index] * step
            model.controls[control_names[control_index]] = BatchParam(values)
        end
        initial = repeat(reshape(Array(closed.v0), 1, 6), state_count, 1)
        batch_closed = solve_batch_closed_orbit(ring, state_count; initial_v0=initial)
        tracked = track_orbits_at_names(ring, batch_closed, bpm_names)
        for (bpm_index, name) in enumerate(bpm_names)
            x = tracked.horizontal[name]
            y = tracked.vertical[name]
            for control_index in eachindex(control_names)
                plus = 2control_index - 1
                minus = 2control_index
                response[bpm_index, control_index] =
                    (1 + bpm_gain_errors[bpm_index, 1]) *
                    (x[plus] - x[minus]) / (2step)
                response[nd + bpm_index, control_index] =
                    (1 + bpm_gain_errors[bpm_index, 2]) *
                    (y[plus] - y[minus]) / (2step)
            end
        end
    finally
        set_commanded_correctors!(
            model, control_names, physical_corrector_gains, commands,
        )
    end
    all(isfinite, response) || error("Non-finite measured orbit-response matrix")
    return response
end

function response_diagnostics(response, relative_cutoff)
    factorization = svd(response; full=false)
    singular = factorization.S
    threshold = relative_cutoff * first(singular)
    rank = count(value -> value >= threshold, singular)
    rank > 0 || error("The orbit-response matrix has zero retained rank")
    condition = first(singular) / singular[rank]
    return factorization, singular, rank, condition
end

function ridge_command(response_factorization, residual, ridge_ratio, relative_cutoff)
    singular = response_factorization.S
    threshold = relative_cutoff * first(singular)
    ridge = ridge_ratio * first(singular)
    factors = [
        value >= threshold ? value / (value^2 + ridge^2) : 0.0
        for value in singular
    ]
    return -response_factorization.V *
           (factors .* (response_factorization.U' * residual))
end

function solve_corrected_machine!(
    model,
    control_names,
    physical_corrector_gains,
    bpm_names,
    target_names,
    bpm_gain_errors,
    reference_measured,
    response,
    initial_closed;
    ridge_ratio,
    relative_cutoff,
    max_iterations,
    tolerance_m,
    max_update,
    line_search_steps,
)
    ring = model.ring
    factorization, singular, retained_rank, condition =
        response_diagnostics(response, relative_cutoff)
    commands = zeros(length(control_names))
    closed = initial_closed
    physical_bpm, physical_target =
        read_machine_orbits(ring, closed, bpm_names, target_names)
    measured = measured_bpm(physical_bpm, bpm_gain_errors)
    residual = flatten_bpm(measured - reference_measured)
    history = fill(NaN, max_iterations + 1)
    history[1] = coordinate_rms(residual)
    accepted_iterations = 0

    for iteration in 1:max_iterations
        history[iteration] <= tolerance_m && break
        update = ridge_command(
            factorization, residual, ridge_ratio, relative_cutoff,
        )
        if max_update > 0 && maximum(abs, update) > max_update
            update .*= max_update / maximum(abs, update)
        end
        update_norm = maximum(abs, update)
        update_norm > eps(Float64) || break

        accepted = false
        alpha = 1.0
        for _ in 0:line_search_steps
            trial_commands = commands .+ alpha .* update
            set_commanded_correctors!(
                model,
                control_names,
                physical_corrector_gains,
                trial_commands,
            )
            trial_closed = scalar_closed_orbit_joint(ring, closed.v0)
            trial_bpm, trial_target =
                read_machine_orbits(ring, trial_closed, bpm_names, target_names)
            trial_measured = measured_bpm(trial_bpm, bpm_gain_errors)
            trial_residual = flatten_bpm(trial_measured - reference_measured)
            trial_rms = coordinate_rms(trial_residual)
            if trial_rms < history[iteration]
                commands = trial_commands
                closed = trial_closed
                physical_bpm = trial_bpm
                physical_target = trial_target
                measured = trial_measured
                residual = trial_residual
                history[iteration + 1] = trial_rms
                accepted_iterations = iteration
                accepted = true
                break
            end
            alpha /= 2
        end
        if !accepted
            set_commanded_correctors!(
                model, control_names, physical_corrector_gains, commands,
            )
            break
        end
    end

    final_rms = coordinate_rms(residual)
    converged = final_rms <= tolerance_m
    return (;
        commands,
        closed,
        physical_bpm,
        physical_target,
        measured_bpm=measured,
        residual,
        history,
        iterations=accepted_iterations,
        converged,
        singular,
        retained_rank,
        condition,
    )
end

function prepare_machine!(
    model,
    case_name,
    machine,
    sextupoles,
    quadrupoles,
    geometry,
    latents,
    quadrupole_alignment_rms,
)
    ring = model.ring
    set_all_controls!(model.controls, 0.0)
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
    apply_joint_machine!(
        ring,
        case_name,
        machine,
        sextupoles,
        quadrupoles,
        geometry,
        latents,
        quadrupole_alignment_rms,
    )
    return nothing
end

function markdown_summary(path, aggregate_rows, response_comparison_rows, metadata)
    open(path, "w") do io
        println(io, "# Paired quadrupole-offset orbit correction")
        println(io)
        println(io, "This latest-lattice SciBmad experiment keeps every maintained static")
        println(io, "random error fixed within each paired machine and switches only the")
        println(io, "quadrupole x/y alignment offsets from zero to 50 micrometers RMS per plane.")
        println(io, "Correction uses BPM observations, corrector commands, and an orbit-response")
        println(io, "matrix; latent quadrupole offsets are never supplied to the solver.")
        println(io)
        println(io, "## Aggregate result")
        println(io)
        println(io, "| response matrix | before BPM RMS [um] | after BPM RMS [um] | reduction | before target 2D RMS [um] | after target 2D RMS [um] | command RMS | max abs command |")
        println(io, "|---|---:|---:|---:|---:|---:|---:|---:|")
        for row in aggregate_rows
            println(
                io,
                "| $(row.method) | $(round(row.before_bpm_rms_um; digits=6)) | " *
                "$(round(row.after_bpm_rms_um; digits=6)) | " *
                "$(round(row.reduction_factor; digits=3))x | " *
                "$(round(row.before_target_2d_rms_um; digits=6)) | " *
                "$(round(row.after_target_2d_rms_um; digits=6)) | " *
                "$(round(row.command_rms; sigdigits=6)) | " *
                "$(round(row.max_abs_command; sigdigits=6)) |",
            )
        end
        println(io)
        println(io, "The reference-ORM row uses the paired zero-quadrupole-offset response")
        println(io, "matrix, representing a stored response measurement from the aligned state.")
        println(io, "The current-ORM row remeasures the response after the offsets are applied.")
        println(io, "Both matrices include the fixed corrector and BPM gains of that machine.")
        if !isempty(response_comparison_rows)
            relative = getproperty.(response_comparison_rows, :relative_l2)
            corrected_delta = getproperty.(
                response_comparison_rows,
                :corrected_bpm_method_difference_rms_um,
            )
            println(io)
            println(
                io,
                "Across paired machines, the current-versus-reference ORM relative-L2 " *
                "difference has median `$(round(median(relative); sigdigits=6))` and " *
                "maximum `$(round(maximum(relative); sigdigits=6))`. The two correction " *
                "methods' final BPM vectors differ by " *
                "`$(round(coordinate_rms(corrected_delta); digits=6))` micrometers RMS.",
            )
        end
        println(io)
        println(io, "## Sextupole scan-range restoration")
        println(io)
        reference_fraction = first(aggregate_rows).reference_truth_outside_scan_radius_fraction
        uncorrected_fraction = first(aggregate_rows).uncorrected_truth_outside_scan_radius_fraction
        scan_radius_mm = 1e3 * getindex(metadata, "scan_radius_m")
        println(
            io,
            "The fraction of beam-relative sextupole centers outside the configured " *
            "$(round(scan_radius_mm; digits=3))-millimeter scan radius is " *
            "`$(round(100 * reference_fraction; digits=3))%` in the zero-offset " *
            "reference and `$(round(100 * uncorrected_fraction; digits=3))%` before correction.",
        )
        println(io)
        for row in aggregate_rows
            println(
                io,
                "- `$(row.method)`: " *
                "`$(round(100 * row.corrected_truth_outside_scan_radius_fraction; digits=3))%` " *
                "outside after correction.",
            )
        end
        println(io)
        println(io, "## Method boundary")
        println(io)
        println(io, "- The reference is the BPM orbit of the paired machine with quadrupole")
        println(io, "  alignment drift disabled; it is not a zero orbit and no oracle target-local")
        println(io, "  coordinates enter correction.")
        println(io, "- BPM readings are noise-free in this bounded test. Stable 1% RMS BPM gains")
        println(io, "  and 1% RMS corrector gains are included in both the observations and the")
        println(io, "  measured response matrices.")
        println(io, "- Corrector commands are the 103 latest-lattice HKICK/VKICK Overlay")
        println(io, "  steering variables in radians. No CESR")
        println(io, "  power-supply or operator limit is asserted.")
        println(io, "- Matching finite BPM readings does not prove exact trajectory equality")
        println(io, "  between BPMs; the target-sextupole orbit residual is reported separately.")
        println(io, "- The latest lattice retains its documented straight-multipole-in-curved-")
        println(io, "  reference qualification. Girder pitch is not varied in this experiment.")
        println(io)
        println(io, "## Reproduction")
        println(io)
        println(io, "```powershell")
        println(io, "julia --project=. sextupole_misalignment/quadrupole_orbit_correction/run_quadrupole_orbit_correction.jl")
        println(io, "python sextupole_misalignment/quadrupole_orbit_correction/validate_orbit_correction.py")
        println(io, "```")
        println(io)
        println(
            io,
            "Machines: `$(getindex(metadata, "machine_count"))`; " *
            "quadrupoles: `$(getindex(metadata, "quadrupole_count"))`; " *
            "BPMs: `$(getindex(metadata, "bpm_count"))`; " *
            "correctors: `$(getindex(metadata, "corrector_count"))`.",
        )
    end
    return path
end

function main_orbit_correction(args=ARGS)
    defaults = Dict(
        "machines" => "16",
        "seed" => "20260823",
        "sextupole-offset-rms-m" => "3.0e-4",
        "corrector-gain-rms" => "0.01",
        "k2-gain-rms" => "0.01",
        "bpm-gain-rms" => "0.01",
        "quadrupole-strength-fraction" => "0.01",
        "quadrupole-roll-rms-rad" => "1.0e-3",
        "quadrupole-alignment-rms-m" => "5.0e-5",
        "response-methods" => join(ORBIT_CORRECTION_METHODS, ','),
        "response-step" => "1.0e-6",
        "ridge-ratio" => "1.0e-2",
        "relative-svd-cutoff" => "1.0e-8",
        "iterations" => "3",
        "tolerance-m" => "1.0e-7",
        "max-update" => "0.0",
        "line-search-steps" => "8",
        "scan-radius-m" => "1.5e-3",
        "output-dir" => joinpath(@__DIR__, "results", "orbit_correction_50um"),
        "overwrite" => "false",
    )
    options = parse_exact11_options(defaults, args)
    machine_count = parse(Int, options["machines"])
    seed = parse(Int, options["seed"])
    methods = parse_method_names(options["response-methods"])
    response_step = parse(Float64, options["response-step"])
    ridge_ratio = parse(Float64, options["ridge-ratio"])
    relative_cutoff = parse(Float64, options["relative-svd-cutoff"])
    max_iterations = parse(Int, options["iterations"])
    tolerance_m = parse(Float64, options["tolerance-m"])
    max_update = parse(Float64, options["max-update"])
    line_search_steps = parse(Int, options["line-search-steps"])
    scan_radius_m = parse(Float64, options["scan-radius-m"])
    quadrupole_alignment_rms = parse(Float64, options["quadrupole-alignment-rms-m"])
    output_dir = abspath(options["output-dir"])
    metadata_path = joinpath(output_dir, "metadata.toml")
    isfile(metadata_path) && lowercase(options["overwrite"]) != "true" &&
        error("Output exists; use --overwrite=true: $metadata_path")

    machine_count > 0 || error("--machines must be positive")
    max_iterations > 0 || error("--iterations must be positive")
    line_search_steps >= 0 || error("--line-search-steps must be nonnegative")
    all(isfinite, (response_step, ridge_ratio, relative_cutoff, tolerance_m, max_update, scan_radius_m)) ||
        error("Correction options must be finite")
    response_step > 0 || error("--response-step must be positive")
    ridge_ratio >= 0 || error("--ridge-ratio must be nonnegative")
    0 < relative_cutoff < 1 || error("--relative-svd-cutoff must lie in (0,1)")
    tolerance_m > 0 || error("--tolerance-m must be positive")
    max_update >= 0 || error("--max-update must be nonnegative")
    scan_radius_m > 0 || error("--scan-radius-m must be positive")
    lowercase(options["overwrite"]) in ("true", "false") ||
        error("--overwrite must be true or false")

    model = load_ring_model(; ring=:latest, zero_value=0.0, rf_on=true)
    ring = model.ring
    sextupoles = active_sextupole_inventory(ring)
    quadrupoles = active_quadrupole_inventory(ring)
    detectors = measurable_bpms(ring)
    control_names = String.(model.metadata.steering_control_names)
    length(control_names) == 103 ||
        error("Expected 103 latest-lattice H/V steering controls, found $(length(control_names))")
    geometry = quadrupole_geometry_joint(ring, quadrupoles)
    bpm_names = String.(base_name.(detectors))
    target_names = String.(getproperty.(sextupoles, :name))
    nominal_v0 = copy(solve_closed_orbit(ring).v0)
    latents = sample_joint_latents(
        machine_count,
        length(sextupoles),
        length(quadrupoles),
        length(control_names),
        length(detectors),
        seed,
        parse(Float64, options["sextupole-offset-rms-m"]),
        parse(Float64, options["corrector-gain-rms"]),
        parse(Float64, options["k2-gain-rms"]),
        parse(Float64, options["bpm-gain-rms"]),
        parse(Float64, options["quadrupole-strength-fraction"]),
        parse(Float64, options["quadrupole-roll-rms-rad"]),
    )

    nm = length(methods)
    nd = length(detectors)
    nt = length(sextupoles)
    nc = length(control_names)
    nq = length(quadrupoles)
    reference_bpm = fill(NaN, machine_count, nd, 2)
    uncorrected_bpm = similar(reference_bpm)
    reference_bpm_measured = similar(reference_bpm)
    uncorrected_bpm_measured = similar(reference_bpm)
    corrected_bpm = fill(NaN, nm, machine_count, nd, 2)
    corrected_bpm_measured = similar(corrected_bpm)
    reference_target = fill(NaN, machine_count, nt, 2)
    uncorrected_target = similar(reference_target)
    corrected_target = fill(NaN, nm, machine_count, nt, 2)
    commands = fill(NaN, nm, machine_count, nc)
    histories = fill(NaN, nm, machine_count, max_iterations + 1)
    singular_values = fill(NaN, nm, machine_count, nc)
    actual_offsets = quadrupole_alignment_rms .* latents.quadrupole_alignment_standard_normals
    per_machine_rows = NamedTuple[]
    response_comparison_rows = NamedTuple[]
    started = time()

    try
        for machine in 1:machine_count
            machine_start = time()
            bpm_gain_errors = Matrix(latents.bpm_gain_errors[machine, :, :])
            physical_corrector_gains = 1 .+ Vector(latents.corrector_gain_errors[machine, :])

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
            reference_closed = scalar_closed_orbit_joint(ring, nominal_v0)
            ref_bpm, ref_target =
                read_machine_orbits(ring, reference_closed, bpm_names, target_names)
            ref_measured = measured_bpm(ref_bpm, bpm_gain_errors)
            reference_bpm[machine, :, :] .= ref_bpm
            reference_bpm_measured[machine, :, :] .= ref_measured
            reference_target[machine, :, :] .= ref_target
            reference_response = nothing
            if "reference_orm" in methods
                reference_response = finite_difference_measured_orm!(
                    model,
                    control_names,
                    physical_corrector_gains,
                    zeros(nc),
                    bpm_names,
                    bpm_gain_errors,
                    reference_closed,
                    response_step,
                )
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
            current_closed = scalar_closed_orbit_joint(ring, nominal_v0)
            cur_bpm, cur_target =
                read_machine_orbits(ring, current_closed, bpm_names, target_names)
            cur_measured = measured_bpm(cur_bpm, bpm_gain_errors)
            uncorrected_bpm[machine, :, :] .= cur_bpm
            uncorrected_bpm_measured[machine, :, :] .= cur_measured
            uncorrected_target[machine, :, :] .= cur_target
            current_response = nothing
            if "current_orm" in methods
                current_response = finite_difference_measured_orm!(
                    model,
                    control_names,
                    physical_corrector_gains,
                    zeros(nc),
                    bpm_names,
                    bpm_gain_errors,
                    current_closed,
                    response_step,
                )
            end

            if !isnothing(reference_response) && !isnothing(current_response)
                response_difference = current_response - reference_response
                push!(response_comparison_rows, (;
                    machine,
                    relative_l2=norm(response_difference) / norm(reference_response),
                    maximum_abs_difference=maximum(abs, response_difference),
                    cosine_similarity=dot(vec(reference_response), vec(current_response)) /
                        (norm(reference_response) * norm(current_response)),
                    corrected_bpm_method_difference_rms_um=NaN,
                    corrector_command_method_difference_rms=NaN,
                ))
            end

            for (method_index, method) in enumerate(methods)
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
                start_closed = scalar_closed_orbit_joint(ring, current_closed.v0)
                response = method == "reference_orm" ? reference_response : current_response
                result = solve_corrected_machine!(
                    model,
                    control_names,
                    physical_corrector_gains,
                    bpm_names,
                    target_names,
                    bpm_gain_errors,
                    ref_measured,
                    response,
                    start_closed;
                    ridge_ratio,
                    relative_cutoff,
                    max_iterations,
                    tolerance_m,
                    max_update,
                    line_search_steps,
                )
                corrected_bpm[method_index, machine, :, :] .= result.physical_bpm
                corrected_bpm_measured[method_index, machine, :, :] .= result.measured_bpm
                corrected_target[method_index, machine, :, :] .= result.physical_target
                commands[method_index, machine, :] .= result.commands
                histories[method_index, machine, :] .= result.history
                singular_values[method_index, machine, :] .= result.singular

                before_bpm = cur_measured - ref_measured
                after_bpm = result.measured_bpm - ref_measured
                before_target = cur_target - ref_target
                after_target = result.physical_target - ref_target
                push!(per_machine_rows, (;
                    method,
                    machine,
                    realized_quad_x_rms_um=coordinate_rms(actual_offsets[machine, :, 1]) * 1e6,
                    realized_quad_y_rms_um=coordinate_rms(actual_offsets[machine, :, 2]) * 1e6,
                    before_bpm_rms_um=coordinate_rms(before_bpm) * 1e6,
                    after_bpm_rms_um=coordinate_rms(after_bpm) * 1e6,
                    bpm_reduction_factor=coordinate_rms(before_bpm) / coordinate_rms(after_bpm),
                    before_bpm_max_abs_um=maximum(abs, before_bpm) * 1e6,
                    after_bpm_max_abs_um=maximum(abs, after_bpm) * 1e6,
                    before_target_2d_rms_um=two_plane_rms(before_target) * 1e6,
                    after_target_2d_rms_um=two_plane_rms(after_target) * 1e6,
                    command_rms=coordinate_rms(result.commands),
                    max_abs_command=maximum(abs, result.commands),
                    iterations=result.iterations,
                    converged=result.converged,
                    retained_response_rank=result.retained_rank,
                    retained_response_condition=result.condition,
                ))
            end
            if !isempty(response_comparison_rows) &&
               response_comparison_rows[end].machine == machine
                reference_index = findfirst(==("reference_orm"), methods)
                current_index = findfirst(==("current_orm"), methods)
                if !isnothing(reference_index) && !isnothing(current_index)
                    row = response_comparison_rows[end]
                    response_comparison_rows[end] = merge(row, (;
                        corrected_bpm_method_difference_rms_um=coordinate_rms(
                            corrected_bpm_measured[reference_index, machine, :, :] -
                            corrected_bpm_measured[current_index, machine, :, :],
                        ) * 1e6,
                        corrector_command_method_difference_rms=coordinate_rms(
                            commands[reference_index, machine, :] -
                            commands[current_index, machine, :],
                        ),
                    ))
                end
            end
            @printf(
                "machine %d/%d complete in %.2f s (elapsed %.2f s)\n",
                machine,
                machine_count,
                time() - machine_start,
                time() - started,
            )
            flush(stdout)
        end
    finally
        set_all_controls!(model.controls, 0.0)
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
    end

    aggregate_rows = NamedTuple[]
    reference_truth = latents.sextupole_offsets - reference_target
    uncorrected_truth = latents.sextupole_offsets - uncorrected_target
    for (method_index, method) in enumerate(methods)
        before_bpm = uncorrected_bpm_measured - reference_bpm_measured
        after_bpm = corrected_bpm_measured[method_index, :, :, :] - reference_bpm_measured
        before_target = uncorrected_target - reference_target
        after_target = corrected_target[method_index, :, :, :] - reference_target
        corrected_truth =
            latents.sextupole_offsets - corrected_target[method_index, :, :, :]
        push!(aggregate_rows, (;
            method,
            before_bpm_rms_um=coordinate_rms(before_bpm) * 1e6,
            after_bpm_rms_um=coordinate_rms(after_bpm) * 1e6,
            reduction_factor=coordinate_rms(before_bpm) / coordinate_rms(after_bpm),
            before_bpm_max_abs_um=maximum(abs, before_bpm) * 1e6,
            after_bpm_max_abs_um=maximum(abs, after_bpm) * 1e6,
            before_target_2d_rms_um=two_plane_rms(reshape(before_target, :, 2)) * 1e6,
            after_target_2d_rms_um=two_plane_rms(reshape(after_target, :, 2)) * 1e6,
            command_rms=coordinate_rms(commands[method_index, :, :]),
            max_abs_command=maximum(abs, commands[method_index, :, :]),
            reference_truth_outside_scan_radius_fraction=mean(
                sqrt.(sum(abs2, reference_truth; dims=3)) .> scan_radius_m,
            ),
            uncorrected_truth_outside_scan_radius_fraction=mean(
                sqrt.(sum(abs2, uncorrected_truth; dims=3)) .> scan_radius_m,
            ),
            corrected_truth_outside_scan_radius_fraction=mean(
                sqrt.(sum(abs2, corrected_truth; dims=3)) .> scan_radius_m,
            ),
            converged_machine_count=count(
                row -> row.method == method && row.converged,
                per_machine_rows,
            ),
        ))
    end

    metadata = Dict(
        "format" => "cesr-paired-quadrupole-orbit-correction-v1",
        "date" => string(Dates.today()),
        "engine" => "SciBmad exact RF-on closed orbit and tracking",
        "lattice" => LATEST_LATTICE,
        "machine_count" => machine_count,
        "quadrupole_count" => nq,
        "bpm_count" => nd,
        "target_count" => nt,
        "corrector_count" => nc,
        "response_methods" => methods,
        "random_seed_base" => seed,
        "quadrupole_alignment_rms_m_per_plane" => quadrupole_alignment_rms,
        "quadrupole_alignment_distribution" => "independent Gaussian x/y per physical quadrupole, coherent across slices",
        "paired_errors" => "all other maintained random draws fixed between zero-offset reference and 50-um/plane offset state",
        "sextupole_offset_rms_m_per_plane" => parse(Float64, options["sextupole-offset-rms-m"]),
        "corrector_gain_rms" => parse(Float64, options["corrector-gain-rms"]),
        "k2_gain_rms" => parse(Float64, options["k2-gain-rms"]),
        "bpm_gain_rms" => parse(Float64, options["bpm-gain-rms"]),
        "quadrupole_strength_halfwidth_fraction" => parse(Float64, options["quadrupole-strength-fraction"]),
        "quadrupole_roll_rms_rad" => parse(Float64, options["quadrupole-roll-rms-rad"]),
        "measurement_noise_rms_m" => 0.0,
        "response_step_corrector_field" => response_step,
        "ridge_ratio_to_largest_singular_value" => ridge_ratio,
        "relative_svd_cutoff" => relative_cutoff,
        "maximum_iterations" => max_iterations,
        "bpm_rms_tolerance_m" => tolerance_m,
        "maximum_update_per_iteration" => max_update,
        "line_search_steps" => line_search_steps,
        "scan_radius_m" => scan_radius_m,
        "calculation_wall_seconds" => time() - started,
        "reference_semantics" => "paired measured BPM closed orbit with quadrupole alignment drift disabled; not a zero-orbit target",
        "solver_inputs" => "reference/current BPM readings and measured corrector-to-BPM response only; latent quadrupole offsets excluded",
        "corrector_registry" => "103 dynamically selected latest-lattice H/V Overlay controls",
        "corrector_units" => "HKICK/VKICK control commands in radians; no CESR hardware limit asserted",
    )

    mkpath(output_dir)
    write_npy(joinpath(output_dir, "reference_bpm_orbits.npy"), reference_bpm)
    write_npy(joinpath(output_dir, "uncorrected_bpm_orbits.npy"), uncorrected_bpm)
    write_npy(joinpath(output_dir, "reference_bpm_readbacks.npy"), reference_bpm_measured)
    write_npy(joinpath(output_dir, "uncorrected_bpm_readbacks.npy"), uncorrected_bpm_measured)
    write_npy(joinpath(output_dir, "corrected_bpm_orbits.npy"), corrected_bpm)
    write_npy(joinpath(output_dir, "corrected_bpm_readbacks.npy"), corrected_bpm_measured)
    write_npy(joinpath(output_dir, "reference_target_orbits.npy"), reference_target)
    write_npy(joinpath(output_dir, "uncorrected_target_orbits.npy"), uncorrected_target)
    write_npy(joinpath(output_dir, "corrected_target_orbits.npy"), corrected_target)
    write_npy(joinpath(output_dir, "corrector_commands.npy"), commands)
    write_npy(joinpath(output_dir, "bpm_rms_history_m.npy"), histories)
    write_npy(joinpath(output_dir, "response_singular_values.npy"), singular_values)
    write_npy(joinpath(output_dir, "quadrupole_offsets_m.npy"), actual_offsets)
    write_npy(joinpath(output_dir, "sextupole_offsets_m.npy"), latents.sextupole_offsets)
    write_lines(joinpath(output_dir, "method_names.txt"), methods)
    write_lines(joinpath(output_dir, "bpm_names.txt"), bpm_names)
    write_lines(joinpath(output_dir, "target_names.txt"), target_names)
    write_lines(joinpath(output_dir, "corrector_names.txt"), control_names)
    write_rows(joinpath(output_dir, "per_machine.csv"), per_machine_rows)
    write_rows(joinpath(output_dir, "aggregate.csv"), aggregate_rows)
    isempty(response_comparison_rows) || write_rows(
        joinpath(output_dir, "response_comparison.csv"),
        response_comparison_rows,
    )
    write_metadata(metadata_path, metadata)
    markdown_summary(
        joinpath(output_dir, "SUMMARY.md"),
        aggregate_rows,
        response_comparison_rows,
        metadata,
    )

    println("Orbit correction complete: $output_dir")
    for row in aggregate_rows
        @printf(
            "  %-13s BPM %.3f -> %.3f um (%.1fx), target %.3f -> %.3f um 2D\n",
            row.method,
            row.before_bpm_rms_um,
            row.after_bpm_rms_um,
            row.reduction_factor,
            row.before_target_2d_rms_um,
            row.after_target_2d_rms_um,
        )
    end
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main_orbit_correction())
end

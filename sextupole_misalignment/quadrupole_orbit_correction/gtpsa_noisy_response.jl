"""GTPSA ORM and noisy-BPM helpers for the paired orbit-correction study."""

module CESROrbitGTPSABackend
include(joinpath(
    @__DIR__, "..", "..", "orbit", "Orbit_Calculation", "benchmark_scibmad.jl",
))
end

const GTPSA_NOISY_CORRECTED_CASE =
    "with_quadrupole_misalignment_gtpsa_noisy_corrected"

const GTPSA_NOMINAL_CORRECTED_CASE =
    "with_all_errors_gtpsa_nominal_corrected"

function measurement_mean_noise_std(noise_rms_m, repeats)
    isfinite(noise_rms_m) && noise_rms_m >= 0 ||
        error("BPM noise RMS must be finite and nonnegative")
    repeats > 0 || error("BPM measurement repeats must be positive")
    return noise_rms_m / sqrt(repeats)
end

function sample_bpm_mean_noise(rng, detector_count, standard_deviation_m)
    standard_deviation_m == 0 && return zeros(detector_count, 2)
    return standard_deviation_m .* randn(rng, detector_count, 2)
end

function noisy_measured_bpm(
    physical,
    gain_errors,
    rng,
    standard_deviation_m,
)
    deterministic = measured_bpm(physical, gain_errors)
    noise = sample_bpm_mean_noise(
        rng, size(deterministic, 1), standard_deviation_m,
    )
    return deterministic + noise, noise
end


"""Return the error-conditioned physical BPM ORM from one first-order GTPSA map."""
function gtpsa_physical_bpm_orm(
    control_names,
    bpm_names,
    nominal_orbit,
    machine,
    case_name,
    latents,
    quadrupole_alignment_rms,
)
    function machine_model_factory(; zero_value=0.0, rf_on=true, kwargs...)
        response_model = load_ring_model(
            ; ring=:latest, zero_value, rf_on,
        )
        response_sextupoles = active_sextupole_inventory(response_model.ring)
        response_quadrupoles = active_quadrupole_inventory(response_model.ring)
        response_geometry = quadrupole_geometry_joint(
            response_model.ring, response_quadrupoles,
        )
        prepare_machine!(
            response_model,
            case_name,
            machine,
            response_sextupoles,
            response_quadrupoles,
            response_geometry,
            latents,
            quadrupole_alignment_rms,
        )
        return response_model
    end

    timed = @timed CESROrbitGTPSABackend.gtpsa_first_order_responses(
        String.(control_names);
        nominal_orbit=vec(Float64.(Array(nominal_orbit))),
        model_factory=machine_model_factory,
    )
    result = timed.value
    detector_index = Dict(
        uppercase(String(name)) => index
        for (index, name) in enumerate(result.detectors)
    )
    selected = Int[]
    for name in bpm_names
        key = uppercase(String(name))
        haskey(detector_index, key) ||
            error("GTPSA response is missing measurable BPM $name")
        push!(selected, detector_index[key])
    end
    all_count = length(result.detectors)
    response = vcat(
        result.detector_response[selected, :],
        result.detector_response[all_count .+ selected, :],
    )
    expected = (2length(bpm_names), length(control_names))
    size(response) == expected ||
        error("Unexpected selected GTPSA ORM shape $(size(response)); expected $expected")
    all(isfinite, response) || error("GTPSA ORM contains a non-finite value")
    return (;
        response,
        seconds=timed.time,
        closure_norm_max=result.closure_norm_max,
        full_detector_count=all_count,
    )
end


"""Return the nominal latest-lattice BPM ORM from one first-order GTPSA map.

Unlike `gtpsa_physical_bpm_orm`, this model factory receives no latent-machine
realization.  It is therefore suitable for the deliberately mismatched
correction protocol in which BPM gains, corrector gains, magnet errors, and
alignment errors remain unknown to the response model.
"""
function nominal_gtpsa_bpm_orm(control_names, bpm_names, nominal_orbit)
    nominal_model_factory(; zero_value=0.0, rf_on=true, kwargs...) =
        load_ring_model(; ring=:latest, zero_value, rf_on)

    timed = @timed CESROrbitGTPSABackend.gtpsa_first_order_responses(
        String.(control_names);
        nominal_orbit=vec(Float64.(Array(nominal_orbit))),
        model_factory=nominal_model_factory,
    )
    result = timed.value
    detector_index = Dict(
        uppercase(String(name)) => index
        for (index, name) in enumerate(result.detectors)
    )
    selected = Int[]
    for name in bpm_names
        key = uppercase(String(name))
        haskey(detector_index, key) ||
            error("Nominal GTPSA response is missing measurable BPM $name")
        push!(selected, detector_index[key])
    end
    all_count = length(result.detectors)
    response = vcat(
        result.detector_response[selected, :],
        result.detector_response[all_count .+ selected, :],
    )
    expected = (2length(bpm_names), length(control_names))
    size(response) == expected ||
        error("Unexpected nominal GTPSA ORM shape $(size(response)); expected $expected")
    all(isfinite, response) || error("Nominal GTPSA ORM contains a non-finite value")
    return (;
        response,
        seconds=timed.time,
        closure_norm_max=result.closure_norm_max,
        full_detector_count=all_count,
    )
end


"""Apply fixed BPM and corrector gains to a physical GTPSA ORM."""
function measured_gtpsa_orm(
    control_names,
    physical_corrector_gains,
    bpm_names,
    bpm_gain_errors,
    nominal_orbit,
    machine,
    case_name,
    latents,
    quadrupole_alignment_rms,
)
    physical = gtpsa_physical_bpm_orm(
        control_names,
        bpm_names,
        nominal_orbit,
        machine,
        case_name,
        latents,
        quadrupole_alignment_rms,
    )
    row_gains = vcat(
        1 .+ bpm_gain_errors[:, 1],
        1 .+ bpm_gain_errors[:, 2],
    )
    response = physical.response .* row_gains .* transpose(physical_corrector_gains)
    all(isfinite, response) || error("Measured GTPSA ORM contains a non-finite value")
    return merge(physical, (; response))
end


"""Iteratively correct to a stored noisy reference using fresh noisy BPM means."""
function solve_noisy_corrected_machine!(
    model,
    control_names,
    physical_corrector_gains,
    bpm_names,
    target_names,
    bpm_gain_errors,
    reference_measured,
    noiseless_reference_measured,
    response,
    initial_closed,
    measurement_rng,
    validation_rng,
    measurement_noise_std_m;
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
    measured, measurement_noise = noisy_measured_bpm(
        physical_bpm,
        bpm_gain_errors,
        measurement_rng,
        measurement_noise_std_m,
    )
    residual = flatten_bpm(measured - reference_measured)
    history = fill(NaN, max_iterations + 1)
    physical_history = fill(NaN, max_iterations + 1)
    history[1] = coordinate_rms(residual)
    physical_history[1] = coordinate_rms(
        measured_bpm(physical_bpm, bpm_gain_errors) - noiseless_reference_measured,
    )
    accepted_iterations = 0

    for iteration in 1:max_iterations
        history[iteration] <= tolerance_m && break
        update = ridge_command(
            factorization, residual, ridge_ratio, relative_cutoff,
        )
        if max_update > 0 && maximum(abs, update) > max_update
            update .*= max_update / maximum(abs, update)
        end
        maximum(abs, update) > eps(Float64) || break

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
            trial_measured, trial_noise = noisy_measured_bpm(
                trial_bpm,
                bpm_gain_errors,
                measurement_rng,
                measurement_noise_std_m,
            )
            trial_residual = flatten_bpm(trial_measured - reference_measured)
            trial_rms = coordinate_rms(trial_residual)
            if trial_rms < history[iteration]
                commands = trial_commands
                closed = trial_closed
                physical_bpm = trial_bpm
                physical_target = trial_target
                measured = trial_measured
                measurement_noise = trial_noise
                residual = trial_residual
                history[iteration + 1] = trial_rms
                physical_history[iteration + 1] = coordinate_rms(
                    measured_bpm(physical_bpm, bpm_gain_errors) -
                    noiseless_reference_measured,
                )
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

    validation_measured, validation_noise = noisy_measured_bpm(
        physical_bpm,
        bpm_gain_errors,
        validation_rng,
        measurement_noise_std_m,
    )
    validation_residual = flatten_bpm(validation_measured - reference_measured)
    final_rms = coordinate_rms(residual)
    return (;
        commands,
        closed,
        physical_bpm,
        physical_target,
        measured_bpm=measured,
        residual,
        history,
        physical_history,
        iterations=accepted_iterations,
        converged=final_rms <= tolerance_m,
        singular,
        retained_rank,
        condition,
        final_measurement_noise=measurement_noise,
        validation_measured_bpm=validation_measured,
        validation_noise,
        validation_residual,
    )
end

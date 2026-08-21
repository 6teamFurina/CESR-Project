#!/usr/bin/env julia

"""Build repaired-lattice two-plane local-orbit knobs for all target sextupoles."""

include(joinpath(@__DIR__, "common.jl"))

function solve_bump_knob(detector_response, target_response, desired, ridge_fraction)
    gram = detector_response' * detector_response
    scale = tr(gram) / size(gram, 1)
    scale > 0 || error("Detector-response Gram matrix has zero scale")
    objective = gram + ridge_fraction * scale * I
    kkt = [objective target_response'; target_response zeros(2, 2)]
    solution = kkt \ vcat(zeros(size(objective, 1)), desired)
    return solution[1:size(objective, 1)]
end

function main(args=ARGS)
    defaults = Dict(
        "ridge-fraction" => "1e-6",
        "design-amplitude-m" => "5e-4",
        "output-dir" => joinpath(@__DIR__, "results", "bump_knobs"),
    )
    options = parse_exact11_options(defaults, args)
    ridge_fraction = parse(Float64, options["ridge-fraction"])
    design_amplitude_m = parse(Float64, options["design-amplitude-m"])
    output_dir = abspath(options["output-dir"])
    all(isfinite.((ridge_fraction, design_amplitude_m))) || error("Options must be finite")
    minimum((ridge_fraction, design_amplitude_m)) > 0 || error("Options must be positive")

    ring = cesr
    controls = independent_corrector_inventory(ring)
    sextupoles = active_sextupole_inventory(ring)
    detectors = measurable_bpms(ring)
    scalar_closed_timed = @timed solve_closed_orbit(ring)
    descriptor = Descriptor(6, 2, length(controls), 1)
    parameters = params(descriptor)
    for (index, control) in enumerate(controls)
        baseline = constant_term(first(control.originals))
        for element_index in control.indices
            if control.axis == :Kn0
                ring.line[element_index].Kn0 = baseline + parameters[index]
            else
                ring.line[element_index].Ks0 = baseline + parameters[index]
            end
        end
    end

    at = vcat(detectors, [ring.line[target.index] for target in sextupoles])
    optics_timed = @timed calculate_twiss(
        ring,
        at,
        scalar_closed_timed.value;
        descriptor,
    )
    optics = optics_timed.value
    index_by_name = table_index_by_name(optics)

    function orbit_parameter_jacobian(values)
        full = Matrix(GTPSA.jacobian(collect(values); include_params=true))
        return full[:, 7:end]
    end

    function response_matrix(names)
        indices = [index_by_name[uppercase(name)] for name in names]
        result = zeros(2 * length(indices), length(controls))
        result[1:2:end, :] .= orbit_parameter_jacobian(optics.table.orbit_x[indices])
        result[2:2:end, :] .= orbit_parameter_jacobian(optics.table.orbit_y[indices])
        return result
    end

    detector_names = String.(base_name.(detectors))
    detector_response = response_matrix(detector_names)
    target_response_all = response_matrix(String.(getproperty.(sextupoles, :name)))
    knob_rows = NamedTuple[]
    summary_rows = NamedTuple[]
    for (target_index, target) in enumerate(sextupoles)
        target_response = target_response_all[2target_index - 1:2target_index, :]
        knob_x = solve_bump_knob(
            detector_response,
            target_response,
            [design_amplitude_m, 0.0],
            ridge_fraction,
        )
        knob_y = solve_bump_knob(
            detector_response,
            target_response,
            [0.0, design_amplitude_m],
            ridge_fraction,
        )
        achieved_x = target_response * knob_x
        achieved_y = target_response * knob_y
        leakage_x = detector_response * knob_x
        leakage_y = detector_response * knob_y
        for index in eachindex(controls)
            push!(knob_rows, (;
                target_sextupole=target.name,
                corrector=controls[index].name,
                field=String(controls[index].axis),
                field_per_x_bump_m=knob_x[index] / design_amplitude_m,
                field_per_y_bump_m=knob_y[index] / design_amplitude_m,
            ))
        end
        push!(summary_rows, (;
            target_sextupole=target.name,
            target_s_m=target.s_m,
            design_amplitude_m,
            achieved_x_from_x_knob_m=achieved_x[1],
            achieved_y_from_x_knob_m=achieved_x[2],
            achieved_x_from_y_knob_m=achieved_y[1],
            achieved_y_from_y_knob_m=achieved_y[2],
            detector_leakage_x_knob_rms_m=sqrt(mean(leakage_x .^ 2)),
            detector_leakage_y_knob_rms_m=sqrt(mean(leakage_y .^ 2)),
            maximum_abs_x_corrector_field=maximum(abs, knob_x),
            maximum_abs_y_corrector_field=maximum(abs, knob_y),
        ))
    end

    knobs_path = write_rows(joinpath(output_dir, "local_bump_knobs.csv"), knob_rows)
    summary_path = write_rows(joinpath(output_dir, "local_bump_knob_summary.csv"), summary_rows)
    metadata = Dict(
        "format" => "cesr-repaired-lattice-local-bump-knobs-v1",
        "date" => string(Dates.today()),
        "engine" => "SciBmad exact RF-on batch closed orbit and tracking",
        "lattice" => LATEST_LATTICE,
        "target_count" => length(sextupoles),
        "bpm_count" => length(detectors),
        "corrector_count" => length(controls),
        "ridge_fraction" => ridge_fraction,
        "design_amplitude_m" => design_amplitude_m,
        "descriptor" => "Descriptor(6, 2, $(length(controls)), 1)",
        "closed_orbit_seconds" => scalar_closed_timed.time,
        "parameterized_twiss_seconds" => optics_timed.time,
        "corrector_policy" => "unambiguous horizontal/vertical Kicker elements; ambiguous bump and pinger elements excluded",
        "interpretation_boundary" => "model-based local-bump design; machine use requires hardware calibration, limits, and operator validation",
    )
    write_metadata(joinpath(output_dir, "local_bump_knob_metadata.toml"), metadata)
    println("Controls: $(length(controls))")
    println("Targets: $(length(sextupoles))")
    println("Knobs: $knobs_path")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end

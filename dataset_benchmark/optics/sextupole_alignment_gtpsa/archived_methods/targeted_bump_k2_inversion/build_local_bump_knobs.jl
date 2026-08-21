#!/usr/bin/env julia

"""Build two nominal corrector-space knobs per target sextupole.

The two hard constraints are unit horizontal or vertical closed-orbit
displacement at the target. The objective minimizes the full-ring detector
orbit response plus a small corrector-norm ridge. The result is a reproducible
model-based local-bump baseline, not yet a machine-approved corrector pattern.
"""

include(joinpath(@__DIR__, "common.jl"))

function parameterized_corrector_twiss(target_names)
    names = corrector_names()
    descriptor = Descriptor(6, 2, length(names), 1)
    parameters = params(descriptor)
    model = load_cesr_model(zero_value=zero(parameters[1]), rf_on=true)
    for index in eachindex(names)
        model.controls[names[index]] = parameters[index]
    end

    scalar_model = load_cesr_model(zero_value=0.0, rf_on=true)
    closed = solve_rf_on_closed_orbit(scalar_model.ring)
    inventory = active_sextupole_inventory(model.ring)
    targets = [find_inventory_entry(inventory, name) for name in target_names]
    at_elements = vcat(
        detector_elements(model.ring),
        [model.ring.line[target.index] for target in targets],
    )
    optics = twiss(
        model.ring;
        GTPSA_descriptor=descriptor,
        at=at_elements,
        v0=closed.v0,
        v0_and_coast=(closed.v0, closed.coasting_beam),
        spin=false,
        RDTs=false,
        normalizing_map=false,
    )
    return (; names, descriptor, optics, targets)
end

function orbit_parameter_jacobian(values)
    full = Matrix(GTPSA.jacobian(collect(values); include_params=true))
    return full[:, 7:end]
end

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
        "targets" => "SEX_08W",
        "ridge-fraction" => "1e-6",
        "design-amplitude-m" => "1e-3",
        "output-dir" => joinpath(@__DIR__, "results", "bump_knobs"),
    )
    options = parse_key_value_args(defaults, args)
    targets = uppercase.(strip.(split(options["targets"], ',')))
    ridge_fraction = parse(Float64, options["ridge-fraction"])
    design_amplitude_m = parse(Float64, options["design-amplitude-m"])
    output_dir = abspath(options["output-dir"])
    ridge_fraction > 0 || error("--ridge-fraction must be positive")
    design_amplitude_m > 0 || error("--design-amplitude-m must be positive")

    timed = @timed parameterized_corrector_twiss(targets)
    setup = timed.value
    table_names = uppercase.(first.(split.(String.(setup.optics.table.name), '!')))
    detector_indices = findall(name -> startswith(name, "DET_"), table_names)
    detector_x = orbit_parameter_jacobian(setup.optics.table.orbit_x[detector_indices])
    detector_y = orbit_parameter_jacobian(setup.optics.table.orbit_y[detector_indices])
    detector_response = vcat(detector_x, detector_y)

    knob_rows = NamedTuple[]
    summary_rows = NamedTuple[]
    for target in setup.targets
        target_index = findfirst(==(target.name), table_names)
        isnothing(target_index) && error("Missing target Twiss row: $(target.name)")
        target_response = vcat(
            orbit_parameter_jacobian([setup.optics.table.orbit_x[target_index]]),
            orbit_parameter_jacobian([setup.optics.table.orbit_y[target_index]]),
        )
        desired_x = [design_amplitude_m, 0.0]
        desired_y = [0.0, design_amplitude_m]
        knob_x = solve_bump_knob(
            detector_response, target_response, desired_x, ridge_fraction,
        )
        knob_y = solve_bump_knob(
            detector_response, target_response, desired_y, ridge_fraction,
        )
        achieved_x = target_response * knob_x
        achieved_y = target_response * knob_y
        leakage_x = detector_response * knob_x
        leakage_y = detector_response * knob_y
        for index in eachindex(setup.names)
            push!(knob_rows, (;
                target_sextupole=target.name,
                corrector=setup.names[index],
                kick_x_rad_per_m=knob_x[index] / design_amplitude_m,
                kick_y_rad_per_m=knob_y[index] / design_amplitude_m,
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
            maximum_abs_x_corrector_rad=maximum(abs, knob_x),
            maximum_abs_y_corrector_rad=maximum(abs, knob_y),
            ridge_fraction,
        ))
    end

    knobs_path = write_rows(joinpath(output_dir, "local_bump_knobs.csv"), knob_rows)
    summary_path = write_rows(joinpath(output_dir, "local_bump_knob_summary.csv"), summary_rows)
    metadata = Dict(
        "format" => "cesr-targeted-sextupole-local-bump-knobs-v1",
        "date" => string(Dates.today()),
        "method" => "corrector-space constrained least squares using one RF-on first-order GTPSA Twiss",
        "targets" => targets,
        "corrector_count" => length(setup.names),
        "detector_count" => length(detector_indices),
        "descriptor" => "Descriptor(6, 2, 119, 1)",
        "design_amplitude_m" => design_amplitude_m,
        "ridge_fraction" => ridge_fraction,
        "parameterized_twiss_seconds" => timed.time,
        "knobs_csv" => knobs_path,
        "summary_csv" => summary_path,
        "interpretation_boundary" => "model-based protocol knob; machine use requires hardware limits and operator validation",
    )
    mkpath(output_dir)
    open(joinpath(output_dir, "local_bump_knob_metadata.toml"), "w") do io
        TOML.print(io, metadata; sorted=true)
    end
    println("Bump knobs: $knobs_path")
    println("Bump summary: $summary_path")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end


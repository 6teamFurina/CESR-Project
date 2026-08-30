#!/usr/bin/env julia

"""Generate nominal latest-lattice models for finite-BPM local-orbit prediction."""

include(joinpath(
    @__DIR__, "..", "quadrupole_affinity", "exact_11_triplet_validation", "common.jl",
))

function map_jacobian(map)
    return [
        Float64(constant_term(GTPSA.deriv(map.v[i], j)))
        for i in 1:6, j in 1:6
    ]
end

function line_index_by_base_name(ring, requested_names)
    wanted = Set(uppercase.(String.(requested_names)))
    result = Dict{String,Int}()
    for (index, element) in enumerate(ring.line)
        name = uppercase(String(base_name(element)))
        name in wanted || continue
        haskey(result, name) || (result[name] = index)
    end
    missing = setdiff(wanted, Set(keys(result)))
    isempty(missing) || error("Missing lattice locations: $(join(sort!(collect(missing)), ", "))")
    return result
end

function main(args=ARGS)
    defaults = Dict(
        "output-dir" => joinpath(@__DIR__, "results", "local_orbit_model"),
    )
    options = parse_exact11_options(defaults, args)
    output_dir = abspath(options["output-dir"])

    ring = cesr
    controls = independent_corrector_inventory(ring)
    sextupoles = active_sextupole_inventory(ring)
    detectors = measurable_bpms(ring)
    bpm_names = String.(base_name.(detectors))
    target_names = String.(getproperty.(sextupoles, :name))

    closed_timed = @timed solve_closed_orbit(ring)
    closed = closed_timed.value

    # Save the nominal closed-orbit anchor used by the two-sided absolute-orbit
    # reconstruction.  These are model quantities, not target-local truth from
    # any latent-machine scan.  A machine-facing estimate transports the
    # measured BPM-minus-nominal residual from the neighboring BPM pair and
    # adds it to this nominal target orbit.
    nominal_tracked = track_orbits_at_names(
        ring,
        closed,
        vcat(bpm_names, target_names),
    )
    nominal_bpm_orbits = zeros(length(detectors), 2)
    nominal_target_orbits = zeros(length(sextupoles), 2)
    for (index, name) in enumerate(bpm_names)
        nominal_bpm_orbits[index, :] .= (
            nominal_tracked.horizontal[name][1],
            nominal_tracked.vertical[name][1],
        )
    end
    for (index, name) in enumerate(target_names)
        nominal_target_orbits[index, :] .= (
            nominal_tracked.horizontal[name][1],
            nominal_tracked.vertical[name][1],
        )
    end

    # Generate ordinary six-dimensional cumulative maps before introducing
    # corrector parameters. This preserves the BeamlineParams attached to the
    # loaded lattice while keeping the order-1 and parameter descriptors apart.
    transport_at = vcat(
        detectors,
        [ring.line[target.index] for target in sextupoles],
        [ring.line[end]],
    )
    transport_timed = @timed calculate_direct_transport(
        ring, transport_at, closed; descriptor=Descriptor(6, 1),
    )
    transport = transport_timed.value
    transport_index = Dict(
        first(split(uppercase(String(name)), '!')) => index
        for (index, name) in enumerate(transport.names)
    )
    bpm_maps = zeros(length(detectors), 6, 6)
    target_maps = zeros(length(sextupoles), 6, 6)
    for (index, name) in enumerate(bpm_names)
        bpm_maps[index, :, :] .= map_jacobian(transport.maps[transport_index[name]])
    end
    for (index, name) in enumerate(target_names)
        target_maps[index, :, :] .= map_jacobian(transport.maps[transport_index[name]])
    end
    one_turn_map = map_jacobian(transport.maps[end])
    location_index = line_index_by_base_name(ring, vcat(bpm_names, target_names))

    # Parameterize every independent physical corrector simultaneously.  The
    # resulting Jacobians use the same control inventory and element-entry
    # coordinate convention as the local-bump generator and saved scan tensor.
    descriptor = Descriptor(6, 2, length(controls), 1)
    parameters = params(descriptor)
    for (index, control) in enumerate(controls)
        baseline = Float64(constant_term(first(control.originals)))
        for element_index in control.indices
            if control.axis == :Kn0
                ring.line[element_index].Kn0 = baseline + parameters[index]
            else
                ring.line[element_index].Ks0 = baseline + parameters[index]
            end
        end
    end

    at = vcat(detectors, [ring.line[target.index] for target in sextupoles])
    response_timed = @timed calculate_twiss(ring, at, closed; descriptor)
    optics = response_timed.value
    optics_index = table_index_by_name(optics)

    function orbit_parameter_jacobian(values)
        full = Matrix(GTPSA.jacobian(collect(values); include_params=true))
        return full[:, 7:end]
    end

    function response_matrix(names)
        indices = [optics_index[uppercase(name)] for name in names]
        result = zeros(2length(indices), length(controls))
        result[1:2:end, :] .= orbit_parameter_jacobian(optics.table.orbit_x[indices])
        result[2:2:end, :] .= orbit_parameter_jacobian(optics.table.orbit_y[indices])
        return result
    end

    bpm_response = response_matrix(bpm_names)
    target_response = response_matrix(target_names)

    bpm_rows = [
        (; bpm_index=index, bpm=name, line_index=location_index[name])
        for (index, name) in enumerate(bpm_names)
    ]
    target_rows = [
        (;
            target_index=index,
            target=name,
            line_index=target.index,
            s_m=target.s_m,
        )
        for (index, (name, target)) in enumerate(zip(target_names, sextupoles))
    ]
    control_rows = [
        (;
            control_index=index,
            corrector=control.name,
            field=String(control.axis),
            total_length_m=control.total_length_m,
        )
        for (index, control) in enumerate(controls)
    ]

    write_npy(joinpath(output_dir, "bpm_control_response.npy"), bpm_response)
    write_npy(joinpath(output_dir, "target_control_response.npy"), target_response)
    write_npy(joinpath(output_dir, "bpm_cumulative_maps.npy"), bpm_maps)
    write_npy(joinpath(output_dir, "target_cumulative_maps.npy"), target_maps)
    write_npy(joinpath(output_dir, "one_turn_map.npy"), one_turn_map)
    write_npy(joinpath(output_dir, "nominal_bpm_orbits.npy"), nominal_bpm_orbits)
    write_npy(joinpath(output_dir, "nominal_target_orbits.npy"), nominal_target_orbits)
    write_rows(joinpath(output_dir, "bpm_locations.csv"), bpm_rows)
    write_rows(joinpath(output_dir, "target_locations.csv"), target_rows)
    write_rows(joinpath(output_dir, "control_inventory.csv"), control_rows)
    metadata = Dict(
        "format" => "cesr-finite-bpm-local-orbit-model-v1",
        "date" => string(Dates.today()),
        "engine" => "SciBmad latest-lattice nominal RF-on response and transport",
        "lattice" => LATEST_LATTICE,
        "bpm_count" => length(detectors),
        "target_count" => length(sextupoles),
        "corrector_count" => length(controls),
        "coordinate_order" => "x, px, y, py, z, pz",
        "location_convention" => "first element-entry occurrence by base name",
        "nominal_orbit_anchor" => "latest-lattice RF-on SciBmad closed orbit at BPM and target entries",
        "response_descriptor" => "Descriptor(6, 2, $(length(controls)), 1)",
        "transport_descriptor" => "Descriptor(6, 1)",
        "closed_orbit_seconds" => closed_timed.time,
        "response_seconds" => response_timed.time,
        "transport_seconds" => transport_timed.time,
    )
    write_metadata(joinpath(output_dir, "model_metadata.toml"), metadata)
    println("Wrote local-orbit models to $output_dir")
    println("BPM / target / corrector counts: $(length(detectors)) / $(length(sextupoles)) / $(length(controls))")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
